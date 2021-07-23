#!/usr/bin/env stack
{- stack
  script
  --optimize
  --resolver lts-14.20
  --package async
  --package bytestring
  --package containers
  --package lens
  --package managed
  --package system-filepath
  --package text
  --package text-show
  --package turtle
  --package foldl
  --package mwc-random
  --package vector
  --package zlib
  --package safe
-}

{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}

import Codec.Compression.GZip (decompress, compress)
import Control.Monad.Managed (MonadManaged)
import Data.Maybe (fromMaybe, fromJust, isJust)
import TextShow (showt)
import Safe (atMay)
import Turtle
import qualified Turtle.Bytes as TB
import qualified Control.Foldl as Fold
import qualified Control.Lens as L
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.Map.Strict as Map
import qualified Data.IntMap.Strict as IMap
import qualified Data.Set as Set
import qualified Data.Text as T
import qualified Data.Text.Read as T
import qualified Data.Text.Encoding as T
import qualified Data.Vector as V
import qualified Filesystem.Path.CurrentOS as FP
import qualified System.Random.MWC as R
import qualified System.Random.MWC.Distributions as R

-- | Safe Shell Text to Shell Line.
toLines :: Shell T.Text -> Shell Line
toLines = fmap unsafeTextToLine . join . reduce (Fold.Fold step initial extract)
  where
    extract (lineAcc, stream) = stream <|> return lineAcc
    initial = (mempty, empty)
    step acc "" = acc
    step (lineAcc, stream) x
      | splitNum < 1 = (lineAcc, stream)
      | splitNum == 1 = (lineAcc <> x, stream)
      | splitNum > 1 = (last split, stream <|> (select . ((lineAcc <> head split) :) . tail . init $ split))
      where
        split = T.splitOn "\n" x
        splitNum = length split

-- | Uncompress file stream.
uncompressFileStream' :: FP.FilePath -> Shell Line
uncompressFileStream' file = do
  temp <- mktempfile "." "temp.tsv"
  _ <- shells ("cat " <> format fp file <> " | gzip -d -c > " <> format fp temp) mempty
  input temp

uncompressFileStream :: T.Text -> FP.FilePath -> Shell Line
uncompressFileStream label file = join
                                . liftIO
                                $ grep (invert (contains "chrY"))
                                . fmap ( unsafeTextToLine
                                       . addBarcode
                                       . T.decodeUtf8
                                       . BL.toStrict
                                       )
                                . select
                                . BL.lines
                                . decompress
                              <$> (BL.readFile . T.unpack . format fp $ file)
  where
    addBarcode = T.intercalate "\t"
               . L.over (L.ix 3) (mappend label . mappend "#")
               . T.splitOn "\t"

-- | Subsample a cellranger matrix.
writePop :: FP.FilePath -> Set.Set T.Text -> (Shell T.Text) -> IO ()
writePop outputPath barcodes =
  TB.output outputPath
    . fmap (BL.toStrict . compress . BL.fromStrict . T.encodeUtf8 . (\x -> T.append x "\n"))
    . mfilter (fromMaybe False . fmap (flip Set.member barcodes) . flip atMay 3 . T.splitOn "\t")

-- Get the intensity of each peak from the peak count matrix.
getPeakCountsFromMat :: FP.FilePath -> FP.FilePath -> Shell (Map.Map T.Text Int)
getPeakCountsFromMat peakPath matPath = do
  let toInt :: T.Text -> Int
      toInt = truncate . either error fst . T.double
      inputGZ :: FP.FilePath -> Shell T.Text
      inputGZ = fmap lineToText
              . Turtle.toLines
              . TB.toUTF8
              . TB.decompress (TB.WindowBits 31)
              . TB.input

  peaks <- reduce Fold.list . inputGZ $ peakPath
  matValsMap <- liftIO . fmap (Map.fromListWith (+) . drop 2)
              . (reduce Fold.list :: Shell (Int, Int) -> IO [(Int, Int)])
              . fmap (((\[x, _, z] -> (toInt x, toInt z)) . T.words) :: T.Text -> (Int, Int))
              . mfilter (not . T.null)
              . inputGZ
              $ matPath

  let peakMap :: IMap.IntMap T.Text
      peakMap = IMap.fromList $ zip [1..] peaks
      res = Map.mapKeys
              (fromMaybe (error "No key.") . flip IMap.lookup peakMap)
              (matValsMap :: Map.Map Int Int)

  return res

-- | Get peaks and matrix files from input peak count matrix.
writeOtherInfo :: FP.FilePath -> FP.FilePath -> Set.Set T.Text -> IO ()
writeOtherInfo outputPath matPath barcodes = sh $ do
  -- Create temporary file with cell whitelist
  tempWhitelistFile <- mktempfile "/tmp" "whitelist.csv"
  output tempWhitelistFile . fmap unsafeTextToLine . select $ barcodes

  -- Run too-many-cells to get filtered matrix and peaks
  procs "too-many-cells" [ "matrix-output"
                         , "-m", format fp matPath
                         , "--cell-whitelist-file", format fp tempWhitelistFile
                         , "--mat-output", format fp outputPath
                         ]
    $ mempty

  -- Format peaks file
  peakCountMap <- getPeakCountsFromMat
                    (outputPath FP.</> "features.tsv.gz")
                    (outputPath FP.</> "matrix.mtx.gz")

  liftIO $ print (outputPath FP.</> "features.tsv.gz")

  output (outputPath FP.</> "peaks.bed")
    . fmap unsafeTextToLine
    . select
    . fmap (\(!x, !y) -> T.intercalate "\t" $ (take 3 $ T.words x >>= T.splitOn "_") `mappend` [showt y])
    . Map.toAscList
    $ peakCountMap

  return ()

-- | Get the input fragment files from a metadata file with a search term.
getInputFiles :: T.Text -> FP.FilePath -> IO [FP.FilePath]
getInputFiles s = reduce Fold.list
                . fmap (FP.fromText . fst . T.breakOn "," . lineToText)
                . grep (contains $ text s)
                . input

-- | Get n random barcodes from label file with a search term and seed.
getBarcodes :: FP.FilePath -> T.Text -> T.Text -> Int -> Int -> IO (Set.Set T.Text)
getBarcodes file label s seed n = do
  barcodes <- reduce Fold.list
            . fmap (fst . T.breakOn "," . lineToText)
            . grep (contains (text s) <|> contains (text label)) -- Should remove header
            . input
            $ file

  g <- R.initialize (V.singleton $ fromIntegral seed)

  subsampled <-
    fmap (Set.fromList . take n . V.toList) . flip R.uniformShuffle g . V.fromList $ barcodes

  return subsampled

main :: IO ()
main = do
  let totalNum = 1000 -- Max sample size.
      percent = 0.005 -- Incremental percent for a rare population.
      maxPercent = 0.05 -- Maximum percent for a rare population.
      minPercent = 0.005 :: Double -- Minimum percent for a rare population.
      runs = 10 -- Number of runs
      seeds = [0..runs - 1]
      inc = round $ fromIntegral totalNum * percent -- Incremental number.
      minNum = round $ fromIntegral totalNum * minPercent -- Minimum number.
      maxNum = round $ fromIntegral totalNum * maxPercent -- Max number.
      outputRoot = "./random_b_cd8_treg/random"
      commonPop = [("B_Cells", "Memory B")]
      rarePop = [("Regulatory_T_Cells", "Treg"), ("Naive_CD8_T_Cells", "Naive CD8 T3")]
      labelFile = "../labels/labels_b_cd8_treg.csv"
      inputs = [
                 ("B_Cells", "../data/satpathy/GSM3722028_B_Cells_fragments.tsv.gz")
               , ("Naive_CD8_T_Cells", "../data/satpathy/GSM3722037_Naive_CD8_T_Cells_fragments.tsv.gz")
               , ("Regulatory_T_Cells", "../data/satpathy/GSM3722030_Regulatory_T_Cells_fragments.tsv.gz")
               ]

  print inputs

  sh $ do
    let stream = fmap lineToText . msum . fmap (uncurry uncompressFileStream) $ inputs

    sh $ do
      rareNum <- select [inc,(inc * 2)..maxNum]

      let commonNum = totalNum - (rareNum * 2)
          pops :: [(Int, (T.Text, T.Text))]
          pops = zip [commonNum, rareNum, rareNum] $ commonPop <> rarePop

      seed <- select seeds
      barcodes <- liftIO
                . fmap mconcat
                . mapM (\(popNum, (popLabel, pop)) -> getBarcodes labelFile popLabel pop seed popNum)
                $ pops

      let outputFragmentsPath = FP.fromText
                    $ format fp outputRoot
                    <> "_seed_"
                    <> showt seed
                    <> "/"
                    <> showt commonNum
                    <> "_fragments.tsv.gz"
          outputOtherInfoPath = FP.fromText
                    $ format fp outputRoot
                    <> "_seed_"
                    <> showt seed
                    <> "/"
                    <> showt commonNum
                    <> "_mat"

      liftIO . print $ [seed, commonNum, rareNum]

      mktree . FP.directory $ outputFragmentsPath
      liftIO $ writePop outputFragmentsPath barcodes stream
      mktree . FP.directory $ outputOtherInfoPath
      liftIO $ writeOtherInfo outputOtherInfoPath "./satpathy/mat" barcodes
