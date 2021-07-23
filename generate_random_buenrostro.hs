#!/usr/bin/env stack
{- stack
  script
  --resolver lts-14.20
  --package async
  --package containers
  --package foldl
  --package system-filepath
  --package text
  --package text-show
  --package turtle
-}

{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}

import Data.Maybe (fromMaybe, fromJust, isJust)
import TextShow (showt)
import Turtle
import qualified Control.Foldl as Fold
import qualified Data.Map.Strict as Map
import qualified Data.IntMap.Strict as IMap
import qualified Data.Text as T
import qualified Data.Text.Read as T
import qualified Filesystem.Path.CurrentOS as FP
import qualified Turtle.Bytes as TB

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
writeOtherInfo :: FP.FilePath -> IO ()
writeOtherInfo matPath = sh $ do
  -- Format peaks file
  peakCountMap <- getPeakCountsFromMat
                    (matPath FP.</> "features.tsv.gz")
                    (matPath FP.</> "matrix.mtx.gz")

  output (matPath FP.</> "peaks.bed")
    . fmap unsafeTextToLine
    . select
    . fmap (\(!x, !y) -> T.intercalate "\t" $ (take 3 $ T.words x >>= T.splitOn "_" >>= T.splitOn ":" >>= T.splitOn "-") `mappend` [showt y])
    . Map.toAscList
    $ peakCountMap

  return ()

-- | Subsample a cellranger matrix.
writePop :: Int -> Int -> T.Text -> T.Text -> IO ()
writePop seed n input output = sh $ do
    procs
        "Rscript"
        [ "/path/to/subsample_cellranger_mat.R"
        , showt seed
        , showt n
        , input
        , output
        ]
        mempty

main :: IO ()
main = sh $ do
    let totalNum = 500 :: Int -- Max sample size.
        percent = 0.01 :: Double -- Incremental percent for a rare population.
        maxPercent = 0.1 :: Double -- Maximum percent for a rare population.
        minPercent = 0.01 :: Double -- Minimum percent for a rare population.
        runs = 10 :: Int -- Number of runs
        seeds = [0..runs - 1] :: [Int]
        inc = round $ fromIntegral totalNum * percent -- Incremental number.
        minNum = round $ fromIntegral totalNum * minPercent -- Minimum number.
        maxNum = round $ fromIntegral totalNum * maxPercent -- Max number.
        inputRoot = "/path/to/buenrostro/"
        outputRoot = "/path/to/buenrostro/CMP_mono_pDC/random"
        commonPop = ["CMP"]
        rarePop = ["mono", "pDC"]

    sh $ do -- Common population.
        dataSet <- select commonPop
        rareNum <- select . fmap (* 2) $ [minNum,(minNum + inc) ..maxNum]
        seed <- select seeds

        let commonNum = totalNum - rareNum
            input = inputRoot <> dataSet
            output = fromText $ outputRoot <> "_seed_" <> showt seed <> "/" <> showt commonNum <> "_" <> dataSet

        liftIO $ print $ unwords ["seed", show seed, "common", show commonNum]
        liftIO $ writePop seed commonNum input $ format fp output
        liftIO $ writeOtherInfo output

    sh $ do -- Rare population.

        dataSet <- select rarePop
        rareNum <- select [minNum,(minNum + inc) ..maxNum]
        seed <- select seeds

        let input = inputRoot <> dataSet
            output = fromText $ outputRoot <> "_seed_" <> showt seed <> "/" <> showt rareNum <> "_" <> dataSet

        liftIO $ writePop seed rareNum input $ format fp output
        liftIO $ writeOtherInfo output
