#!/usr/bin/env stack
{- stack
  script
  --resolver lts-15.3
  --package async
  --package lens
  --package exceptions
  --package system-filepath
  --package safe
  --package text
  --package text-show
  --package turtle
  --package foldl
  --package criterion
  --package criterion-measurement
-}

{-# LANGUAGE OverloadedStrings #-}

import Control.Concurrent.Async (mapConcurrently_)
import Data.Function (on)
import Data.Monoid ((<>))
import Safe
import System.Environment (getArgs)
import TextShow (showt)
import Turtle
import Turtle.Line
import Control.Monad.Catch
import qualified Control.Foldl as Fold
import qualified Control.Lens as L
import qualified Data.List as List
import qualified Data.Text as T
import qualified Filesystem.Path.CurrentOS as FP
import qualified Criterion.Main as C
import qualified Criterion.Types as C
import qualified Criterion.Measurement as CM
import qualified Criterion.Measurement.Types as CM hiding (measure)

data Clustering = TooLSA
                | Snap
                | CisTopic
                | Cicero
                | Signac
                | Apec
                | EpiScanpy
                | Cusanovich2018
                | CisTopicStandard
                deriving (Eq, Ord, Read, Show)

-- | Shell with failed clustering
handleSh :: (MonadCatch m, MonadIO m) => T.Text -> m () -> m ()
handleSh c = handle (\e -> append "log" $ select [unsafeTextToLine $ c <> " failed: " <> T.pack (show (e :: ExitCode))])

-- | Join clustering results with labels.
assignLabels :: FP.FilePath -> Shell Line -> Shell Line
assignLabels labelPath = inproc "csvcut" ["-c", "label,cluster,item"]
                       . inproc "csvjoin" ["-c", "item", format fp labelPath, "-"]

clusterTooLSA :: T.Text -> [T.Text] -> IO ()
clusterTooLSA labelPath files = sh . handleSh "TooLSA" $ do
    let csvCommands = concat $ zipWith (\x y -> [x, y]) (repeat "--matrix-path") files
        outputFile = "too-many-cells-lsa_clustering.csv"
        args = ["make-tree"]
            <> csvCommands
            <> [ "-o", "too-many-cells-lsa_clustering"
               , "-l", labelPath
               , "--draw-mark", "MarkModularity"
               , "--filter-thresholds", "(1000, 1)"
               , "--binwidth", "5000"
               , "--normalization", "NoneNorm"
               , "--lsa", "50"
               , "--blacklist-regions-file", "/home/gw/research/genomes/Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz" -- Depending
               , "+RTS", "-N6"
               ]

    liftIO . print $ csvCommands

    tmp <- using $ mktempfile "/tmp" "tooOutput.csv"

    output tmp
        . sed ("cell,cluster,path" *> pure "item,cluster,path")
        . inproc "too-many-cells" args
        $ mempty

    output (FP.fromText outputFile)
      . assignLabels (FP.fromText labelPath)
      . input
      $ tmp

clusterGenericR :: T.Text -> T.Text -> T.Text -> [T.Text] -> IO ()
clusterGenericR command clustering labelPath fullPaths = sh . handleSh clustering $ do
    let files "snap" = fullPaths
        files _ = fmap (T.replace "_fragments.tsv.gz" "_mat") fullPaths
        outputFile = clustering <> "_clustering.csv"
        args "cicero" = command
                      : labelPath
                      : "0"
                      : (clustering <> "_clustering")
                      : files clustering
        args "cistopic" = command
                        : labelPath
                        : "0"
                        : (clustering <> "_clustering")
                        : files clustering
        args "cusanovich2018" = command
                        : labelPath
                        : "0"
                        : (clustering <> "_clustering")
                        : files clustering
        args "cistopic_standard" = command
                        : labelPath
                        : "0"
                        : (clustering <> "_clustering")
                        : files clustering
        args "signac" = command
                        : labelPath
                        : "0"
                        : (clustering <> "_clustering")
                        : files clustering
        args _ = command
               : labelPath
               : (clustering <> "_clustering")
               : files clustering
        upper (Header x) = x
        upper (Row _ x)  = unsafeTextToLine
                         . T.intercalate ","
                         . L.over (L.ix 0) T.toUpper
                         . T.splitOn ","
                         . lineToText
                         $ x
        postProcess "snap" = fmap upper . header
        postProcess _ = id

    liftIO . print $ args clustering

    tmp <- mktempfile "/tmp" $ clustering <> "Output.csv"
    tmpLabel <- mktempfile "/tmp" $ clustering <> "Label.csv"

    output tmp
        . sed (",NA" *> pure "")
        . inproc "Rscript" (args clustering)
        $ mempty

    output tmpLabel
        . postProcess clustering -- Need uppercase for snapatac.
        . input
        $ fromText labelPath

    output (FP.fromText outputFile)
      . inproc "csvcut" ["-c", "label,cluster,item"]
      . inproc "csvjoin" ["-c", "item", format fp tmpLabel, format fp tmp]
      $ mempty

clusterGenericPy :: T.Text -> T.Text -> T.Text -> [T.Text] -> IO ()
clusterGenericPy command clustering labelPath fullPaths = sh . handleSh clustering $ do
    let files = fmap (T.replace "_fragments.tsv.gz" "_mat") fullPaths
        outputFile = clustering <> "_clustering.csv"
        args = command
             : labelPath
             : (clustering <> "_clustering")
             : files
        pyCommand "apec" = "/path/to/apec/python3"
        pyCommand _ = "/path/to/python3"

    liftIO . print $ args

    tmp <- mktempfile "/tmp" $ clustering <> "Output.csv"
    tmpLabel <- mktempfile "/tmp" $ clustering <> "Label.csv"

    output tmp
        . sed (",NA" *> pure "")
        . inproc (pyCommand clustering) args
        $ mempty

    output tmpLabel
        . input
        $ fromText labelPath

    output (FP.fromText outputFile)
      . inproc "csvcut" ["-c", "label,cluster,item"]
      . inproc "csvjoin" ["-c", "item", format fp tmpLabel, format fp tmp]
      . sed ("#" <> plus (notChar ',') *> "")  -- To remove samples (unfortunately)
      $ mempty

mainGo :: Clustering -> IO ()
mainGo clusterType = do
    let root = "/path/to/buenrostro/"
        labelPath = "/path/to/labels.csv" :: T.Text
        cluster TooLSA      = clusterTooLSA
        cluster TooLSAPeaks = clusterTooLSAPeaks
        cluster Snap        = clusterGenericR "snapatac.R" "snap"
        cluster CisTopic    = clusterGenericR "cistopic.R" "cistopic"
        cluster CisTopicStandard = clusterGenericR "cistopic_louvain.R" "cistopicLouvain"
        cluster Cicero      = clusterGenericR "cicero.R" "cicero"
        cluster EpiScanpy   = clusterGenericPy "run_episcanpy.py" "episcanpy"
        cluster Signac      = clusterGenericR "signac.R" "signac"
        mat = [root <> "mat"]

    cluster clusterType labelPath mat
    return ()

main :: IO ()
main = do
  let clusterTypes = [ ("Snap", Snap)
                     , ("TooManyCellsLSA", TooLSA)
                     , ("Cicero", Cicero)
                     , ("EpiScanpy", EpiScanpy)
                     , ("CisTopic", CisTopic)
                     , ("CisTopicLouvain", CisTopicLouvain)
                     , ("Signac", Signac)
                     , ("Apec", Apec)
                     , ("Cusanovich2018", Cusanovich2018)
                     ]
  print "clustering,time"
  mapM_ (\(title, c) -> putStrLn . (\x -> title <> "," <> show x) . CM.secs . CM.measTime . fst =<< (flip CM.measure 1 . C.nfIO $ mainGo c))
    . mconcat
    . replicate 3
    $ clusterTypes
