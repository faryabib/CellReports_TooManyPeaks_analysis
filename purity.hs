#!/usr/bin/env stack
{- stack
  runghc
  --resolver lts-10.0
  --package async
  --package async-pool-0.9.1
  --package lens
  --package system-filepath
  --package safe
  --package text
  --package text-show
  --package turtle
  --package foldl
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
import qualified Control.Concurrent.Async.Pool as Async
import qualified Control.Foldl as Fold
import qualified Control.Lens as L
import qualified Data.List as List
import qualified Data.Text as T
import qualified Filesystem.Path.CurrentOS as FP

data Clustering = TooLSA
                | Snap
                | CisTopic
                | CisTopicLouvain
                | Cicero
                | EpiScanpy
                | Cusanovich2018
                | Signac
                | Apec
                | PAGA
                deriving (Eq, Ord, Read, Show)

-- | Join clustering results with labels.
assignLabels :: FP.FilePath -> Shell Line -> Shell Line
assignLabels labelPath = inproc "csvcut" ["-c", "label,cluster,item"]
                       . inproc "csvjoin" ["-c", "item", format fp labelPath, "-"]

clusterTooLSA :: T.Text -> [T.Text] -> IO ()
clusterTooLSA labelPath files = sh $ do
    let csvCommands = concat $ zipWith (\x y -> [x, y]) (repeat "--matrix-path") files
        outputFile = "too-many-cells-lsa_clustering.csv"
        args = ["make-tree"]
            <> csvCommands
            <> [ "-o", "too-many-cells-lsa_clustering"
               , "-l", labelPath
               , "--cell-whitelist-file", "/path/to/whitelist_groups_cd34_b.csv"
               , "--draw-mark", "MarkModularity"
               , "--filter-thresholds", "(1000, 1)"
               , "--binwidth", "5000"
               , "--normalization", "NoneNorm"
               , "--lsa", "50"
               , "--blacklist-regions-file", "/path/to/Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz" -- Depending on data set
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
clusterGenericR command clustering labelPath fullPaths = sh $ do
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
        args _ = command
               : labelPath
               : "0"
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
clusterGenericPy command clustering labelPath fullPaths = sh $ do
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

main :: IO ()
main = do
    args <- liftIO getArgs

    let root = "/path/to/satpathy/"
        labelPath = "/path/to/labels_groups_cd34_b.csv" :: T.Text
        clusterType    = maybe TooLSA read . headMay $ args
        cluster TooLSA      = clusterTooLSA
        cluster Snap        = clusterGenericR "snapatac.R" "snap"
        cluster CisTopic    = clusterGenericR "cistopic.R" "cistopic"
        cluster CisTopicLouvain = clusterGenericR "cistopic_louvain.R" "cistopicLouvain"
        cluster Cicero      = clusterGenericR "cicero.R" "cicero"
        cluster Cusanovich2018 = clusterGenericR "cusanovich2018.R" "cusanovich2018"
        cluster Signac      = clusterGenericR "signac.R" "signac"
        cluster Apec        = clusterGenericPy "run_apec.py" "apec"
        cluster EpiScanpy   = clusterGenericPy "run_episcanpy.py" "episcanpy"
        cluster PAGA   = clusterGenericPy "run_episcanpy_paga.py" "paga"
        algorithms = [ TooLSA
                     , Snap
                     , CisTopic
                     , Cicero
                     , EpiScanpy
                     , CisTopicLouvain
                     , PAGA
                     , Apec
                     , Signac
                     ]

    inFiles <- fold ( fmap lineToText
             . grep (contains $ "CD34")
             . fmap (unsafeTextToLine . format fp)
             $ ls (fromText root)
             )
             Fold.list

    mapConcurrently_ (\ a -> (cluster a) labelPath ["/path/to/CD34_mat"])
      $ algorithms
