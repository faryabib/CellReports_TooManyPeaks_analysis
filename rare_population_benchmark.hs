#!/usr/bin/env stack
{- stack
  script
  --resolver lts-14.16
  --package async
  --package system-filepath
  --package lens
  --package safe
  --package text
  --package text-show
  --package turtle
  --package foldl
-}

{-# LANGUAGE OverloadedStrings #-}

import Data.Monoid ((<>))
import Safe
import System.Environment (getArgs)
import TextShow (showt)
import Turtle
import Turtle.Line
import qualified Control.Foldl as Fold
import qualified Control.Lens as L
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
                deriving (Eq, Ord, Read, Show)

clusterTooLSA :: T.Text -> T.Text -> [T.Text] -> IO ()
clusterTooLSA outputLabel labelPath files = sh $ do
    let csvCommands = concat $ zipWith (\x y -> [x, y]) (repeat "--matrix-path") files
        outputFile = outputLabel <> "_too-many-cells-lsa_clustering.csv"
        args = ["make-tree"]
            <> csvCommands
            <> [ "-l", labelPath
               , "-o", outputLabel <> "_too-many-cells-lsa_clustering"
               , "--draw-mark", "MarkModularity"
               , "--filter-thresholds", "(1, 1)"
               , "--normalization", "NoneNorm"
               , "--lsa", "50"
               , "--binwidth", "5000"
               , "--blacklist-regions-file", "/path/to/Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz" -- Depending on data set
               , "--dendrogram-output", "dendrogram.svg"
               , "--draw-collection", "PieRing"
               , "+RTS", "-N2"
               ]

    liftIO . print $ csvCommands

    tmp <- using $ mktempfile "/tmp" "tooOutput.csv"

    output tmp
        . sed ("cell,cluster,path" *> pure "item,cluster,path")
        . inproc "too-many-cells" args
        $ mempty

    output (FP.fromText outputFile)
      . inproc "csvcut" ["-c", "label,cluster,item"]
      . inproc "csvjoin" ["-c", "item", labelPath, format fp tmp]
      $ mempty

clusterGenericR :: T.Text -> T.Text -> T.Text -> T.Text -> [T.Text] -> IO ()
clusterGenericR command clustering outputLabel labelPath files = sh $ do
    let outputFile = outputLabel <> "_" <> clustering <> "_clustering.csv"
        args _ = command
             : labelPath
             : "0"
             : (outputLabel <> "_" <> clustering <> "_clustering")
             : files

    liftIO . print $ args clustering

    tmp <- mktempfile "/tmp" $ clustering <> "Output.csv"
    tmpLabel <- mktempfile "/tmp" $ clustering <> "Label.csv"

    output tmp
        . sed (",NA" *> pure "")
        . inproc "Rscript" (args clustering)
        $ mempty

    output tmpLabel
        . input
        $ fromText labelPath

    output (FP.fromText outputFile)
      . inproc "csvcut" ["-c", "label,cluster,item"]
      . inproc "csvjoin" ["-c", "item", format fp tmpLabel, format fp tmp]
      $ mempty

clusterGenericPy :: T.Text -> T.Text -> T.Text -> T.Text -> [T.Text] -> IO ()
clusterGenericPy command clustering outputLabel labelPath fullPaths = sh $ do
    let outputFile = outputLabel <> "_" <> clustering <> "_clustering.csv"
        args = command
             : labelPath
             : (outputLabel <> "_" <> clustering <> "_clustering")
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
      $ mempty

main :: IO ()
main = sh $ do
    args <- liftIO getArgs

    let totalNum       = 1000 -- Max sample size.
        percent        = 0.005 -- Incremental percent for a rare population.
        maxPercent     = 0.05 -- Maximum percent for a rare population.
        runs = 10 :: Int -- Number of runs
        seeds = [0..runs - 1]
        inc            = round $ fromIntegral totalNum * percent -- Incremental number.
        maxNum         = round $ fromIntegral totalNum * maxPercent -- Max number.
        outputRoot     = "/path/to/random_b_cd8_treg/random" :: T.Text
        labelPath = "/path/to/labels_b_cd8_treg.csv"
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

    sh $ do -- Common population.
        rareNum <- select [inc,(inc * 2)..maxNum]
        seed <- select seeds

        let commonNum = totalNum - (2 * rareNum) :: Int
            dataSet = outputRoot <> "_seed_" <> showt seed <> "/" <> showt commonNum <> "_mat"
            outputLabel = showt commonNum <> "_seed_" <> showt seed

        liftIO . (cluster clusterType) outputLabel labelPath $ [dataSet]
