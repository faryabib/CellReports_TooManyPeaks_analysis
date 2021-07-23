{- Types
Gregory W. Schwartz

Collects the types used in the program
-}

{-# LANGUAGE StrictData #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE FlexibleContexts #-}

module ClusteringContingency.Types where

-- Remote
import Frames
import qualified Data.Map.Strict as Map
import qualified Data.Text as T

-- Local

declareColumn "cluster" ''T.Text
declareColumn "label" ''T.Text
declareColumn "label1" ''T.Text
declareColumn "label2" ''T.Text
declareColumn "count" ''Int
declareColumn "percent" ''Double

-- Basic
type Row = Record ["label" :-> T.Text, "cluster" :-> T.Text]
type OutputRow =
    Record ["label1" :-> T.Text, "label2" :-> T.Text, "count" :-> Int, "percent" :-> Double]
newtype ClusterN = ClusterN T.Text deriving (Eq, Ord)
newtype LabelN = LabelN {unLabelN :: T.Text } deriving (Eq, Ord)
newtype ClusterMap = ClusterMap
    { unClusterMap :: Map.Map ClusterN (Map.Map LabelN Int)
    }

-- Advanced
