{- ClusteringContingency.Contingency
Gregory W. Schwartz

Collects the functions pertaining to the counting of frequencies for the
contingency table.
-}

{-# LANGUAGE DataKinds #-}

module ClusteringContingency.Contingency
    ( getContingency
    , getClusterMap
    ) where

-- Remote
import Frames
import Control.Lens
import qualified Control.Foldl as L
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set

-- Local
import ClusteringContingency.Types

-- | Count the number of each label in each cluster.
countLabel :: L.Fold Row ClusterMap
countLabel =
    L.Fold
        (\ (ClusterMap acc) x ->
            ClusterMap
                $ Map.unionWith
                    (Map.unionWith (+))
                    acc
                . Map.singleton (ClusterN $ view cluster x)
                . Map.singleton (LabelN $ view label x)
                $ 1
        )
        (ClusterMap Map.empty)
        id

-- | Get the ClusterMap.
getClusterMap :: Frame Row -> ClusterMap
getClusterMap = L.fold countLabel

-- | Get the labels from a ClusterMap.
getLabels :: ClusterMap -> Set.Set LabelN
getLabels = Set.fromList . concat . Map.elems . Map.map Map.keys . unClusterMap

-- | Get the count for a comparison.
getCount :: ClusterMap -> LabelN -> LabelN -> Int
getCount cm l1 l2 =
    Map.foldl' (+) 0 . Map.map getClusterCount . unClusterMap $ cm
  where
    getClusterCount :: Map.Map LabelN Int -> Int
    getClusterCount lm =
        if l1 == l2
            then round
               $ ( (fromIntegral $ Map.findWithDefault 0 l1 lm)
                 * (fromIntegral $ Map.findWithDefault 0 l2 lm - 1)
                 )
               / 2
            else Map.findWithDefault 0 l1 lm * Map.findWithDefault 0 l2 lm

-- | Get the contingency table from the ClusterMap.
getContingency :: ClusterMap -> Frame OutputRow
getContingency cm = fmap (over rsubset applyPercentage) baseFrame
  where
    applyPercentage :: Record '[Percent] -> Record '[Percent]
    applyPercentage = mapMono (/ fromIntegral maxCount)
    maxCount = L.fold L.sum (fmap (view count) baseFrame)
    baseFrame = toFrame getCounts :: Frame OutputRow
    labels = Set.toList . getLabels $ cm
    getCounts =
        (\l1 l2 -> unLabelN l1 &: unLabelN l2 &: getCount cm l1 l2 &: fromIntegral (getCount cm l1 l2) &:Nil)
            <$> labels
            <*> labels
