{- clustering-contingency
Gregory W. Schwartz

Blah's the blah in the blah.
-}

{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeOperators     #-}

module Main where

-- Remote
import Options.Generic
import Frames
import Frames.CSV
import Pipes ((>->), runEffect)
import Pipes.Prelude (stdoutLn)
import Pipes.Prelude.Text (fromHandleLn)
import System.IO (stdin, openFile, IOMode (..))
import qualified Data.Text as T
import qualified Pipes.Prelude as P

-- Local
import ClusteringContingency.Types
import ClusteringContingency.Contingency

-- | Command line arguments
data Options = Options { input  :: Maybe String
                               <?> "(FILE) The input file in the format \"label,cluster\"."
                       }
               deriving (Generic)

modifiers :: Modifiers
modifiers = lispCaseModifiers { shortNameModifier = firstLetter }

instance ParseRecord Options where
    parseRecord = parseRecordWithModifiers modifiers

main :: IO ()
main = do
    opts <- getRecord "clustering-contingency, Gregory W. Schwartz.\
                      \ Gets the contingency table for clustered data."

    hIn <- maybe (return stdin) (flip openFile ReadMode) . unHelpful . input $ opts

    frame <- runSafeT $ inCoreAoS $ fromHandleLn hIn >-> pipeTable

    let clusterMap = getClusterMap frame

    runEffect $ produceCSV (getContingency clusterMap) >-> stdoutLn

    return ()
