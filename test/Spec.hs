module Main where

import qualified Data.Graph.MarkovSpec
import qualified Data.Graph.SparseSpec
import qualified Data.GraphSpec
import Test.Hspec

main :: IO ()
main = hspec spec

spec :: Spec
spec = do
    describe "Data.Graph" Data.GraphSpec.spec
    describe "Data.Graph.Sparse" Data.Graph.SparseSpec.spec
    describe "Data.Graph.Markov" Data.Graph.MarkovSpec.spec
