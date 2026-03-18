{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Data.GraphSpec (spec) where

import Data.Graph.Base
import qualified Data.HashMap.Strict as HM
import qualified Data.HashSet as HS
import Data.Hashable (Hashable)
import Test.Hspec
import Test.QuickCheck

instance (Arbitrary a, Arbitrary b, Hashable a, Eq a) => Arbitrary (Graph a b) where
    arbitrary = do
        nodes <- arbitrary
        edges <- arbitrary
        return $ mkDiGraph nodes edges

spec :: Spec
spec = do
    describe "mkUniGraph" $ do
        it "creates a symmetric graph" $ do
            let g = mkUniGraph [] [((1, 2), 1.0)] :: Graph Int Double
            hasEdge g 1 2 `shouldBe` True
            hasEdge g 2 1 `shouldBe` True

    describe "mkDiGraph" $ do
        it "creates a directed graph" $ do
            let g = mkDiGraph [] [((1, 2), 1.0)] :: Graph Int Double
            hasEdge g 1 2 `shouldBe` True
            hasEdge g 2 1 `shouldBe` False

    describe "invertGraph" $ do
        it "is its own inverse" $
            property $
                \(g :: Graph Int Int) -> invertGraph (invertGraph g) == g

        it "swaps edges" $
            property $
                \(g :: Graph Int Int) ->
                    let gt = invertGraph g
                        edges = graphToList g
                     in all (\((u, v), w) -> getGraphEdge gt (v, u) == Just w) edges

    describe "getSubGraph" $ do
        it "contains only requested nodes" $
            property $
                \(g :: Graph Int Int) (ns :: [Int]) ->
                    let sg = getSubGraph g ns
                        nodes = HM.keys (graph sg)
                     in all (`elem` ns) nodes

    describe "connComp" $ do
        it "partitions all nodes in an undirected graph" $
            property $
                \(nodes :: [Int]) (edges :: [((Int, Int), Int)]) ->
                    let g = mkUniGraph nodes edges
                        comps = connComp g
                        allCompNodes = HS.unions comps
                        allGraphNodes = HS.fromList $ HM.keys (graph g)
                     in allCompNodes == allGraphNodes
