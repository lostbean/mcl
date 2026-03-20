{-# LANGUAGE ScopedTypeVariables #-}

module Data.Graph.MarkovSpec (spec) where

import Data.Graph.Base
import Data.Graph.Markov
import Data.Graph.Sparse
import qualified Data.List as L
import Data.TestOrphans ()
import qualified Data.Vector.Unboxed as U
import Test.Hspec
import Test.QuickCheck

spec :: Spec
spec = do
    let myCfg = defaultMCL{selfLoop = 1.0, inflation = 2.0}
    describe "runMCL" $ do
        it "clusters a simple disconnected graph" $ do
            let g = mkUniGraph [] [((1, 2), 1.0), ((3, 4), 1.0)] :: Graph Int Double
                clusters = runMCL myCfg g
            length clusters `shouldBe` 2
            let sortedClusters = L.sort (map L.sort clusters)
            sortedClusters `shouldBe` [[1, 2], [3, 4]]

        it "clusters a simple graph with two dense components" $ do
            -- Two triangles connected by a single bridge
            let g =
                    mkUniGraph
                        []
                        [ ((1, 2), 1.0)
                        , ((2, 3), 1.0)
                        , ((3, 1), 1.0)
                        , ((4, 5), 1.0)
                        , ((5, 6), 1.0)
                        , ((6, 4), 1.0)
                        , ((3, 4), 0.1) -- bridge
                        ] ::
                        Graph Int Double
                clusters = runMCL myCfg g
            length clusters `shouldBe` 2
            let sortedClusters = L.sort (map L.sort clusters)
            sortedClusters `shouldBe` [[1, 2, 3], [4, 5, 6]]

    describe "normalizeVecS" $ do
        it "produces a vector that sums to 1 (if not empty)" $
            property $
                \(v :: VecS Double) ->
                    let v' = normalizeVecS v
                        sumV = U.foldl' (\acc (_, x) -> acc + x) 0 (vecS v')
                     in if U.null (vecS v) then sumV == 0 else abs (sumV - 1.0) < 1e-9

    describe "pruneVecS" $ do
        it "removes elements below threshold in FixPrune" $
            property $
                \(v :: VecS Double) (Positive (t :: Double)) ->
                    let v' = pruneVecS (FixPrune t) v
                        elems = U.toList (vecS v')
                     in all (\(_, x) -> x > t) elems
