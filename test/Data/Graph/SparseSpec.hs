{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Data.Graph.SparseSpec (spec) where

import Data.Graph.Base
import Data.Graph.Sparse
import qualified Data.HashMap.Strict as HM
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Test.Hspec
import Test.QuickCheck

instance Arbitrary (VecS Double) where
    arbitrary = do
        pairs <- arbitrary
        return $ mkVecS pairs

instance Arbitrary (CRS Double) where
    arbitrary = do
        rows <- arbitrary
        return $ CRS $ V.fromList rows

spec :: Spec
spec = do
    describe "VecS" $ do
        it "adds two sparse vectors correctly" $ do
            let v1 = mkVecS [(1, 1.0), (3, 2.0)]
                v2 = mkVecS [(1, 0.5), (2, 3.0)]
                res = addVV v1 v2
            vecS res `shouldBe` U.fromList [(1, 1.5), (2, 3.0), (3, 2.0)]

        it "multiplies two sparse vectors (inner product) correctly" $ do
            let v1 = mkVecS [(1, 1.0), (3, 2.0)]
                v2 = mkVecS [(1, 0.5), (3, 1.5), (4, 4.0)]
            multVV v1 v2 `shouldBe` (1.0 * 0.5 + 2.0 * 1.5)

        it "multiplies a sparse vector by a constant correctly" $ do
            let v = mkVecS [(1, 1.0), (3, 2.0)]
                res = multKV 2.0 v
            vecS res `shouldBe` U.fromList [(1, 2.0), (3, 4.0)]

        it "addVV is commutative" $
            property $
                \(v1 :: VecS Double) (v2 :: VecS Double) ->
                    addVV v1 v2 == addVV v2 v1

        it "multVV is commutative" $
            property $
                \(v1 :: VecS Double) (v2 :: VecS Double) ->
                    multVV v1 v2 == multVV v2 v1

    describe "CRS" $ do
        it "converts Graph to CRS and back" $ do
            let g = mkDiGraph [] [((1, 2), 1.0), ((2, 3), 0.5)] :: Graph Int Double
                m = graphToCRS g
                g' = crsToGraph m
            getGraphEdge g (1, 2) `shouldBe` Just 1.0
            getGraphEdge g (2, 3) `shouldBe` Just 0.5
            getGraphEdge g' (1, 2) `shouldBe` Just 1.0
            getGraphEdge g' (2, 3) `shouldBe` Just 0.5

        it "transposes a sparse matrix correctly" $ do
            let g = mkDiGraph [] [((1, 2), 1.0), ((2, 3), 0.5)] :: Graph Int Double
                m = graphToCRS g
                mt = transpose m
                gt = crsToGraph mt
            getGraphEdge gt (2, 1) `shouldBe` Just 1.0
            getGraphEdge gt (3, 2) `shouldBe` Just 0.5

        it "transpose is its own inverse (non-zero entries)" $
            property $
                \(m :: CRS Double) ->
                    let mt = transpose m
                        mtt = transpose mt
                        nonZeros mat = HM.fromList $ graphToList (crsToGraph mat)
                     in nonZeros mtt == nonZeros m

        it "multiplies two sparse matrices (multMM) correctly" $ do
            -- [ 0 1 ]   [ 0 2 ]   [ 0 3 ]
            -- [ 0 0 ] * [ 3 0 ] = [ 0 0 ]
            let g1 = mkDiGraph [] [((0, 1), 1.0)] :: Graph Int Double
                g2 = mkDiGraph [] [((1, 0), 3.0), ((0, 1), 2.0)] :: Graph Int Double
                m1 = graphToCRS g1
                m2 = graphToCRS g2
                m3 = multMM m1 m2
                g3 = crsToGraph m3
            getGraphEdge g3 (0, 0) `shouldBe` Just 3.0
            getGraphEdge g3 (0, 1) `shouldBe` Nothing

        it "multiplies two sparse matrices (multMMsmrt) correctly" $ do
            let g1 = mkDiGraph [] [((0, 1), 1.0)] :: Graph Int Double
                g2 = mkDiGraph [] [((1, 0), 3.0), ((0, 1), 2.0)] :: Graph Int Double
                m1 = graphToCRS g1
                m2 = graphToCRS g2
                m3 = multMMsmrt m1 m2
                g3 = crsToGraph m3
            getGraphEdge g3 (0, 0) `shouldBe` Just 3.0
