{-# LANGUAGE FlexibleInstances #-}
{-# OPTIONS_GHC -Wno-orphans #-}

module Data.TestOrphans () where

import Data.Graph.Base
import Data.Graph.Sparse
import Data.Hashable (Hashable)
import qualified Data.Vector as V
import Test.QuickCheck

instance Arbitrary (VecS Double) where
    arbitrary = do
        pairs <- listOf $ do
            k <- arbitrary
            v <- abs <$> arbitrary `suchThat` (/= 0)
            return (k, v)
        return $ mkVecS pairs

instance Arbitrary (CRS Double) where
    arbitrary = do
        rows <- arbitrary
        return $ CRS $ V.fromList rows

instance (Arbitrary a, Arbitrary b, Hashable a, Eq a) => Arbitrary (Graph a b) where
    arbitrary = do
        nodes <- arbitrary
        edges <- arbitrary
        return $ mkDiGraph nodes edges
