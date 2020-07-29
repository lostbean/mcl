{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE RecordWildCards  #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies     #-}

module Data.Graph.Sparse
       ( -- * Sparse Vector
         VecS (vecS)
       , unsafeMkVecS
       , mkVecS
       , addVV
       , multVV
       , multKV
         -- * Sparse Vector
       , CRS (CRS, crs)
       , multMV
       , multMM
       , multMMsmrt
       , transpose
       , graphToCRS
       , crsToGraph
       ) where

import Control.DeepSeq
import Control.Monad.ST
import Data.Function
import Data.Vector.Unboxed (Vector)
import qualified Data.HashMap.Strict          as HM
import qualified Data.Vector                  as V
import qualified Data.Vector.Mutable          as MV
import qualified Data.Vector.Unboxed          as U
import qualified Data.Vector.Generic          as G
import qualified Data.Vector.Unboxed.Mutable  as MU
import qualified Data.Vector.Algorithms.Intro as S

import Data.Graph

-- ================================== Sparse Vector ======================================

-- | Sparse vector with ordered indexes.
newtype VecS a = VecS { vecS :: Vector (Int, a) } deriving (Show, Eq)

instance (NFData a)=> NFData (VecS a) where
  rnf (VecS{..}) = rnf vecS

-- | Constructs a sparse vectors.
mkVecS :: [(Int, Double)] -> VecS Double
mkVecS = VecS . unsafeSort . U.fromList

-- | Unsafely constructs a sparse vectors without sorting. Therefore, sorted input is
-- expected.
unsafeMkVecS :: Vector (Int, Double) -> VecS Double
unsafeMkVecS = VecS

-- | Sparse vector addition.
addVV :: VecS Double -> VecS Double -> VecS Double
addVV v1 v2 = VecS $ U.reverse $ U.fromList $ go 0 0 []
  where
    l1 = U.length $ vecS v1
    l2 = U.length $ vecS v2
    go !i1 !i2 !acc
      | i1 >= l1 && i2 >= l2  = acc
      | k1 >  k2  = go i1 (i2+1) ((k2, x2) : acc)
      | k2 >  k1  = go (i1+1) i2 ((k1, x1) : acc)
      | otherwise = go (i1+1) (i2+1) ((k1, x1 + x2) : acc)
      where
        (k1, x1) = U.unsafeIndex (vecS v1) i1
        (k2, x2) = U.unsafeIndex (vecS v2) i2
{-# NOINLINE addVV #-}

-- | Sparse vector inner product (dot product).
multVV :: VecS Double -> VecS Double -> Double
multVV v1 v2 = go 0 0 0
  where
    l1 = U.length $ vecS v1
    l2 = U.length $ vecS v2
    go !i1 !i2 !acc
      | i1 >= l1  = acc
      | i2 >= l2  = acc
      | k1 >  k2  = acc `deepseq` go i1 (i2+1) acc
      | k2 >  k1  = acc `deepseq` go (i1+1) i2 acc
      | otherwise = acc `deepseq` go (i1+1) (i2+1) (acc + x1 * x2)
      where
        (k1, x1) = U.unsafeIndex (vecS v1) i1
        (k2, x2) = U.unsafeIndex (vecS v2) i2
{-# INLINE multVV #-}

-- | Multiplies a sparse vector by a constant value.
multKV :: Double -> VecS Double -> VecS Double
multKV k = VecS . G.map (\(i, x) -> (i, k * x)) . vecS
{-# INLINE multKV #-}

-- ================================== Sparse Matrix ======================================

-- | Sparse matrix representation by rows. All rows are stored (including the empty ones).
-- Each row is represented by @VecS@.
-- Ex.: [[],[(1,0.3), (3,0.1)], [], [(5,1.0)]]
newtype CRS a = CRS { crs :: V.Vector (VecS a) } deriving (Show, Eq)

instance (NFData a)=> NFData (CRS a) where
  rnf (CRS{..}) = rnf crs

graphToCRS :: Graph Int Double -> CRS Double
graphToCRS = CRS . func . HM.map (mkVecS . HM.toList) . graph
  where
    func xs = let
      n     = HM.foldlWithKey' (\acc k _ -> max k acc) 0 xs
      vinit = G.replicate (n+1) (VecS G.empty)
      in vinit G.// (HM.toList xs)

crsToGraph :: CRS Double -> Graph Int Double
crsToGraph = Graph . HM.fromList . V.toList . G.filter (HM.null . snd) . G.imap func . crs
  where func i x = (i, HM.fromList $ U.toList $ vecS x)

-- | Transposes a sparse matrix.
transpose :: CRS Double -> CRS Double
transpose CRS{..} = CRS $ V.map mkVecS mat
  where
    n :: Int
    n   = G.maximum $ G.map (G.foldl' (\acc x -> max acc (fst x)) 0 . vecS) crs
    mat = V.create $ do
      m <- MV.replicate (n+1) []
      let func i = do
            let (VecS row) = G.unsafeIndex crs i
            U.mapM_ (\(j, v) -> MV.read m j >>= (MV.write m j . ((i,v):))) row
      mapM_ func [0 .. G.length crs - 1]
      return m

-- | Multiply a sparse matrix (row based) by a columnar vector.
multMV :: CRS Double -> VecS Double -> VecS Double
multMV CRS{..} v = VecS . U.convert $ V.imap func crs
  where func i x = (i, multVV v x)

-- | Sparse matrix (row based) multiplication.
multMM :: CRS Double -> CRS Double -> CRS Double
multMM m1 m2 = CRS $ G.map func (crs m1)
  where
    cols = crs $ transpose m2
    func row1 = VecS $ foo cols
      where
        foo = G.filter ((/= 0) . snd) . G.convert .
              G.imap (\j col2 -> (j, multVV row1 col2))

-- | Sparse matrix (row based) multiplication.
multMMsmrt :: CRS Double -> CRS Double -> CRS Double
multMMsmrt m1 m2 = CRS $ G.map func (crs m1)
  where func = flatter . smartRowsGetter m2

-- | Super efficient sparse multiplication based on the SpGEMM algorithm as nicely explained
-- on: "An Efﬁcient GPU General Sparse Matrix-Matrix Multiplication for Irregular Data"
smartRowsGetter :: CRS Double -> VecS Double -> V.Vector (Double, VecS Double)
smartRowsGetter CRS{..} row = V.map func $ G.convert (vecS row)
  where func (i1, x) = (x, crs G.! i1)
{-# INLINE smartRowsGetter #-}

-- | Super efficient sparse multiplication based on the SpGEMM algorithm as nicely explained
-- on: "An Efﬁcient GPU General Sparse Matrix-Matrix Multiplication for Irregular Data"
flatter :: V.Vector (Double, VecS Double) -> VecS Double
flatter vs = VecS $ U.filter ((/= 0) . snd) $ U.imap (,) $ U.create $ do
  let
    imax :: Int
    imax = G.foldl' (\acc (_, VecS x) -> if G.null x
                                         then acc
                                         else max acc (fst $ G.last x)
                    ) 0 vs
  vec <- MU.replicate (imax+1) 0
  let func (k, VecS v) = U.mapM_ (\(i,x) ->  MU.unsafeRead vec i >>=
                                            (MU.unsafeWrite vec i . (+ x * k))
                                 ) v
  V.mapM_ func vs
  return vec
{-# INLINE flatter #-}

-- Slow attempts of the flatter function
{--
flatter2 :: V.Vector (VecS Double) -> VecS Double
flatter2 vs = VecS $ G.filter ((/= 0) . snd) $ G.imap (,) vec
  where
    imax = G.foldl' (\acc (VecS x) -> if G.null x then acc else max acc (fst $ G.last x)) 0 vs
    vec0 = G.replicate (imax+1) 0
    vec  = G.foldl' (\acc (VecS x) -> G.unsafeAccumulate (+) acc x) vec0 vs

flatter_slow :: V.Vector (VecS Double) -> VecS Double
flatter_slow = f2 . f1
  where
    f1 = HM.fromListWith (+) . G.toList . G.concatMap (G.convert . vecS)
    f2 = VecS . G.fromList . HM.toList
--}

-- ===================================== Tools ===========================================

-- | Unsafe sorting. Sort in-place by first element of a vectors of pairs.
unsafeSort :: (G.Vector v1 (a, b), G.Vector v (a, b), Ord a, G.Mutable v1 ~ G.Mutable v) =>
              v (a, b) -> v1 (a, b)
unsafeSort v = runST $ do
  mv <- G.unsafeThaw v
  S.sortBy (compare `on` fst) mv
  G.unsafeFreeze mv