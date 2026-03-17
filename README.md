# mcl -- Markov Cluster Algorithm

A Haskell library implementing the Markov Cluster Algorithm (MCL), a fast and scalable unsupervised clustering algorithm for graphs. MCL was originally conceived by Stijn van Dongen and simulates stochastic flow on a graph to identify natural clusters by alternating expansion and inflation operations on a column-stochastic (Markov) matrix.

## Algorithm Overview

The MCL algorithm discovers cluster structure in graphs by simulating random walks. It operates on the column-stochastic matrix derived from the graph's adjacency matrix through three core operations applied iteratively until convergence:

### Expansion

Expansion computes the square of the current stochastic matrix (`M * M`), simulating longer random walks. Nodes within the same dense cluster will strengthen their mutual transition probabilities, while nodes in different clusters will see theirs weaken. The implementation uses an efficient sparse matrix-matrix multiplication (SpGEMM) based on the algorithm described in _"An Efficient GPU General Sparse Matrix-Matrix Multiplication for Irregular Data"_.

### Inflation

Inflation raises every element of the matrix to a power `r` (the inflation parameter) and then re-normalizes each column so that it sums to 1. This has the effect of boosting strong connections and demoting weak ones, sharpening the distinction between intra-cluster and inter-cluster flow. Higher inflation values produce finer-grained (more numerous, smaller) clusters; lower values produce coarser clusters.

### Pruning

After inflation, many matrix entries become negligibly small. Pruning removes these entries to maintain sparsity and keep computation tractable. The library supports four pruning strategies:

| Strategy | Type | Description |
|---|---|---|
| `FixPrune t` | Fixed threshold | Removes all entries below a fixed value `t`. Simple and predictable. |
| `VarPrune a` | Variable threshold | Computes a threshold from the statistical distribution of the vector (based on mean and max), scaled by factor `a`. Adapts to the data. |
| `DynPrune k` | Dynamic budget | Ensures the total mass of pruned entries does not exceed a fraction `k` of the column sum. Provides a pruning budget guarantee. |
| `NoPrune` | No pruning | Retains all entries. Useful for debugging or very small graphs. |

### Convergence

The algorithm tracks a "chaos" metric: for each column, it compares the maximum entry to the sum of squared entries. When these are close (i.e., each column has converged to a near-idempotent state), the matrix has stabilized and the clusters can be read off. Iteration stops when the chaos falls below `minError` or the maximum number of iterations (`maxIter`) is reached.

## Architecture

The library is organized into three layers, each building on the one below:

```
Data.Graph.Markov   (MCL algorithm: configuration, iteration loop, cluster extraction)
       |
Data.Graph.Sparse   (CRS sparse matrix and sparse vector operations)
       |
Data.Graph          (HashMap-based graph: construction, traversal, algebraic ops)
```

### Layer 1: Graph (`Data.Graph`) -- internal module

A general-purpose graph representation built on `HashMap`. Supports both undirected and directed graphs with arbitrary node and edge types.

**Core data type:**

```haskell
newtype Graph a b = Graph { graph :: HashMap a (HashMap a b) }
```

**Construction:**
- `mkUniGraph :: [a] -> [((a, a), b)] -> Graph a b` -- Build an undirected graph from a node list and edge list. Each edge `(a, b)` is stored in both directions.
- `mkDiGraph :: [a] -> [((a, a), b)] -> Graph a b` -- Build a directed graph. Each edge is stored only in its given direction.

**Query:**
- `getChilderen :: Graph a b -> a -> HashMap a b` -- Get all neighbors (children) of a node.
- `hasEdge :: Graph a b -> a -> a -> Bool` -- Check whether an edge exists.
- `getGraphEdge :: Graph a b -> (a, a) -> Maybe b` -- Look up the weight of a specific edge.
- `getSubGraph :: Graph a b -> [a] -> Graph a b` -- Extract the induced subgraph over a set of nodes.
- `graphToList :: Graph a b -> [((a, a), b)]` -- Flatten the graph to an edge list.

**Graph search:**
- `search :: Graph a b -> a -> HashSet a` -- BFS-like traversal returning all reachable nodes.
- `dfs :: Graph a b -> a -> HashSet a` -- Depth-first search returning all reachable nodes.
- `connComp :: Graph a b -> [HashSet a]` -- Find all connected components.

**Algebraic operations (for `Graph Int b` used as sparse matrices):**
- `addGraph :: Graph Int b -> Graph Int b -> Graph Int b` -- Element-wise matrix addition.
- `multGraph :: Graph Int b -> Graph Int b -> Graph Int b` -- Matrix-matrix multiplication.
- `invertGraph :: Graph a b -> Graph a b` -- Transpose (invert edge directions).

### Layer 2: Sparse Matrix (`Data.Graph.Sparse`) -- exposed module

High-performance sparse linear algebra using unboxed vectors. This layer provides the Compressed Row Storage (CRS) representation that the MCL algorithm operates on directly.

**Sparse vector (`VecS`):**

```haskell
newtype VecS a = VecS { vecS :: Vector (Int, a) }
```

A sparse vector stored as an unboxed vector of `(index, value)` pairs, sorted by index.

- `mkVecS :: [(Int, Double)] -> VecS Double` -- Construct a sparse vector (sorts input).
- `unsafeMkVecS :: Vector (Int, Double) -> VecS Double` -- Construct from pre-sorted input (no sort).
- `addVV :: VecS Double -> VecS Double -> VecS Double` -- Sparse vector addition via merge.
- `multVV :: VecS Double -> VecS Double -> Double` -- Sparse dot product (inner product).
- `multKV :: Double -> VecS Double -> VecS Double` -- Scalar-vector multiplication.

**Sparse matrix (`CRS`):**

```haskell
newtype CRS a = CRS { crs :: V.Vector (VecS a) }
```

A row-major sparse matrix where each row is a `VecS`. All rows are stored, including empty ones (for index consistency).

- `graphToCRS :: Graph Int Double -> CRS Double` -- Convert a `Graph` to CRS format.
- `crsToGraph :: CRS Double -> Graph Int Double` -- Convert CRS back to a `Graph`.
- `transpose :: CRS Double -> CRS Double` -- Matrix transpose using mutable vectors.
- `multMV :: CRS Double -> VecS Double -> VecS Double` -- Matrix-vector multiplication.
- `multMM :: CRS Double -> CRS Double -> CRS Double` -- Matrix-matrix multiplication (transpose-based).
- `multMMsmrt :: CRS Double -> CRS Double -> CRS Double` -- Matrix-matrix multiplication using the SpGEMM scatter approach. This is the method used by the MCL iteration.

**Parallelism:**

Matrix multiplication (`multMM`, `multMMsmrt`) is parallelized using GHC's `Control.Parallel.Strategies`. The work is partitioned into chunks sized proportionally to the number of GHC capabilities (`-Nx` RTS option).

### Layer 3: Markov Cluster (`Data.Graph.Markov`) -- exposed module

The top-level module that implements the MCL algorithm and re-exports `Data.Graph`.

**Configuration:**

```haskell
data MCLCfg = MCLCfg
  { inflation :: Double   -- Inflation exponent (controls cluster granularity)
  , pruneMode :: Pruning  -- Pruning strategy
  , maxIter   :: Int      -- Maximum iterations
  , minError  :: Double   -- Convergence threshold (chaos metric)
  , selfLoop  :: Double   -- Self-loop weight added to each node (0 = none)
  }
```

**Default configuration:**

```haskell
defaultMCL :: MCLCfg
defaultMCL = MCLCfg
  { inflation = 1.2
  , pruneMode = FixPrune 0.0001
  , maxIter   = 200
  , minError  = 0.0001
  , selfLoop  = 0
  }
```

**Entry point:**

```haskell
runMCL :: MCLCfg -> Graph Int Double -> [[Int]]
```

Takes a configuration and a weighted graph, returns a list of clusters where each cluster is a list of node IDs.

**Internal pipeline (inside `runMCL`):**

1. `addSelfNodes` -- Optionally adds self-loops to every node.
2. `graphToCRS` -- Converts the graph to CRS sparse matrix format.
3. `normalizeCRS` -- Normalizes each column to sum to 1 (making it column-stochastic).
4. Iterative loop applying `stepMCL`:
   - `expandCRS` -- Squares the matrix (`multMMsmrt x x`).
   - `inflateCRS` -- Applies element-wise power, pruning, and re-normalization.
5. `getChaos` -- Measures convergence; stops when below `minError`.
6. `getMClClusters` -- Extracts clusters from the converged matrix.

## Usage

### Basic example

```haskell
import Data.Graph.Markov

main :: IO ()
main = do
  -- Define a graph with two triangles connected by one edge
  let g = mkUniGraph []
            [ ((1, 2), 1), ((2, 3), 1), ((3, 1), 1)   -- triangle 1
            , ((4, 5), 1), ((5, 6), 1), ((6, 4), 1)   -- triangle 2
            , ((3, 4), 1)                               -- bridge
            ]

  -- Run MCL with default configuration
  let clusters = runMCL defaultMCL g
  print clusters
  -- Expected: two clusters, e.g. [[1,2,3],[4,5,6]]
```

### Custom configuration

```haskell
import Data.Graph.Markov

main :: IO ()
main = do
  let g = mkUniGraph [] [((0,1), 1), ((1,2), 1), ((2,0), 1), ((3,4), 1), ((4,3), 1)]

  let cfg = defaultMCL
        { inflation = 2.0         -- higher inflation = smaller clusters
        , pruneMode = DynPrune 0.05  -- dynamic pruning with 5% budget
        , selfLoop  = 1.0         -- add self-loops with weight 1
        , maxIter   = 100
        , minError  = 1e-5
        }

  print $ runMCL cfg g
```

### Using the sparse matrix layer directly

```haskell
import Data.Graph.Markov
import Data.Graph.Sparse

main :: IO ()
main = do
  let g = mkDiGraph []
            [ ((1,1), 2), ((1,2), 5), ((1,3), 1)
            , ((2,1), 4), ((2,2), 3), ((2,3), 1)
            ]

  -- Convert to CRS and perform matrix operations
  let m = graphToCRS g
  let mt = transpose m
  let product = multMMsmrt m mt
  print product
```

## Configuration Reference

| Parameter | Type | Default | Description |
|---|---|---|---|
| `inflation` | `Double` | `1.2` | Inflation exponent. Higher values produce more, smaller clusters. Typical range: 1.2 -- 5.0. |
| `pruneMode` | `Pruning` | `FixPrune 0.0001` | Strategy for removing negligible entries after inflation. |
| `maxIter` | `Int` | `200` | Hard upper bound on the number of expansion/inflation cycles. |
| `minError` | `Double` | `0.0001` | Convergence threshold. Iteration stops when chaos drops below this value. |
| `selfLoop` | `Double` | `0` | Weight of self-loops added to each node before clustering. Set to a positive value (e.g., 1.0) to bias towards keeping nodes in their own cluster. |

## Dependencies

| Package | Version | Purpose |
|---|---|---|
| `base` | >= 4, < 5 | Standard library |
| `containers` | >= 0.4.2.1 | Standard container types |
| `deepseq` | >= 1.2 | Deep evaluation for parallelism |
| `hashable` | >= 1.2 | Hashing for HashMap/HashSet |
| `parallel` | >= 3.0 | Parallel evaluation strategies |
| `random` | >= 1.0 | Random number generation |
| `unordered-containers` | >= 0.2 | HashMap and HashSet (graph representation) |
| `vector` | >= 0.10 | Boxed vectors (CRS rows, parallel strategies) |
| `vector-algorithms` | >= 0.5 | In-place sorting for sparse vector construction |

## Building

### With Cabal

```bash
cabal build
```

### With Nix (development shell)

```bash
nix develop
cabal build --allow-newer
```

The Nix flake provides a development shell with GHC, cabal-install, HLS, and formatters (fourmolu, cabal-fmt, nixpkgs-fmt).

### Compiler flags

The library is compiled with aggressive optimizations by default:

- `-O2` -- Optimization level 2
- `-threaded` -- Enable threaded runtime (required for parallel matrix multiplication)
- `-funbox-strict-fields` -- Unbox strict fields for better performance
- `-rtsopts` -- Enable RTS options (use `+RTS -N` to set number of capabilities)

### Running with parallelism

To benefit from parallel matrix multiplication, run the program with multiple capabilities:

```bash
./my-program +RTS -N4    # Use 4 cores
./my-program +RTS -N     # Use all available cores
```

## Module Summary

| Module | Visibility | Description |
|---|---|---|
| `Data.Graph` | Internal | HashMap-based graph with construction, traversal (BFS, DFS), connected components, and sparse-matrix-style algebraic operations (add, multiply, transpose). |
| `Data.Graph.Sparse` | Exposed | Compressed Row Storage (CRS) sparse matrix and sparse vector (`VecS`) with parallel matrix multiplication, transpose, and scalar operations. |
| `Data.Graph.Markov` | Exposed | MCL algorithm entry point (`runMCL`), configuration (`MCLCfg`), pruning strategies, and convergence detection. Re-exports `Data.Graph`. |

## References

- [Markov Cluster Algorithm (Wikipedia)](https://en.wikipedia.org/wiki/Markov_clustering)

## Author

Edgar Gomes de Araujo (<talktoedgar@gmail.com>)

## License

MIT -- see [LICENSE](./LICENSE).
