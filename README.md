# mcl -- Markov Cluster Algorithm

A Haskell library implementing the Markov Cluster Algorithm (MCL), a fast and scalable unsupervised clustering algorithm for graphs.

## What is this?

MCL was originally conceived by Stijn van Dongen. It discovers cluster structure in graphs by simulating random walks (stochastic flow). The core idea: if you drop a random walker on a graph, it tends to stay within densely connected regions. MCL exploits this by repeatedly strengthening in-cluster flow and weakening between-cluster flow until the graph naturally separates into clusters.

The algorithm operates on a column-stochastic (Markov) matrix derived from the graph's adjacency matrix. It is a general-purpose algorithm applicable wherever you need to find natural groupings in a graph: social networks, protein interaction networks, image segmentation, community detection, etc.

## How the algorithm works

MCL alternates three operations until convergence:

### Expansion

Squares the current stochastic matrix (M x M), which is equivalent to simulating random walks of increasing length. After expansion, two nodes in the same dense cluster will have high mutual transition probabilities, while nodes in different clusters will have low ones. The implementation uses an efficient Sparse General Matrix-Matrix Multiplication (SpGEMM) with parallel evaluation across cores.

### Inflation

Raises every matrix entry to a power *r* (the inflation parameter) and re-normalizes each column to sum to 1. This is the key step that creates separation: strong connections get stronger, weak connections get weaker. Higher inflation values produce more, smaller clusters; lower values produce fewer, larger clusters. Typical range: 1.2 to 5.0.

### Pruning

After inflation, many entries become negligibly small. Pruning removes them to maintain sparsity and keep computation tractable. Four strategies are available: fixed threshold, variable threshold (adapts to data distribution), dynamic budget (limits total pruned mass per column), or no pruning.

### Convergence

The algorithm tracks a "chaos" metric comparing the maximum entry to the sum of squared entries in each column. When each column has converged to a near-idempotent state (one dominant entry), the matrix has stabilized and clusters can be read off.

## Example

```haskell
import Data.Graph.Markov

-- Two triangles connected by a single bridge edge
let g = mkUniGraph []
          [ ((1, 2), 1), ((2, 3), 1), ((3, 1), 1)   -- triangle 1
          , ((4, 5), 1), ((5, 6), 1), ((6, 4), 1)   -- triangle 2
          , ((3, 4), 1)                               -- bridge
          ]

let clusters = runMCL defaultMCL g
-- Expected: [[1,2,3],[4,5,6]]
```

## Where is it used?

- **hammer** -- uses MCL for graph-based analysis of microstructure connectivity in polycrystalline materials

## How to build

```bash
# With Nix (recommended)
nix develop
cabal build --allow-newer

# With Cabal
cabal build

# Run tests
cabal test
```

For parallel matrix multiplication, run with multiple cores:

```bash
./program +RTS -N4    # Use 4 cores
./program +RTS -N     # Use all available cores
```

## References

- [Markov Cluster Algorithm (Wikipedia)](https://en.wikipedia.org/wiki/Markov_clustering)

## License

MIT -- see [LICENSE](./LICENSE).
