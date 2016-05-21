# Description

This is an Haskell implementation of the [MCL](https://en.wikipedia.org/wiki/K-d_tree) (Markov Chain Clustering) algorithm. The MCL find cluster of strong connected
nodes in graph without a predifined and fixed number of cluster (like K-means).

## Usage

````

inputGraph :: Graph Int Double
inputGraph = mkUniGraph [20, 2] ([((1,2), 1), ((2,3), 1), ((5,2), 1), ((5,3), 1), ((60,1), 1) , ((7,4), 1)]::[((Int,Int), Double)])

main = runMCL defaultMCL {inflation = 2} inputGraph >>= mapM_ print

````

The granularity of the cluster depend on the value of the 'inflation' factor. The minimum value in 1.2 and increasing it results in smaller clusters.
