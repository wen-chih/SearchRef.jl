# Performance comparison with BuildHull
 
 
-**BuildHull.jl** implements the algorithm BuildHull ([Dula, 2011](https://doi.org/10.1287/ijoc.1100.0400)). BuildHull is the state-of-the-art algorithm for DEA computation, and is used as the benchmark for SearchRef.jl. Here, currently, it solves the input-oriented VRS models. 

-**BuildHull.jl** implements

-**SearchRef_batch.jl** will call SeachRef.jl to compute DEA efficiencies for performance comparison. The settings for BuildHull_bacth.jl and SearchRef_batch.jl are idential. The settings and results can be found in [A reference-searching–based algorithm for large-scale data envelopment analysis computation](https://arxiv.org/abs/1710.10482/).   
