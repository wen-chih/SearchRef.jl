# Performance comparison with BuildHull
  
-**BuildHull.jl** implements the algorithm BuildHull ([Dula, 2011](https://doi.org/10.1287/ijoc.1100.0400)). BuildHull is the state-of-the-art algorithm for DEA computation, and is used as the benchmark for SearchRef.jl. Here, currently, it solves the input-oriented VRS models. 

-**Khezrimotlagh.jl** implements the algorithm proposed by [Khezrimotlagh et al. (2019)](https://doi.org/10.1016/j.ejor.2018.10.044). 

-**fullsize.jl** implements the naive method, solving full-size LPs, to determine the efficiencies of all DMUs.

-**SearchRef_batch.jl** will call SeachRef.jl to compute DEA efficiencies for performance comparison. 

Please refer ([Dula, 2011](https://doi.org/10.1287/ijoc.1100.0400)) for the full data sets. The settings for implementing BuildHull and SearchRef are idential. The detailed settings and results can be found in [A reference-searching–based algorithm for large-scale data envelopment analysis computation](https://arxiv.org/abs/1710.10482/).   
