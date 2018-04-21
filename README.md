# SearchRef.jl
#### Copyright © 2018 by Wen-Chih Chen.  Released under the MIT License.

SearchRef.jl is a Julia package for Frontier Efficiency Analysis (aka Data Envelopment Analysis, DEA) computation. In particular, it determines a firm’s relative efficiency. SearchRef.jl is desiged to enhance the large-scale DEA computation. The theoretical foundation can be found in the paper: [A reference-searching–based algorithm for large-scale data envelopment analysis computation](https://arxiv.org/abs/1710.10482/). Currently, it solves the input-oriented VRS model. 

## Requirements
SearchRef.jl requires Gurobi solver.

## Example

```julia
# The example determines the efficiencies for all DMUs based on the VRS input-oriented model (VRS model)
using Gurobi
include("SearchRef.jl")

dataMatrix = readcsv("data.csv") # reed input-ouput data
numInput = 2 # the first two columns are inputs and the others are outputs
scale, dimension = size(dataMatrix)

GurobiTol = 10.0^-6
efficiency = 0.0
 
env = Gurobi.Env()
setparams!(env; OutputFlag=0, OptimalityTol=GurobiTol) 

for k = 1:scale
    efficiency, extremeValueSet = SearchRef(k, dataMatrix, numInput, env)
    # k: DMU (row) k of the dataMatrix
end

# The following can give better performance.
#eSet = Array{Int64}(0)
#for k = 1:scale
#    efficiency, eSet = SearchRef(k, dataMatrix, numInput, env, extremeValueSet=eSet, initialSampling = "extreme value")
#    # k: DMU (row) k of the dataMatrix
#end
```


## Parameters

>
**incrementSize** : the incremental size to expand the sample (default value: 150).

	SearchRef(k, dataMatrix, numInput, env, incrementSize = 200) # set the incremental size to 200

>
**tol** : the solution tolerance for solving DEA problem (default value: 1e-6). It also resets the dual feasibility tolerance in the solver to the given value.
<br>

	SearchRef(k, dataMatrix, numInput, env, tol = 10.0^-4) # set the solution tolerance to 1e-4

>
**lpUB** : the size limit of the LP, i.e. the limitation of number of variables in the LP (default value: Inf).
<br>

	SearchRef(k, dataMatrix, numInput, env, lpUB = 300) # set the LP size limitation to 300 variables

>
**initialSize** : to size of the initial sample (default value: dimension+1).
<br>

	SearchRef(k, dataMatrix, numInput, env, initialSize = 300) # set the size of the initial sample to be 300 

>
**M** : a posive-value parameter to avoid the degenercy (default value: 10).
<br>

	SearchRef(k, dataMatrix, numInput, env, M = 5) # set M = 5

>
**maxIteration** : the maximum number of iteration of the computation (default value: 300).
<br>

	SearchRef(k, dataMatrix, numInput, env, maxIteration = 100) # set maxIteration = 100

>    
**initialSampling** : the initial sampling policy: can be "random" or "extreme value"), which will select DMUs with extreme values on each dimension (default value: "random").
<br>

	SearchRef(k, dataMatrix, numInput, env, initialSampling = "extreme value") 

>
**extremeValueSet** : the index set of DMUs with extreme vale for each dimension (default is empty).
<br>

	SearchRef(k, dataMatrix, numInput, env, extremeValueSet = aSet) 
	# aSet is Vector{Int} containing indices of DMUs with extreme values 



## Citation
If you find FrontierEfficiencyAnalysis useful in your work, we kindly request that you cite the following papers

	@misc{chen2017,
	Author = {Wen-Chih Chen},
	Title = {A reference-searching–based algorithm for large-scale data envelopment analysis computation},
	Year = {2017},
	Eprint = {https://arxiv.org/abs/1710.10482},
	}

## Acknowledgements
SearchRef.jl has been developed under the financial support of the Ministry of Science and Technology, Taiwan (Grant No. 104-2410-H-009-026-MY2). The contributors include Yueh-Shan Chung and Hao-Yun Chen.
