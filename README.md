# SearchRef.jl
#### Copyright © 2018 by Wen-Chih Chen.  Released under the MIT License.

SearchRef.jl is a Julia package for Frontier Efficiency Analysis (aka Data Envelopment Analysis, DEA) computation. It is desiged to enhance the large-scale DEA computation. It is the implementation of the algorithm presented in "[A reference-searching–based algorithm for large-scale data envelopment analysis computation](https://arxiv.org/abs/1710.10482/)". Currently, it solves the input-oriented VRS model. 

## Usage
DEA is a linear program (LP)-based method used to determine a firm’s relative efficiency.

## Example

```julia
# The example determine the efficiencies for all DMUs based on the VRS input-oriented model (VRS model)
using Gurobi
include("SearchRef.jl")

dataPath = "data.csv" # input-output data source

dataMatrix = readcsv(dataPath)
scale, dimension = size(dataMatrix)
numInput = 2 # the first two columns are inputs and the others are outputs

GurobiTol = 10.0^-6
efficiency = 0.0
iterNum = 0
 
env = Gurobi.Env()
setparams!(env; OutputFlag=0, OptimalityTol=GurobiTol) 

for k = 1:scale
    efficiency, iterNum, extremeValueSet = SearchRef(k, dataMatrix, numInput, env)
end
```

<br>

## Parameters

>
**incrementSize** : the incremental size to expand the sample ( default value: 100 ).

	solveDEA(model, incrementSize = 200) # set the incremental size to 200

>
**tol** : the solution tolerance for solving DEA problem (default value: 1e-6). It also resets the dual feasibility tolerance in the solver to the given value.
<br>

	solveDEA(model, tol = 10^-4) # set the solution tolerance to 1e-4

>
**lpUB** : the size limit of the LP, i.e. the limitation of number of variables in the LP (default value: Inf).
<br>

	solveDEA(model, lpUB = 300) # set the LP size limitation to 300 variables

>
**extremeValueSetFlag** : to enable (=1) or disable (=0) performing initial sampling by selecing extreme value in each input/output dimension (default value: 0).
<br>


	solveDEA(model, extremeValueSetFlag = 1) # enable




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
