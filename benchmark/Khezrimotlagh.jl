# the object of this code is to record the time spent to solve each problem by algorithm
# model: standard VRS
# algorithm: Khezrimotlagh et al. (2019)
# solver: "Gurobi" with method "primal simplex"

using Gurobi

function run_and_output(dataFile::AbstractString, inputSize::Int64, rootPath::AbstractString)
    resultFilePath = "$rootPath/results$dataFile.csv"
    efficiency, runTime = runEJOR2019(dataFile, inputSize, rootPath)

    # println("tol time = $(runTime[end]), 6.1 = $(runTime[end-2]), 6.2 = $(runTime[end-1])")
    # my_log("tol time = $(runTime[end]), 6.1 = $(runTime[end-2]), 6.2 = $(runTime[end-1])")

    writecsv(resultFilePath, efficiency)
    f = open(resultFilePath, "a")
    write(f, "\n$(join(runTime, ','))\r\n")
    close(f)
end

function my_log(msg::AbstractString)
    resultFilePath = "$rootPath/Khezrimotlagh/ejor_log.txt"
    f = open(resultFilePath, "a")
    write(f, "\n$msg")
    close(f)
end

function runEJOR2019(dataFile::AbstractString, inputSize::Int64, rootPath::AbstractString)

    println("\n$dataFile begins @ $(Dates.hour(now())):$(Dates.minute(now())):$(Dates.second(now()))")
    dataPath = "$rootPath/$dataFile.csv" # input-output data source

    denseFullsizeConstMatrix = transpose(readcsv(dataPath)) # input User's (.csv) data path

    # initial varibles related with data model
    iterations = 0
    dimension, scale = size(denseFullsizeConstMatrix)
    lpSize = 0

    # initialization
    outputSize = dimension - inputSize
    tol = 10.0^-5

    # paper related parameters
    k = 0.25  # percentile threadhold (top 100k%)
    # subsample_size = 30 # p in paper
    subsample_size = convert(Int64, round(sqrt(scale), 0)) # p in paper
    # subsample_size = convert(Int64, round(min(sqrt(inputSize*scale), scale/2), 0)) # p in paper
    # subsample_size = convert(Int64, round(min(sqrt(dimension*scale), scale/2), 0)) # p in paper

    # create dictionary for fulldatamatrix mapping
    fulldatamatrix_dict = Dict([denseFullsizeConstMatrix[1:end,i] => i for i in 1:length(denseFullsizeConstMatrix[1,:])])
    dataConstMatrix = denseFullsizeConstMatrix[1:end, 1:scale]

    tic()
    # Step 3: select a subsample of D

    extremeValueSet = zeros(Int64, dimension)
    for i in 1 : dimension
      if i <= inputSize
        extremeValueSet[i] = findmin(dataConstMatrix[i,:])[2]  # may be duplicated
      else
        extremeValueSet[i] = findmax(dataConstMatrix[i,:])[2]  # may be duplicated
      end
    end
    extremeValueSet =  sort(unique(extremeValueSet))  # may be duplicated
    smallLP = dataConstMatrix[:,extremeValueSet[:]]
    slotleft = subsample_size - length(extremeValueSet)# dimension

    # score unselected DMU  (pseudocode 4)
    total_length = length(dataConstMatrix[1,:])
    score_record = zeros(total_length)
    for i in 1:length(dataConstMatrix[:,1])
      sorted = sortperm(dataConstMatrix[i,:])
      top_percentage = Int(floor((1.0-k)*total_length))
      score_record[sorted[top_percentage:length(sorted)]] += 1
    end

    smallLP = [smallLP dataConstMatrix[:, sortperm(score_record, rev=true)[1:slotleft]] ]  # only take top subsample_size - smallLP

    runTime3 = toq()
    tic()

    # ready for computing scores
    results = -1.0 * ones(scale) # vector to store final results (efficiency)

    env = Gurobi.Env()
    setparams!(env; Method=0, OutputFlag=0, OptimalityTol=tol) #using primal simplex method and not displaying solver's status
    # setparams!(env; OutputFlag=0) #using non-deterministic concurrent (default) method and not displaying solver's status


    # step 4
    subscores = computingScores(smallLP, smallLP, inputSize, env)

    for i = 1:length(subscores)
      results[fulldatamatrix_dict[smallLP[:, i]]] = subscores[i]
    end

    subsample = find(results.>=0)  # need to consider the tol
    bestpractice = find(results.>=1-tol) # B^s in step 4  # need to consider the tol
    not_subsample = find(results.==-1)

    runTime4 = toq()

    tic()
    # step 5
    # now I have bestpractice wrt dataMatrix
    subscores5 = computingScores(denseFullsizeConstMatrix[:,bestpractice], denseFullsizeConstMatrix[:,not_subsample], inputSize, env)
    subscores5_tmp = zeros(Int64, scale)
    dataMatrix5 = denseFullsizeConstMatrix[:,not_subsample]
    for i = 1:length(subscores5)
      results[fulldatamatrix_dict[dataMatrix5[:, i]]] = subscores5[i]
	  if subscores5[i] >= 1.0 - tol#> 1 # >1+tol
          subscores5_tmp[fulldatamatrix_dict[dataMatrix5[:, i]]] = 1
      end
    end

    exterior = find(subscores5_tmp.!=0)  # score > 1
# print(exterior)
    # step 6
    if isempty(exterior)
      runTime5n6 = toq()
      return results, [runTime3, runTime4 , runTime5n6, 0.0, 0.0, runTime3+runTime4+runTime5n6]
    end

    runTime5n6 = toq()

    tic()

    # step 6.1
    subscores6_1 = computingScores(denseFullsizeConstMatrix[:, union(subsample, exterior)], denseFullsizeConstMatrix[:, union(subsample, exterior)], inputSize, env)

    dataMatrix6_1 = denseFullsizeConstMatrix[:, union(subsample, exterior)]

    subscores6_1_tmp = zeros(Int64, scale)
    for i = 1:length(subscores6_1)
      results[fulldatamatrix_dict[dataMatrix6_1[:, i]]] = subscores6_1[i]
      if subscores6_1[i] >= 1-tol
          subscores6_1_tmp[fulldatamatrix_dict[dataMatrix6_1[:, i]]] = 1
      end
    end

    bestpractice2 = find(subscores6_1_tmp.!=0)

    not_union_subsample_exterior = setdiff(collect(1:1:scale), union(subsample, exterior)) # ?????

    runTime6_1 = toq()
    tic()
    # step 6.2  (2nd stage in BuildHull)
    subscores6_2 = computingScores(denseFullsizeConstMatrix[:, bestpractice2], denseFullsizeConstMatrix[:, not_union_subsample_exterior], inputSize, env)

    dataMatrix6_2 = denseFullsizeConstMatrix[:, not_union_subsample_exterior]
    for i = 1:length(subscores6_2)
      results[fulldatamatrix_dict[dataMatrix6_2[:, i]]] = subscores6_2[i]
    end
    runTime6_2 = toq()

    return results, [runTime3, runTime4 , runTime5n6, runTime6_1, runTime6_2, runTime3+runTime4+runTime5n6+runTime6_1+runTime6_2]
end


# determining efficiency of each entry in rhsDataMatrix w.r.t. peerMatrix
function computingScores(peerMatrix::Matrix{Float64}, rhsDataMatrix::Matrix{Float64}, numInput::Int64, env::Gurobi.Env)
  # parameters
  # peerMatrix: the peer data
  # rhsDataMatrix: the data to be evaluated. Each column is a record to be evaluated
  # numInput: number of input factors
  # env: Gurobi environment


  # data preprocess, initialization
  dimension, scale = size(rhsDataMatrix)
  peerdimension, peerscale = size(peerMatrix)
  if dimension != peerdimension
    println("error")  # throw error
  end

  numOutput = dimension - numInput
  lambdaZero = zeros(peerscale) # the vector all zero with length equal to scale
  lambdaOne = ones(peerscale) # the vector all one with length equal to scale
  inputZero = zeros(numInput)
  outputZero = zeros(numOutput)

  objValue = 0.0
  scores = zeros(scale) # vector storing the efficiency scores for rhsDataMatrix

  # building model
  model = Gurobi.Model(env, "lp", :minimize)
  # defining variables
  add_cvars!(model, lambdaZero, lambdaZero, Inf)  # defining lambda(r), r from 1 to scale
  add_cvar!(model, 1.0, -Inf, Inf) # defining theta
  update_model!(model)

  # defining constraints
  # input constraints
  add_constrs!(model, [peerMatrix[1:numInput, :] -rhsDataMatrix[1:numInput, 1]], '<', inputZero)
  # output constraints
  add_constrs!(model, [peerMatrix[numInput+1:numInput+numOutput, :] outputZero], '>', rhsDataMatrix[numInput+1:numInput+numOutput, 1])
  # convexity constraint
  add_constr!(model, [lambdaOne; 0], '=', 1.0)
  update_model!(model)
  # end building model

  # solving model
  optimize(model)

  if get_status(model) == :optimal # if optimal
    scores[1] = get_objval(model)
  else
    scores[1] = 98765.4321  # infeasible or unbounded
  end


  for k = 2:scale  # modifying the coefficient with warm start
    # update RHS and solve LP again
    scores[k] = updateAndSolve4EJOR(model, peerscale, numInput, numOutput, rhsDataMatrix[:, k])
  end

  return scores
end

# update RHS and solve LP
function updateAndSolve4EJOR(model::Gurobi.Model, scale::Int64, numInput::Int64, numOutput::Int64, dmuValue::Vector{Float64})
  # model: model that variables are added into
  # scale: # of DMUs in the peer
  # numInput: number of inputs
  # numOutput: number of outputs
  # dmuValue: vector of the new DMU to be evaluated

  # changing RHS of model
  Gurobi.set_dblattrarray!(model, "RHS", numInput + 1, numOutput, dmuValue[numInput+1:numInput+numOutput])
  # changing coefficients

  # cind = collect(1:numInput) #find(transDataMatrix[1:numInput, k])
  vind = fill(scale+1, numInput)
  chgCoeffs!(model, numInput, collect(1:numInput), vind, -dmuValue[1:numInput])
  # chgCoeffs!(model, numInput, collect(2:3), [scale+1], dmuValue[1:numInput])
  # chgCoeffs!(model, number of coefficients to be changed, index of constraint, index of variable, new coefficients' value)

  update_model!(model)

  # solving model
  optimize(model)
  if get_status(model) == :optimal # if optimal
    return get_objval(model)
  else
    return 98765.4321  # infeasible or unbounded
  end
end

# a macro of gurobi c call function
macro grb_ccall(func, args...)
    f = "GRB$(func)"
    args = map(esc,args)

    # default use gurobi80

    is_unix() && return quote
        ccall(($f,:gurobi80), $(args...))
    end
    is_windows() && VERSION < v"0.6-" && return quote
        ccall(($f,:gurobi80), stdcall, $(args...))
    end
    is_windows() && VERSION >= v"0.6-" && return quote
        ccall(($f,:gurobi80), $(esc(:stdcall)), $(args...))
    end
end

# chgCoeffs!: changing the coefficients of variables on a constraint
function chgCoeffs!(model::Gurobi.Model, n::Int64, cind::Vector{Int64}, vind::Vector{Int64}, val::Vector{Float64})
  # parameters
  # input
  # model: model that variables are added into
  # n: number of coefficients to be changed
  # cind: index of constraint
  # vind: index of variable
  # val: new coefficients' value

  # main call
  ret = @grb_ccall(chgcoeffs, Cint, (
                   Ptr{Void},
                   Cint,
                   Ptr{Cint},
                   Ptr{Cint},
                   Ptr{Float64},
                   ),
                   model,
                   n,
                   Gurobi.ivec(cind.-1),
                   Gurobi.ivec(vind.-1),
                   Gurobi.fvec(val))

  if ret != 0
    throw(Gurobi.GurobiError(model.env, ret))
  end
  nothing
end
