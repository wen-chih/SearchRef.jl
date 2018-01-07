# for Chen, W.C. "A reference-searching–based algorithm for large-scale data envelopment analysis computation".
# solving input-oriented VRS model

using Gurobi
# include("functions.jl")

function SearchRef(k::Int64, dataMatrix::Matrix{Float64}, numInput::Int, env::Gurobi.Env; incrementSize = 150, M = 10.0, tol = 10.0^-6, lpUB = 0, initialSize = 0, maxIteration = 300, initialSampling = "random", extremeValueSet=Array{Int64}(0))

    # k: index of the DMU to be evaluated
    # dataMatrix: full-size data (scale by dimension) matrix
    # numInput: number of inputs
    # env: Gurobi setting
    # incrementSize (default: 150): incrementation size to expand the sample
    # M (default: 10.0): the multiplier of algorithm "anti-degeneracy" that make DMU worse
    # tol (10.0^-6): optimal tolerance for the searchRef()  #10.0^-6 value of OptimalityTol in Gurobi
    # lpUB (default: scale): the size limit for each LP solving
    # initialSize (default: dimension+1): size for the initial sample
    # maxIteration (default: 300): iteration limit
    # initialSampling (default: "random"): can be "random" or "extreme value"), which select DMU with extreme value on each dimension
    # extremeValueSet (default is empty): index set of DMUs with extreme vale for each dimension

    scale, dimension = size(dataMatrix)
    numOutput = dimension - numInput
    if lpUB <= 0 # using the real default value = scale
      lpUB = scale # the limit size of LP
    end
    if initialSize <= 0 # using the real default value = dimension + 1
      initialSize = dimension + 1
    end
    increment = [incrementSize, incrementSize] # incremental plan

    # -------------------------------------------------------------------
    # initialization of the computation (to optimize the Julia performance)
    lambdaZero = zeros(initialSize) # same size of lambda in DEA model
    lambdaOne = ones(initialSize) # same size of lambda in DEA model
    inputZero = zeros(numInput) # size of input no.
    outputZero = zeros(numOutput) # size of output no.
    dual = Array{Float64}(dimension+1)  # the dual values of constraints
    iterations = 0  #the iterations that algorithm uses to solve model
    objValue = 0.0
    realObjValue = 0.0
    deleteIdx = 2  # index of DMUs which is selected to be deleted from set (k is the first)
    optStatus = false # the optimality status of model
    additiveSize = 0
    lpSize = initialSize + 1 # size of lp:  lambdas + theta

    # --------------------------------------------------------
    # algorithm begins here to solve efficiency of k

    # making k worse, the RHS of DEA model
    originK = copy(dataMatrix[k, :])
    dataMatrix[k, 1:numInput] = dataMatrix[k, 1:numInput]*M

    # --------------------------------------------------------
    # initial sampling
    if !(initialSampling == "random" || initialSampling == "extreme value")
      error("The initial sampling policy is not supported")
    else
      if initialSampling == "random" # select sample randomly
        # do nothing
      else # initialSampling == "extreme value"
        # get extreme value
        if isempty(extremeValueSet)
          extremeValueSet = getExtremeValueIdx(dataMatrix, numInput)
        end
      end
    end

    set, notset = get2Sets(scale, initialSize, k, sort(unique([k; extremeValueSet])))
    # set: sample; notset: nonsample
    # end of initial sampling

    # ------------------------------------
    # solving LP by Gurobi
    model = Gurobi.Model(env, "lp", :minimize) # building model

    # defining variables
    add_cvar!(model, 1.0, -Inf, Inf) # defining theta
    add_cvars!(model, lambdaZero, lambdaZero, Inf)  # defining lambda(r), r from 1 to initialSize-1
    update_model!(model)
    # input constraints
    # add_constrs_t!(model, At, rel, rhs)  # here At can be dense or sparse
    add_constrs_t!(model, [-reshape(dataMatrix[k, 1:numInput], 1, numInput); dataMatrix[set, 1:numInput]], '<', inputZero)

    # output constraints
    add_constrs_t!(model, [zeros(1, numOutput); dataMatrix[set, numInput+1:dimension]], '>', dataMatrix[k, numInput+1:dimension])

    # convexity constraint
    add_constr!(model, [0.0; lambdaOne], '=', 1.0)

    update_model!(model)
    # end building model

    # ----------------------------------------
    # loop until model is optimal or algorithm processing over maxIterations
    # search for proper sample
    while optStatus == false && iterations <= maxIteration
      iterations = iterations + 1

      optimize(model) # solving the LP model
      dual = Gurobi.get_dblattrarray(model, "Pi", 1, num_constrs(model))
      objValue = Gurobi.get_objval(model)

      # modifying additiveSize, let lpSize + additiveSize can't be more than lpUB
      additiveSize = min(max(lpUB-lpSize,1), increment[min(iterations, length(increment))])
      # if lpSize < lpUB # there's room to increase sample size
      #   additiveSize = increment[min(iterations, length(increment))]
      #   if lpSize + additiveSize > lpUB
      #     additiveSize = lpUB - lpSize
      #   end
      # else
      #   # if the lpSize is already equal to lpUB, then let additiveSize to be 1
      #   additiveSize = 1
      # end

      # checking optimality and getting the most infeasible DMU which is added to set
      optStatus, overIdx = checkOpt(dual, dataMatrix, notset, additiveSize, tol)
      # overIdx is the set to be added
      # end checking optimality

      # resampling
      if optStatus==false # redefining sample if the model is not optimal

        if additiveSize > 1 # there's room to increase sample size, adding without dropping
          # add_vars!(model, n, A, c, lb, ub)
          add_vars!(model, additiveSize, [dataMatrix[notset[overIdx],:]'; ones(1, additiveSize)], zeros(additiveSize), zeros(additiveSize), fill(Inf, additiveSize))
          update_model!(model)

          set, notset = regruping(overIdx, set, notset) # regrouping (redefining the DMUs in or nt in the sample)
          lpSize = lpSize + additiveSize
        else # additiveSize = 1, adding and droping
          # deleteIdx wrt set
          deleteIdx = rand(2:length(set)-1)

          # get the index of constraint to be updated
          cind = collect(1:1:dimension)
          # get the index of variables to be updated
          vind = fill(deleteIdx + 1, dimension) # index of variables to be changed, the 1st is theta
          # update the coefficieint values
          # vval = vec(dataMatrix[notset[overIdx[1]], :]') # new values
          vval = dataMatrix[notset[overIdx[1]], :] # new values
          chgCoeffs!(model, length(cind), cind, vind, vval)
          # chgCoeffs!(model, # of constraints to be changed, ...)
          update_model!(model)

          # adding the most infeasible DMU (overIdx[1]) from notset to set and
          # move deleteIdx from set to noset (remove one frm sample)
          swapSample!(set, notset, deleteIdx, overIdx[1])
        end # if additiveSize > 1

      else # optStatus== true
        # objValue = get_objval(model)
        # if original k passes the kkt conditions, k is inefficient, if not, k is efficient
        if (originK'*dual[1:dimension])[1] + dual[dimension+1] <= tol
          # Result 1 in the paper
          realObjValue = objValue*M
        else # original k is infeasible to KKT
          # Result 2 in the paper
          if abs(objValue - 1) <= tol # degenerate, indeed its real effcieincy is one.
            realObjValue = 10.0 # degenerate, indeed its real effcieincy is one.
          else
            realObjValue = 1.0
          end
        end

      end # optStatus == false
      # end redefining set
    end # end loop while optStatus == false && iterations <= maxIteration

    dataMatrix[k, :] = originK

    return  realObjValue, extremeValueSet
end

# ==============================================================================
# a macro of gurobi c call function
macro grb_ccall(func, args...)
    f = "GRB$(func)"
    args = map(esc,args)

    is_windows() && VERSION < v"0.6-" && return quote
        ccall(($f,:gurobi70), stdcall, $(args...))
    end
    is_windows() && VERSION >= v"0.6-" && return quote
        ccall(($f,:gurobi70), $(esc(:stdcall)), $(args...))
    end
    is_unix() && return quote
        ccall(($f,libgurobi), $(args...))
    end
end

# add_vars!: adding variables to exist constraints
# call c function
function add_vars!(model::Gurobi.Model, n::Int64, matrix::Matrix{Float64}, c::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64})
  # parameters
  # input
  # model: model that variables are added into
  # n: number of variables to be added
  # matrix: the coefficients of added variables on each constraint, number of collumns = number of variable, number of rows = number of constraints
  # c: the coefficients of added variables on objective function
  # lb: the lower bounds of added variables
  # ub: the upper bounds of added variables

  sparseMatrix = sparse(matrix)
  vbeg = sparseMatrix.colptr[1:sparseMatrix.n]
  vind = sparseMatrix.rowval
  vval = sparseMatrix.nzval
  nnz = length(vval)

  # main call
  ret = @grb_ccall(addvars, Cint, (
        Ptr{Void},  # model
        Cint,       # numvars
        Cint,       # numnz
        Ptr{Cint},  # vbeg
        Ptr{Cint},  # vind
        Ptr{Float64}, # vval
        Ptr{Float64}, # obj
        Ptr{Float64}, # lb
        Ptr{Float64}, # ub
        Ptr{Cchar},   # vtypes
        Ptr{Ptr{UInt8}}, # varnames
        ),
        model,
        n,
        nnz,
        Gurobi.ivec(vbeg.-1),
        Gurobi.ivec(vind.-1),
        Gurobi.fvec(vval),
        Gurobi.fvec(c),
        Gurobi.fvecx(lb, n),
        Gurobi.fvecx(ub, n),
        Gurobi.cvecx(GRB_CONTINUOUS, n),
        C_NULL)

  if ret != 0
    throw(Gurobi.GurobiError(model.env, ret))
  end
  nothing
end

# checkOpt: checking whether the solution of small model is also the optimal solution of full-size model or not
function checkOpt(duals::Vector{Float64}, dataMatrix::Matrix{Float64}, notset::Vector{Int64}, additiveSize::Int64, tol::Float64)
  # parameter
  # input
  # duals: the dual values of small model
  # dataMatrix: data matrix of full-size problem
  # notset: the unselected DMUs, which are the indices w.r.t the variable vector
  # additiveSize: number of DMUs which are added to set and deleted from notset
  # output
  # optimal: boolean value, which represent to whether the solution of small model is also the optimal solution of full-size model or not
  # infeasibleIdx: top n most infeasible DMUs, which are the indices w.r.t. notset, n is "additiveSize"

  # initialization
  optimal = true

  infeasibleIdx = Array{Int64}(additiveSize)
  infeasibleValue = Array{Float64}(additiveSize)

  # result of a part of kkt condition, it's equal to the product of dataMatrix and duals
  # for more details, please refer to Chen and Lai (2017), "Determining radial efficiency with a large data set by solving small-size linear programs"
  kkt = dataMatrix*duals[1:length(duals)-1] + duals[length(duals)]

  # if additiveSize == 1 it represent to that the size of small model built by set is equal to the limit size
  # so it has to delete a DMU from set, then there is a space to add a new DMU into set
  if additiveSize == 1
    # finding the most infeasible DMU (with the minimun value of kkt), which is the index w.r.t notset
    infeasibleValue, infeasibleIdx = findmax(kkt[notset])

    # cheking whether the solution is optimal or not
    # if the minimun value of kkt is more than tol, then the solution is optimal
    # if it is smaller less than tol, then the solution is not optimal
    optimal = infeasibleValue[1] <= tol #boolean value

  # the set having enough space to add DMUs, the number of DMUs to be added is additiveSize
  else
    # finding top n most infeasible DMUs (with the minimun value of kkt), which is the index w.r.t notset
    infeasibleValue, infeasibleIdx = largestn(kkt[notset], additiveSize)

    # cheking whether the solution is optimal or not
    # if the minimun value of kkt is more than tol, then the solution is optimal
    # if it is smaller less than tol, then the solution is not optimal
    optimal = findmax(infeasibleValue)[1] <= tol # boolean value

    # if !optimal
    #   #not optimal, only keep those violating the constraints (> tol)
    #   a = find(infeasibleValue .> 0)
    #   infeasibleIdx = infeasibleIdx[a]
    #
    # end
  end

  return optimal, infeasibleIdx
end

# chgCoeffs!: changing the coefficients of variables on a constraint
# call c function
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
                   Gurobi.ivec(vind),
                   Gurobi.fvec(val))

  if ret != 0
    throw(Gurobi.GurobiError(model.env, ret))
  end
  nothing
end

# return sample and nonsample
function get2Sets(scale::Int64, sampleSize::Int64, k::Int64, aSet::Vector{Int64})
  # aSet is the default set including target MDU should be sorted
  # k: the target DMU index
  # if length(sample) < sampleSize, rest sample elements are randomly generated
  # the first sampleSize elements are sample; the others are nonsample
  output = collect(1:scale)
  temp = 0

  for pt = 1:sampleSize
    if pt <= length(aSet) # get element from sample
      next = aSet[pt]
    else # get element randomly
      next = rand(pt:scale)
    end
    # swap (noting in the sample)
    temp = output[pt]
    output[pt] = output[next]
    output[next] = temp
  end

  # move k to the first
  a = findfirst(output[1:length(aSet)] .== k)
  temp = output[a]  # a is a Vector
  output[a] = output[1]
  output[1] = temp

  return output[1:sampleSize], output[sampleSize+1:end]
end

# extremeValueSampling: get min of an input coumn and max of an output column
function getExtremeValueIdx(A::Matrix{Float64}, numInput::Int64)
  # function getExtremeValueIdx(scale::Int64, dataMatrix::Matrix{Float64}, numInput::Int64, numOutput::Int64)
  #input:
  #dataMatrix: the full data
  #numInput: no of inputs
  #output:
  #set: the selected DMUs, which are the indices w.r.t the variable vector
  scale, dimension = size(A)
  numOutput = dimension - numInput
  set = zeros(Int64, numInput+numOutput)
  #get idx of min value in each input column
  set[1:numInput] = findmin(A[:, 1:numInput], 1)[2]
  #get idx of max value in each output column
  set[numInput+1: end] = findmax(A[:, numInput+1: end], 1)[2]
  for i = 1:numInput+numOutput
    set[i] = ind2sub((scale,dimension), set[i])[1]
  end
  return sort(unique(set))
end

# dynamicReSampling: adding top n most infeasible DMUs into set without deleting DMUs, n is additiveSize
# move overIdx from notset to set
function regruping(toMove::Vector{Int64}, set::Vector{Int64}, notset::Vector{Int64})
  # toMove: move from notset to set, wrt to the index in notset
  set = [set; notset[toMove]]
  get2Sets!(notset, length(toMove), sort(toMove)) # so that first length(toMove) elements are to be removed
  notset = notset[length(toMove)+1:end]
  # notset[overIdx] = notset[length(notset)-length(overIdx)+1:length(notset)]
  # notset = notset[1:length(notset)-length(overIdx)]

  return set, notset
end

# swap p[i] and q[j]
function swapSample!(p::Vector{Int64}, q::Vector{Int64}, i::Int64, j::Int64)
  temp = p[i]
  p[i] = q[j]
  q[j] = temp
end
