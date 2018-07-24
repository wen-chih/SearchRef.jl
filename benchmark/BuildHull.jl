using Gurobi

function runBuildHull(fileName::AbstractString, inputSize::Int64, outputSize::Int64, root::AbstractString)
  println("\n$fileName begins @ $(Dates.hour(now())):$(Dates.minute(now())):$(Dates.second(now()))")

  # loading data
  data = readcsv("$root/RichmondDataSet/$fileName.csv") # loading input and output data
  # phase 1

  theFrame, tP1 = findFrame(data, inputSize, outputSize, "$root/BuildHull/phase1/$(fileName)P1.txt") #phase 1
  writecsv("$root/BuildHull/phase1/dulaFrame$fileName.csv", theFrame)
  # phase 2
  tP2 = runComputationGurobi(fileName, inputSize, outputSize, root) #phase 2

  # output
  resultPath = "$root/BuildHull/phase2/BuildHullP2_$fileName.csv" # output path
  resultFile = open(resultPath, "a")
  write(resultFile, "$(tP1+tP2)\r\n")
  close(resultFile)

end

# finding frame in Phase 1
function findFrame(data::Array{Float64, 2}, numOfIn::Int64, numOfOut::Int64, filePath::AbstractString)

  tol = 10.0^-5

  # ==========================================
  # begining of data preprocess
  tic() # beginning of tDataProcess
  # data preprocess
  fullsizeDataSet = dataPreprocess(data, numOfIn, numOfOut)

  # data information
  dimension, scale = size(fullsizeDataSet)
  dimensionOne = ones(dimension)

  # initialization
  temp = collect(1:scale)      # temporary index of full size data
  deaFrameIdx = findmax(fullsizeDataSet[1, :])[2]     # select a extreme point to initialize the frame element
  count = 1
  temp[deaFrameIdx] = temp[length(temp)]
  temp = temp[1:scale-count]

  tDataProcess = toq() # data preprocessing time
  # end of data preprocess

  # ==========================================
  # begining of building model
  tic() # beginning of tModelBuilding
  #set solver
  env = Gurobi.Env()
  setparams!(env; Method=0, OutputFlag=0, OptimalityTol=tol)
  # setparams!(env; Method=0, OutputFlag=0, OptimalityTol=tol)
  # setparams!(env; OutputFlag=0)

  # initialization
  bIdx = temp[1]
  frameNum = length(deaFrameIdx)
  frameZero = zeros(frameNum)
  frameOne = ones(frameNum)

  # building model
  model = Gurobi.Model(env, "lp", :minimize)

  # defining variables
  add_cvar!(model, 1, -Inf, Inf)
  add_cvars!(model, frameZero, frameZero, Inf)
  update_model!(model)

  # defining constraints
  add_constrs!(model, [dimensionOne fullsizeDataSet[:, deaFrameIdx]], '>', fullsizeDataSet[:, bIdx])
  add_constr!(model, [0; frameOne], '=', 1.0)
  update_model!(model)

  tModelBuilding = toq()
  # end of building model

  tSolve = 0.0
  tFindEff = 0.0

  #begin loop
  while length(temp) > 0

    # solving LP
    tic() # gegin of tSolve
    optimize(model)
    result = get_objval(model)
    tSolve += toq()

    tic()  # begin of tFindEff
    if result > 0

      dual = Gurobi.get_dblattrarray(model, "Pi", 1, num_constrs(model))
      beta = dual[length(dual)]

      #find
      resulFor4A = (dual[1:dimension].'*fullsizeDataSet[:, temp]) + beta
      aValue, aStar = findmax(resulFor4A)

      #update dea frame
      deaFrameIdx = [deaFrameIdx; temp[aStar]]

      #update work-space
      count = count + 1
      temp[aStar] = temp[length(temp)]
      temp = temp[1:scale-count]

      if length(temp) > 0
        bIdx = temp[1]
      end

      cind = [find( x -> x != 0.0, fullsizeDataSet[:, deaFrameIdx[length(deaFrameIdx)]]); num_constrs(model)]
      nnz = length(cind)
      vval = [fullsizeDataSet[cind[1:length(cind)-1], deaFrameIdx[length(deaFrameIdx)]]; 1.0]

      # adding new variable
      add_var!(model, nnz, cind, vval, 0.0, 0.0, Inf, GRB_CONTINUOUS)
      update_model!(model)

      # changing RHS of model
      Gurobi.set_dblattrarray!(model, "RHS", 1, num_constrs(model)-1, fullsizeDataSet[:, bIdx])
      update_model!(model)
    else

      #update work-space
      count = count + 1
      temp[1] = temp[length(temp)]
      temp = temp[1:scale-count]
      if length(temp) > 0
        bIdx = temp[1]
      end

      # changing RHS of model
      Gurobi.set_dblattrarray!(model, "RHS", 1, num_constrs(model)-1, fullsizeDataSet[:, bIdx])
      update_model!(model)
    end

    tFindEff += toq()

  end

  println("phase 1 done! $(tDataProcess+tModelBuilding+tSolve+tFindEff)")

  f = open(filePath, "w")
  write(f, "Time of data preprocess: $tDataProcess\r\n")
  write(f, "Time of model building: $tModelBuilding\r\n")
  write(f, "Time of model solving: $tSolve\r\n")
  write(f, "Time of finding efficient points: $tFindEff\r\n")
  write(f, "Total time: $(tDataProcess+tModelBuilding+tSolve+tFindEff)\r\n")
  close(f)

  return deaFrameIdx, tDataProcess+tModelBuilding+tSolve+tFindEff
end

#data preprocess in Phase 1
function dataPreprocess(data::Array{Float64, 2}, numOfIn::Int64, numOfOut::Int64)

  transposeData = transpose(data)
  return [-transposeData[1:numOfIn, :]; transposeData[numOfIn + 1:numOfIn + numOfOut, :]]

end


# Phase 2
# the object of this code is to record the time spent to solve each problem by algorithm
# the second phase of BuildHull
# model: standard input-oriented VRS
# algorithm: full-size
# solver: "Gurobi" with method "primal simplex"
# runComputationGurobi: input file name and the number of inputs and outputs to use algorithm to run all computation by Gurobi
function runComputationGurobi(fileName::AbstractString, numInput::Int64, numOutput::Int64, root::AbstractString)
  # parameters
  # input
  # fileName: name of the file run
  # numInput: number of inputs of each DMU
  # numOutput: number of outputs of each DMU

  println("\n$fileName Phase 2 begins @ $(Dates.hour(now())):$(Dates.minute(now())):$(Dates.second(now()))")

  tol = 10.0^-5

  # loading data
  # dataMatrix = readcsv("C:/Users/mb516/Documents/Wen-Chih/SearchRef-2018/RichmondDataSet/$fileName.csv") # loading input and output data
  dataMatrix = readcsv("$root/RichmondDataSet/$fileName.csv") # loading input and output data
  deaFrame = convert(Array{Int64,2}, readcsv("$root/BuildHull/phase1/dulaFrame$fileName.csv"))[:,1]
  frameNum = length(deaFrame)

  # path of output file recording information of time and solution
  resultPath = "$root/BuildHull/phase2/BuildHullP2_$fileName.csv" # output path

  # computation begins

  tic()
  # setting gurobi parameter
  env = Gurobi.Env()
  setparams!(env; Method=0, OutputFlag=0, OptimalityTol=tol) #using primal simplex method and not displaying solver's status
  # setparams!(env; OutputFlag=0) #using non-deterministic concurrent (default) method and not displaying solver's status

  # data preprocess
  transDataMatrix = transpose(dataMatrix)
  dimension, scale = size(transDataMatrix)
  # frameMatrix = transDataMatrix[:, deaFrame] # only takeing frame, not fall data, for computation

  frameZero = zeros(frameNum) # the vector all zero with length equal to scale
  frameOne = ones(frameNum) # the vector all one with length equal to scale
  inputZero = zeros(numInput)
  outputZero = zeros(numOutput)

  # recoding time initialization
  tSolving = 0.0
  tChgRHS = 0.0
  tIdentifyEff = 0.0
  objValue = 0.0
  tChgRHS4k = 0.0
  tSolving4k = 0.0

  # building model
  model = Gurobi.Model(env, "lp", :minimize)

  # defining variables
  add_cvars!(model, frameZero, frameZero, Inf)  # defining lambda(r), r from 1 to scale
  add_cvar!(model, 1.0, -Inf, Inf) # defining theta
  update_model!(model)

  # defining constraints
  # input constraints
  # add_constrs!(model, [transDataMatrix[1:numInput, deaFrame] -transDataMatrix[1:numInput, 1]], '<', inputZero)
  # add_constrs_t!(model, At, rel, rhs)  # here At can be dense or sparse
  add_constrs_t!(model, [dataMatrix[deaFrame, 1:numInput]; -reshape(dataMatrix[1, 1:numInput], 1, numInput)], '<', inputZero)

  # output constraints
  # add_constrs!(model, [transDataMatrix[numInput+1:numInput+numOutput, deaFrame] outputZero], '>', transDataMatrix[numInput+1:numInput+numOutput, 1])
  add_constrs_t!(model, [dataMatrix[deaFrame, numInput+1:numInput+numOutput]; zeros(1, numOutput)], '>', dataMatrix[1, numInput+1:dimension])

  # convexity constraint
  add_constr!(model, [frameOne; 0], '=', 1.0)
  update_model!(model)

  tBuilding = toq()
  # end of initialization

  # solving 1st LP model
  tic()
  optimize(model)
  tSolving = toq()

  resultFile = open(resultPath, "w")
  write(resultFile, "1,$(get_objval(model))\r\n")
  close(resultFile)


  for k = 2:scale

    if length(find(x -> x == k, deaFrame)) != 0 # DMU in the frame, efficiency = 1
      objValue = 1.0
      resultFile = open(resultPath, "a")
      write(resultFile, "$k,$objValue\r\n")
      close(resultFile)

    else
      # update RHS and solve LP again
      objValue, tChgRHS4k, tSolving4k = updateAndSolve(model, k, frameNum, numInput, numOutput, transDataMatrix)
      tChgRHS += tChgRHS4k
      tSolving += tSolving4k

      resultFile = open(resultPath, "a")
      write(resultFile, "$k,$objValue\r\n")
      close(resultFile)
    end
  end

  println("phase 2 done! $(tBuilding+tSolving+tChgRHS)")
  resultFile = open(resultPath, "a")
  write(resultFile, "$(tBuilding+tSolving+tChgRHS),$tBuilding,$tSolving,$tChgRHS\r\n")
  close(resultFile)

  return tBuilding+tSolving+tChgRHS
end


# update RHS and solve LP
function updateAndSolve(model::Gurobi.Model, k::Int64, scale::Int64, numInput::Int64, numOutput::Int64, transDataMatrix::Matrix{Float64})

  tic()
  # changing RHS of model
  Gurobi.set_dblattrarray!(model, "RHS", numInput+1, numOutput, transDataMatrix[numInput+1:numInput+numOutput, k])
  # changing coefficients
  cind = find(transDataMatrix[1:numInput, k])
  vind = fill(scale+1, length(cind))
  vval = -transDataMatrix[cind, k]

  chgCoeffs!(model, length(cind), cind,  vind, vval)
  update_model!(model)
  tChgRHS = toq()

  tic()
  # solving model
  optimize(model)
  tSolving = toq()

  return get_objval(model), tChgRHS, tSolving
end

# ========================================================================
# a macro of gurobi c call function
macro grb_ccall(func, args...)
    f = "GRB$(func)"
    args = map(esc,args)

    is_windows() && VERSION < v"0.6-" && return quote
        ccall(($f,:gurobi80), stdcall, $(args...))
    end
    is_windows() && VERSION >= v"0.6-" && return quote
        ccall(($f,:gurobi80), $(esc(:stdcall)), $(args...))
    end
    is_unix() && return quote
        ccall(($f,libgurobi), $(args...))
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
