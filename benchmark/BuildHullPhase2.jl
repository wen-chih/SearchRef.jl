# the object of this code is to record the time spent to solve each problem by algorithm
# the second phase of BuildHull
# model: standard input-oriented VRS
# algorithm: full-size
# solver: "Gurobi" with method "primal simplex"
using Gurobi

# runComputationGurobi: input file name and the number of inputs and outputs to use algorithm to run all computation by Gurobi
function runComputationGurobi(fileName::AbstractString, numInput::Int, numOutput::Int)
  # parameters
  # input
  # fileName: name of the file run
  # numInput: number of inputs of each DMU
  # numOutput: number of outputs of each DMU

  println(fileName)

  tol = 10.0^-5

  # loading data
  dataMatrix = readcsv("$fileName.csv") # loading input and output data
  deaFrame = readcsv("dulaFrame$fileName.csv")
  frameNum = length(deaFrame)

  # opening the file to record the result
  # recording information of time and solution
  filePath = "Phase2_$fileName.csv"
  # filePath = "results/phase2/Phase2_$fileName.csv"
  resultFile = open(filePath, "w")

  # resultFile = open(join(["versionSpeedUp/dula/phase2/standard/result/standard_VRS_dulaPhase2_", fileName, ".csv"], ""), "w")

  beginProcess = time()
  # setting gurobi parameter
  env = Gurobi.Env()
  setparams!(env; Method=0, OutputFlag=0, OptimalityTol=tol) #using primal simplex method and not displaying solver's status
  # setparams!(env; OutputFlag=0) #using non-deterministic concurrent (default) method and not displaying solver's status

  # data preprocess
  transDataMatrix = transpose(dataMatrix)
  dimension, scale = size(transDataMatrix)
  frameMatrix = transDataMatrix[:, deaFrame] # only takeing frame, not fall data, for computation

  frameZero = zeros(frameNum) # the vector all zero with length equal to scale
  frameOne = ones(frameNum) # the vector all one with length equal to scale
  inputZero = zeros(numInput)
  outputZero = zeros(numOutput)

  # recoding time
  tSolving = 0.0
  tChgRHS = 0.0
  tIdentifyEff = 0.0

  # building model
  model = Gurobi.Model(env, "lp", :minimize)

  # defining variables
  add_cvars!(model, frameZero, frameZero, Inf)  # defining lambda(r), r from 1 to scale
  add_cvar!(model, 1.0, -Inf, Inf) # defining theta
  update_model!(model)

  # defining constraints
  # input constraints
  add_constrs!(model, [frameMatrix[1:numInput, :] -transDataMatrix[1:numInput, 1]], '<', inputZero)
  # output constraints
  add_constrs!(model, [frameMatrix[numInput+1:numInput+numOutput, :] outputZero], '>', transDataMatrix[numInput+1:numInput+numOutput, 1])
  # convexity constraint
  add_constr!(model, [frameOne; 0], '=', 1.0)
  update_model!(model)
  # end building model
  tBuilding = time() - beginProcess

  beginSolving = time()
  # solving model
  optimize(model)
  tSolving = tSolving + time() - beginSolving

  write(resultFile, join((1, get_objval(model)), ", "), "\r\n")

  for k = 2:scale
    objValue = 1.0

    if length(find(x -> x == k, deaFrame)) != 0 # DMU in the frame, efficiency = 1
      objValue = 1.0
      write(resultFile, join((k, objValue), ", "), "\r\n")
    else
      beginChgRHS = time()
      # changing RHS of model
      Gurobi.set_dblattrarray!(model, "RHS", numInput + 1, numOutput, transDataMatrix[numInput+1:numInput+numOutput, k])
      # changing coefficients
      cind = find(transDataMatrix[1:numInput, k])
      vind = fill(frameNum+1, length(cind))
      vval = -transDataMatrix[cind, k]

      chgCoeffs!(model, length(cind), cind,  vind, vval)
      update_model!(model)
      tChgRHS = tChgRHS + time() - beginChgRHS

      beginSolving = time()
      # solving model
      optimize(model)
      objValue = get_objval(model)
      tSolving = tSolving + time() - beginSolving

      write(resultFile, join((k, objValue), ", "), "\r\n")
    end
  end

  write(resultFile, join((tBuilding + tSolving + tChgRHS, tBuilding, tSolving, tChgRHS), ", "), "\r\n")

  close(resultFile)
end

========================================================================
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

# chgCoeffs!: changing the coefficients of variables on a constraint
function chgCoeffs!(model::Gurobi.Model, n::Int, cind::Vector{Int}, vind::Vector{Int}, val::Vector{Float64})
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
