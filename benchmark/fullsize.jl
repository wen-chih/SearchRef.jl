# the object of this code is to record the time spent to solve each problem by algorithm
# model: standard VRS
# algorithm: full-size
# solver: "Gurobi" with method "primal simplex"

using Gurobi


# runComputationGurobi: input file name and the number of inputs and outputs to use algorithm to run all computation by Gurobi
function runComputationGurobi(fileName::AbstractString, numInput::Int64, numOutput::Int64, root::AbstractString)
  # parameters
  # input
  # fileName: name of the file run
  # numInput: number of inputs of each DMU
  # numOutput: number of outputs of each DMU

  println("$fileName>> $(Dates.hour(now())):$(Dates.minute(now())):$(Dates.second(now()))")
  tol = 10.0^-5
  # loading data
  dataMatrix = readcsv("$root/RichmondDataSet/$fileName.csv") # loading input and output data
  # recording information of time and solution
  filePath = "$root/fullsize/fullsize$fileName.csv"


  # computation begins

  tic()
  # setting gurobi parameter
  env = Gurobi.Env()
  setparams!(env; Method=0, OutputFlag=0, OptimalityTol=tol) #using primal simplex method and not displaying solver's status
  # setparams!(env; OutputFlag=0) #using non-deterministic concurrent (default) method and not displaying solver's status

  # data preprocess
  transDataMatrix = transpose(dataMatrix)
  dimension, scale = size(transDataMatrix)

  lambdaZero = zeros(scale) # the vector all zero with length equal to scale
  lambdaOne = ones(scale) # the vector all one with length equal to scale
  inputZero = zeros(numInput)
  outputZero = zeros(numOutput)

  # initialization
  tChgRHS4k = 0.0
  tSolving4k = 0.0
  objValue = 0.0

  # recoding time
  tSolving = 0.0
  tChgRHS = 0.0
  tIdentifyEff = 0.0

  # building model
  model = Gurobi.Model(env, "lp", :minimize)

  # defining variables
  add_cvars!(model, lambdaZero, lambdaZero, Inf)  # defining lambda(r), r from 1 to scale
  add_cvar!(model, 1.0, -Inf, Inf) # defining theta
  update_model!(model)

  # defining constraints
  # input constraints
  add_constrs!(model, [transDataMatrix[1:numInput, :] -transDataMatrix[1:numInput, 1]], '<', inputZero)
  # output constraints
  add_constrs!(model, [transDataMatrix[numInput+1:numInput+numOutput, :] outputZero], '>', transDataMatrix[numInput+1:numInput+numOutput, 1])
  # convexity constraint
  add_constr!(model, [lambdaOne; 0], '=', 1.0)
  update_model!(model)
  # end building model
  tBuilding = toq()

  tic()
  # solving model
  optimize(model)
  tSolving += toq()


  resultFile = open(filePath, "w")
  write(resultFile, "1,$(get_objval(model))\r\n")
  close(resultFile)

  for k = 2:scale

    # update RHS and solve LP again
    objValue, tChgRHS4k, tSolving4k = updateAndSolve(model, k, scale, numInput, numOutput, transDataMatrix)
    tChgRHS += tChgRHS4k
    tSolving += tSolving4k

    resultFile = open(filePath, "a")
    write(resultFile, "$k,$objValue\r\n")
    close(resultFile)

  end

  println("done! $(tBuilding+tSolving+tChgRHS)")
  resultFile = open(filePath, "a")
  # write(resultFile, join((tBuilding + tSolving + tChgRHS, tBuilding, tSolving, tChgRHS), ", "), "\r\n")
  write(resultFile, "$(tBuilding+tSolving+tChgRHS),$tBuilding,$tSolving,$tChgRHS\r\n")
  close(resultFile)
end

# update RHS and solve LP
function updateAndSolve(model::Gurobi.Model, k::Int64, scale::Int64, numInput::Int64, numOutput::Int64, transDataMatrix::Matrix{Float64})

  tic()
  # changing RHS of model
  Gurobi.set_dblattrarray!(model, "RHS", numInput + 1, numOutput, transDataMatrix[numInput+1:numInput+numOutput, k])
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
