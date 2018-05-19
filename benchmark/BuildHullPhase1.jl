using Gurobi

#data preprocess
function dataPreprocess(data::Array{Float64, 2}, numOfIn::Int, numOfOut::Int)

  transposeData = transpose(data)

  return [-transposeData[1:numOfIn, :]; transposeData[numOfIn + 1:numOfIn + numOfOut, :]]

end

#main function
function runDula(data::Array{Float64, 2}, numOfIn::Int, numOfOut::Int, dataName::AbstractString)

  println(dataName)
  filePath = "phase1_$dataName.txt"
  # filePath = "results/phase1/dula_$dataName.txt"

  tol = 10.0^-5

  # begining of data preprocess
  begin_TDPP = time()
  #############################################
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
  #############################################
  TDPP = time()- begin_TDPP

  # begining of building model
  begin_TModelBuilding = time()
  #############################################
  #set solver
  env = Gurobi.Env()
  setparams!(env; Method=0, OutputFlag=0, OptimalityTol=tol)
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
  #############################################
  TModelBuilding = time() - begin_TModelBuilding

  Tsolve = 0
  findEff = 0

  #begin loop
  while length(temp) > 0

    begin_Tsolve = time()
    #################################
    optimize(model)
    result = get_objval(model)
    #println(result)
    #################################
    Tsolve = Tsolve + time() - begin_Tsolve

    begin_findEff = time()
    #################################
    if result > 0

      dual = Gurobi.get_dblattrarray(model, "Pi", 1, num_constrs(model))
      beta = dual[length(dual)]

      #find
      resulFor4A = (dual[1:dimension].'*fullsizeDataSet[:, temp]) + beta
      #biggerThanOne = find(x-> x>0, resulFor4A)
      #aStar = findmax(resulFor4A[biggerThanOne])[2]
      aValue, aStar = findmax(resulFor4A)
      #maxSet = find(x -> x==aValue, resulFor4A)

      #if length(maxSet) != 1
      #  println("Warning!!!!!")
        #dont know how to solve this problem
      #  aStar = temp[findmax(fullsizeDataSet[1, temp[maxSet]])[2]]
      #  println(temp[aStar])
      #end

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
    ###################################
    findEff = findEff + time() - begin_findEff

    # println(length(temp))
  end

  f = open(filePath, "w")
  write(f, join(("Time of data preprocess: ", TDPP), " "), "\r\n")
  write(f, join(("Time of model building: ", TModelBuilding), " "), "\r\n")
  write(f, join(("Time of model solving: ", Tsolve), " "), "\r\n")
  write(f, join(("Time of finding efficient points: ", findEff), " "), "\r\n")
  write(f, join(("Total time: ", TDPP+TModelBuilding+Tsolve+findEff), ""), "\r\n")
  # write(f, join(deaFrameIdx, " "), "\r\n")

  close(f)
  return deaFrameIdx
end
