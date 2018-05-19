
using Gurobi
include("SearchRef.jl")

function batch()

  resultFilePath = "result.csv"
  dataPath = "5-2-100k_01.csv" # input-output data source

  dataMatrix = readcsv(dataPath)
  scale, dimension = size(dataMatrix)

  GurobiTol = 10.0^-6
  efficiency = 0.0
  iterNum = 0
  numInput = 2

  env = Gurobi.Env()
  setparams!(env; OutputFlag=0, OptimalityTol=GurobiTol) #using non-deterministic concurrent (default) method and not displaying solver's status

  SearchRef(1, dataMatrix, 2, env) # "warm up" for Julia, performance not recorded

  aSet = Array{Int64}(0)

  startingTime = time()

  # computing efficiencies for all DMUs
  @time for k = 1:scale

    efficiency, iterNum, aSet = SearchRef(k, dataMatrix, numInput, env, incrementSize = 100, extremeValueSet=aSet, initialSampling = "extreme value")

    # resultFile = open(resultFilePath, "a")
    # write(resultFile, "$k,$efficiency,$iterNum\r\n")
    # close(resultFile)

  end
  b = time()-startingTime

  println("done!  $b")
  # # recording time
  # resultFile = open(resultFilePath, "a")
  # write(resultFile, "$runTime\r\n")
  # close(resultFile)
  # #--------------------------

end


batch()
