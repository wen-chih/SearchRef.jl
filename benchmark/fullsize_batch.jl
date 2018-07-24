# it only uses Gurobi 8.0

include("fullsize.jl")

rootPath = "your file path"
# warm up
runComputationGurobi("test", 2, 3, rootPath)
runComputationGurobi("test", 2, 3, rootPath)

# computation begins
runComputationGurobi("5-2-100k_01", 2, 3, rootPath)
runComputationGurobi("5-2-100k_10", 2, 3, rootPath)
runComputationGurobi("5-2-100k_25", 2, 3, rootPath)

runComputationGurobi("10-5-100k_01", 5, 5, rootPath)
runComputationGurobi("10-5-100k_10", 5, 5, rootPath)
runComputationGurobi("10-5-100k_25", 5, 5, rootPath)

runComputationGurobi("15-7-100k_01", 7, 8, rootPath)
runComputationGurobi("15-7-100k_10", 7, 8, rootPath)
runComputationGurobi("15-7-100k_25", 7, 8, rootPath)

runComputationGurobi("20-10-100k_01", 10, 10, rootPath)
runComputationGurobi("20-10-100k_10", 10, 10, rootPath)
runComputationGurobi("20-10-100k_25", 10, 10, rootPath)
