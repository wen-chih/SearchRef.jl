# it only uses Gurobi 8.0

#loading function
include("BuildHull.jl")

rootPath = "your file path"
# warm up
runBuildHull("test", 2, 3, rootPath)
runBuildHull("test", 2, 3, rootPath)

# computation begins
runBuildHull("5-2-100k_01", 2, 3, rootPath)
runBuildHull("5-2-100k_10", 2, 3, rootPath)
runBuildHull("5-2-100k_25", 2, 3, rootPath)

runBuildHull("10-5-100k_01", 5, 5, rootPath)
runBuildHull("10-5-100k_10", 5, 5, rootPath)
runBuildHull("10-5-100k_25", 5, 5, rootPath)

runBuildHull("15-7-100k_01", 7, 8, rootPath)
runBuildHull("15-7-100k_10", 7, 8, rootPath)
runBuildHull("15-7-100k_25", 7, 8, rootPath)

runBuildHull("20-10-100k_01", 10, 10, rootPath)
runBuildHull("20-10-100k_10", 10, 10, rootPath)
runBuildHull("20-10-100k_25", 10, 10, rootPath)
