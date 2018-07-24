# it only uses Gurobi 8.0

#loading function
include("SearchRef.jl")

rootPath = "your file path"

runSearchRef("test", 2, rootPath)
runSearchRef("test", 2, rootPath)

# computation begins
runSearchRef("5-2-100k_01", 2, rootPath)
runSearchRef("5-2-100k_10", 2, rootPath)
runSearchRef("5-2-100k_25", 2, rootPath)

runSearchRef("10-5-100k_01", 5, rootPath)
runSearchRef("10-5-100k_10", 5, rootPath)
runSearchRef("10-5-100k_25", 5, rootPath)

runSearchRef("15-7-100k_01", 7, rootPath)
runSearchRef("15-7-100k_10", 7, rootPath)
runSearchRef("15-7-100k_25", 7, rootPath)

runSearchRef("20-10-100k_01", 10, rootPath)
runSearchRef("20-10-100k_10", 10, rootPath)
runSearchRef("20-10-100k_25", 10, rootPath)
