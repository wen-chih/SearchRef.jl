# loading package

#loading function
include("BuildHullPhase1.jl")
include("BuildHullPhase2.jl")

fileName = "5-2-100_01"
# for a filename a-b-c-d.csv, a is the dimension, b is the number of inpput, c is the scale and d is the density
data = readcsv("$fileName.csv") # data source.
theFrame = runDula(data, 2, 3, fileName) # "warm up", the performance is not recorded
writecsv("dulaFrame$fileName.csv")
runComputationGurobi(fileName, 2, 3) # "warm up", the performance is not recorded

# formally begins
fileName = "5-2-100k_01"
data = readcsv("$fileName.csv") # data source.
theFrame = runDula(data, 2, 3, fileName) #phase 1
writecsv("dulaFrame$fileName.csv")
runComputationGurobi(fileName, 2, 3) #phase 2

fileName = "20-10-100k_25"
data = readcsv("$fileName.csv") # data source.
theFrame = runDula(data, 10, 10, fileName) #phase 1
writecsv("dulaFrame$fileName.csv")
runComputationGurobi(fileName, 10, 10) #phase 2
