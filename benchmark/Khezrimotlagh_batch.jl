include("Khezrimotlagh.jl")

rootPath = "C:/your_path"
run_and_output("test", 2, rootPath)
run_and_output("test", 2, rootPath)

# run_and_output("5-2-100k_01", 2, rootPath)
# run_and_output("5-2-100k_10", 2, rootPath)
# run_and_output("5-2-100k_25", 2, rootPath)
#
# run_and_output("10-5-100k_01", 5, rootPath)
# run_and_output("10-5-100k_10", 5, rootPath)
# run_and_output("10-5-100k_25", 5, rootPath)
#
# run_and_output("15-7-100k_01", 7, rootPath)
# run_and_output("15-7-100k_10", 7, rootPath)
# run_and_output("15-7-100k_25", 7, rootPath)

run_and_output("20-10-100k_01", 10, rootPath)
run_and_output("20-10-100k_10", 10, rootPath)
run_and_output("20-10-100k_25", 10, rootPath)
#
run_and_output("11-5-19939_DulaBanking", 6, rootPath)
println("done!")
