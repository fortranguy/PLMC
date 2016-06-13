if size(ARGS, 1) == 0
    error("Please provide a file with observables.")
end
data = readdlm(ARGS[1])[find(x -> x >=0, readdlm(ARGS[1])[:, 1]), 2:end]
println("mean:", mean(data, 1))
println("standard deviation:", sqrt(var(data, 1)))
