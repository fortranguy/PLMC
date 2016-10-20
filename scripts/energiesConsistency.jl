if size(ARGS, 1) == 0
    error("Please provide a energies.out files.")
end

relativeError(expr, ref) = abs((ref - expr) ./ ref)

for iFile in 1:size(ARGS, 1)
    println(ARGS[iFile])
    file = readdlm(ARGS[iFile])
    println(relativeError(file[end-1, 2:end], file[end, 2:end]))
end
