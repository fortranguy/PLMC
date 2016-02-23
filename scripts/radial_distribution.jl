if size(ARGS, 1) == 0
    error("Please provide a configuration file.")
end
files = ARGS

#todo: json input
max_distance = 20.
delta = 0.1
distrib_i = zeros(Float64, round(Int, max_distance/delta))
distrib = zeros(distrib_i)

for file in files
    pos = readdlm(file)[:, 1:3]
    distrib_i = zeros(distrib_i)
    for i_part=1:size(pos, 1)
        for j_part=i_part+1:size(pos, 1)
            dist_ij = norm(pos[j_part, :] - pos[i_part, :]) # todo: real PBC
            i_dist = round(Int, dist_ij/delta) + 1
            distrib_i[i_dist] = distrib_i[i_dist] + 1
        end
    end
    distrib_i = distrib_i / size(pos, 1)
    distrib = distrib + distrib_i
end

writedlm("test.out", distrib)
