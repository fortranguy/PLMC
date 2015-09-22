import JSON
json = JSON

if size(ARGS, 1) == 0
    error("Please provide a .json file.")
end
input_data = json.parsefile(ARGS[1]; ordered=false, use_mmap=true)
if !input_data["Particles"]["are dipolar"]
    error("Particles aren't dipolar: orientations are useless (yet).")
end
num_particles = input_data["Particles"]["number"]
orientations = randn(num_particles, 3)
for i=1:size(orientations, 1)
    orientations[i, :] = orientations[i, :] / norm(orientations[i, :])
end

outputFile = "orientations.in"
writedlm(outputFile, orientations)
println("Positions written in ", outputFile)
