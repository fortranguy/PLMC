module PLMC

    type Component
        num :: Int64
        positions :: Array{Float64, 2}
        orientations :: Array{Float64, 2}
    end

    function folded(x::Array{Float64, 1}, box_size::Array{Float64, 1})
        v = mod(x, box_size)
        for i=1:3
            if (v[i] > box_size[i]/2)
                v[i] -= box_size[i]
            end
        end
        return v
    end

    function distance(x::Array{Float64, 1}, y::Array{Float64, 1}, box_size::Array{Float64, 1})
        return norm(folded(y - x, box_size))
    end

end
