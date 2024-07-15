module PyMYAMI

using JSON

project_path(parts...) = normpath(joinpath(@__DIR__, "..", parts...))
pymyami_coefs = JSON.parsefile(project_path("..", "coefficients", "polynomial_coefficients.json"))

function generate_polynomial_features(
    TempC::Union{AbstractFloat, Integer, AbstractArray{<:AbstractFloat}}=25.0,
    Sal::Union{AbstractFloat, Integer, AbstractArray{<:AbstractFloat}}=35.0,
    Mg::Union{AbstractFloat, Integer, AbstractArray{<:AbstractFloat}}=0.0528171,
    Ca::Union{AbstractFloat, Integer, AbstractArray{<:AbstractFloat}}=0.0102821)
    """
    Generate polynomial features for the given input data.
    """
    TempC = TempC isa AbstractArray ? TempC : fill(TempC, 1)
    Sal = Sal isa AbstractArray ? Sal : fill(Sal, 1)
    Mg = Mg isa AbstractArray ? Mg : fill(Mg, 1)
    Ca = Ca isa AbstractArray ? Ca : fill(Ca, 1)

    max_length = maximum([length(TempC), length(Sal), length(Mg), length(Ca)])
    
    TempC = length(TempC) == max_length ? TempC : repeat(TempC, outer=ceil(Int, max_length / length(TempC)))
    Sal = length(Sal) == max_length ? Sal : repeat(Sal, outer=ceil(Int, max_length / length(Sal)))
    Mg = length(Mg) == max_length ? Mg : repeat(Mg, outer=ceil(Int, max_length / length(Mg)))
    Ca = length(Ca) == max_length ? Ca : repeat(Ca, outer=ceil(Int, max_length / length(Ca)))

    TempK = TempC .+ 273.15
    X = hcat(TempK, log.(TempK), Sal, Mg, Ca)
    
    out = hcat(
        fill(1.0, max_length),
        X[:,1],
        X[:,2],
        X[:,3],
        X[:,4],
        X[:,5],
        X[:,1].^2,
        X[:,1].*X[:,2],
        X[:,1].*X[:,3],
        X[:,1].*X[:,4],
        X[:,1].*X[:,5],
        X[:,2].^2,
        X[:,2].*X[:,3],
        X[:,2].*X[:,4],
        X[:,2].*X[:,5],
        X[:,3].^2,
        X[:,3].*X[:,4],
        X[:,3].*X[:,5],
        X[:,4].^2,
        X[:,4].*X[:,5],
        X[:,5].^2,
        X[:,1].^3,
        X[:,1].^2 .*X[:,2],
        X[:,1].^2 .*X[:,3],
        X[:,1].^2 .*X[:,4],
        X[:,1].^2 .*X[:,5],
        X[:,1].*X[:,2].^2,
        X[:,1].*X[:,2].*X[:,3],
        X[:,1].*X[:,2].*X[:,4],
        X[:,1].*X[:,2].*X[:,5],
        X[:,1].*X[:,3].^2,
        X[:,1].*X[:,3].*X[:,4],
        X[:,1].*X[:,3].*X[:,5],
        X[:,1].*X[:,4].^2,
        X[:,1].*X[:,4].*X[:,5],
        X[:,1].*X[:,5].^2,
        X[:,2].^3,
        X[:,2].^2 .*X[:,3],
        X[:,2].^2 .*X[:,4],
        X[:,2].^2 .*X[:,5],
        X[:,2].*X[:,3].^2,
        X[:,2].*X[:,3].*X[:,4],
        X[:,2].*X[:,3].*X[:,5],
        X[:,2].*X[:,4].^2,
        X[:,2].*X[:,4].*X[:,5],
        X[:,2].*X[:,5].^2,
        X[:,3].^3,
        X[:,3].^2 .*X[:,4],
        X[:,3].^2 .*X[:,5],
        X[:,3].*X[:,4].^2,
        X[:,3].*X[:,4].*X[:,5],
        X[:,3].*X[:,5].^2,
        X[:,4].^3,
        X[:,4].^2 .*X[:,5],
        X[:,4].*X[:,5].^2,
        X[:,5].^3
    )
    return out
end

function approximate_Fcorr(
    TempC::Union{AbstractFloat, Integer, AbstractArray{<:AbstractFloat}}=25.0,
    Sal::Union{AbstractFloat, Integer, AbstractArray{<:AbstractFloat}}=35.0,
    Mg::Union{AbstractFloat, Integer, AbstractArray{<:AbstractFloat}}=0.0528171,
    Ca::Union{AbstractFloat, Integer, AbstractArray{<:AbstractFloat}}=0.0102821)
    """
    Approximate the correction factor for the given input data.
    """
    X = generate_polynomial_features(TempC, Sal, Mg, Ca)

    Fcorr = Dict()
    for (k, coefs) in pymyami_coefs
        Fcorr[k] = sum(X .* reshape(coefs, 1, :), dims=2)
    end

    return Fcorr

end

end