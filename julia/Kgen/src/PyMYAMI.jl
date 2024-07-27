module PyMYAMI

export approximate_Fcorr

using JSON

project_path(parts...) = normpath(joinpath(@__DIR__, "..", parts...))
pymyami_coefs = JSON.parsefile(project_path("..", "coefficients", "polynomial_coefficients.json"))

function generate_polynomial_features(
    TempC::AbstractFloat=25.0,
    Sal::AbstractFloat=35.0,
    Mg::AbstractFloat=0.0528171,
    Ca::AbstractFloat=0.0102821)
    """
    Generate polynomial features for the given input data.
    """

    TempK = TempC + 273.15
    X = hcat(TempK, log(TempK), Sal, Mg, Ca)
    
    out = hcat(
        1.0,
        X[1],
        X[2],
        X[3],
        X[4],
        X[5],
        X[1]^2,
        X[1]*X[2],
        X[1]*X[3],
        X[1]*X[4],
        X[1]*X[5],
        X[2]^2,
        X[2]*X[3],
        X[2]*X[4],
        X[2]*X[5],
        X[3]^2,
        X[3]*X[4],
        X[3]*X[5],
        X[4]^2,
        X[4]*X[5],
        X[5]^2,
        X[1]^3,
        X[1]^2 *X[2],
        X[1]^2 *X[3],
        X[1]^2 *X[4],
        X[1]^2 *X[5],
        X[1]*X[2]^2,
        X[1]*X[2]*X[3],
        X[1]*X[2]*X[4],
        X[1]*X[2]*X[5],
        X[1]*X[3]^2,
        X[1]*X[3]*X[4],
        X[1]*X[3]*X[5],
        X[1]*X[4]^2,
        X[1]*X[4]*X[5],
        X[1]*X[5]^2,
        X[2]^3,
        X[2]^2 *X[3],
        X[2]^2 *X[4],
        X[2]^2 *X[5],
        X[2]*X[3]^2,
        X[2]*X[3]*X[4],
        X[2]*X[3]*X[5],
        X[2]*X[4]^2,
        X[2]*X[4]*X[5],
        X[2]*X[5]^2,
        X[3]^3,
        X[3]^2 *X[4],
        X[3]^2 *X[5],
        X[3]*X[4]^2,
        X[3]*X[4]*X[5],
        X[3]*X[5]^2,
        X[4]^3,
        X[4]^2 *X[5],
        X[4]*X[5]^2,
        X[5]^3
    )
    return out
end

function approximate_Fcorr(;
    K::String,
    TempC::AbstractFloat=25.0,
    Sal::AbstractFloat=35.0,
    Mg::AbstractFloat=0.0528171,
    Ca::AbstractFloat=0.0102821)
    """
    Approximate the correction factor for the given input data.
    """
    X = generate_polynomial_features(TempC, Sal, Mg, Ca)

    (X * pymyami_coefs[K])[1]
    
end

function approximate_Fcorrs(;
    TempC::AbstractFloat=25.0,
    Sal::AbstractFloat=35.0,
    Mg::AbstractFloat=0.0528171,
    Ca::AbstractFloat=0.0102821)
    """
    Approximate the correction factor for the given input data.
    """
    X = generate_polynomial_features(TempC, Sal, Mg, Ca)

    Fcorr = Dict()
    for (k, coefs) in pymyami_coefs
        Fcorr[k] = (X * coefs)[1]
    end

    return Fcorr

end

end