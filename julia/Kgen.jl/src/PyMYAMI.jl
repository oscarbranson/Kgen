"""
    PyMYAMI

Polynomial approximation of the MyAMI seawater-composition correction factors.

The coefficients in `coefficients/polynomial_coefficients.json` are exported from
[pymyami](https://github.com/PalaeoCarb/pymyami) and are shared verbatim with the R and
MATLAB implementations. They were fitted by `sklearn.preprocessing.PolynomialFeatures`
with `degree=3`, so the 56 terms must be generated in that function's exact order —
see [`generate_polynomial_features`](@ref).
"""
module PyMYAMI

using JSON

export approximate_seawater_correction, approximate_seawater_corrections

"Ks for which pymyami provides a polynomial approximation."
const POLY_NAMES = (:K0, :K1, :K2, :KW, :KB, :KS, :KspA, :KspC)

const N_TERMS = 56

function _load_polynomial_coefficients()
    path = joinpath(@__DIR__, "coefficients", "polynomial_coefficients.json")
    # See the note in Kgen._load_coefficients: without this, refreshing the coefficients via
    # update_pymyami.py would not take effect until the precompile cache happened to be cleared.
    include_dependency(path)
    raw = JSON.parsefile(path)
    NamedTuple{POLY_NAMES}(
        Tuple(NTuple{N_TERMS,Float64}(Float64.(raw[String(name)])) for name in POLY_NAMES)
    )
end

const POLY_COEFS = _load_polynomial_coefficients()

"""
    generate_polynomial_features(temp_c, sal, magnesium, calcium) -> NTuple{56,Float64}

Build the degree-3 polynomial feature vector in `(temp_k, log(temp_k), sal, magnesium,
calcium)` that the pymyami coefficients were fitted against.

The term order reproduces `sklearn.preprocessing.PolynomialFeatures(degree=3)`: the bias
term, then the five linear terms, then all degree-2 combinations *with replacement* in
index order, then all degree-3 combinations. Changing this order silently invalidates
every correction factor, so it is written out longhand rather than generated.
"""
function generate_polynomial_features(temp_c::Real, sal::Real, magnesium::Real, calcium::Real)
    temp_k = temp_c + 273.15
    a, b, c, d, e = float(temp_k), log(temp_k), float(sal), float(magnesium), float(calcium)

    return (
        1.0,
        a, b, c, d, e,
        a*a, a*b, a*c, a*d, a*e,
        b*b, b*c, b*d, b*e,
        c*c, c*d, c*e,
        d*d, d*e,
        e*e,
        a*a*a, a*a*b, a*a*c, a*a*d, a*a*e,
        a*b*b, a*b*c, a*b*d, a*b*e,
        a*c*c, a*c*d, a*c*e,
        a*d*d, a*d*e,
        a*e*e,
        b*b*b, b*b*c, b*b*d, b*b*e,
        b*c*c, b*c*d, b*c*e,
        b*d*d, b*d*e,
        b*e*e,
        c*c*c, c*c*d, c*c*e,
        c*d*d, c*d*e,
        c*e*e,
        d*d*d, d*d*e,
        d*e*e,
        e*e*e,
    )
end

"Dot product of the feature and coefficient tuples. `mapreduce` over tuples unrolls, so
this compiles to a straight-line multiply-add chain with no intermediate allocation."
_evaluate(features::NTuple{N_TERMS,Float64}, coefficients::NTuple{N_TERMS,Float64}) =
    mapreduce(*, +, features, coefficients)

"""
    approximate_seawater_correction(K; temp_c=25.0, sal=35.0, magnesium=0.0528171, calcium=0.0102821)

Approximate the MyAMI correction factor for a single equilibrium constant `K`.

Returns `1.0` for Ks that pymyami does not correct (see [`POLY_NAMES`](@ref)), so the
result can always be multiplied through unconditionally.
"""
function approximate_seawater_correction(K::Symbol;
                                         temp_c::Real=25.0,
                                         sal::Real=35.0,
                                         magnesium::Real=0.0528171,
                                         calcium::Real=0.0102821)
    K in POLY_NAMES || return 1.0
    features = generate_polynomial_features(temp_c, sal, magnesium, calcium)
    return _evaluate(features, POLY_COEFS[K])
end

"""
    approximate_seawater_corrections(; temp_c=25.0, sal=35.0, magnesium=0.0528171, calcium=0.0102821)

Approximate the MyAMI correction factors for all Ks pymyami covers, as a `NamedTuple`
keyed by [`POLY_NAMES`](@ref).
"""
function approximate_seawater_corrections(; temp_c::Real=25.0,
                                            sal::Real=35.0,
                                            magnesium::Real=0.0528171,
                                            calcium::Real=0.0102821)
    features = generate_polynomial_features(temp_c, sal, magnesium, calcium)
    # Written out rather than `map`ped over POLY_COEFS: mapping a closure over a NamedTuple
    # of 56-element tuples defeats inlining and heap-allocates the result.
    return (
        K0   = _evaluate(features, POLY_COEFS.K0),
        K1   = _evaluate(features, POLY_COEFS.K1),
        K2   = _evaluate(features, POLY_COEFS.K2),
        KW   = _evaluate(features, POLY_COEFS.KW),
        KB   = _evaluate(features, POLY_COEFS.KB),
        KS   = _evaluate(features, POLY_COEFS.KS),
        KspA = _evaluate(features, POLY_COEFS.KspA),
        KspC = _evaluate(features, POLY_COEFS.KspC),
    )
end

end # module