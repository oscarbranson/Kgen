"""
    Kgen

Calculate the stoichiometric equilibrium constants of the marine carbonate system at a
given temperature, salinity and pressure, optionally corrected for non-modern seawater
Mg/Ca composition.

This is the Julia implementation of [Kgen](https://github.com/PalaeoCarb/Kgen); it is kept
numerically identical to the Python, R and MATLAB implementations by the cross-language
test harness in `crosscheck/`.

The core functions are scalar and allocation-free. To evaluate many conditions, broadcast:

```julia
calc_K.(:K1, temp_c=[0.0, 10.0, 25.0], sal=35.0)
```
"""
module Kgen

include("PyMYAMI.jl")

using JSON

export calc_K, calc_Ks

"""
Ks calculated by this package, in the order they are returned by [`calc_Ks`](@ref).
"""
const K_NAMES = (:K0, :K1, :K2, :KW, :KB, :KS, :KspA, :KspC, :KP1, :KP2, :KP3, :KSi, :KF)

const MODERN_MAGNESIUM = 0.0528171  # mol/kgsw
const MODERN_CALCIUM = 0.0102821    # mol/kgsw

function _load_coefficients(filename)
    path = joinpath(@__DIR__, "coefficients", filename)
    # The coefficients are baked into a `const` at precompile time, so Julia must be told to
    # invalidate the cache when the JSON changes - otherwise edits are silently ignored.
    include_dependency(path)
    raw = JSON.parsefile(path)["coefficients"]
    NamedTuple{K_NAMES}(Tuple(Tuple(Float64.(raw[String(name)])) for name in K_NAMES))
end

function _load_fundamental_constants()
    path = joinpath(@__DIR__, "coefficients", "fundamental_constants.json")
    include_dependency(path)
    return JSON.parsefile(path)["coefficients"]
end

"""
Gas constant in cm³ bar / (mol K), read from `fundamental_constants.json`.

Note that Python and R hardcode 83.1451 here, the value inherited from CO2SYS. 83.14472 is
the exact conversion of R = 8.314472 J/(mol K), and is what MATLAB uses. The difference
shifts pressure correction factors by ~2e-6 - well inside the 1e-4 cross-language tolerance.
"""
const GAS_CONSTANT = Float64(_load_fundamental_constants()["R_P"])

const K_COEFS = _load_coefficients("K_calculation.json")

# Every pressure correction takes exactly 5 coefficients, so unlike K_COEFS this NamedTuple
# is homogeneous and can be indexed by a runtime Symbol without losing type stability.
const K_PRESSCORR_COEFS = _load_coefficients("K_pressure_correction.json")

#####################################################################################
# Composition of seawater
#####################################################################################

"""
    calc_ionic_strength(sal)

Ionic strength of seawater at salinity `sal`, after Dickson (1990a); see Dickson et al. (2007).
"""
calc_ionic_strength(sal::Real) = 19.924 * sal / (1000 - 1.005 * sal)

"""
    calc_sulphate(sal=35.0)

Total sulphate in mol/kg-SW, from Dickson et al. (2007), Table 2.

Note that `sal / 1.80655` is chlorinity.
"""
calc_sulphate(sal::Real=35.0) = 0.14 * sal / 1.80655 / 96.062

"""
    calc_fluorine(sal=35.0)

Total fluorine in mol/kg-SW, from Dickson et al. (2007), Table 2.

Note that `sal / 1.80655` is chlorinity.
"""
calc_fluorine(sal::Real=35.0) = 6.7e-5 * sal / 1.80655 / 18.9984

#####################################################################################
# Individual K functions
#
# Each takes its coefficient tuple, temperature in Celsius and salinity, matching the
# signatures used by the Python, R and MATLAB implementations.
#####################################################################################

"""
    calc_K1K2(coefficients, temp_c, sal)

K1 or K2 — the first or second dissociation constant of carbonic acid.
"""
function calc_K1K2(coefficients::NTuple{5,Float64}, temp_c::Real, sal::Real)
    temp_k = temp_c + 273.15
    p = coefficients
    return 10^(
        p[1] +
        p[2] / temp_k +
        p[3] * log(temp_k) +
        p[4] * sal +
        p[5] * sal * sal
    )
end

"""
    calc_KW(coefficients, temp_c, sal)

KW — the ion product of water.
"""
function calc_KW(coefficients::NTuple{7,Float64}, temp_c::Real, sal::Real)
    temp_k = temp_c + 273.15
    log_temp_k = log(temp_k)
    p = coefficients
    return exp(
        p[1] +
        p[2] / temp_k +
        p[3] * log_temp_k +
        (p[4] / temp_k + p[5] + p[6] * log_temp_k) * sqrt(sal) +
        p[7] * sal
    )
end

"""
    calc_KB(coefficients, temp_c, sal)

KB — the dissociation constant of boric acid.
"""
function calc_KB(coefficients::NTuple{12,Float64}, temp_c::Real, sal::Real)
    temp_k = temp_c + 273.15
    sqrt_sal = sqrt(sal)
    p = coefficients
    return exp(
        (p[1] + p[2] * sqrt_sal + p[3] * sal) +
        (
            p[4] +
            p[5] * sqrt_sal +
            p[6] * sal +
            p[7] * sal * sqrt_sal +
            p[8] * sal * sal
        ) / temp_k +
        (p[9] + p[10] * sqrt_sal + p[11] * sal) * log(temp_k) +
        p[12] * sqrt_sal * temp_k
    )
end

"""
    calc_K0(coefficients, temp_c, sal)

K0 — the solubility of CO2 in seawater.
"""
function calc_K0(coefficients::NTuple{6,Float64}, temp_c::Real, sal::Real)
    temp_k = temp_c + 273.15
    p = coefficients
    return exp(
        p[1] +
        p[2] * 100 / temp_k +
        p[3] * log(temp_k / 100) +
        sal * (p[4] + p[5] * temp_k / 100 + p[6] * (temp_k / 100) * (temp_k / 100))
    )
end

"""
    calc_KS(coefficients, temp_c, sal)

KS — the dissociation constant of bisulphate.
"""
function calc_KS(coefficients::NTuple{11,Float64}, temp_c::Real, sal::Real)
    temp_k = temp_c + 273.15
    log_temp_k = log(temp_k)
    ionic_strength = calc_ionic_strength(sal)
    sqrt_istr = sqrt(ionic_strength)
    p = coefficients
    return exp(
        p[1] +
        p[2] / temp_k +
        p[3] * log_temp_k +
        sqrt_istr * (p[4] / temp_k + p[5] + p[6] * log_temp_k) +
        ionic_strength * (p[7] / temp_k + p[8] + p[9] * log_temp_k) +
        p[10] / temp_k * ionic_strength * sqrt_istr +
        p[11] / temp_k * ionic_strength * ionic_strength +
        log(1 - 0.001005 * sal)
    )
end

"""
    calc_Ksp(coefficients, temp_c, sal)

KspA or KspC — the solubility product of aragonite or calcite.
"""
function calc_Ksp(coefficients::NTuple{9,Float64}, temp_c::Real, sal::Real)
    temp_k = temp_c + 273.15
    sqrt_sal = sqrt(sal)
    p = coefficients
    return 10^(
        p[1] +
        p[2] * temp_k +
        p[3] / temp_k +
        p[4] * log10(temp_k) +
        (p[5] + p[6] * temp_k + p[7] / temp_k) * sqrt_sal +
        p[8] * sal +
        p[9] * sal * sqrt_sal
    )
end

"""
    calc_KP(coefficients, temp_c, sal)

KP1 or KP2 — the first or second dissociation constant of phosphoric acid.
"""
function calc_KP(coefficients::NTuple{7,Float64}, temp_c::Real, sal::Real)
    temp_k = temp_c + 273.15
    p = coefficients
    return exp(
        p[1] / temp_k +
        p[2] +
        p[3] * log(temp_k) +
        (p[4] / temp_k + p[5]) * sqrt(sal) +
        (p[6] / temp_k + p[7]) * sal
    )
end

"""
    calc_KP3(coefficients, temp_c, sal)

KP3 — the third dissociation constant of phosphoric acid.

Unlike KP1 and KP2 this has no `log(temp_k)` term, so it takes 6 coefficients rather
than 7 and cannot share [`calc_KP`](@ref).
"""
function calc_KP3(coefficients::NTuple{6,Float64}, temp_c::Real, sal::Real)
    temp_k = temp_c + 273.15
    p = coefficients
    return exp(
        p[1] / temp_k +
        p[2] +
        (p[3] / temp_k + p[4]) * sqrt(sal) +
        (p[5] / temp_k + p[6]) * sal
    )
end

"""
    calc_KSi(coefficients, temp_c, sal)

KSi — the dissociation constant of silicic acid.
"""
function calc_KSi(coefficients::NTuple{9,Float64}, temp_c::Real, sal::Real)
    temp_k = temp_c + 273.15
    ionic_strength = calc_ionic_strength(sal)
    p = coefficients
    return exp(
        p[1] / temp_k +
        p[2] +
        p[3] * log(temp_k) +
        (p[4] / temp_k + p[5]) * sqrt(ionic_strength) +
        (p[6] / temp_k + p[7]) * ionic_strength +
        (p[8] / temp_k + p[9]) * ionic_strength * ionic_strength
    ) * (1 - 0.001005 * sal)
end

"""
    calc_KF(coefficients, temp_c, sal)

KF — the dissociation constant of hydrogen fluoride.
"""
function calc_KF(coefficients::NTuple{3,Float64}, temp_c::Real, sal::Real)
    temp_k = temp_c + 273.15
    p = coefficients
    return exp(p[1] / temp_k + p[2] + p[3] * sqrt(sal))
end

#####################################################################################
# Pressure correction
#####################################################################################

"""
    calc_pressure_correction(coefficients, p_bar, temp_c)

Pressure correction factor for a thermodynamic K, from Millero et al. (2007,
doi:10.1021/cr0503557) eqns 38–40.

```
K_at_pressure = calc_pressure_correction(...) * K_at_surface
```
"""
function calc_pressure_correction(coefficients::NTuple{5,Float64}, p_bar::Real, temp_c::Real)
    p = coefficients
    delta_volume = p[1] + p[2] * temp_c + p[3] * temp_c^2
    # NB: CO2SYS carries a factor of 1000 here, which is folded into the coefficients.
    compressibility = p[4] + p[5] * temp_c
    RT = GAS_CONSTANT * (temp_c + 273.15)
    return exp((-delta_volume + 0.5 * compressibility * p_bar) * p_bar / RT)
end

"""
    calc_pH_scale_conversions(temp_c, sal, p_bar, sulphate, fluorine)

pH scale conversion factors that bracket the pressure correction, returned as
`(tot_to_sws_surface, sws_to_tot_deep)`.

The Millero pressure corrections are defined on the seawater scale, but Kgen reports Ks on
the total scale, so a K must be converted to SWS at surface pressure, corrected, then
converted back to total at depth.
"""
function calc_pH_scale_conversions(temp_c::Real, sal::Real, p_bar::Real,
                                   sulphate::Real, fluorine::Real)
    KS_surf = calc_KS(K_COEFS.KS, temp_c, sal)
    KS_deep = KS_surf * calc_pressure_correction(K_PRESSCORR_COEFS.KS, p_bar, temp_c)
    KF_surf = calc_KF(K_COEFS.KF, temp_c, sal)
    KF_deep = KF_surf * calc_pressure_correction(K_PRESSCORR_COEFS.KF, p_bar, temp_c)

    tot_to_sws_surface = (1 + sulphate / KS_surf + fluorine / KF_surf) / (1 + sulphate / KS_surf)
    sws_to_tot_deep = (1 + sulphate / KS_deep) / (1 + sulphate / KS_deep + fluorine / KF_deep)

    return tot_to_sws_surface, sws_to_tot_deep
end

#####################################################################################
# Public API
#####################################################################################

# Resolving coefficients by a runtime Symbol would be type-unstable, because K_COEFS holds
# tuples of differing length. Branching on the name instead keeps every path Float64.
function _calc_surface_K(K::Symbol, temp_c::Real, sal::Real)
    K === :K0   ? calc_K0(K_COEFS.K0, temp_c, sal)     :
    K === :K1   ? calc_K1K2(K_COEFS.K1, temp_c, sal)   :
    K === :K2   ? calc_K1K2(K_COEFS.K2, temp_c, sal)   :
    K === :KW   ? calc_KW(K_COEFS.KW, temp_c, sal)     :
    K === :KB   ? calc_KB(K_COEFS.KB, temp_c, sal)     :
    K === :KS   ? calc_KS(K_COEFS.KS, temp_c, sal)     :
    K === :KspA ? calc_Ksp(K_COEFS.KspA, temp_c, sal)  :
    K === :KspC ? calc_Ksp(K_COEFS.KspC, temp_c, sal)  :
    K === :KP1  ? calc_KP(K_COEFS.KP1, temp_c, sal)    :
    K === :KP2  ? calc_KP(K_COEFS.KP2, temp_c, sal)    :
    K === :KP3  ? calc_KP3(K_COEFS.KP3, temp_c, sal)   :
    K === :KSi  ? calc_KSi(K_COEFS.KSi, temp_c, sal)   :
    K === :KF   ? calc_KF(K_COEFS.KF, temp_c, sal)     :
    throw(ArgumentError("$K is not a valid K. Should be one of $(K_NAMES)."))
end

function _check_MyAMI_mode(MyAMI_mode::Symbol)
    if MyAMI_mode === :approximate
        return nothing
    elseif MyAMI_mode === :calculate
        throw(ArgumentError(
            "MyAMI_mode=:calculate runs the full MyAMI model, which is not implemented in " *
            "the Julia version of Kgen. Use MyAMI_mode=:approximate, or the Python " *
            "implementation if you need the full model."
        ))
    else
        throw(ArgumentError("Unknown MyAMI_mode $MyAMI_mode - must be :approximate or :calculate."))
    end
end

"""
    calc_K(K, temp_c, sal, p_bar=0.0, magnesium=0.0528171, calcium=0.0102821; kwargs...)
    calc_K(K; temp_c=25.0, sal=35.0, p_bar=0.0, magnesium=0.0528171, calcium=0.0102821, kwargs...)

Calculate a single stoichiometric equilibrium constant on the total pH scale.

The keyword form mirrors the Python, R and MATLAB implementations. The positional form
exists because Julia does not broadcast keyword arguments — use it to evaluate many
conditions at once (see below).

# Arguments
- `K`: name of the constant, e.g. `:K1`. One of `K_NAMES`. A `String` is also accepted.
- `temp_c`: temperature in Celsius.
- `sal`: salinity in PSU.
- `p_bar`: pressure in bar.
- `magnesium`: *average* seawater magnesium in mol/kgsw. Used to correct the Ks via MyAMI.
- `calcium`: *average* seawater calcium in mol/kgsw. Used to correct the Ks via MyAMI.

# Keyword arguments
- `sulphate`: total sulphate in mol/kgsw. Calculated from `sal` if not given.
- `fluorine`: total fluorine in mol/kgsw. Calculated from `sal` if not given.
- `MyAMI_mode`: only `:approximate` is available in Julia; see [`PyMYAMI`](@ref).

# Examples
```julia
julia> calc_K(:K1, temp_c=25.0, sal=35.0)
1.4212669153166358e-6
```

Broadcast over the positional form to evaluate many conditions at once:

```julia
julia> calc_K.(:K1, [0.0, 25.0], 35.0)
2-element Vector{Float64}:
 8.379e-7
 1.4213e-6
```
"""
function calc_K(K::Symbol,
                temp_c::Real,
                sal::Real,
                p_bar::Real=0.0,
                magnesium::Real=MODERN_MAGNESIUM,
                calcium::Real=MODERN_CALCIUM;
                sulphate::Union{Nothing,Real}=nothing,
                fluorine::Union{Nothing,Real}=nothing,
                MyAMI_mode::Symbol=:approximate)
    _check_MyAMI_mode(MyAMI_mode)

    result = _calc_surface_K(K, temp_c, sal)

    if p_bar != 0
        total_sulphate = isnothing(sulphate) ? calc_sulphate(sal) : float(sulphate)
        total_fluorine = isnothing(fluorine) ? calc_fluorine(sal) : float(fluorine)
        tot_to_sws, sws_to_tot = calc_pH_scale_conversions(
            temp_c, sal, p_bar, total_sulphate, total_fluorine
        )
        result *= tot_to_sws *
                  calc_pressure_correction(K_PRESSCORR_COEFS[K], p_bar, temp_c) *
                  sws_to_tot
    end

    # Applied unconditionally, matching R and MATLAB. The factor is ~1 at modern composition
    # but not exactly 1 (up to 9e-4 away over the usual T/S range), so skipping it when
    # magnesium and calcium are modern would put Julia outside the cross-language tolerance.
    result *= PyMYAMI.approximate_seawater_correction(
        K; temp_c=temp_c, sal=sal, magnesium=magnesium, calcium=calcium
    )

    return result
end

function calc_K(K::Symbol;
                temp_c::Real=25.0,
                sal::Real=35.0,
                p_bar::Real=0.0,
                magnesium::Real=MODERN_MAGNESIUM,
                calcium::Real=MODERN_CALCIUM,
                kwargs...)
    calc_K(K, temp_c, sal, p_bar, magnesium, calcium; kwargs...)
end

calc_K(K::AbstractString, args::Real...; kwargs...) = calc_K(Symbol(K), args...; kwargs...)
calc_K(K::AbstractString; kwargs...) = calc_K(Symbol(K); kwargs...)

"""
    calc_Ks(temp_c, sal, p_bar=0.0, magnesium=0.0528171, calcium=0.0102821; kwargs...)
    calc_Ks(; temp_c=25.0, sal=35.0, p_bar=0.0, magnesium=0.0528171, calcium=0.0102821, kwargs...)

Calculate all 13 stoichiometric equilibrium constants on the total pH scale, returned as a
`NamedTuple` keyed by `K_NAMES`.

Arguments are as for [`calc_K`](@ref).

# Examples
```julia
julia> ks = calc_Ks(temp_c=25.0, sal=35.0);

julia> ks.K1
1.4212669153166358e-6
```

Broadcasting the positional form gives one `NamedTuple` per condition, which can be passed
straight to a `DataFrame`:

```julia
calc_Ks.([0.0, 25.0], 35.0)
```
"""
function calc_Ks(temp_c::Real,
                 sal::Real,
                 p_bar::Real=0.0,
                 magnesium::Real=MODERN_MAGNESIUM,
                 calcium::Real=MODERN_CALCIUM;
                 sulphate::Union{Nothing,Real}=nothing,
                 fluorine::Union{Nothing,Real}=nothing,
                 MyAMI_mode::Symbol=:approximate)
    _check_MyAMI_mode(MyAMI_mode)

    ks = (
        K0   = calc_K0(K_COEFS.K0, temp_c, sal),
        K1   = calc_K1K2(K_COEFS.K1, temp_c, sal),
        K2   = calc_K1K2(K_COEFS.K2, temp_c, sal),
        KW   = calc_KW(K_COEFS.KW, temp_c, sal),
        KB   = calc_KB(K_COEFS.KB, temp_c, sal),
        KS   = calc_KS(K_COEFS.KS, temp_c, sal),
        KspA = calc_Ksp(K_COEFS.KspA, temp_c, sal),
        KspC = calc_Ksp(K_COEFS.KspC, temp_c, sal),
        KP1  = calc_KP(K_COEFS.KP1, temp_c, sal),
        KP2  = calc_KP(K_COEFS.KP2, temp_c, sal),
        KP3  = calc_KP3(K_COEFS.KP3, temp_c, sal),
        KSi  = calc_KSi(K_COEFS.KSi, temp_c, sal),
        KF   = calc_KF(K_COEFS.KF, temp_c, sal),
    )

    if p_bar != 0
        total_sulphate = isnothing(sulphate) ? calc_sulphate(sal) : float(sulphate)
        total_fluorine = isnothing(fluorine) ? calc_fluorine(sal) : float(fluorine)
        # The scale conversions and the KS/KF pressure corrections they depend on are the
        # same for every K, so they are computed once here rather than inside the loop.
        tot_to_sws, sws_to_tot = calc_pH_scale_conversions(
            temp_c, sal, p_bar, total_sulphate, total_fluorine
        )
        scale = tot_to_sws * sws_to_tot
        # Spelled out rather than `map`ped over K_PRESSCORR_COEFS: mapping a closure over a
        # NamedTuple inlines unreliably here, and heap-allocates when it does not.
        ks = (
            K0   = ks.K0   * scale * calc_pressure_correction(K_PRESSCORR_COEFS.K0, p_bar, temp_c),
            K1   = ks.K1   * scale * calc_pressure_correction(K_PRESSCORR_COEFS.K1, p_bar, temp_c),
            K2   = ks.K2   * scale * calc_pressure_correction(K_PRESSCORR_COEFS.K2, p_bar, temp_c),
            KW   = ks.KW   * scale * calc_pressure_correction(K_PRESSCORR_COEFS.KW, p_bar, temp_c),
            KB   = ks.KB   * scale * calc_pressure_correction(K_PRESSCORR_COEFS.KB, p_bar, temp_c),
            KS   = ks.KS   * scale * calc_pressure_correction(K_PRESSCORR_COEFS.KS, p_bar, temp_c),
            KspA = ks.KspA * scale * calc_pressure_correction(K_PRESSCORR_COEFS.KspA, p_bar, temp_c),
            KspC = ks.KspC * scale * calc_pressure_correction(K_PRESSCORR_COEFS.KspC, p_bar, temp_c),
            KP1  = ks.KP1  * scale * calc_pressure_correction(K_PRESSCORR_COEFS.KP1, p_bar, temp_c),
            KP2  = ks.KP2  * scale * calc_pressure_correction(K_PRESSCORR_COEFS.KP2, p_bar, temp_c),
            KP3  = ks.KP3  * scale * calc_pressure_correction(K_PRESSCORR_COEFS.KP3, p_bar, temp_c),
            KSi  = ks.KSi  * scale * calc_pressure_correction(K_PRESSCORR_COEFS.KSi, p_bar, temp_c),
            KF   = ks.KF   * scale * calc_pressure_correction(K_PRESSCORR_COEFS.KF, p_bar, temp_c),
        )
    end

    # See the note in calc_K: applied unconditionally to stay consistent with R and MATLAB.
    corrections = PyMYAMI.approximate_seawater_corrections(
        temp_c=temp_c, sal=sal, magnesium=magnesium, calcium=calcium
    )
    # Only some Ks are corrected; merge puts the corrected subset back without disturbing
    # the field order of `ks`.
    ks = merge(ks, map(*, ks[PyMYAMI.POLY_NAMES], corrections))

    return ks
end

function calc_Ks(; temp_c::Real=25.0,
                   sal::Real=35.0,
                   p_bar::Real=0.0,
                   magnesium::Real=MODERN_MAGNESIUM,
                   calcium::Real=MODERN_CALCIUM,
                   kwargs...)
    calc_Ks(temp_c, sal, p_bar, magnesium, calcium; kwargs...)
end

end # module
