module Kgen

include("PyMYAMI.jl")

using JSON

project_path(parts...) = normpath(joinpath(@__DIR__, "..", parts...))

K_coefs = JSON.parsefile(project_path("..", "coefficients", "K_calculation.json"))["coefficients"]
K_presscorr_coefs = JSON.parsefile(project_path("..", "coefficients", "K_pressure_correction.json"))["coefficients"]

function calc_Istr(Sal::AbstractFloat)
    return (
        19.924 * Sal / (1000 - 1.005 * Sal)
    )
end

function fn_K1K2(;
    p::Vector,
    TK::AbstractFloat,
    lnTK::AbstractFloat,
    S::AbstractFloat,
    sqrtS::AbstractFloat)
    """
    Return K1 or K2 at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::AbstractFloat`: Temperature in Kelvin
    - `lnTK::AbstractFloat`: natural log of temperature in kelvin
    - `S::AbstractFloat`: Salinity
    - `sqrtS::AbstractFloat`: square root of salinity
    ...
    """
    10^(
        p[1] +
        p[2] / TK +
        p[3] * lnTK +
        p[4] * S +
        p[5] * S * S
    ) 
end


function fn_KW(;
    p::Vector, 
    TK::AbstractFloat,
    lnTK::AbstractFloat,
    S::AbstractFloat,
    sqrtS::AbstractFloat)
    """
    Return K1 or K2 at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::AbstractFloat`: Temperature in Kelvin
    - `lnTK::AbstractFloat`: natural log of temperature in kelvin
    - `S::AbstractFloat`: Salinity
    - `sqrtS::AbstractFloat`: square root of salinity
    ...
    """
    exp(
        p[1] +
        p[2] / TK +
        p[3] * lnTK +
        (p[4] / TK + p[5] + p[6] * lnTK) * sqrtS +
        p[7] * S
    )
end

function fn_KB(;
    p::Vector, 
    TK::AbstractFloat,
    lnTK::AbstractFloat,
    S::AbstractFloat,
    sqrtS::AbstractFloat)
    """
    Return K1 or K2 at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::AbstractFloat`: Temperature in Kelvin
    - `lnTK::AbstractFloat`: natural log of temperature in kelvin
    - `S::AbstractFloat`: Salinity
    - `sqrtS::AbstractFloat`: square root of salinity
    ...
    """
    exp(
        (p[1] + p[2] * sqrtS + p[3] * S) + 
        (
            p[4] +
            p[5] * sqrtS +
            p[6] * S +
            p[7] * S * sqrtS +
            p[8] * S * S
        ) / TK +
        (p[9] + p[10] * sqrtS + p[11] * S) * lnTK +
        p[12] * sqrtS * TK
    )
end

function fn_K0(;
    p::Vector, 
    TK::AbstractFloat,
    lnTK::AbstractFloat,
    S::AbstractFloat,
    sqrtS::AbstractFloat)
    """
    Return K1 or K2 at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::AbstractFloat`: Temperature in Kelvin
    - `lnTK::AbstractFloat`: natural log of temperature in kelvin
    - `S::AbstractFloat`: Salinity
    - `sqrtS::AbstractFloat`: square root of salinity
    ...
    """
    exp(
        p[1] +
        p[2] * 100 / TK +
        p[3] * log(TK / 100) +
        S * (p[4] + p[5] * TK / 100 + p[6] * (TK / 100) * (TK / 100))
    )
end

function fn_KS(;
    p::Vector, 
    TK::AbstractFloat,
    lnTK::AbstractFloat,
    S::AbstractFloat,
    sqrtS::AbstractFloat)
    """
    Return KS at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::AbstractFloat`: Temperature in Kelvin
    - `lnTK::AbstractFloat`: natural log of temperature in kelvin
    - `S::AbstractFloat`: Salinity
    - `sqrtS::AbstractFloat`: square root of salinity
    ...
    """
    
    Istr = calc_Istr(S)
    # Ionic strength after Dickson 1990a; see Dickson et al 2007

    exp(
        p[1] +
        p[2] / TK +
        p[3] * lnTK +
        Istr^0.5 * (p[4] / TK + p[5] + p[6] * lnTK) +
        Istr * (p[7] / TK + p[8] + p[9] * lnTK) +
        p[10] / TK * Istr * Istr^0.5 +
        p[11] / TK * Istr^2 +
        log(1 - 0.001005 * S)
    )
end

function fn_Ksp(;
    p::Vector, 
    TK::AbstractFloat,
    lnTK::AbstractFloat,
    S::AbstractFloat,
    sqrtS::AbstractFloat)
    """
    Return KspA or KspC at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::AbstractFloat`: Temperature in Kelvin
    - `lnTK::AbstractFloat`: natural log of temperature in kelvin
    - `S::AbstractFloat`: Salinity
    - `sqrtS::AbstractFloat`: square root of salinity
    ...
    """

    10^(
        p[1] + 
        p[2] * TK +
        p[3] / TK +
        p[4] * log10(TK) +
        (p[5] + p[6] * TK + p[7] / TK) * sqrtS +
        p[8] * S +
        p[9] * S * sqrtS
    )
end

function fn_KP(;
    p::Vector, 
    TK::AbstractFloat,
    lnTK::AbstractFloat,
    S::AbstractFloat,
    sqrtS::AbstractFloat)
    """
    Return KP1 or KP2 at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::AbstractFloat`: Temperature in Kelvin
    - `lnTK::AbstractFloat`: natural log of temperature in kelvin
    - `S::AbstractFloat`: Salinity
    - `sqrtS::AbstractFloat`: square root of salinity
    ...
    """

    exp(
        p[1] / TK +
        p[2] +
        p[3] * lnTK +
        (p[4] / TK + p[5]) * sqrtS +
        (p[6] / TK + p[7]) * S
    )
end

function fn_KP3(;
    p::Vector, 
    TK::AbstractFloat,
    lnTK::AbstractFloat,
    S::AbstractFloat,
    sqrtS::AbstractFloat)
    """
    Return KP3 at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::AbstractFloat`: Temperature in Kelvin
    - `lnTK::AbstractFloat`: natural log of temperature in kelvin
    - `S::AbstractFloat`: Salinity
    - `sqrtS::AbstractFloat`: square root of salinity
    ...
    """

    exp(
        p[1] / TK +
        p[2] +
        (p[3] / TK + p[4]) * sqrtS +
        (p[6] / TK + p[6]) * S
    )
end

function fn_KSi(;
    p::Vector, 
    TK::AbstractFloat,
    lnTK::AbstractFloat,
    S::AbstractFloat,
    sqrtS::AbstractFloat)
    """
    Return KSi at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::AbstractFloat`: Temperature in Kelvin
    - `lnTK::AbstractFloat`: natural log of temperature in kelvin
    - `S::AbstractFloat`: Salinity
    - `sqrtS::AbstractFloat`: square root of salinity
    ...
    """

    Istr = calc_Istr(S)

    exp(
        p[1] / TK + 
        p[2] +
        p[3] * log(TK) +
        (p[4] / TK + p[5]) * Istr^0.5 +
        (p[6] / TK + p[7]) * Istr +
        (p[8] / TK + p[9]) * Istr^2
    ) * (1 - 0.001005 * S)
end

function fn_KF(;
    p::Vector, 
    TK::AbstractFloat,
    lnTK::AbstractFloat,
    S::AbstractFloat,
    sqrtS::AbstractFloat)
    """
    Return KF at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::AbstractFloat`: Temperature in Kelvin
    - `lnTK::AbstractFloat`: natural log of temperature in kelvin
    - `S::AbstractFloat`: Salinity
    - `sqrtS::AbstractFloat`: square root of salinity
    ...
    """

    exp(
        p[1] / TK + 
        p[2] + 
        p[3] * sqrtS
    )
end

K_fns = Dict(
    "K0" => fn_K0,
    "K1" => fn_K1K2,
    "K2" => fn_K1K2,
    "KW" => fn_KW,
    "KB" => fn_KB,
    "KS" => fn_KS,
    "KspA" => fn_Ksp,
    "KspC" => fn_Ksp,
    "KP1" => fn_KP,
    "KP2" => fn_KP,
    "KP3" => fn_KP3,
    "KSi" => fn_KSi,
    "KF" => fn_KF
)

function calc_TS(Sal::AbstractFloat=35.0)
    """
    Calculate total Sulphur in mol/kg-SW- lifted directly from CO2SYS.m

    From Dickson et al., 2007, Table 2
    Note: Sal / 1.80655 = Chlorinity
    """
    return 0.14 * Sal / 1.80655 / 96.062 # mol/kg-SW
end

function calc_TF(Sal::AbstractFloat=35.0)
    """
    Calculate total Fluorine in mol/kg-SW

    From Dickson et al., 2007, Table 2
    Note: Sal / 1.80655 = Chlorinity
    """
    return 6.7e-5 * Sal / 1.80655 / 18.9984 # mol/kg-SW
end

function prescorr(;
    p::Vector, 
    P::AbstractFloat,
    TC::AbstractFloat)
    """
    Return pressure correction factor for thermodynamic Ks.

    From Millero et al (2007, doi:10.1021/cr0503557)
    Eqns 38-40

    Usage:
    K_at_pressure / K_orig = [output]
    K_at_pressure = [output] * K_orig

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `P::AbstractFloat`: Pressure in bar
    - `TC::AbstractFloat`: Temperature in Celsius
    ...
    """
    dV = p[1] + p[2] * TC + p[3] * TC^2
    dk = (p[4] + p[5] * TC)  # NB: there is a factor of 1000 in CO2sys, which has been incorporated into the coefficients for the function.    
    RT = 83.1451 * (TC + 273.15)
    exp((-dV + 0.5 * dk * P) * P / RT)
end

function calc_K(k::String;
    TempC::AbstractFloat=25.0,
    Sal::AbstractFloat=35.0,
    Pres::AbstractFloat=0.0,
    Mg::AbstractFloat=0.0528171,
    Ca::AbstractFloat=0.0102821,
    MyAMI_mode::String="approximate"
    )

    TK = TempC + 273.15
    lnTK = log(TK)
    sqrtS = Sal^0.5

    K = K_fns[k](p=K_coefs[k], TK=TK, lnTK=lnTK, S=Sal, sqrtS=sqrtS)

    if Pres > 0
        if haskey(K_presscorr_coefs, k)
            TF = calc_TF(Sal)
            TS = calc_TS(Sal)

            KS_surf = K_fns["KS"](p=K_coefs["KS"], TK=TK, lnTK=lnTK, S=Sal, sqrtS=sqrtS)
            KS_deep = KS_surf * prescorr(p=K_presscorr_coefs["KS"], P=Pres, TC=TempC)
            KF_surf = K_fns["KF"](p=K_coefs["KF"], TK=TK, lnTK=lnTK, S=Sal, sqrtS=sqrtS)
            KF_deep = KF_surf * prescorr(p=K_presscorr_coefs["KF"], P=Pres, TC=TempC)

            tot_to_sws_surface = (1 + TS / KS_surf) / (1 + TS / KS_surf + TF / KF_surf)  # convert from TOT to SWS before pressure correction
            sws_to_tot_deep = (1 + TS / KS_deep + TF / KF_deep) / (1 + TS / KS_deep)  # convert from SWS to TOT after pressure correction
            
            K *= tot_to_sws_surface * prescorr(p=K_presscorr_coefs[k], P=Pres, TC=TempC) * sws_to_tot_deep
        end
    end

    if MyAMI_mode == "approximate"
        if any(Mg != 0.0528171) || any(Ca != 0.0102821)
            Fcorr = PyMYAMI.approximate_Fcorr(K=k, TempC=TempC, Sal=Sal, Mg=Mg, Ca=Ca)
            K *= Fcorr
        end
    end

    return K

end

mutable struct Ks
    K0::AbstractFloat
    K1::AbstractFloat
    K2::AbstractFloat
    KW::AbstractFloat
    KB::AbstractFloat
    KS::AbstractFloat
    KspA::AbstractFloat
    KspC::AbstractFloat
    KP1::AbstractFloat
    KP2::AbstractFloat
    KP3::AbstractFloat
    KSi::AbstractFloat
    KF::AbstractFloat

    # default constructor
    Ks() = new(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end

function calc_Ks(;
    TempC::AbstractFloat=25.0,
    Sal::AbstractFloat=35.0,
    Pres::AbstractFloat=0.0,
    Mg::AbstractFloat=0.0528171,
    Ca::AbstractFloat=0.0102821,
    MyAMI_mode::String="approximate"
    )

    if MyAMI_mode == "approximate"
        if (Mg != 0.0528171) || (Ca != 0.0102821)
            Fcorr = PyMYAMI.approximate_Fcorrs(TempC=TempC, Sal=Sal, Mg=Mg, Ca=Ca)
            composition_correct = true
        else
            composition_correct = false
        end
    end

    TK = TempC + 273.15
    lnTK = log(TK)
    sqrtS = Sal^0.5

    result = Ks()
    for (k, fn) in K_fns
        K_tmp = fn(p=K_coefs[k], TK=TK, lnTK=lnTK, S=Sal, sqrtS=sqrtS)

        if Pres > 0
            if haskey(K_presscorr_coefs, k)
                TF = calc_TF(Sal)
                TS = calc_TS(Sal)
    
                KS_surf = K_fns["KS"](p=K_coefs["KS"], TK=TK, lnTK=lnTK, S=Sal, sqrtS=sqrtS)
                KS_deep = KS_surf * prescorr(p=K_presscorr_coefs["KS"], P=Pres, TC=TempC)
                KF_surf = K_fns["KF"](p=K_coefs["KF"], TK=TK, lnTK=lnTK, S=Sal, sqrtS=sqrtS)
                KF_deep = KF_surf * prescorr(p=K_presscorr_coefs["KF"], P=Pres, TC=TempC)
    
                tot_to_sws_surface = (1 + TS / KS_surf) / (1 + TS / KS_surf + TF / KF_surf)  # convert from TOT to SWS before pressure correction
                sws_to_tot_deep = (1 + TS / KS_deep + TF / KF_deep) / (1 + TS / KS_deep)  # convert from SWS to TOT after pressure correction
                
                K_tmp *= tot_to_sws_surface * prescorr(p=K_presscorr_coefs[k], P=Pres, TC=TempC) * sws_to_tot_deep
            end
        end
        if composition_correct
            if haskey(Fcorr, k)
                K_tmp *= Fcorr[k]
            end
        end

        setproperty!(result, Symbol(k), K_tmp)
    end
    
    return result
end

end
