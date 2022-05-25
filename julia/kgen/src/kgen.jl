module kgen

function fn_K1K2(p::Vector, 
    TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}})
    """
    Return K1 or K2 at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Temperature in Kelvin
    - `lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: natural log of temperature in kelvin
    - `S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Salinity
    - `sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: square root of salinity
    ...
    """
    10 .^ (
        p[1] .+
        p[2] ./ TK .+
        p[3] .* lnTK .+
        p[4] .* S .+
        p[5] .* S .* S
    ) 
end


function fn_KW(p::Vector, 
    TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}})
    """
    Return K1 or K2 at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Temperature in Kelvin
    - `lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: natural log of temperature in kelvin
    - `S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Salinity
    - `sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: square root of salinity
    ...
    """
    exp.(
        p[1] .+
        p[2] ./ TK .+
        p[3] .* lnTK .+
        (p[4] ./ TK .+ p[5] .+ p[6] .* lnTK) .* sqrtS .+
        p[7] .* S
    )
end

function fn_KB(p::Vector, 
    TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}})
    """
    Return K1 or K2 at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Temperature in Kelvin
    - `lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: natural log of temperature in kelvin
    - `S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Salinity
    - `sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: square root of salinity
    ...
    """
    exp.(
        (p[1] .+ p[2] .* sqrtS .+ p[3] .* S) .+ 
        (
            p[4] .+
            p[5] .* sqrtS .+
            p[6] .* S .+
            p[7] .* S .* sqrtS .+
            p[8] .* S .* S
        ) ./ TK .+
        (p[9] .+ p[10] .* sqrtS .+ p[11] .* S) .* lnTK .+
        p[12] .* sqrtS .* TK
    )
end

function fn_K0(p::Vector, 
    TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}})
    """
    Return K1 or K2 at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Temperature in Kelvin
    - `lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: natural log of temperature in kelvin
    - `S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Salinity
    - `sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: square root of salinity
    ...
    """
    exp.(
        p[1] .+
        p[2] .* 100 ./ TK .+
        p[3] .* log.(TK ./ 100) .+
        S .* (p[4] .+ p[5] .* TK ./ 100 .+ p[6] .* (TK ./ 100) .* (TK ./ 100))
    )
end

function fn_KS(p::Vector, 
    TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}})
    """
    Return KS at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Temperature in Kelvin
    - `lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: natural log of temperature in kelvin
    - `S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Salinity
    - `sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: square root of salinity
    ...
    """
    
    Istr = (
        19.924 * S / (1000 - 1.005 * S)
    )
    # Ionic strength after Dickson 1990a; see Dickson et al 2007

    exp.(
        p[1] .+
        p[2] / TK .+
        p[3] * lnTK .+
        Istr .^ 0.5 .* (p[4] ./ TK .+ p[5] .+ p[6] .* lnTK) .+
        Istr .* (p[7] ./ TK .+ p[8] .+ p[9] .* lnTK) .+
        p[10] ./ TK .* Istr .* Istr .^ 0.5 .+
        p[11] ./ TK .* Istr .^ 2 .+
        log.(1 .- 0.001005 .* S)
    )
end

function fn_Ksp(p::Vector, 
    TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}})
    """
    Return KspA or KspC at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Temperature in Kelvin
    - `lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: natural log of temperature in kelvin
    - `S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Salinity
    - `sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: square root of salinity
    ...
    """

    10 .^ (
        p[1] .+ 
        p[2] .* TK .+
        p[3] ./ TK .+
        p[4] .* log10.(TK) .+
        (p[5] .+ p[6] .* TK .+ p[7] ./ TK) .* sqrtS .+
        p[8] .* S .+
        p[9] .* S .* sqrtS
    )
end

function fn_KP(p::Vector, 
    TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}})
    """
    Return KP1 or KP2 at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Temperature in Kelvin
    - `lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: natural log of temperature in kelvin
    - `S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Salinity
    - `sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: square root of salinity
    ...
    """

    exp.(
        p[1] ./ TK .+
        p[2] .+
        p[3] .* lnTK .+
        (p[4] ./ TK .+ p[5]) .* sqrtS .+
        (p[6] ./ TK .+ p[7]) .* S
    )
end

function fn_KP3(p::Vector, 
    TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}})
    """
    Return KP3 at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Temperature in Kelvin
    - `lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: natural log of temperature in kelvin
    - `S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Salinity
    - `sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: square root of salinity
    ...
    """

    exp.(
        p[1] ./ TK .+
        p[2] .+
        (p[3] ./ TK .+ p[4]) .* sqrtS .+
        (p[6] ./ TK .+ p[6]) .* S
    )
end

function fn_KSi(p::Vector, 
    TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}})
    """
    Return KSi at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Temperature in Kelvin
    - `lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: natural log of temperature in kelvin
    - `S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Salinity
    - `sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: square root of salinity
    ...
    """

    Istr = 19.924 .* S ./ (1000 .- 1.005 .* S)

    exp.(
        p[1] ./ TK .+ 
        p[2] .+
        p[3] .* log.(TK) .+
        (p[4] ./ TK .+ p[5]) .* Istr .^ 0.5 .+
        (p[6] ./ TK .+ p[7]) .* Istr .+
        (p[8] ./ TK .+ p[9]) .* Istr .^ 2
    ) .* (1 .- 0.001005 .* S)
end

function fn_KF(p::Vector, 
    TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}},
    sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}})
    """
    Return KF at given T and S.

    ...
    # Arguments
    ----------
    - `p::Vector`: parameters to be used for for K calculation
    - `TK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Temperature in Kelvin
    - `lnTK::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: natural log of temperature in kelvin
    - `S::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: Salinity
    - `sqrtS::Union{AbstractFloat, AbstractArray{<:AbstractFloat}}`: square root of salinity
    ...
    """

    exp.(
        p[1] ./ TK .+ 
        p[2] .+ 
        p[3] .* sqrtS
    )
end

K_fns = {
    "K0": fn_K0,
    "K1": fn_K1K2,
    "K2": fn_K1K2,
    "KW": fn_KW,
    "KB": fn_KB,
    "KS": fn_KS,
    "KspA": fn_Ksp,
    "KspC": fn_Ksp,
    "KP1": fn_KP,
    "KP2": fn_KP,
    "KP3": fn_KP3,
    "KSi": fn_KSi,
    "KF": fn_KF
}

end

