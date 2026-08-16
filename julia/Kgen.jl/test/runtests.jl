using Kgen
using Test
using JSON

# The check values are shared with the Python, R and MATLAB test suites, so they are read
# from the repository root rather than duplicated here. This means the tests must be run
# from a checkout of the Kgen repository, not from an installed copy of the package.
const CHECK_VALUE_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "check_values"))

isdir(CHECK_VALUE_DIR) || error(
    "Cannot find $CHECK_VALUE_DIR. The Kgen.jl tests must be run from a checkout of the " *
    "Kgen repository, because they share check values with the other language implementations."
)

read_check_values(name) = JSON.parsefile(joinpath(CHECK_VALUE_DIR, name))

"""
Number of decimal places a check value is quoted to, so that each K is compared at the
precision it was actually published to. Mirrors the `sigfig` logic in `python/test.py`.
"""
function quoted_decimal_places(value)
    decimals = split(rstrip(string(float(value)), '0'), '.')[2]
    return length(decimals)
end

@testset "Kgen" begin

    @testset "K values against published check values" begin
        check = read_check_values("check_Ks.json")
        sal = check["input_conditions"]["S"]
        temp_c = check["input_conditions"]["TC"]

        ks = calc_Ks(temp_c=temp_c, sal=sal)

        for (name, expected) in check["check_values"]
            calculated = log(ks[Symbol(name)])
            tolerance = 0.5 * 10.0^-quoted_decimal_places(expected)
            @test isapprox(calculated, expected; atol=tolerance)
        end
    end

    @testset "pressure correction factors" begin
        check = read_check_values("check_presscorr.json")
        temp_c = check["input_conditions"]["TC"]
        p_bar = check["input_conditions"]["P"]

        for (name, expected) in check["check_values"]
            factor = Kgen.calc_pressure_correction(
                Kgen.K_PRESSCORR_COEFS[Symbol(name)], p_bar, temp_c
            )
            @test isapprox(factor, expected; atol=1e-5)
        end
    end

    @testset "pH scale conversions bracket the pressure correction" begin
        # Guards against silently inverting the conversions, which shifts every
        # pressure-corrected K by ~0.37% - far more than the 0.01% cross-language tolerance.
        tot_to_sws, sws_to_tot = Kgen.calc_pH_scale_conversions(
            25.0, 35.0, 300.0, Kgen.calc_sulphate(35.0), Kgen.calc_fluorine(35.0)
        )
        @test tot_to_sws > 1
        @test sws_to_tot < 1

        # Value computed from the reference formulation at temp_c=25, sal=35, p_bar=300, using
        # R_P from fundamental_constants.json and before the MyAMI correction. Divide the
        # correction back out so this pins the pressure and pH-scale path alone.
        # (The equivalent figure with Python and R's hardcoded 83.1451 is 1.8632444248743345e-6.)
        seawater_correction = Kgen.PyMYAMI.approximate_seawater_correction(
            :K1, temp_c=25.0, sal=35.0
        )
        @test calc_K(:K1, temp_c=25.0, sal=35.0, p_bar=300.0) / seawater_correction ≈
              1.8632467272350019e-6
    end

    @testset "gas constant comes from fundamental_constants.json" begin
        constants = JSON.parsefile(
            joinpath(@__DIR__, "..", "src", "coefficients", "fundamental_constants.json")
        )["coefficients"]
        @test Kgen.GAS_CONSTANT == constants["R_P"]
        # Guards against silently reverting to CO2SYS's 83.1451, which Python and R still use.
        @test Kgen.GAS_CONSTANT != 83.1451
    end

    @testset "seawater composition correction" begin
        modern = calc_Ks(temp_c=25.0, sal=35.0)
        # Passing modern Mg/Ca explicitly must be identical to omitting them.
        @test calc_Ks(temp_c=25.0, sal=35.0, magnesium=0.0528171, calcium=0.0102821) == modern

        # The correction is applied unconditionally, matching R and MATLAB, so even at modern
        # composition it shifts the corrected Ks slightly away from the uncorrected value.
        # Skipping it instead would push Julia outside the cross-language tolerance.
        @test modern.K1 != Kgen._calc_surface_K(:K1, 25.0, 35.0)
        @test modern.KP1 == Kgen._calc_surface_K(:KP1, 25.0, 35.0)  # no polynomial for KP1

        altered = calc_Ks(temp_c=25.0, sal=35.0, magnesium=0.03, calcium=0.02)
        @test altered.K1 != modern.K1
        # KP1 has no MyAMI polynomial, so it must be left untouched.
        @test altered.KP1 == modern.KP1

        # Correction factors are ~1 at modern composition.
        for factor in Kgen.PyMYAMI.approximate_seawater_corrections()
            @test isapprox(factor, 1.0; atol=2e-5)
        end
    end

    @testset "argument handling" begin
        # Integers must work - python/test.py exercises the equivalent case.
        @test calc_K(:K1, temp_c=25, sal=35) isa Float64
        @test calc_Ks(temp_c=30, sal=36, p_bar=2).K1 isa Float64

        # Strings are accepted for parity with the other implementations.
        @test calc_K("K1", temp_c=25.0, sal=35.0) == calc_K(:K1, temp_c=25.0, sal=35.0)

        # Positional and keyword forms must agree.
        @test calc_K(:K1, 25.0, 35.0, 300.0) == calc_K(:K1, temp_c=25.0, sal=35.0, p_bar=300.0)
        @test calc_Ks(25.0, 35.0, 300.0) == calc_Ks(temp_c=25.0, sal=35.0, p_bar=300.0)

        # Explicit sulphate/fluorine override the values derived from salinity.
        @test calc_K(:K1, temp_c=25.0, sal=35.0, p_bar=300.0, sulphate=0.03) !=
              calc_K(:K1, temp_c=25.0, sal=35.0, p_bar=300.0)

        @test_throws ArgumentError calc_K(:K9)
        @test_throws ArgumentError calc_K(:K1, MyAMI_mode=:nonsense)
        # The full MyAMI model is not implemented in Julia; it must fail loudly.
        @test_throws ArgumentError calc_K(:K1, MyAMI_mode=:calculate)
        @test_throws ArgumentError calc_Ks(MyAMI_mode=:calculate)
    end

    @testset "type stability and allocations" begin
        @test Base.return_types(calc_K, (Symbol, Float64, Float64))[1] === Float64
        @test all(t -> t === Float64, Base.return_types(calc_Ks, (Float64, Float64))[1].types)

        # The scalar path is the fast path for broadcasting, so it must not allocate.
        calc_Ks(25.0, 35.0, 300.0, 0.03, 0.02)
        @test @allocated(calc_Ks(25.0, 35.0, 300.0, 0.03, 0.02)) == 0
    end

    @testset "broadcasting" begin
        temps = [0.0, 12.5, 25.0]
        broadcasted = calc_K.(:K1, temps, 35.0)
        @test broadcasted == [calc_K(:K1, t, 35.0) for t in temps]

        all_ks = calc_Ks.(temps, 35.0)
        @test length(all_ks) == 3
        @test all_ks[end].K1 == calc_K(:K1, 25.0, 35.0)
    end

    @testset "K_NAMES covers every coefficient set" begin
        coefficients = JSON.parsefile(
            joinpath(@__DIR__, "..", "src", "coefficients", "K_calculation.json")
        )["coefficients"]
        @test Set(Symbol.(keys(coefficients))) == Set(Kgen.K_NAMES)
        @test length(calc_Ks()) == length(Kgen.K_NAMES)
    end

end
