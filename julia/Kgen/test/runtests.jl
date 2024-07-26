using Kgen
using Test
using CSV

project_path(parts...) = normpath(joinpath(@__DIR__, "..", parts...))

check_Ks = JSON.parsefile(project_path("..", "check_values", "check_Ks.json"))

check_presscorr = JSON.parsefile(project_path("..", "check_values", "check_presscorr.json"))

@testset "K values" begin
    calced_Ks = Kgen.calc_Ks(TempC=check_Ks["input_conditions"]["TC"], Sal=check_Ks["T"])

    for (k, v) in calced_Ks
        # assert almost equal calced_k to check_k
    end

end

@testset "Pressure correction factors" begin
    pressure_factors = Kgen.prescorr(TC=check_presscorr["input_conditions"]["TC"], P=check_presscorr["input_conditions"]["P"])
    
    for (k, f) in pressure_factors
        # assert f almost equal to reference values
    end
end