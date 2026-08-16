using DelimitedFiles

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "julia", "Kgen.jl"))
Pkg.instantiate()

using Kgen

# 1. Load test_conditions.csv
raw, header = readdlm(joinpath(@__DIR__, "test_conditions.csv"), ',', Float64, '\n'; header=true)
columns = Dict(strip(name) => raw[:, i] for (i, name) in enumerate(vec(header)))

temp_c = columns["temp_c"]
sal = columns["sal"]
p_bar = columns["p_bar"]
magnesium = columns["magnesium"]
calcium = columns["calcium"]

# 2. Run Kgen over those inputs. The scalar functions are broadcast rather than vectorised;
#    see julia/Kgen.jl/README.md.
Ks = calc_Ks.(temp_c, sal, p_bar, magnesium, calcium)

# 3. Save to ./generated_Ks as julia_approximated.csv.
#    Julia only implements the polynomial approximation of MyAMI, so it produces no
#    julia_calculated.csv - crosscheck.py discovers whatever files are present.
output = joinpath(@__DIR__, "generated_Ks", "julia_approximated.csv")
open(output, "w") do io
    println(io, join(String.(Kgen.K_NAMES), ','))
    for row in Ks
        println(io, join(row, ','))
    end
end

println("Wrote $output")
