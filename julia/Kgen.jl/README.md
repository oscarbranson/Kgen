# Kgen.jl

Julia implementation of [Kgen](https://github.com/PalaeoCarb/Kgen) — stoichiometric
equilibrium constants for the marine carbonate system, calculated at a given temperature,
salinity and pressure, and optionally corrected for non-modern seawater Mg/Ca composition
using a polynomial approximation of [MyAMI](https://github.com/PalaeoCarb/pymyami).

Kgen is maintained in Python, R, MATLAB and Julia. All four are held to agreement within
0.01% by the cross-language test harness in `crosscheck/`.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/PalaeoCarb/Kgen", subdir="julia/Kgen.jl")
```

## Usage

`calc_Ks` returns all 13 constants as a `NamedTuple`:

```julia
using Kgen

ks = calc_Ks(temp_c=25.0, sal=35.0, p_bar=300.0)
ks.K1     # 1.8632444248743343e-6
```

`calc_K` returns a single constant:

```julia
calc_K(:K1, temp_c=25.0, sal=35.0)          # 1.4212669153166358e-6
calc_K("K1", temp_c=25.0, sal=35.0)         # strings work too
```

Non-modern seawater composition, in mol/kgsw:

```julia
calc_Ks(temp_c=25.0, sal=35.0, magnesium=0.03, calcium=0.02)
```

### Constants

`:K0 :K1 :K2 :KW :KB :KS :KspA :KspC :KP1 :KP2 :KP3 :KSi :KF`, all on the **total** pH scale.

### Arguments

| argument | meaning | default |
|---|---|---|
| `temp_c` | temperature in Celsius | `25.0` |
| `sal` | salinity in PSU | `35.0` |
| `p_bar` | pressure in bar | `0.0` |
| `magnesium` | *average* seawater Mg in mol/kgsw | `0.0528171` |
| `calcium` | *average* seawater Ca in mol/kgsw | `0.0102821` |
| `sulphate` | total sulphate in mol/kgsw | calculated from `sal` |
| `fluorine` | total fluorine in mol/kgsw | calculated from `sal` |
| `MyAMI_mode` | `:approximate` only (see below) | `:approximate` |

### Many conditions at once

Julia does not broadcast keyword arguments, so both functions also take their numeric
arguments positionally, in the order `temp_c, sal, p_bar, magnesium, calcium`:

```julia
calc_K.(:K1, [0.0, 10.0, 25.0], 35.0)        # Vector{Float64}
calc_Ks.([0.0, 10.0, 25.0], 35.0)            # Vector{NamedTuple}
```

The scalar functions are type-stable and allocation-free, so broadcasting is the fast path —
there is no separate vectorised implementation.

## Differences from the other implementations

`MyAMI_mode=:calculate`, which runs the full MyAMI model rather than the polynomial
approximation, is **not implemented in Julia**. Passing it raises an informative error.
Use the Python implementation if you need the full model.

The pressure correction takes its gas constant from `coefficients/fundamental_constants.json`
(`R_P = 83.14472`, the exact conversion of R = 8.314472 J/(mol K)), as MATLAB does. Python and
R instead hardcode 83.1451, the value inherited from CO2SYS. This shifts pressure-corrected Ks
by ~2e-6 relative to those two — twenty times inside the 1e-4 cross-language tolerance, and
zero at the surface.

## Development

```julia
Pkg.develop(path="julia/Kgen.jl")
Pkg.test("Kgen")
```

The coefficient files in `src/coefficients/` are copies of those in the repository root, so
that the package is installable on its own. `polynomial_coefficients.json` is refreshed by
`update_pymyami.py` at the repository root, alongside the R and MATLAB copies.
