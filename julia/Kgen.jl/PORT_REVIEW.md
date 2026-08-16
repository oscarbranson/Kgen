# Kgen.jl — critical review and remediation plan

> **Working document.** This records why the Julia port is being rewritten and tracks the work.
> Decide before merging whether to keep it (as a record) or delete it.
>
> **To resume an interrupted session:** read the Progress tracker below and continue from the
> first unticked box.

## Progress

- [x] Phase 0  plan persisted to `julia/Kgen.jl/PORT_REVIEW.md`
- [x] Phase 1a vendor coefficients into `src/coefficients/`; add Julia to `update_pymyami.py` + `RELEASE.md`
- [x] Phase 1b `const` NTuple coefficient tables; concrete `Real` signatures
- [x] Phase 1c fix KP3 coefficient index (bug 1)
- [x] Phase 1d fix inverted pH-scale conversions in BOTH call sites (bug 2)
- [x] Phase 1e factor shared pressure-correction block; hoist out of `calc_Ks` loop
- [x] Phase 1f `calc_Ks` returns NamedTuple; delete `mutable struct Ks`
- [x] Phase 1g `calc_K` dispatcher + `ArgumentError`; `sulphate`/`fluorine` kwargs; `MyAMI_mode` validation
- [x] Phase 1h rename `fn_*` → `calc_*` per 0.3.0 harmonisation; move docstrings above functions; export
- [x] Phase 1i **added during implementation:** positional methods for `calc_K`/`calc_Ks`, because
      Julia does not broadcast keyword arguments — the keyword-only API could not be vectorised at all
- [x] Phase 2  rewrite `PyMYAMI.jl` (NTuple, no `hcat`, keyword args)
- [x] Phase 3  `Project.toml` compat; drop `Manifest.toml` + nested `.github`; write README
- [x] Phase 4a rewrite `test/runtests.jl` with real assertions (incl. integer-input test) — 58/58 pass
- [x] Phase 4b root `.github/workflows/julia-tests.yml` (+ `make test-julia`)
- [x] Phase 4c `crosscheck/gen_julia.jl` — no registration needed, `crosscheck.py` globs
      `generated_Ks/*.csv`. Julia added to `test_polynomial_coefficients`; makefile updated.
- [x] Phase 4d **found during testing:** `include_dependency` on the coefficient JSON — the
      `const` tables are built at precompile time, so without it a coefficient refresh was
      silently ignored until the cache happened to be cleared
- [x] Resolved: MyAMI correction now applied **unconditionally**, matching R/MATLAB — see below
- [x] Verification steps 1–7 all pass

## Resolved — MyAMI gating semantics

Python gates the correction on `np.any(magnesium != 0.0528171)` (array-global: if any row is
non-modern, every row is corrected); R and MATLAB apply it unconditionally. The first draft of
this rewrite decided per row, which left the 125 modern-composition rows of
`test_conditions.csv` uncorrected and put Julia 9.0e-4 outside the 1e-4 crosscheck tolerance —
because the polynomial returns up to 1.0009 at modern composition where it should return 1.

**Decision (user, 2026-08-16): apply unconditionally**, matching R and MATLAB. The residual is
well inside the repo's own `APPROX_TOLERANCE` of 4.1e-3. Consequence to be aware of: Julia's
`calc_K(:K1)` at default modern composition differs from Python's *scalar* `calc_K('K1')` by
~7.6e-6, because Python's `np.any` skips the correction for scalar inputs. Python is
internally inconsistent here (scalar skips, array corrects); no single Julia behaviour matches
both.

Crosscheck result with both sides on the same pymyami build: **max relative difference
5.6e-10** across all 3125 conditions, ~180,000× inside tolerance.

## Gas constant now read from fundamental_constants.json

`coefficients/fundamental_constants.json` was not initially vendored, because nothing in the
Julia code read it — `calc_pressure_correction` hardcoded `83.1451`, exactly as Python
(`K_functions.py:333`) and R (`K_functions.R:229`) still do. R ships the file unused; MATLAB
is the only implementation that actually reads it.

Reviewing this surfaced a latent cross-language inconsistency: the file says
`R_P = 83.14472` — the exact conversion of R = 8.314472 J/(mol K) to cm³ bar / (mol K) —
while three of the four implementations hardcode CO2SYS's legacy `83.1451`. The effect is
~2e-6 on pressure correction factors, i.e. **50× below the crosscheck tolerance, so the
harness structurally cannot detect it.**

**Decision (user, 2026-08-16):** Julia now reads `R_P` from the vendored file, siding with
MATLAB and with the physically correct value. Consequences, all verified:

- Crosscheck vs Python: worst case 4.6e-6, still 20× inside the 1e-4 tolerance. Surface rows
  are unaffected (5.6e-10), confirming the divergence is exactly and only the gas constant.
- `check_presscorr.json` still passes at `atol=1e-5`; the shift is ~1.6e-6.
- Test suite pins the value and asserts `GAS_CONSTANT != 83.1451` so a silent revert fails.

**Still open at repo level:** Python and R should arguably be moved onto `R_P` too, so all
four agree and the JSON is genuinely the single source of truth. That changes their published
outputs slightly and was left out of scope here.

## Note on pymyami versions

The repo pins `pymyami==2.1.0` in `requirements.txt`, `setup.cfg`, `update_pymyami.py` and
both CI workflows, and the vendored coefficient files in R, MATLAB and now Julia are all
byte-identical 2.1.0. A local environment running a different build (e.g. 2.1.3 in
`~/.venvs/py3`) makes `test_polynomial_coefficients` fail for **all three** languages
uniformly — that is environment drift, not a code defect.

**Decisions taken** (confirmed with the user before starting):

1. Adopt the cross-language harmonised API (`calc_K(K; temp_c, sal, p_bar, magnesium, calcium)`),
   matching Python/R/MATLAB. Breaking, but the package is unreleased (v0.1.0, unregistered).
2. Keep the core **scalar and type-stable**; broadcasting (`calc_K.(...)`) is the vectorisation
   story rather than explicit array methods.
3. Edits are manually approved as they are made.

---

## Context

**What Kgen.jl is for.** Kgen is a multi-language toolkit (Python, R, MATLAB, Julia) that
calculates the 13 stoichiometric equilibrium constants of the marine carbonate system —
`K0 K1 K2 KW KB KS KspA KspC KP1 KP2 KP3 KSi KF` — from temperature, salinity and pressure,
with an optional correction for non-modern seawater Mg/Ca composition via MyAMI. The whole
point of the repo is that the *four* implementations agree to within 0.01%
(`crosscheck/crosscheck.py`, `RDIFF_TOLERANCE = 0.0001`).

`julia/Kgen.jl/` is the Julia port. The only change in the working tree at the time of review
was the directory rename `julia/Kgen` → `julia/Kgen.jl` (correct Julia convention). The code
itself was unchanged since the 2024-07 port and **had never been run by CI or the crosscheck**.

**Why this change is needed.** The review below found that the port is numerically wrong in
two places, cannot be installed as a package, has no working tests, and is ~130× slower than
it should be. It is not currently fit to publish. The intended outcome is a Julia package
that is registrable, numerically identical to the Python reference, joined to the crosscheck
harness, and idiomatically fast.

Everything below was verified by running the code, not read off the source.

---

## Review findings

### P0 — Correctness (the port produces wrong numbers)

**1. `fn_KP3` indexes the wrong coefficient — KP3 is wrong by a factor of 194.**
`src/Kgen.jl:254` reads `(p[6] / TK + p[6]) * S`; it should be `(p[5] / TK + p[6]) * S`.
Python's `calc_KP3` (`python/kgen/K_functions.py:230-236`) uses `coefficients[4]`/`coefficients[5]`.
Measured:

```
Kgen.calc_K("KP3", TempC=25.0, Sal=35.0)  →  ln K = -14.968
check_values/check_Ks.json                →  ln K = -20.24
ratio = 194.5×
```

**2. Both pH-scale conversion factors are inverted — every pressure-corrected K is off by ~0.37%.**
`src/Kgen.jl:406-407`, duplicated at `:479-480`. Julia has the reciprocals of what Python
(`K_functions.py:457-458`), R (`r/R/main.R:109,112`) and MATLAB (`kgen_static.m:244-245`) use.
This is exactly the pre-fix code from commit `5ea188f` (2023-07, *"Fixed the tot_to_sws and
sws_to_tot conversions"*), which touched matlab/python/r but not julia — and the Julia
pH-scale block was then written a year later reproducing the bug. Measured at
T = 25 °C, S = 35, P = 300 bar:

```
reference K1 = 1.8632444e-06
Kgen.jl  K1 = 1.8563834e-06     →  -0.368%,  i.e. 37× the crosscheck tolerance
```

**3. `MyAMI_mode="calculate"` throws `UndefVarError`.** `composition_correct` is assigned only
inside the `if MyAMI_mode == "approximate"` branch (`:452-459`) but read unconditionally at
`:485`. Confirmed: `UndefVarError(:composition_correct, …, :local)`. There is no full-MyAMI
path in Julia at all, and no validation of the mode string — Python raises a clear
`ValueError`. Julia's default is also `"approximate"` where Python's is `"calculate"`.

**4. The MyAMI polynomial coefficients Julia reads are stale.** `update_pymyami.py` writes
`matlab/polynomial_coefficients.json` and `r/inst/coefficients/polynomial_coefficients.json`
only. The root `coefficients/polynomial_coefficients.json` that Julia alone reads was added by
the Julia author and is never resynced — `md5sum` confirms it differs from the R/MATLAB pair
(which match each other). `RELEASE.md`'s pymyami checklist has no Julia entry.

### P1 — The package is not installable and not tested

**5. Coefficient JSON lives outside the package.** `project_path("..", "coefficients", …)`
(`src/Kgen.jl:9-10`, `src/PyMYAMI.jl:8`) resolves through the `julia/coefficients →
../coefficients` symlink to the repo root. `Pkg.add`/`Pkg.develop` copies only
`julia/Kgen.jl/`, so the path won't exist and — because these are top-level calls — the module
fails at **load** time. Python, R and MATLAB all ship the data inside the package.

**6. `test/runtests.jl` cannot run.** Every line is broken: `using CSV` (not in `[deps]` or
`[extras]`), `JSON.parsefile` without `using JSON`, `check_Ks["T"]` (no such key), integer
`TC=25` against an `::AbstractFloat` signature, iterating a non-iterable `Ks` struct,
`prescorr` called without its required `p` kwarg — and **both `@testset` bodies contain zero
assertions**, just `# assert almost equal` comments. Bugs 1 and 2 would have been caught
instantly by a real test.

**7. No CI, no crosscheck, no docs.** `julia/Kgen.jl/.github/workflows/` is dead — GitHub only
runs workflows at the repo root, and root `.github/workflows/` has python/r/matlab/crosscheck
jobs but nothing for Julia. There is no `crosscheck/gen_julia.jl` beside `gen_python.py`/
`gen_r.r`/`gen_matlab.m`. `grep -ril julia` over `.github/`, `crosscheck/`, `makefile`,
`README.md`, `HISTORY.md`, `RELEASE.md` returns nothing. README.md is the 7-byte string
`# kgen`. A stale `Manifest.toml` (resolved under Julia 1.10.4, listing Revise and friends
that aren't in `Project.toml`) is committed, which a library should not do. `[compat]` has no
bound for JSON, which blocks General-registry registration.

**8. Docstrings are inside function bodies.** All 20-odd `"""…"""` blocks sit *after* the
`function` line, making them inert string literals rather than docstrings. Confirmed with
`length(Docs.meta(Kgen))`: the original module registers **0** documented bindings; the
rewrite registers 19. (An earlier draft of this review cited `Base.Docs.doc(Kgen.calc_K)`
throwing a `MethodError` as the evidence — that was wrong, as it throws regardless of whether
docs are attached. The finding itself stands on the `Docs.meta` count.) Three docstrings
(`fn_KW`, `fn_KB`, `fn_K0`) are also copy-paste errors reading *"Return K1 or K2…"*, and
`prescorr`'s `@param` list documents Python's `TC`/`sal` rather than its actual `P`/`TC`.
Nothing is `export`ed either.

### P2 — Not idiomatic, and slow

Measured on Julia 1.12.6:

| call | time | allocations |
|---|---|---|
| `calc_K("K1")` | 496 ns | 416 B |
| `calc_Ks()` | 15.8 µs | 9.3 kB |
| `calc_Ks(Pres=300.0)` | 59.8 µs | 46.3 kB |
| `calc_Ks(Mg=0.03, Ca=0.02)` | 35.7 µs | 34.2 kB |

`Base.return_types(Kgen.calc_K, (String,))` is **`Any`** — the whole package is type-unstable.
Causes, in order of impact:

- **Non-`const` globals holding `Dict{String,Any}`** whose values are `Vector{Any}`. Every
  coefficient access is a boxed, dynamically-typed load.
- **`mutable struct Ks` with 13 `::AbstractFloat` fields** — abstract fields mean pointer
  indirection and heap allocation per field; `setproperty!(result, Symbol(k), …)` (`:491`) is
  dynamic field assignment by runtime Symbol.
- **`::AbstractFloat` in every signature** is simultaneously too loose to help the compiler and
  too tight to accept users' input: `calc_K("K1", TempC=25, Sal=35)` is a `TypeError`.
  `python/test.py::test_calls` specifically tests integer inputs. `Real` is the right bound.
- **The 25-line pressure-correction block is copy-pasted verbatim** between `calc_K`
  (`:396-411`) and `calc_Ks` (`:469-484`) — and inside `calc_Ks` it sits *inside* the per-K
  loop, so `TS`, `TF`, `KS_surf`, `KF_surf` and two `prescorr` calls are recomputed 13×.
- **Keyword-only `p=`, `TK=`, `lnTK=`, `S=`, `sqrtS=` plumbing.** Passing pre-computed `lnTK`
  and `sqrtS` is a numpy-era micro-optimisation that Python itself deleted in `d06f00c`; in
  Julia it buys nothing and costs an unsplattable kwarg call at every site.
- `hcat` in `generate_polynomial_features` (`src/PyMYAMI.jl:20-79`) allocates a 1×5 matrix just
  to index it, then a 1×56 matrix, then does `Matrix * Vector{Any}`.
- `any(Mg != 0.0528171)` (`:414`) — vestigial numpy; a no-op on scalars. The same test is
  written *without* `any` at `:453`.
- `calc_K` never validates `k`; an unknown name gives a raw `KeyError` rather than Python's
  `ValueError(f'{K} is not valid. Should be one of …')`.
- `fn_*` / `calc_Istr` / `calc_TS` / `calc_TF` / `prescorr` are all names the other three
  languages deliberately migrated *away* from in the 0.3.0 harmonisation (`HISTORY.md`).
  Julia is the only implementation still on the old vocabulary.

**A prototype confirms the headroom is real.** Rewriting just the `K1` path with `const`
`NTuple{N,Float64}` coefficients, concrete `Float64` signatures and no kwarg plumbing:

```
prototype calc_K1(25.0, 35.0, 300.0) = 1.8632444248743343e-06   (matches reference)
  no pressure : 113.8 ns, 0 bytes
  Pres=300    : 118.4 ns, 0 bytes
return type   : Float64
```

That is **~4× faster than the current single-K call and zero-allocation**; extrapolated across
`calc_Ks` it is roughly 130× (15.8 µs → ~0.12 µs).

---

## Plan

### Phase 1 — Restructure `src/Kgen.jl`

Replace the module wholesale. Target shape:

- **Vendor the data.** Copy `coefficients/*.json` to `julia/Kgen.jl/src/coefficients/` (or add
  an `Artifacts.toml`; plain files are simpler and match R/MATLAB). Read them once at
  precompile time into `const` tables. Add a Julia section to `update_pymyami.py` alongside the
  existing matlab/r blocks (lines 41, 48) and a Julia line to `RELEASE.md`'s pymyami checklist.
  This also fixes finding 4 — the stale root `coefficients/polynomial_coefficients.json` should
  be removed once nothing reads it (note it is currently being bundled into the *Python* wheel
  by `package_data={'kgen': ['coefficients/*.json']}`, where it is unused).

- **`const` concrete coefficient tables.** A `NamedTuple` keyed by the 13 K names whose values
  are `NTuple{N,Float64}` (N = 5,6,7,9,11,12 per K — see `K_calculation.json`). `const`,
  concrete, immutable, stack-resident. Same for the 5-element pressure-correction tuples and
  the 56-element polynomial tuples (8 Ks only: `K0 K1 K2 KB KS KW KspA KspC`).

- **Signatures.** `calc_K1K2(coefficients, temp_c, sal)` style, matching Python: positional,
  `Real` bounds, deriving `temp_k`/`log`/`sqrt` internally. Rename per the 0.3.0 harmonisation:
  `fn_*` → `calc_*`, `calc_Istr` → `calc_ionic_strength`, `calc_TS`/`calc_TF` →
  `calc_sulphate`/`calc_fluorine`, `prescorr` → `calc_pressure_correction`.

- **Fix the two numerical bugs** (findings 1 and 2) while rewriting.

- **Factor the duplicated pressure block** into one `calc_pressure_correction_factors(temp_c,
  sal, p_bar, sulphate, fluorine)` returning the `tot_to_sws_surface`/`sws_to_tot_deep` pair,
  computed **once** outside the loop in `calc_Ks`. This removes the copy-paste *and* the 13×
  redundant recomputation in a single change.

- **`calc_Ks` returns a `NamedTuple`**, built by unrolled `map` over the coefficient
  `NamedTuple`. Deletes `mutable struct Ks`, the abstract fields and the dynamic
  `setproperty!` in one stroke; `ks.K1` still works and `ks[:K1]` does too.

- **`calc_K(K::Symbol; …)`** as a thin, explicit `if/elseif` dispatcher over the 13 names with
  an informative `ArgumentError` on an unknown name (finding: raw `KeyError`). Accept `String`
  too for parity with the other languages.

- **Add `sulphate`/`fluorine` keyword overrides** (Python has them, Julia doesn't), validate
  `MyAMI_mode` with a clear error, and either implement or explicitly `error()` on
  `"calculate"` rather than throwing `UndefVarError` (finding 3). Default should match Python's.

- **`export calc_K, calc_Ks`**, and move every docstring **above** its `function` line.
  Fix the three "Return K1 or K2" copy-paste docstrings and `calc_pressure_correction`'s
  wrong `@param` list.

### Phase 2 — `src/PyMYAMI.jl`

Rewrite `generate_polynomial_features` to return an `NTuple{56,Float64}` built directly from
`temp_k, log(temp_k), sal, mg, ca` with no `hcat`, and evaluate the dot product against the
`const NTuple{56,Float64}` coefficients. Zero allocation. **Keep the existing term ordering** —
it was independently verified correct against sklearn/pymyami and must not change. Switch to
keyword args for consistency with the rest of the package.

### Phase 3 — Package metadata

- Add `[compat] JSON = "0.21"` (registration blocker). Add `Test` to `[extras]`/`targets`
  plus whatever the tests actually need.
- Delete the committed `Manifest.toml` and add it to `.gitignore` (library convention).
- Delete `julia/Kgen.jl/.github/` — those workflows never execute from a subdirectory.
- Write a real `README.md`: install, `calc_K`/`calc_Ks` usage, broadcasting example, pointer to
  the parent repo.

### Phase 4 — Tests and CI

- Rewrite `test/runtests.jl` against `check_values/check_Ks.json` and
  `check_values/check_presscorr.json`, mirroring `python/test.py`: real `@test isapprox`
  assertions with sig-fig-derived tolerances, plus a `test_calls`-equivalent that passes
  **integers** (`temp_c=30, sal=36, p_bar=2`) to lock in the `Real` signatures. Check values
  must be vendored or read via a path that works from a test-time package install.
- Add `.github/workflows/julia-tests.yml` at the **repo root**, beside `python-tests.yml`.
- Add `crosscheck/gen_julia.jl` mirroring `gen_r.r` — read `test_conditions.csv`, broadcast
  `calc_Ks` over the rows, write `generated_Ks/julia_{calculated,approximated}.csv`. Register
  Julia in `crosscheck/crosscheck.py` and add it to the `test-crosscheck` makefile target.

### Ordering note

Do Phase 4's test rewrite *early* — ideally right after Phase 1's bug fixes — so the KP3 and
pH-scale corrections are demonstrated against `check_values/` rather than asserted.

---

## Verification

1. **Numerical correctness.** `cd julia/Kgen.jl && julia --project=. -e 'using Pkg; Pkg.test()'`.
   Specifically confirm `log(calc_K(:KP3, temp_c=25.0, sal=35.0)) ≈ -20.24` and that
   `calc_K(:K1, temp_c=25.0, sal=35.0, p_bar=300.0)` returns `1.8632444e-06`, not the current
   `1.8563834e-06`.
2. **Cross-language agreement.** `make test-crosscheck` with the Julia generator added — Julia
   must land within `RDIFF_TOLERANCE = 0.0001` of Python/R/MATLAB across all rows of
   `crosscheck/test_conditions.csv`.
3. **Installability.** From a clean temp depot, `Pkg.develop(path="julia/Kgen.jl")` in a
   directory *outside* the repo, then `using Kgen` — must load without a coefficient-path error.
4. **Type stability.** `Base.return_types(calc_K, (Symbol,))[1] === Float64`, and
   `@code_warntype calc_Ks(temp_c=25.0, sal=35.0)` shows no red. Also verify
   `calc_K(:K1, temp_c=25, sal=35)` (integers) no longer raises `TypeError`.
5. **Performance.** `@benchmark` (or an `@elapsed`/`@allocated` harness)
   `calc_Ks(temp_c=25.0, sal=35.0, p_bar=300.0)` — expect sub-microsecond and **0 bytes**
   allocated, against today's 59.8 µs / 46.3 kB.
6. **Docs attach.** `Base.Docs.doc(Kgen.calc_K)` returns the docstring instead of erroring.
7. **Broadcasting.** `calc_K.(:K1, temp_c=[0.0, 25.0], sal=[30.0, 35.0])` — sanity-check the
   documented vectorisation story actually works.