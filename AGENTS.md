# Working on storm-analysis

Notes for automated coding tools, and for anyone touching the code for the
first time. `README.md` is the entry point for *using* this package; this file
is about changing it.

Everything below is checked against the tree rather than remembered, so if
something here disagrees with the code, the code is right and this file is
stale. Please fix it.


## Build

The C libraries are built by SCons, not by setuptools, and `pip install -e .`
runs SCons for you through `_custom_build.py`:

```
pip install -e .
```

FFTW3 and LAPACK need to be present first. `doc/source/install.rst` has the
per-platform details, and the CI workflow in `.github/workflows` runs exactly
those steps on Linux, macOS and Windows.

**After editing any `.c` or `.h` file, rebuild before running anything:**

```
scons -Q
```

Nothing reminds you. The Python side loads whatever shared library is already
sitting in `storm_analysis/c_libraries/`, so a test run against a stale library
looks like a normal result and quietly answers a question you did not ask.

On Windows SCons has to be told which compiler to use; `_custom_build.py`
defaults to `compiler=mingw` and reads `SA_SCONS_COMPILER` if you want
something else.


## Tests

From the repository root:

```
pytest
```

280 tests, around 45 seconds. Configuration lives in the `[tool.pytest]` table
of `pyproject.toml` — that is pytest's TOML mode, which is why `minversion` is
9.0. The older `[tool.pytest.ini_options]` spelling means ini mode; both are
valid and neither is a mistake to be corrected into the other.

Two things about reading the result:

- **Trust the exit status, not the last line printed.** Several of the C
  libraries write progress to stdout through the C runtime. Those buffers are
  flushed when the interpreter exits, which is after pytest has printed its
  summary, so the final lines on screen are frequently library chatter and the
  `N passed` line is nowhere in sight. Use `--junitxml=results.xml` when you
  need counts rather than a yes or no.

- **The `__main__` blocks are not covered.** Most test modules end with a block
  that calls the tests so the file can be run directly. pytest never executes
  it. Rename a test function without updating that block and the suite stays
  green while the file is broken for anyone running it by hand — this has
  happened, and went unnoticed across six pull requests.


## Diagnostics

`storm_analysis/diagnostics/` holds end-to-end runs against simulated data.
They are the check that matters after changing a fitter or anything in the
analysis pipeline: the unit tests confirm the code runs, the diagnostics
confirm it still recovers the right positions. Read
`storm_analysis/diagnostics/README.txt` first — each one is a
`settings.py` / `configure.py` / `make_data.py` / `analyze_data.py` /
`collate.py` sequence run from a working directory.

Two things worth knowing before you start one:

- Set `STORM_ANALYSIS_HEADLESS=1` unless you are sitting in front of the
  machine. Without it several of these open a plot window and block until you
  close it. Figures written to disk are unaffected.
- `multiplane/configure.py` requires `--psf-model`, one of `psf_fft`,
  `pupilfn` or `spline`. It has no default and exits 2 without it. This is the
  only required argument in any of them. `fista_decon`, `spliner`,
  `spliner_2d` and `multiplane` also accept an optional `--no-splines`;
  everything else runs bare.


## Layout

| Path | What's in it |
| --- | --- |
| `sa_library/` | Shared machinery: fitters (`dao_fit.c`, `multi_fit.c`), readers and writers, the HDF5 localization format |
| `sa_utilities/` | Pipeline stages: drift correction, tracking, format conversion |
| `spliner/` | Cubic-spline PSF fitting, and the tools that measure a spline from bead data |
| `multi_plane/`, `pupilfn/`, `psf_fft/` | Multi-plane and pupil-function fitting |
| `simulator/` | Synthetic movie generation, used heavily by the tests and diagnostics |
| `diagnostics/` | End-to-end accuracy runs (above) |
| `test/` | The pytest suite |
| `c_libraries/` | Build output. Not source — SCons writes here |

Python talks to C through `ctypes`, and the wrappers are named for what they
wrap: `sa_library/dao_fit_c.py` is the interface to `sa_library/dao_fit.c`. A
change to a C signature has to be made in both places, and see the ctypes note
under Traps.


## Conventions

These are consistent across the tree and are deliberate. Matching them matters
more than matching any external style guide.

- **`import numpy`, never `import numpy as np`.** All 242 modules that import
  numpy do it this way, and `numpy.` appears in full at every use.
- **Class methods are camelCase** — `getPeaks()`, `cleanUp()`. 704 of them,
  against 6 in snake_case. Module-level functions are more mixed, and test
  functions are snake_case because pytest collects on the `test_` prefix.
- **Line endings are LF**, enforced by `.gitattributes`.

Where the tree is genuinely inconsistent, follow the file you are editing
rather than imposing a rule. Spacing around keyword arguments is the main one:
`f(x = 1)` and `f(x=1)` are close to evenly split overall, and 108 of 254 files
contain both.


## Traps

**`decon_storm/` is not built.** It contains `.c` files, but `SConstruct` has
no entry for them and nothing in the package imports it — it is mostly MATLAB
kept for reference. Two consequences: changes there have no effect, and a
function definition found there is not the one that runs. A search that keys C
definitions by name and stops at the first match will find the dead copy and
hide the live one.

**Stale bytecode can survive an edit.** Python decides whether a `.pyc` is
current from the source file's size and modification time, and mtime has
one-second resolution. An edit that leaves the size unchanged and lands within
the same second as the previous one keeps the old bytecode. This is a real
problem when comparing two versions of a small change: clear `__pycache__`
between the two runs or you will measure the same code twice.

**The ctypes boundary is checked on one side only.** A wrapper declaring
`ndpointer(dtype=numpy.int32)` guarantees what Python passes; it says nothing
about what the C function expects. If the two disagree the result is a wrong
stride through valid memory, not an exception. C functions on this boundary
take `int32_t` rather than `int` for that reason. Scalars passed as
`ctypes.c_int` are fine as they are.
