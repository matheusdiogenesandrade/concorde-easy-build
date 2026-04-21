# Concorde

[Concorde](http://www.math.uwaterloo.ca/tsp/concorde/index.html) is likely the most famous solver for the [Travelling Salesman Problem](http://www.math.uwaterloo.ca/tsp/concorde/index.html).
It was developed by four legends in the Operational Research community: David Applegate, Robert Bixby, Vasek Chvatal, and William Cook.
The solver is blazing fast, and is able to solve enormous instances.

However, as most of its users know, building Concorde can be painful.
So painful that people have written
[technical reports](https://www.researchgate.net/publication/324485167_Concorde_solver_installation_and_use),
blog posts
[[1](https://www.leandro-coelho.com/install-and-run-concorde-with-cplex/),
[2](https://www.leandro-coelho.com/installing-concorde-tsp-with-cplex-linux/),
[3](https://hackaday.io/project/158802-improve-tool-path-planning-in-cura/log/147747-using-concorde-tsp-solver)],
stackoverflow questions
[[4](https://stackoverflow.com/questions/48284456/concorde-installation-need-to-link-an-lp-solver-to-use-this-function),
[5](https://stackoverflow.com/questions/29056498/cant-build-concorde-tsp-solver-on-mac-yosemite)],
(to link a few) dedicated to this topic.

![](build-concorde.gif)

# Build Concorde like it's no longer 1997

This repository allows you to build concorde using modern toolchains.

## Standalone build (binary in `./build`)

Configure and compile **from this directory** so all artifacts stay under **`build/`** next to `CMakeLists.txt` (self-contained; no parent project required):

```bash
export CPLEX_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio2211   # Studio root: must contain cplex/ and concert/
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build --target concorde_cli -j$(nproc)
```

After a successful build you will have at least:

| Output | Path |
| :--- | :--- |
| Upstream TSP driver | **`build/concorde_cli`** |
| Static library merge | `build/libconcorde_full.a` |
| Shared library | `build/libconcorde_full.so` |
| Amalgamated header | `build/concorde.h` |

Run the driver on a TSPLIB file (symmetric `.tsp` only):

```bash
./build/concorde_cli -s 1 data/tsp/att48.tsp
```

Optional smoke (same `build/` tree; requires **`bash`** and **`python3`**):

The default smoke runs **`concorde_cli`** on all **111** symmetric `.tsp` fixtures under `data/tsp`, in the same **discover → sort by `DIMENSION` → 10-way partition** order as [`test_concorde_symmetric_suite.cpp`](../../tests/test_concorde_symmetric_suite.cpp) in `concorde_wrapper` (bucket 1/10 through 10/10). Only **process exit code** is checked; the C++ suite proves optimality against `solutions.txt`.

This can take **substantial wall time** (many exact solves). `ctest` sets a **7200 s** timeout on `concorde_cli_tsp_smoke`.

```bash
export CONCORDE_BIN="$(pwd)/build/concorde_cli"
export TSP_DIR="$(pwd)/data/tsp"
bash scripts/run_concorde_cli_tsp_smoke.sh
```

- **`QUICK=1`**: run only `att48.tsp`, `eil51.tsp`, `pr76.tsp` (fast local check).
- **`ALL=1`**: same as default (full bucket-aligned run); kept as a no-op alias for older docs.

Print bucket layout without invoking Concorde:

```bash
python3 scripts/run_concorde_cli_tsp_buckets.py --tsp-dir data/tsp --list-only
```

CTest (after configure; CMake runs **`find_package(Python3 COMPONENTS Interpreter REQUIRED)`** when `CONCORDE_BUILD_CLI` is on): `cd build && ctest -R concorde_cli_tsp_smoke` (or from anywhere: `ctest --test-dir build -R concorde_cli_tsp_smoke`).

- Turn off the CLI target with **`-DCONCORDE_BUILD_CLI=OFF`**.
- **`build/`** is a normal CMake out-of-tree directory; keep it out of version control if your workflow ignores local build trees.

### Built via `concorde_wrapper` / COPAlgorithms

When this folder is pulled in with `add_subdirectory` from [`concorde_wrapper`](../CMakeLists.txt), CMake may place the **same** targets under the **parent** build directory (e.g. `…/build/Release/external/concorde_wrapper/concorde-easy-build/`) instead of `./build` here. To always use **`./build/concorde_cli`**, configure this directory standalone as above.

## License

In this repository, I store a copy of the 2003 version of Concorde, together with a minimal build system.
Let me stress that this is a totally pirate version of the software.
In no way any of the authors have authorised me to fork their solver and create this repository.
Therefore, I have in no way any claim on authorship of this software.
If you are using Concorde from this repository you have to agree to exactly the same license as if you had downloaded it from the [official website](http://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm).
In short, this means that Concorde is free for you if and only if you use it for academic purposes.
In all other cases you should contact [William Cook](mailto:bico@uwaterloo.ca) to discuss licensing options.

For what concern my contribution, which is basically a `CMakeLists.txt` file, I release it to public domain under the [unlicense](https://unlicense.org/).

## Assumptions

I assumed that:

* You are on Linux, with a recent version of GCC.
* You use CMake as your build system.
* You use a recent version of CPLEX as your LP solver.

Users not corresponding to this identikit might have to adapt this solution.

I also assume that you are interested in having the following at the end of a **standalone** build (under `build/`):

* The upstream CLI **`concorde_cli`** (from `src/TSP/concorde.c`).
* **`libconcorde_full.so`**, **`libconcorde_full.a`**, and **`concorde.h`** for embedding Concorde in your own code (you still need to link CPLEX and the usual `pthread` / `m` / `dl`).

Older docs referred to `concorde-bin`; that target was renamed to **`concorde_cli`** to avoid clashing with the `concorde` **library** target in CMake.
