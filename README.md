# Cgl 0.60.9

[![A COIN-OR Project](https://coin-or.github.io/coin-or-badge.png)](https://www.coin-or.org)

Projects such as this one are maintained by a small group of volunteers under
the auspices of the non-profit [COIN-OR Foundation](https://www.coin-or.org)
and we need your help! Please consider [sponsoring our
activities](https://github.com/sponsors/coin-or) or [volunteering](mailto:volunteer@coin-or.org) to help!

[![Latest Release](https://img.shields.io/github/v/release/coin-or/Cgl?sort=semver)](https://github.com/coin-or/Cgl/releases)

_This file is auto-generated from [config.yml](.coin-or/config.yml) using the 
[generate_readme](.coin-or/generate_readme) script.
To make changes, please edit [config.yml](.coin-or/config.yml) or the generation scripts
[here](.coin-or/generate_readme) and [here](https://github.com/coin-or/coinbrew/blob/master/scripts/generate_readme)._

The COIN-OR Cut Generation Library (`Cgl`) is a collection of cut generators that can be 
used with other COIN-OR packages that make use of cuts, such as, among others, the linear solver 
[Clp](https://github.com/coin-or/Clp) or the mixed integer linear programming solvers 
[Cbc](https://github.com/coin-or/Cbc) or [BCP](https://github.com/coin-or/Bcp). 
`Cgl` uses the abstract class `OsiSolverInterface` 
(see [Osi](https://github.com/coin-or/Osi)) to use or communicate with a solver.
It does not directly call a solver.

Each cut generator is in a separate directory with its own maintainer.
All generators are combined in one library when `Cgl` is compiled.

Available cut generators are: 

 * Combinatorial cuts:
   * [CglAllDifferent](https://github.com/coin-or/Cgl/wiki/CglAllDifferent)
   * CglBKClique
   * [CglClique](https://github.com/coin-or/Cgl/wiki/CglClique)
   * [CglKnapsackCover](https://github.com/coin-or/Cgl/wiki/CglKnapsackCover)
   * [CglOddHole](https://github.com/coin-or/Cgl/wiki/CglOddHole)
   * CglOddWheel
   * CglZeroHalf

 * Flow cover cuts:
   * [CglFlowCover](https://github.com/coin-or/Cgl/wiki/CglFlowCover)

 * Gomory cuts and variants:
   * [CglGomory](https://github.com/coin-or/Cgl/wiki/CglGomory)
   * CglGMI
   * [CglRedSplit](https://github.com/coin-or/Cgl/wiki/CglRedSplit)
   * CglRedSplit2

 * Lift-and-project cuts:
   * [CglLiftAndProject](https://github.com/coin-or/Cgl/wiki/CglLiftAndProject)
   * [CglLandP](https://github.com/coin-or/Cgl/wiki/CglLandP)

 * Mixed integer rounding cuts and variants:
   * [CglMixedIntegerRounding](https://github.com/coin-or/Cgl/wiki/CglMixedIntegerRounding)
   * [CglMixedIntegerRounding2](https://github.com/coin-or/Cgl/wiki/CglMixedIntegerRounding2)
   * [CglTwomir](https://github.com/coin-or/Cgl/wiki/CglTwomir)
   * [CglResidualCapacity](https://github.com/coin-or/Cgl/wiki/CglResidualCapacity)

 * Strengthening:
   * CglCliqueStrengthening
   * [CglDuplicateRow](https://github.com/coin-or/Cgl/wiki/CglDuplicateRow)
   * [CglPreprocess](https://github.com/coin-or/Cgl/wiki/CglPreprocess)
   * [CglProbing](https://github.com/coin-or/Cgl/wiki/CglProbing)
   * [CglSimpleRounding](https://github.com/coin-or/Cgl/wiki/CglSimpleRounding)

CoinUtils is an open-source collection of classes and helper functions
that are generally useful to multiple COIN-OR projects.
These utilities include:
 * classes for storing and manipulating sparse matrices and vectors,
 * performing matrix factorization,
 * parsing input files in standard formats, e.g. MPS,
 * building representations of mathematical programs,
 * performing simple presolve operations,
 * warm starting algorithms for mathematical programs,
 * comparing floating point numbers with a tolerance
 * classes for storing and manipulating conflict graphs, and
 * classes for searching and storing cliques and odd cycles in conflict graphs, among others.

The project managers of Cgl are Robin Lougee (@rlougee) and Francois Margot.


Cgl is written in C++ and is released as open source under the [Eclipse Public License 2.0](http://www.opensource.org/licenses/EPL-2.0).

It is distributed under the auspices of the [COIN-OR Foundation](https://www.coin-or.org).

The Cgl development site is https://github.com/coin-or/Cgl.

## CITE

Code: [![DOI](https://zenodo.org/badge/173502902.svg)](https://zenodo.org/badge/latestdoi/173502902)

## CURRENT BUILD STATUS

[![Windows Builds](https://github.com/coin-or/Cgl/actions/workflows/windows-ci.yml/badge.svg?branch=releases/0.60.9)](https://github.com/coin-or/Cgl/actions/workflows/windows-ci.yml?query=branch%3Areleases/0.60.9)

[![Linux and MacOS Builds](https://github.com/coin-or/Cgl/actions/workflows/linux-ci.yml/badge.svg?branch=releases/0.60.9)](https://github.com/coin-or/Cgl/actions/workflows/linux-ci.yml?query=branch%3Areleases/0.60.9)

## DOWNLOAD

What follows is a quick start guide for obtaining or building
Cgl on common platforms. More detailed information is
available [here](https://coin-or.github.io/user_introduction.html).

### Docker image

There is a Docker image that provides Cgl, as well as other projects
in the [COIN-OR Optimization
Suite](https://github.com/coin-or/COIN-OR-OptimizationSuite) [here](https://hub.docker.com/repository/docker/coinor/coin-or-optimization-suite)

### Binaries

For newer releases, binaries will be made available as assets attached to
releases in Github
[here](https://github.com/coin-or/Cgl/releases). Older binaries
are archived as part of Cbc
[here](https://www.coin-or.org/download/binary/Cbc).

 * *Linux* (see https://repology.org/project/coin-or-cgl/versions for a complete listing): 
   * arch:
     ```
     $ sudo pacman -S  coin-or-cgl
     ```
   * Debian/Ubuntu:
     ```
     $ sudo apt-get install  coinor-cgl coinor-libcgl-dev
     ```
   * Fedora/Redhat/CentOS:
     ```
     $ sudo yum install  coin-or-Cgl coin-or-Cgl-devel
     ```
   * freebsd:
     ```
     $ sudo pkg install math/cgl
     ```
   * linuxbrew:
     ```
     $ brew install cgl
     ```
 * *Windows*: The easiest way to get Cgl on Windows is to download an archive as described above.
 * *Mac OS X*: The easiest way to get Cgl on Mac OS X is through [Homebrew](https://brew.sh).
     ```
     $ brew tap coin-or-tools/coinor
     $ brew install coin-or-tools/coinor/cgl
     ```

* *conda* (cross-platform, no Windows for now):
     ```
     $ conda install coin-or-cgl
     ```

Due to license incompatibilities, pre-compiled binaries lack some 
functionality. If binaries are not available for your platform for the latest 
version and you would like to request them to be built and posted, feel free 
to let us know on the mailing list. 

### Source

Source code can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of Cgl from the
 [releases](https://github.com/coin-or/Cgl/releases) page,
 * Cloning this repository from [Github](https://github.com/coin-or/Cgl), or 
 * Using the [coinbrew](https://github.com/coin-or/coinbrew) script to get the project and all dependencies (recommended, see below).   

### Dependencies

Cgl has a number of dependencies, which are detailed in
[config.yml](.coin-or/config.yml). Dependencies on other COIN-OR projects are
automatically downloaded when obtaining the source with `coinbrew`. For some
of the remaining third-party dependencies, automatic download scripts and
build wrappers are provided (and will also be automatically run for required
and recommended dependencies), while other libraries that are aeasy to obtain
must be installed using an appropriate package manager (or may come with your
OS by default). 

## BUILDING from source

These quick start instructions assume you are in a bash shell. 

### Using `coinbrew`

To download and build Cgl from source, execute the 
following on the command line. 
```
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew fetch Cgl@0.60.9
./coinbrew build Cgl
```
For more detailed instructions on coinbrew, see https://coin-or.github.io/coinbrew.
The `coinbrew` script will fetch the additional projects specified in the Dependencies section of [config.yml](.coin-or/config.yml).

### Without `coinbrew` (Expert users)

 * Download the source code, e.g., by cloning the git repo https://github.com/coin-or/Cgl
 * Download and install the source code for the dependencies listed in [config.yml](.coin-or/config.yml)
 * Build the code as follows (make sure to set PKG_CONFIG_PTH to install directory for dependencies).

```
./configure -C
make
make test
make install
```

## Doxygen Documentation

If you have `Doxygen` available, you can build a HTML documentation by typing

`make doxydoc` 

in the build directory. If Cgl was built via `coinbrew`, then the build
directory will be `./build/Cgl/0.60.9` by default. The doxygen documentation main file
is found at `<build-dir>/doxydoc/html/index.html`.

If you don't have `doxygen` installed locally, you can use also find the
documentation [here](http://coin-or.github.io/Cgl/Doxygen).


## Project Links

 * [Code of Conduct](https://www.coin-or.org/code-of-conduct/)
 * [COIN-OR Web Site](http://www.coin-or.org/)
 * [COIN-OR general discussion forum](https://github.com/orgs/coin-or/discussions)
 * [Cgl Discussion forum](https://github.com/coin-or/Cgl/discussions)
 * [Report a bug](https://github.com/coin-or/Cgl/issues/new)
 * [Doxygen generated documentation](http://coin-or.github.io/Cgl/Doxygen)

---------

## Information for Subproject Managers

A cut generator in `Cgl` must conform to the following:

 * Its main class `CglCutGeneratorDeriv` is derived from the class `CglCutGenerator`.
 * It has three related classes used for data, parameters and information with respect to the enumeration tree:
   * A class `CglDataDeriv` derived from `CglData`; it should contain pointers on all data used by the generator that might be obtained from an `OsiSolverInterface` object when calling `generateCuts()` with an `OsiSolverInterface` object as parameter. The class `CglDataDeriv` might be `CglData` if the latter is sufficient. An exception is made for generators needing information deemed too expensive to collect from the solver (for example the optimal Simplex tableau); in this case `CglDataDeriv` might still contain a pointer on the `OsiSolverInterface` object, but its use should be limited to obtaining the "expensive" information from the solver.
   * A class `CglParamDeriv` derived from `CglParam`. It should contain parameters of the generator that can be set by the user. The parameters in the class `CglParamDeriv` must be taken into account during the cut generation. The class `CglParamDeriv` might be `CglParam` if the latter is sufficient.
   * A class `CglTreeInfoDeriv` derived from `CglTreeInfo`. The class `CglTreeInfoDeriv` might be `CglTreeInfo` if the latter is sufficient.
 * The class `CglCutGeneratorDeriv` must have 
   * A member of type `CglParamDeriv` used to store the current parameters.
   * A method `getParam()` that returns the object storing the current parameters.
   * A method `generateCuts(const OsiSolverInterface & si, OsiCuts & cs, const CglTreeInfoDeriv info)`
   * A method `generateCuts(const CglDataDeriv &data, OsiCuts & cs, const CglTreeInfoDeriv info)`
 * The data class `CglDataDeriv` must have methods `getMember()` and `setMember()` for each `member` of the class. Data members in `CglData` irrelevant for a generator are completely ignored. If a data member that is used by a generator is not available when `generateCuts(const CglDataDeriv &data, OsiCuts & cs, const CglTreeInfoDeriv info)` is called, the call is aborted, as if no cuts were found. A warning message might be printed.
 * The class `CglParamDeriv` must have methods `getMember()` and `setMember()` for each `member` of the class. All parameters must have default values. Each cut generator with a derived class is free to change the default values for all the members of `CglParamDeriv`, including those from `CglParam`.
 * Once an object of the cut generator class is created, it should be possible to call generateCuts() several times in a row without having to destroy and re-create the object.
 * By default, a successful call to `generateCuts()` should not generate any output. If an error occurs, a message might be printed.

