# AFLOW

[![Paper](https://img.shields.io/badge/HighEntropyAlloysMater-s44210--025--00058--2-0B51B1)](https://doi.org/10.1007/s44210-025-00058-2)
[![Paper](https://img.shields.io/badge/ComMatSci-2022.111889-0B51B1)](https://doi.org/10.1016/j.commatsci.2022.111889)
[![Paper](https://img.shields.io/badge/ComMatSci-2015.07.019-0B51B1)](https://doi.org/10.1016/j.commatsci.2015.07.019)

The Automatic-Flow (__AFLOW__) framework is a class of interconnected software tools designed for high-throughput materials
discovery, leveraging both first-principles quantum mechanical calculations using density functional theory (DFT), and 
data informatics as a screening strategy.

In tandem, [aflow.org](https://aflow.org), a web ecosystem of applications mirroring the functionality of __AFLOW__,
and databases following the FAIR (Findable, Accessible, Interoperable, and Reusable) principles provide both graphical 
and programmatic methods to access over 3.5 million material entries. The databases also contain over 2000 crystallographic 
prototypes, used to uniquely classify crystal structures. For users who require only a subset of __AFLOW__ features, 
such as to create convex hulls for thermodynamic analysis, or who are unfamiliar with command line interfaces, 
[aflow.org](https://aflow.org) is an ideal choice to find the necessary information for their research.

--- 

## Quick Start

The easiest way to install `aflow` is to utilize our binary releases via the different package managers or as direct download.

### Linux

#### snap

[Snap packages](https://snapcraft.io/) are a modern package format designed for Linux, developed by Canonical. They
encapsulate software and all its dependencies, making it easy to install applications without worrying about system-wide
conflicts or missing libraries. It ensures also that `aflow` is kept up to date. Snaps are secure by default, running in
isolated environments (sandboxes), and they work across various Linux distributions. The drawback is that `aflow` can
only interact with files in your home directory.

```shell
snap install aflow
```

#### .deb (Ubuntu & Debian)

`.deb` packages are the standard software package format for Debian-based Linux distributions. They contain precompiled
binaries. We pre-build `.deb` packages for different systems, and they can be downloaded from the [GitHub Release]() page. 

They can either be installed by double-clicking on an Ubuntu desktop system or using `dpkg` on the commandline.
```shell
curl -O https://github.com/aflow-org/AFLOW/archive/refs/tags/aflow_4.0.4.deb #TODO update
sudo dpkg -i ./aflow_4.0.4.deb
```

#### .rpm (RedHat & CentOS) 

#### .sh (Self-extracting installer script)

This self-extracting installer allows for an interactive installation experience by simply running the following:
```shell
sh aflow-4.0.4-Linux.sh # TODO update
```
Use the `--help` option to see available installation options.

### MacOS

[Homebrew](https://brew.sh) is a widely-used and reliable package manager for macOS that simplifies installing software.
It allows you to easily install tools and applications without manually downloading or compiling them.

```shell
brew tap aflow-org/aflow
brew install aflow
```

## Installing from source

`aflow` can also be built directly from source. This method allows users to modify the source code and is recommended 
for the best possible performance.

### Get the source code

If `git` is available it is recommended to download the source code directly from GitHub.
```
git clone --recurse-submodules https://github.com/aflow-org/AFLOW.git
```
Otherwise, the full source can be downloaded on the [GitHub Release]() page.
```shell
curl -O https://github.com/aflow-org/AFLOW/archive/refs/tags/aflow_4.0.4.orig.tar.xz #TODO update
tar -xzmf aflow_4.0.4.orig.tar.xz
```

### Build tools
To build `aflow` you need [cmake](https://cmake.org/), [ninja](https://ninja-build.org/), and a C/C++ compiler.

On Debian based systems they can be installed with:
```shell
sudo apt update
sudo apt install build-essential cmake ninja-build pkg-config
```
In case the available versions of the build tools provided by your operating system are too old, install the tools locally 
from the [cmake](https://github.com/Kitware/CMake/releases) or [ninjia](https://github.com/ninja-build/ninja/releases) release pages. 


On macOS, the C/C++ compiler is provided by the _Xcode Command Line Tools_, which are installed with [Homebrew](https://brew.sh).
[Homebrew](https://brew.sh) is used to install the other build dependencies.
```shell
brew update
brew install cmake ninja
```

### Dependencies
`aflow` utilizes a few external libraries that can either be provided system-wide as shared libraries or built alongside 
`aflow` using the [vcpkg](https://vcpkg.io/en/) package manager. Using system-wide libraries reduces the build times and potentially reduces
the memory footprint of aflow instances. Utilizing [vcpkg](https://vcpkg.io/en/) will reduce the dependency on the host system which is 
especially helpful on HPC systems.

The following libraries are used:
- [libarchive](https://www.libarchive.org/)
  - for compression like [xz](https://xz.tukaani.org/) and [tar](https://manned.org/tar.1)
- [OpenSSL](https://www.openssl.org)
  - for hashing function
  - for SSL & TLS connections
- [libcurl](https://curl.se/libcurl/)
  - for querying AFLOW API

Install dependencies system-wide on Debian based systems:
```shell
sudo apt update
sudo apt install libarchive libssl-dev libcurl4-openssl-dev
```
Install dependencies system-wide on macOS:
```shell
brew update
brew install libarchive openssl@3 curl
```

For the static approach [vcpkg](https://vcpkg.io/en/) needs to be available. It can be [set up locally](https://learn.microsoft.com/en-us/vcpkg/get_started/get-started?pivots=shell-bash). 
It is important that [vcpkg](https://vcpkg.io/en/) is initialized and the `VCPKG_ROOT` environment variable is set. 

### Building and installing aflow
First the buildsystem needs to prepare the system. Different presets are available that already cover most needed cases.
- all available presets can be listed with `cmake -S src --list-presets`
- presets including `release` are optimized builds intended for production use, while `debug` should be used during development
- presets with `vcpkg` will build a static linked version of `aflow` that can be helpful for clustered systems, or where libraries can't be easy installed or requested by users
  - for all `vcpkg` presets the environment variable `VCPKG_ROOT` needs to be set
- customized presets can be added to a `CMakeUserPresets.json` (see [CMakeUserPresets.example.json](src/CMakeUserPresets.example.json) for an example)

To generate `aflow` using shared libraries follow these steps:
- create and set up the build folder with
  - `cmake -S src --preset release`
  - this step configures the build environment with the available tools (compilers, linkers ...)
  - it will inform you if there are missing dependencies
- build the `aflow` executable in parallel with 8 threads
  - `cmake --build build/release --parallel 8`
- `aflow` can then be executed
  - `build/release/aflow --version`
  - `ctest --test-dir build/release/ -L quick --parallel 8`
- if the binary works as expected it can be installed with:
  - `sudo cmake --install build/release/` 

To generate `aflow` alongside the needed libraries use the `release_vcpkg` preset:
```shell
cmake -S src --preset release_vcpkg
cmake --build build/release_vcpkg --parallel 8
build/release_vcpkg/aflow --version
ctest --test-dir build/release_vcpkg/ -L quick --parallel 8
sudo cmake --install build/release_vcpkg/
```

### Troubleshooting
- if the vcpkg build fails due to problems with an external library:
  - remove the binary cache of vcpkg `rm -rf ~/.cache/vcpkg/archives/`
  - remove the build folder `rm -rf build/release_vcpkg`
  - rerun the cmake preparation step `cmake --preset release_vcpkg`
  - note: this problem can happen if the compiler or base libraries on an HPC machine changes
- Make sure you have no spaces in the path to vcpkg or the project. Vcpkg will complain and not work.



## Testing
Multiple tests are defined in [CMakeLists.txt](src/CMakeLists.txt) to check the integrity of the build `aflow` binary. When `aflow` is 
built from source the full test suite can be run in the build directory. For a quick check a special test group is 
available: `ctest -L quick --parallel 8`, while `ctest -E 'ut.structure_gen.proto|ut.aflowlib.lib2raw' --parallel 8`
should be used for a detailed test run. If a test failed, `ctest --rerun-failed --parallel 8 --extra-verbose` reruns only 
the failed test with more details.

If aflow is installed from a binary release, a subset of tests can be run by using  `aflow --unit_test`. A common quick 
test would be `aflow --unit_test aurostd`.

## Documentation
The usage of `aflow` is described in multiple peer-reviewed publications. These should also be cited if `aflow` is used
in any scientific context.

> **AFLOW4: Heading Toward Disorder**  
> S. Divilov, H. Eckert, S.D. Thiel, S.D. Griesemer, R. Friedrich, N.H. Anderson, M.J. Mehl, D. Hicks, M. Esters, 
> N. Hotz, X. Campilongo, A. Calzolari, and S. Curtarolo  
> *High Entropy Alloys & Materials*  (2025)  
> [DOI: 10.1007/s44210-025-00058-2](https://doi.org/10.1007/s44210-025-00058-2)

> **aflow++: A C++ framework for autonomous materials design**  
> C. Oses, M. Esters, D. Hicks, S. Divilov, H. Eckert, R. Friedrich, M.J. Mehl, A. Smolyanyuk, X. Campilongo,
> A. van de Walle, J. Schroers, A.G. Kusne, I. Takeuchi, E. Zurek, M.B. Nardelli, M. Fornari, Y. Lederer, O. Levy, 
> C. Toher, and S. Curtarolo  
> *Computational Materials Science* **217**, 111889 (2023).  
> [DOI: 10.1016/j.commatsci.2022.111889](https://doi.org/10.1016/j.commatsci.2022.111889)
> | [BibTeX](http://materials.duke.edu/auro/AUROARTICULA/BIBITEMS/curtarolo:art191.txt) 

> **The AFLOW standard for high-throughput materials science calculations**  
> C.E. Calderon, J.J. Plata, C. Toher, C. Oses, O. Levy, M. Fornari, A. Natan, M. J. Mehl,
> G. Hart, M.B. Nardelli, and S. Curtarolo  
> *Computational Materials Science* **108A**, 233-238 (2015).  
> [DOI: 10.1016/j.commatsci.2015.07.019](https://doi.org/10.1016/j.commatsci.2015.07.019)
> | [BibTeX](http://materials.duke.edu/auro/AUROARTICULA/BIBITEMS/curtarolo:art104.txt)

To document the inner working of `aflow` [doxygen](https://doxygen.nl/) is used. Details on how to build the code documentation and 
style are described in [docs/README.md](docs/README.md).
