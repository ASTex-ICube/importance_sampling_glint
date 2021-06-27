Importance Sampling of Glittering BSDFs based on Finite Mixture Distributions
=============================================================================

Code associated to the article *Importance Sampling of Glittering BSDFs based on Finite Mixture Distributions*, by [Xavier Chermain](http://igg.unistra.fr/People/chermain/) ([ICUBE](https://icube.unistra.fr/en/)), 
[Basile Sauvage](https://igg.icube.unistra.fr/index.php/Basile_Sauvage)
([ICUBE](https://icube.unistra.fr/en/)), 
[Jean-Michel Dishler](https://dpt-info.u-strasbg.fr/~dischler/)
([ICUBE](https://icube.unistra.fr/en/)) and 
[Carsten Dachsbacher](https://cg.ivd.kit.edu/english/dachsbacher/)
([KIT](https://www.kit.edu/english/index.php)).

Accepted for [EGSR 2021](https://egsr.eu/2021/).

* [Project page](http://igg.unistra.fr/People/chermain/importance_sampling_glint/)
* [Paper](http://igg.unistra.fr/People/chermain/assets/pdf/Chermain2021ImportanceSampling.pdf)
* [Video](http://igg.unistra.fr/People/chermain/assets/mp4/Chermain2021ImportanceSampling.mp4)
* [Bibtex](http://igg.unistra.fr/People/chermain/assets/Chermain2021ImportanceSampling.txt)

## Outline

* [File organisation](#file_organisation)
* [Building instructions](#building)
* [Render scenes and material usage](#render)
* [Tool commands](#tools)
* [pbrt-v3 readme](#pbrt-v3)

## <a name="file_organisation"></a>File organisation

We integrate our glittering materials in the [`pbrt-v3`](https://github.com/mmp/pbrt-v3) renderer:
* The implementation of the glittering *conductor* material is located in the files `src/materials/glitteringconductor.*`.
* The implementation of the glittering *dielectric* material is located in the files `src/materials/glitteringdielectric.*`.
* The implementation of the glint tools (convergence comparisons and chi square tests) is located in the file `src/tools/glinttools.cpp`.

## <a name="building"></a>Building instructions

The following commands build the project
```
git clone --recursive https://github.com/ASTex-ICube/importance_sampling_glint.git
cd importance_sampling_glint
mkdir build-release
cd build-release
cmake ../
make -j<number of threads>
```
Add `pathtopbrt/build-release` to your `PATH` environment variable to use the
`pbrt` and ` glinttool` commands anywhere. See the `readme.md` of
[`pbrt-v3`](https://github.com/mmp/pbrt-v3) for more information
concerning the building.

## <a name="render"></a>Render scenes and material usage

Scenes using our glittering materials are available in the `pbrt-v3-scenes` folder. Here an example usage of one of them:
```
pbrt fig1_left.pbrt
```
where the command line is launched from the `pbrt-v3-scenes` folder and when the command `pbrt` is included in the `PATH` environment variable.

The following python script
```
pbrt-v3-scene/launchallscenes.py
```
launches all the renderings of the paper. Python <img src="https://render.githubusercontent.com/render/math?math=\ge 3.5"> with `subprocess` dependency is required to use the python script.

We use a modified version of the `pbrt-v3` path tracing algorithm for rendering.
In the original path tracer, the algorithm uniformly samples *one* light for each ray -- scene intersection.
In our version, the algorithm uniformly samples *all* lights for each intersection. See the source files `src/integrators/path.*` for more details.

In the following, we use the same parameter presentation as [`pbrt-v3`](https://github.com/mmp/pbrt-v3).

The `glitteringconductor` material models reflection from glittering conductors. Its parameters are:
| **Type** | **Name** | **Default Value** | **Description** | 
| -------- | -------- | ----------------- | --------------- |
| spectrum texture | eta | (copper) | Index of refraction to use in computing the material's reflectance.  |
|spectrum texture | k | (copper) | Absorption coefficient to use in computing the material's reflectance.  |
|float texture | alphax | 0.5 | Beckmann roughness in the x direction.  |
|float texture | alphay | 0.5 | Beckmann roughness in the y direction.  |
|float texture | rho | 0. | Slope correlation factor.  |
|float texture | logmicrofacetdensity | 20. | The logarithm of the microfacet density, without the microfacet relative area parameter applied. Set to a high value (e.g. 40) to have a glossy material (without glints). |
|float texture | microfacetrelativearea | 1. | Percentage of the surface without microfacets. Note: *Effective* microfacet density = exp(logmicrofacetdensity) * microfacetrelativearea  |
|float texture | alphaxbasematerial | 0.01 | Roughness in the x direction for the base material[^1]. |
|float texture | alphaybasematerial | 0.01 | Roughness in the y direction for the base material[^1]. |
|float texture | rhobasematerial | 0. | Slope correlation factor for the base material[^1]. |
|float | densityrandomisation | 2. | Randomly changes the density of microfacets per cell. More precisely, this parameter is the standard deviation of a normal distribution sampled to randomise the microfacet density. |
|bool | fresnelnoop | false | If true, the Fresnel term is always 1 (useful for white furnace test).|
|bool | samplevisiblearea | true | If true, samples the visible area of the normal distribution function.|
|bool | sampleapproximation | false | If true, samples the Gaussian approximation of the normal distribution function.|
|spectrum texture | dictionary | n/a | Dictionary of multi-scale, piecewise linear 1D distributions. |
|integer | nlevels | n/a | Number of levels of detail of the dictionary. |
|integer | N | n/a | Number of multi-scale 1D distributions in the dictionary. |
|float | alpha_dict | n/a | Roughness used to generate the dictionary. |

[^1]: These parameters are only used when the microfacetrelativearea < 1, i.e., when (1 - microfacetrelativearea) of the surface is covered by a base material.

The `glitteringdielectric` material models reflection and transmission from glittering dielectrics. This material has the same parameters as the `glitteringconductor` material, without `eta`, `k` and `fresnelnoop`, and with the following specific parameters.

| **Type** | **Name** | **Default Value** | **Description** |
| -- | -- | -- | -- |
|spectrum texture | Kr | 1 | The reflectivity of the surface.  | 
|spectrum texture | Kt | 1 | The transmissivity of the surface.  |
|float texture | index | 1.5 | The index of refraction of the inside of the object. (`pbrt` implicitly assumes that the exterior of objects is a vacuum, with IOR of 1.)  |

## <a name="tools"></a>Tool commands

We provide the following `C++` commands with the `glinttool`:
* `plotglitteringndf`
* `chisquaretestglitteringvndf`
* `chisquaretestglitteringbrdf`
* `chisquaretestglitteringbsdf`
* `convergencecomparisons`

The source code of these commands can be found in `src/tools/glinttool.cpp`.
We also provide python scripts which call these `C++` commands. Python <img src="https://render.githubusercontent.com/render/math?math=\ge 3.5"> with `matplotlib`, `numpy` and `subprocess` dependencies is required to use the python scripts (`glinttool` must also be included in the `PATH` environment variable).

These commands generate python scripts using the `matplotlib` library. Execute the generated script to visualise the result. For example:
```
python generatedfile.py
```
For Chi square tests, the commands also display the result test result (success or failure).

### Plot glittering NDF

Manual of the command `plotglitteringndf`:
```
plotglitteringndf: Plots the glittering ndf of Chermain et al. 2020. Filename 1:
path to the dictionary. Filename 2: output matplotlib filename.
  options:
    --imagesize        Size of the output image. Default: 256
    --alphax           Alpha roughness of the surface in the s direction.
                       Default: 0.5
    --alphay           Alpha roughness of the surface in the t direction.
                       Default: 0.5
    --dsdx             Partial derivative of s (first component of the surface
                       position) with respect to x (first component of the
                       pixel coordinate). Default: 0.0005
    --dtdx             Partial derivative of t (second component of the surface
                       position) with respect to x (first component of the
                       pixel coordinate). Default: 0.0
    --dsdy             Partial derivative of s (first component of the surface
                       position) with respect to y (second component of the
                       pixel coordinate). Default: 0.0
    --dtdy             Partial derivative of t (second component of the surface
                       position) with respect to y (second component of the
                       pixel coordinate). Default: 0.0005
```

Example:
```
glinttool plotglitteringndf -imagesize 256 -alphax 0.6 -alphay 0.6 -dsdx
0.00052 -dtdy 0.00052 dict_N768_nLevels8.exr plotglitteringndf.py
```
where `dict_N768_nLevels8.exr` is located in the current directory. See also the python script using this command:
```
command_plotglitteringndf/command_plotglitteringndf.py
```
We use this command to plot the NDFs in the paper.

### Chi square tests

The correctness of our sampling algorithms is verified with chi square tests. We validate the sampling procedure of the glittering VNDF, BRDF and BSDF with the commands `chisquaretestglitteringvndf`, `chisquaretestglitteringbrdf` and `chisquaretestglitteringbsdf`, respectively.

#### Glittering VNDF

Manual of the command `chisquaretestglitteringvndf`:
```
chisquaretestglitteringvndf: Validates our sampling of the vndf of Chermain et
al. 2020 with a chi square test. Filename 1: path to the dictionary.
Filename 2: output matplotlib filename. The plot shows differences between the
analytic PDF and the histogram built by sampling the PDF.
  options:
    --nsamplehisto     Number of samples to compute the histogram.
                       Default: 1000000
    --res              Size of integration grid. Default: 512
    --alphax           Alpha roughness of the surface in the s direction.
                       Default: 0.3
    --alphay           Alpha roughness of the surface in the t direction.
                       Default: 0.3
    --rho              Slope correlation factor. Default: 0
    --stx              X component of the pixel footprint center. Default: 0.
    --dstdxx           Partial derivative of s (first component of the surface
                       position) with respect to x (first component of the
                       pixel coordinate). Default: 0.001
    --dstdyy           Partial derivative of t (second component of the surface
                       position) with respect to y (second component of the
                       pixel coordinate). Default: 0.001
    --thetao           Polar angle of the observation direction. Default: 1.5
    --phio             Azimuthal angle of the observation direction.
                       Default: 0.
```

Example:
```
glinttool chisquaretestglitteringvndf -nsamplehisto 1000000 -res 512 -alphax 0.3
-alphay 0.3 -rho 0. -stx 0. -dstdxx 0.001 -dstdyy 0.001 -thetao 1.5 -phio 0.
dict_N768_nLevels8.exr plotchisquaretestglitteringvndf.py
```
where `dict_N768_nLevels8.exr` is located in the current directory. See also the python script using this command:
```
command_chisquaretestglitteringvndf/command_chisquaretestglitteringvndf.py
```

#### Glittering BRDF

Manual of the command `chisquaretestglitteringbrdf`:
```
chisquaretestglitteringbrdf: Validates our sampling of the glittering BRDF of
Chermain et al. 2020 with a chi square test. Filename 1: path to the dictionary.
Filename 2: output matplotlib filename. The plot shows differences between the
analytic PDF and the histogram built by sampling the PDF.
  options:
    --nsamplehisto     Number of samples to compute the histogram.
                       Default: 1000000
    --res              Size of integration grid. Default: 512
    --alphax           Alpha roughness of the surface in the s direction.
                       Default: 0.3
    --alphay           Alpha roughness of the surface in the t direction.
                       Default: 0.3
    --rho              Slope correlation factor. Default: 0
    --mra              Microfacet relative area. Default: 1.
    --stx              X component of the pixel footprint center. Default: 0.
    --dstdxx           Partial derivative of s (first component of the surface
                       position) with respect to x (first component of the
                       pixel coordinate). Default: 0.001
    --dstdyy           Partial derivative of t (second component of the surface
                       position) with respect to y (second component of the
                       pixel coordinate). Default: 0.001
    --thetao           Polar angle of the observation direction. Default: 0.2
    --phio             Azimuthal angle of the observation direction.
                       Default: 0.
```

Example:
```
glinttool chisquaretestglitteringbrdf -nsamplehisto 1000000 -res 512 -alphax 0.3
-alphay 0.3 -rho 0. -stx 0. -dstdxx 0.001 -dstdyy 0.001 -thetao 0.2 -phio 0. -mra 1.
dict_N768_nLevels8.exr plotchisquaretestglitteringbrdf.py
```
where `dict_N768_nLevels8.exr` is located in the current directory. See also the python script using this command:
```
command_chisquaretestglitteringbrdf/command_chisquaretestglitteringbrdf.py
```

#### Glittering BSDF

Manual of the command `chisquaretestglitteringbsdf`:
```
chisquaretestglitteringbsdf: Validates our sampling of the glittering BSDF with
a chi square test. Filename 1: path to the dictionary. Filename 2: output
matplotlib filename. The plot shows differences between the analytic PDF and
the histogram built by sampling the PDF.
  options:
    --nsamplehisto     Number of samples to compute the histogram.
                       Default: 256000000
    --res              Size of integration grid. Default: 4096
    --alphax           Alpha roughness of the surface in the s direction.
                       Default: 0.25
    --alphay           Alpha roughness of the surface in the t direction.
                       Default: 0.25
    --rho              Slope correlation factor. Default: 0
    --mra              Microfacet relative area. Default: 1.
    --stx              X component of the pixel footprint center. Default: 0.
    --dstdxx           Partial derivative of s (first component of the surface
                       position) with respect to x (first component of the
                       pixel coordinate). Default: 0.001
    --dstdyy           Partial derivative of t (second component of the surface
                       position) with respect to y (second component of the
                       pixel coordinate). Default: 0.001
    --thetao           Polar angle of the observation direction. Default: 0.2
    --phio             Azimuthal angle of the observation direction.
                       Default: 0.
```

Example:
```
glinttool chisquaretestglitteringbsdf -nsamplehisto 256000000 -res 4096 -alphax 0.25
-alphay 0.25 -rho 0. -stx 0. -dstdxx 0.001 -dstdyy 0.001 -thetao 0.2 -phio 0. -mra 1.
dict_N768_nLevels8.exr plotchisquaretestglitteringbsdf.py
```
where `dict_N768_nLevels8.exr` is located in the current directory. See also the python script using this command:
```
command_chisquaretestglitteringbsdf/command_chisquaretestglitteringbsdf.py
```

#### Convergence comparisons

Manual of the command `convergencecomparisons`:
```
convergencecomparisons: Compares convergences of two importance sampling
schemes. The first (our) uses sampling of the multi-lobe component of the
glittering BSDF, and the second (previous) uses sampling of the gaussian
approximation of the multi-lobe component. Filename 1: path to the dictionary.
Filename 2: output matplotlib filename. Output: Two matplotlib files. The first
named <filename 2>_radiance_n.py shows eight convergence curves for each
sampling strategy. The second, named <filename 2>_pointwise_boxplots.py shows
the pointwise boxplot of the <nruns> for each sampling strategy. 
  options:
    --nsamples         Number of samples to compute the integral.
                       Default: 10000
    --nruns            Number of runs. Default: 8
    --pgf              If 1, matplotlib will use pgf exporter. Default: 1
    --alphax           Alpha roughness of the surface in the s direction.
                       Default: 0.3
    --alphay           Alpha roughness of the surface in the t direction.
                       Default: 0.3
    --rho              Slope correlation factor. Default: 0
    --mra              Microfacet relative area. Default: 1.
    --stx              X component of the pixel footprint center. Default: 0.
    --dstdxx           Partial derivative of s (first component of the surface
                       position) with respect to x (first component of the
                       pixel coordinate). Default: 0.001
    --dstdyy           Partial derivative of t (second component of the surface
                       position) with respect to y (second component of the
                       pixel coordinate). Default: 0.001
    --thetao           Polar angle of the observation direction. Default: 0.2
    --phio             Azimuthal angle of the observation direction.
                       Default: 0.
```

Example:
```
glinttool convergencecomparisons -stx 0. -thetao 1.5 -alphax 0.6 -alphay 0.6
-dstdxx 0.0001 -dstdyy 0.0001 -pgf 0 -nruns 8 dict_N768_nLevels8.exr plot.py
```
where `dict_N768_nLevels8.exr` is located in the current directory. See also the python scripts using this command:
```
command_convergencecomparisons/convergencecomparisons.py
command_convergencecomparisons/allconvergencecomparisons.py
```
The last python script generates the pointwise boxplots of the [supplemental material 1](http://igg.unistra.fr/People/chermain/assets/pdf/Chermain2021ImportanceSamplingSupplemental1.pdf).

<a name="pbrt-v3"></a> pbrt, Version 3
===============

[![Build Status](https://travis-ci.org/mmp/pbrt-v3.svg?branch=master)](https://travis-ci.org/mmp/pbrt-v3)
[![Build status](https://ci.appveyor.com/api/projects/status/mlm9g91ejxlcn67s/branch/master?svg=true)](https://ci.appveyor.com/project/mmp/pbrt-v3/branch/master)

This repository holds the source code to the version of pbrt that is
described in the third edition of *Physically Based Rendering: From
Theory to Implementation*, by [Matt Pharr](http://pharr.org/matt), [Wenzel
Jakob](http://www.mitsuba-renderer.org/~wenzel/), and Greg Humphreys.  As
before, the code is available under the BSD license.

The [pbrt website](http://pbrt.org) has general information about both the
*Physically Based Rendering* book as well as many other resources for pbrt.
As of October 2018, the full [text of the book](http://www.pbr-book.org) is
now available online, for free.

Example scenes
--------------

Over 8GB of example scenes are available for download. (Many are new and
weren't available with previous versions of pbrt.)  See the [pbrt-v3 scenes
page](http://pbrt.org/scenes-v3.html) on the pbrt website for information
about how to download them.

After downloading them, see the `README.md.html` file in the scene
distribution for more information about the scenes and preview images.

Additional resources
--------------------

* There is a [pbrt Google
  Groups](https://groups.google.com/forum/#!forum/pbrt) mailing list that can
  be a helpful resource.
* Please see the [User's Guide](http://pbrt.org/users-guide.html) for more
  information about how to check out and build the system as well as various
  additional information about working with pbrt.
* Should you find a bug in pbrt, please report it in the [bug
  tracker](https://github.com/mmp/pbrt-v3/issues).
* Please report any errors you find in the *Physically Based Rendering*
  book to authors@pbrt.org.

Note: we tend to let bug reports and book errata emails pile up for a few
months for processing them in batches. Don't think we don't appreciate
them. :-)

Building pbrt
-------------

To check out pbrt together with all dependencies, be sure to use the
`--recursive` flag when cloning the repository, i.e.
```bash
$ git clone --recursive https://github.com/mmp/pbrt-v3/
```
If you accidentally already cloned pbrt without this flag (or to update an
pbrt source tree after a new submodule has been added, run the following
command to also fetch the dependencies:
```bash
$ git submodule update --init --recursive
```

pbrt uses [cmake](http://www.cmake.org/) for its build system.  On Linux
and OS X, cmake is available via most package management systems.  To get
cmake for Windows, or to build it from source, see the [cmake downloads
page](http://www.cmake.org/download/).  Once you have cmake, the next step
depends on your operating system.

### Makefile builds (Linux, other Unixes, and Mac) ###

Create a new directory for the build, change to that directory, and run
`cmake [path to pbrt-v3]`. A Makefile will be created in the current
directory.  Next, run `make` to build pbrt, the obj2pbrt and imgtool
utilities, and an executable that runs pbrt's unit tests.  Depending on the
number of cores in your system, you will probably want to supply make with
the `-j` parameter to specify the number of compilation jobs to run in
parallel (e.g. `make -j8`).

By default, the makefiles that are created that will compile an optimized
release build of pbrt. These builds give the highest performance when
rendering, but many runtime checks are disabled in these builds and
optimized builds are generally difficult to trace in a debugger.

To build a debug version of pbrt, set the `CMAKE_BUILD_TYPE` flag to
`Debug` when you run cmake to create build files to make a debug build.  To
do so, provide cmake with the argument `-DCMAKE_BUILD_TYPE=Debug` and build
pbrt using the resulting makefiles. (You may want to keep two build
directories, one for release builds and one for debug builds, so that you
don't need to switch back and forth.)

Debug versions of the system run much more slowly than release
builds. Therefore, in order to avoid surprisingly slow renders when
debugging support isn't desired, debug versions of pbrt print a banner
message indicating that they were built for debugging at startup time.

### Xcode ###

To make an Xcode project on OS X, run `cmake -G Xcode [path to pbrt-v3]`.
A `PBRT-V3.xcodeproj` project file that can be opened in Xcode.  Note that
the default build settings have an optimization level of "None"; you'll
almost certainly want to choose "Faster" or "Fastest".

### MSVC on Windows ###

On Windows, first point the cmake GUI at the directory with pbrt's source
code.  Create a separate directory to hold the result of the build
(potentially just a directory named "build" inside the pbrt-v3 directory)
and set that for "Where to build the binaries" in the GUI.

Next, click "Configure".  Note that you will want to choose the "Win64"
generator for your MSVC installation unless you have a clear reason to need
a 32-bit build of pbrt.  Once cmake has finished the configuration step,
click "Generate"; when that's done, there will be a "PBRT-V3.sln" file in
the build directory you specified. Open that up in MSVC and you're ready to
go.

### Build Configurations ###

There are two configuration settings that must be set when configuring the
build. The first controls whether pbrt uses 32-bit or 64-bit values for
floating-point computation, and the second controls whether tristimulus RGB
values or sampled spectral values are used for rendering.  (Both of these
aren't amenable to being chosen at runtime, but must be determined at
compile time for efficiency).  The cmake configuration variables
`PBRT_FLOAT_AS_DOUBLE` and `PBRT_SAMPLED_SPECTRUM` configure them,
respectively.

If you're using a GUI version of cmake, those settings should be available
in the list of configuration variables; set them as desired before choosing
'Generate'.

With command-line cmake, their values can be specified when you cmake via
`-DPBRT_FLOAT_AS_DOUBLE=1`, for example.
