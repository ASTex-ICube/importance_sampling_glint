//
// glinttool.cpp
//
// Various useful operations on glittery BSDFs.
//

#include <core/api.h>
#include <ctype.h>
#include <glog/logging.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <fstream>

#include "fileutil.h"
#include "imageio.h"
#include "materials/glitteringconductor.h"
#include "materials/glitteringdielectric.h"
#include "microfacet.h"
#include "parallel.h"
#include "pbrt.h"
#include "rng.h"
#include "spectrum.h"

using namespace pbrt;

#define CHI2_RUNS 1
#define CHI2_MINFREQ 5
#define CHI2_SLEVEL 0.01

static void usage(const char *msg = nullptr, ...) {
    if (msg) {
        va_list args;
        va_start(args, msg);
        fprintf(stderr, "glinttool: ");
        vfprintf(stderr, msg, args);
        fprintf(stderr, "\n");
    }
    fprintf(stderr, R"(usage: glinttool <command> [options] <filenames...>

commands: plotglitteringndf, chisquaretestglitteringvndf,
          chisquaretestglitteringbrdf, chisquaretestglitteringbsdf,
          convergencecomparisons

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

)");
    exit(1);
}

BxDF *createGlitteringBRDF(
    MemoryArena &arena, Float alphaX, Float alphaY, Float rho,
    Float alphaXBaseMaterial, Float alphaYBaseMaterial, Float rhoBaseMaterial,
    Point2f st, Vector2f &dstdx, Vector2f &dstdy, int NLevels, int N,
    Float logMicrofacetDensity, Float alphaDict,
    std::vector<PiecewiseLinearDistribution1D> *distributions,
    Float microfacetRelativeArea, Float densityRandomization,
    bool sampleVisibleArea, bool sampleApproximation) {
    Fresnel *fresnel = ARENA_ALLOC(arena, FresnelNoOp)();
    RNG *rngEval = ARENA_ALLOC(arena, RNG)();
    RNG *rngSample = ARENA_ALLOC(arena, RNG)();

    std::pair<Distribution2D *, Distribution2D *> adjacentWP;
    std::pair<Bounds2i, Bounds2i> adjacentFootprintAABB;
    ClampRayFootprint(st, dstdx, dstdy);

    return ARENA_ALLOC(arena, GlitteringConductorReflection)(
        1, fresnel, rngEval, rngSample, st, dstdx, dstdy, alphaX, alphaY, rho,
        alphaXBaseMaterial, alphaYBaseMaterial, rhoBaseMaterial,
        distributions->at(0).Count() + 1, distributions, NLevels, N, alphaDict,
        logMicrofacetDensity, microfacetRelativeArea, densityRandomization,
        sampleVisibleArea, sampleApproximation);
}

BxDF *createGlitteringBSDF(
    MemoryArena &arena, Float alphaX, Float alphaY, Float rho,
    Float alphaXBaseMaterial, Float alphaYBaseMaterial, Float rhoBaseMaterial,
    Point2f st, Vector2f &dstdx, Vector2f &dstdy, int NLevels, int N,
    Float logMicrofacetDensity, Float alphaDict,
    std::vector<PiecewiseLinearDistribution1D> *distributions,
    Float microfacetRelativeArea, Float densityRandomization,
    bool sampleVisibleArea, bool sampleApproximation) {
    RNG *rngEval = ARENA_ALLOC(arena, RNG)();
    RNG *rngSample = ARENA_ALLOC(arena, RNG)();

    std::pair<Distribution2D *, Distribution2D *> adjacentWP;
    std::pair<Bounds2i, Bounds2i> adjacentFootprintAABB;
    ClampRayFootprint(st, dstdx, dstdy);

    return ARENA_ALLOC(arena, GlitteringDielectricScattering)(
        1., 1., 1., 1.5, TransportMode::Importance, rngEval, rngSample, st,
        dstdx, dstdy, alphaX, alphaY, rho, alphaXBaseMaterial,
        alphaYBaseMaterial, rhoBaseMaterial, distributions->at(0).Count() + 1,
        distributions, NLevels, N, alphaDict, logMicrofacetDensity,
        microfacetRelativeArea, densityRandomization, sampleVisibleArea,
        sampleApproximation);
}

int plotglitteringndf(int argc, char *argv[]) {
    MemoryArena arena;
    Options opt;
    pbrtInit(opt);

    int imageSize = 256;
    Float alphaX = 0.6;
    Float alphaY = 0.6;
    Float rho = 0.;
    Point2f st(0, 0);
    Vector2f dstdx(0.0005, 0.);
    Vector2f dstdy(0., 0.0005);

    int i;
    auto parseArg = [&]() -> std::pair<std::string, std::string> {
        const char *ptr = argv[i];
        // Skip over a leading dash or two.
        CHECK_EQ(*ptr, '-');
        ++ptr;
        if (*ptr == '-') ++ptr;

        // Copy the flag name to the string.
        std::string flag;
        while (*ptr && *ptr != '=') flag += *ptr++;
        if (!*ptr && i + 1 == argc)
            usage("missing value after %s flag", argv[i]);
        const char *value = (*ptr == '=') ? (ptr + 1) : argv[++i];
        return {flag, value};
    };

    std::pair<std::string, double> arg;
    for (i = 0; i < argc; ++i) {
        if (argv[i][0] != '-')
            break;
        else {
            std::pair<std::string, std::string> arg = parseArg();
            if (std::get<0>(arg) == "imagesize") {
                imageSize = int(atof(std::get<1>(arg).c_str()));
            } else if (std::get<0>(arg) == "alphax") {
                alphaX = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "alphay") {
                alphaY = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "dsdx") {
                dstdx[0] = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "dtdx") {
                dstdx[1] = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "dsdy") {
                dstdy[0] = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "dtdy") {
                dstdy[1] = atof(std::get<1>(arg).c_str());
            } else
                usage();
        }
    }

    if (i + 1 >= argc)
        usage("missing second filename for \"plotglitteringndf\"");
    else if (i >= argc)
        usage("missing filenames for \"plotglitteringndf\"");

    const char *dictFilename = argv[i], *outFilename = argv[i + 1];
    std::string dataPath(dictFilename);

    BxDF *bxdf = nullptr;

    std::shared_ptr<std::vector<PiecewiseLinearDistribution1D>> distributions =
        std::make_shared<std::vector<PiecewiseLinearDistribution1D>>();
    std::shared_ptr<ImageTexture<RGBSpectrum, Spectrum>> distributionsDict;

    int nLevels = 8;
    int N = 768;
    Float alphaXBaseMaterial = 0.1;
    Float alphaYBaseMaterial = 0.1;
    Float rhoBaseMaterial = 0.;
    Float logMicrofacetDensity = 20.;
    Float microfacetRelativeArea = 1.;
    Float densityRandomization = 0.01;
    bool sampleVisibleArea = true;
    bool sampleApproximation = false;

    // Load the texture
    std::unique_ptr<TextureMapping2D> map;
    Float su = 1.;
    Float sv = 1.;
    Float du = 0.;
    Float dv = 0.;
    map.reset(new UVMapping2D(su, sv, du, dv));
    distributionsDict = std::make_shared<ImageTexture<RGBSpectrum, Spectrum>>(
        std::move(map), dataPath, false, 8.f, ImageWrap::Clamp, 1.f, false);

    int distResolution = distributionsDict->Width() / (N / 3);
    std::vector<Float> dist1D(distResolution);
    for (int lvl = 0; lvl < nLevels; lvl++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < distResolution; j++) {
                dist1D[j] = distributionsDict->Texel(
                    0, i / 3 * distResolution + j, lvl)[i % 3];
            }
            distributions->push_back(
                PiecewiseLinearDistribution1D(dist1D.data(), distResolution));
        }
    }

    bxdf = createGlitteringBRDF(
        arena, alphaX, alphaY, rho, alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial, st, dstdx, dstdy, nLevels, N, logMicrofacetDensity,
        0.5, distributions.get(), microfacetRelativeArea, densityRandomization,
        sampleVisibleArea, sampleApproximation);

    GlitteringConductorReflection *glitteringConductorReflection =
        dynamic_cast<GlitteringConductorReflection *>(bxdf);

    std::cout << "Number of microfacets within pixel footprint: "
              << glitteringConductorReflection->GetNMicrofacetWithinP(st, dstdx,
                                                                      dstdy)
              << std::endl;

    std::vector<Float> result(imageSize * imageSize, 0.);
    Float integral = 0;

    for (int ix = 0; ix < imageSize; ++ix) {
        for (int iy = 0; iy < imageSize; ++iy) {
            Vector2f projectedNormal(Float(ix) / Float(imageSize),
                                     Float(iy) / Float(imageSize));
            projectedNormal = projectedNormal * 2. - Vector2f(1., 1.);
            float microNormalZSqr = 1. - projectedNormal.x * projectedNormal.x -
                                    projectedNormal.y * projectedNormal.y;
            if (microNormalZSqr <= 0.0001) continue;
            Vector3f microNormal(projectedNormal.x, projectedNormal.y,
                                 std::sqrt(microNormalZSqr));

            Float NDFValue = 0.f;
            GlitteringConductorReflection *glitteringConductorReflection =
                dynamic_cast<GlitteringConductorReflection *>(bxdf);

            NDFValue = glitteringConductorReflection->D_P(microNormal, st,
                                                          dstdx, dstdy);

            result[ix + iy * imageSize] = NDFValue;

            // integral
            integral += NDFValue * (2. / imageSize) * (2. / imageSize);
        }
    }

    std::cout << "Integral: " << integral << std::endl;

    // Extract the output file name without the extension
    std::string outputfilename(outFilename);
    std::size_t found = outputfilename.rfind(".");
    CHECK(found != std::string::npos);
    std::string outFilenameWithoutExtension =
        outputfilename.substr(std::size_t(0), found);

    // Create python file containing the plot
    std::ofstream pythonFile(outFilename);

    pythonFile << "import matplotlib.pyplot as plt" << std::endl;
    pythonFile << "import numpy as np" << std::endl;
    pythonFile << "pdf = [";
    // Write PDF
    for (int iy = 0; iy < imageSize; ++iy) {
        pythonFile << "[";
        for (int ix = 0; ix < imageSize; ++ix) {
            pythonFile << result.at(iy * imageSize + ix);
            if (ix != imageSize - 1)
                pythonFile << ",";
            else
                pythonFile << "]" << std::endl;
        }
        if (iy != imageSize - 1)
            pythonFile << ",";
        else
            pythonFile << "]" << std::endl;
    }
    pythonFile << "x = np.linspace(" << -1 << ", " << 1 << ", " << imageSize
               << ")" << std::endl;
    pythonFile << "y = np.linspace(" << -1 << ", " << 1 << ", " << imageSize
               << ")" << std::endl;
    pythonFile << "fig, axs = plt.subplots(1,1, figsize=(5, 5))" << std::endl;
    pythonFile << "pdf = np.array(pdf)" << std::endl;
    pythonFile << "vmaxval = pdf.max()" << std::endl;
    pythonFile << "pdf_plot = axs.imshow(pdf, "
                  "interpolation=\'sinc\', vmin=0, vmax=vmaxval, extent=[-"
               << 1 << ", " << 1 << ", -" << 1 << ", " << 1
               << "], origin=\'lower\')" << std::endl;
    pythonFile << "axs.set_axis_off()" << std::endl;
    // pythonFile << "axs.title.set_text(\'NDF\')" << std::endl;
    // pythonFile << "props = dict(fraction=0.046, pad=0.04)" << std::endl;
    // pythonFile << "fig.colorbar(pdf_plot, ax=axs, **props)" << std::endl;
    pythonFile << "plt.tight_layout()" << std::endl;
    pythonFile << "plt.savefig('" << outFilenameWithoutExtension
               << ".png', bbox_inches='tight')" << std::endl;
    // pythonFile << "plt.show()" << std::endl;

    pbrtCleanup();
    return 0;
}

int plotglitteringsdf(int argc, char *argv[]) {
    MemoryArena arena;
    Options opt;
    pbrtInit(opt);

    int imageSize = 256;
    Float alphaX = 0.6;
    Float alphaY = 0.6;
    Float rho = 0.;
    Point2f st(0, 0);
    Vector2f dstdx(0.0005, 0.);
    Vector2f dstdy(0., 0.0005);

    int i;
    auto parseArg = [&]() -> std::pair<std::string, std::string> {
        const char *ptr = argv[i];
        // Skip over a leading dash or two.
        CHECK_EQ(*ptr, '-');
        ++ptr;
        if (*ptr == '-') ++ptr;

        // Copy the flag name to the string.
        std::string flag;
        while (*ptr && *ptr != '=') flag += *ptr++;
        if (!*ptr && i + 1 == argc)
            usage("missing value after %s flag", argv[i]);
        const char *value = (*ptr == '=') ? (ptr + 1) : argv[++i];
        return {flag, value};
    };

    std::pair<std::string, double> arg;
    for (i = 0; i < argc; ++i) {
        if (argv[i][0] != '-')
            break;
        else {
            std::pair<std::string, std::string> arg = parseArg();
            if (std::get<0>(arg) == "imagesize") {
                imageSize = int(atof(std::get<1>(arg).c_str()));
            } else if (std::get<0>(arg) == "alphax") {
                alphaX = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "alphay") {
                alphaY = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "dsdx") {
                dstdx[0] = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "dtdx") {
                dstdx[1] = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "dsdy") {
                dstdy[0] = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "dtdy") {
                dstdy[1] = atof(std::get<1>(arg).c_str());
            } else
                usage();
        }
    }

    if (i + 1 >= argc)
        usage("missing second filename for \"plotglitteringndf\"");
    else if (i >= argc)
        usage("missing filenames for \"plotglitteringndf\"");

    const char *dictFilename = argv[i], *outFilename = argv[i + 1];
    std::string dataPath(dictFilename);

    BxDF *bxdf = nullptr;

    std::shared_ptr<std::vector<PiecewiseLinearDistribution1D>> distributions =
        std::make_shared<std::vector<PiecewiseLinearDistribution1D>>();
    std::shared_ptr<ImageTexture<RGBSpectrum, Spectrum>> distributionsDict;

    int nLevels = 8;
    int N = 768;
    Float alphaXBaseMaterial = 0.1;
    Float alphaYBaseMaterial = 0.1;
    Float rhoBaseMaterial = 0.;
    Float logMicrofacetDensity = 20.;
    Float microfacetRelativeArea = 1.;
    Float densityRandomization = 0.01;
    bool sampleVisibleArea = true;
    bool sampleApproximation = false;

    // Load the texture
    std::unique_ptr<TextureMapping2D> map;
    Float su = 1.;
    Float sv = 1.;
    Float du = 0.;
    Float dv = 0.;
    map.reset(new UVMapping2D(su, sv, du, dv));
    distributionsDict = std::make_shared<ImageTexture<RGBSpectrum, Spectrum>>(
        std::move(map), dataPath, false, 8.f, ImageWrap::Clamp, 1.f, false);

    int distResolution = distributionsDict->Width() / (N / 3);
    std::vector<Float> dist1D(distResolution);
    for (int lvl = 0; lvl < nLevels; lvl++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < distResolution; j++) {
                dist1D[j] = distributionsDict->Texel(
                    0, i / 3 * distResolution + j, lvl)[i % 3];
            }
            distributions->push_back(
                PiecewiseLinearDistribution1D(dist1D.data(), distResolution));
        }
    }

    bxdf = createGlitteringBRDF(
        arena, alphaX, alphaY, rho, alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial, st, dstdx, dstdy, nLevels, N, logMicrofacetDensity,
        0.5, distributions.get(), microfacetRelativeArea, densityRandomization,
        sampleVisibleArea, sampleApproximation);

    GlitteringConductorReflection *glitteringConductorReflection =
        dynamic_cast<GlitteringConductorReflection *>(bxdf);

    std::cout << "Number of microfacets within pixel footprint: "
              << glitteringConductorReflection->GetNMicrofacetWithinP(st, dstdx,
                                                                      dstdy)
              << std::endl;

    std::vector<Float> result(imageSize * imageSize, 0.);
    Float integral = 0;

    Float minMaxSlope = std::max(alphaX, alphaY) * 3.;

    for (int ix = 0; ix < imageSize; ++ix) {
        for (int iy = 0; iy < imageSize; ++iy) {
            Vector2f microSlope(Float(ix) / Float(imageSize),
                                Float(iy) / Float(imageSize));  // 0 - 1
            microSlope = microSlope * 2. * minMaxSlope -
                         Vector2f(minMaxSlope, minMaxSlope);

            Float SDFValue = 0.f;
            GlitteringConductorReflection *glitteringConductorReflection =
                dynamic_cast<GlitteringConductorReflection *>(bxdf);

            SDFValue = glitteringConductorReflection->P22_P(microSlope, st,
                                                            dstdx, dstdy);

            result[ix + iy * imageSize] = SDFValue;

            // integral
            integral += SDFValue * (2. * minMaxSlope / imageSize) * (2. * minMaxSlope / imageSize);
        }
    }

    std::cout << "Integral: " << integral << std::endl;

    // Extract the output file name without the extension
    std::string outputfilename(outFilename);
    std::size_t found = outputfilename.rfind(".");
    CHECK(found != std::string::npos);
    std::string outFilenameWithoutExtension =
        outputfilename.substr(std::size_t(0), found);

    // Create python file containing the plot
    std::ofstream pythonFile(outFilename);

    pythonFile << "import matplotlib.pyplot as plt" << std::endl;
    pythonFile << "import numpy as np" << std::endl;
    pythonFile << "pdf = [";
    // Write PDF
    for (int iy = 0; iy < imageSize; ++iy) {
        pythonFile << "[";
        for (int ix = 0; ix < imageSize; ++ix) {
            pythonFile << result.at(iy * imageSize + ix);
            if (ix != imageSize - 1)
                pythonFile << ",";
            else
                pythonFile << "]" << std::endl;
        }
        if (iy != imageSize - 1)
            pythonFile << ",";
        else
            pythonFile << "]" << std::endl;
    }
    pythonFile << "x = np.linspace(" << -1 << ", " << 1 << ", " << imageSize
               << ")" << std::endl;
    pythonFile << "y = np.linspace(" << -1 << ", " << 1 << ", " << imageSize
               << ")" << std::endl;
    pythonFile << "fig, axs = plt.subplots(1,1, figsize=(5, 5))" << std::endl;
    pythonFile << "pdf = np.array(pdf)" << std::endl;
    pythonFile << "vmaxval = pdf.max()" << std::endl;
    pythonFile << "pdf_plot = axs.imshow(pdf, "
                  "interpolation=\'sinc\', vmin=0, vmax=vmaxval, extent=[-"
               << 1 << ", " << 1 << ", -" << 1 << ", " << 1
               << "], origin=\'lower\')" << std::endl;
    pythonFile << "axs.set_axis_off()" << std::endl;
    // pythonFile << "axs.title.set_text(\'NDF\')" << std::endl;
    // pythonFile << "props = dict(fraction=0.046, pad=0.04)" << std::endl;
    // pythonFile << "fig.colorbar(pdf_plot, ax=axs, **props)" << std::endl;
    pythonFile << "plt.tight_layout()" << std::endl;
    pythonFile << "plt.savefig('" << outFilenameWithoutExtension
               << ".png', bbox_inches='tight')" << std::endl;
    // pythonFile << "plt.show()" << std::endl;

    pbrtCleanup();
    return 0;
}

int plotglitteringp22pdiscretelod(int argc, char *argv[]) {
    MemoryArena arena;
    Options opt;
    pbrtInit(opt);

    int imageSize = 256;
    Float alphaX = 0.6;
    Float alphaY = 0.6;
    Float rho = 0.;
    Point2f st(0, 0);
    Vector2f dstdx(0.0005, 0.);
    Vector2f dstdy(0., 0.0005);
    int l = 0, ldist = 0;

    int i;
    auto parseArg = [&]() -> std::pair<std::string, std::string> {
        const char *ptr = argv[i];
        // Skip over a leading dash or two.
        CHECK_EQ(*ptr, '-');
        ++ptr;
        if (*ptr == '-') ++ptr;

        // Copy the flag name to the string.
        std::string flag;
        while (*ptr && *ptr != '=') flag += *ptr++;
        if (!*ptr && i + 1 == argc)
            usage("missing value after %s flag", argv[i]);
        const char *value = (*ptr == '=') ? (ptr + 1) : argv[++i];
        return {flag, value};
    };

    std::pair<std::string, double> arg;
    for (i = 0; i < argc; ++i) {
        if (argv[i][0] != '-')
            break;
        else {
            std::pair<std::string, std::string> arg = parseArg();
            if (std::get<0>(arg) == "imagesize") {
                imageSize = int(atof(std::get<1>(arg).c_str()));
            } else if (std::get<0>(arg) == "alphax") {
                alphaX = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "alphay") {
                alphaY = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "dsdx") {
                dstdx[0] = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "dtdx") {
                dstdx[1] = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "dsdy") {
                dstdy[0] = atof(std::get<1>(arg).c_str());
            } else if (std::get<0>(arg) == "dtdy") {
                dstdy[1] = atof(std::get<1>(arg).c_str());
            } else
                usage();
        }
    }

    if (i + 1 >= argc)
        usage("missing second filename for \"plotglitteringndf\"");
    else if (i >= argc)
        usage("missing filenames for \"plotglitteringndf\"");

    const char *dictFilename = argv[i], *outFilename = argv[i + 1];
    std::string dataPath(dictFilename);

    BxDF *bxdf = nullptr;

    std::shared_ptr<std::vector<PiecewiseLinearDistribution1D>> distributions =
        std::make_shared<std::vector<PiecewiseLinearDistribution1D>>();
    std::shared_ptr<ImageTexture<RGBSpectrum, Spectrum>> distributionsDict;

    int nLevels = 8;
    int N = 768;
    Float alphaXBaseMaterial = 0.1;
    Float alphaYBaseMaterial = 0.1;
    Float rhoBaseMaterial = 0.;
    Float logMicrofacetDensity = 20.;
    Float microfacetRelativeArea = 1.;
    Float densityRandomization = 0.01;
    bool sampleVisibleArea = true;
    bool sampleApproximation = false;

    // Load the texture
    std::unique_ptr<TextureMapping2D> map;
    Float su = 1.;
    Float sv = 1.;
    Float du = 0.;
    Float dv = 0.;
    map.reset(new UVMapping2D(su, sv, du, dv));
    distributionsDict = std::make_shared<ImageTexture<RGBSpectrum, Spectrum>>(
        std::move(map), dataPath, false, 8.f, ImageWrap::Clamp, 1.f, false);

    int distResolution = distributionsDict->Width() / (N / 3);
    std::vector<Float> dist1D(distResolution);
    for (int lvl = 0; lvl < nLevels; lvl++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < distResolution; j++) {
                dist1D[j] = distributionsDict->Texel(
                    0, i / 3 * distResolution + j, lvl)[i % 3];
            }
            distributions->push_back(
                PiecewiseLinearDistribution1D(dist1D.data(), distResolution));
        }
    }

    bxdf = createGlitteringBRDF(
        arena, alphaX, alphaY, rho, alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial, st, dstdx, dstdy, nLevels, N, logMicrofacetDensity,
        0.5, distributions.get(), microfacetRelativeArea, densityRandomization,
        sampleVisibleArea, sampleApproximation);

    GlitteringConductorReflection *glitteringConductorReflection =
        dynamic_cast<GlitteringConductorReflection *>(bxdf);

    std::cout << "Number of microfacets within pixel footprint: "
              << glitteringConductorReflection->GetNMicrofacetWithinP(st, dstdx,
                                                                      dstdy)
              << std::endl;

    std::vector<Float> result(imageSize * imageSize, 0.);
    Float integral = 0;

    Float minMaxSlope = std::max(alphaX, alphaY) * 3.;

    for (int ix = 0; ix < imageSize; ++ix) {
        for (int iy = 0; iy < imageSize; ++iy) {
            Vector2f microSlope(Float(ix) / Float(imageSize),
                                Float(iy) / Float(imageSize));  // 0 - 1
            microSlope = microSlope * 2. * minMaxSlope -
                         Vector2f(minMaxSlope, minMaxSlope);

            Float SDFValue = 0.f;
            GlitteringConductorReflection *glitteringConductorReflection =
                dynamic_cast<GlitteringConductorReflection *>(bxdf);

            SDFValue = glitteringConductorReflection->P22_P(microSlope, st,
                                                            dstdx, dstdy);

            result[ix + iy * imageSize] = SDFValue;

            // integral
            integral += SDFValue * (2. * minMaxSlope / imageSize) * (2. * minMaxSlope / imageSize);
        }
    }

    std::cout << "Integral: " << integral << std::endl;

    // Extract the output file name without the extension
    std::string outputfilename(outFilename);
    std::size_t found = outputfilename.rfind(".");
    CHECK(found != std::string::npos);
    std::string outFilenameWithoutExtension =
        outputfilename.substr(std::size_t(0), found);

    // Create python file containing the plot
    std::ofstream pythonFile(outFilename);

    pythonFile << "import matplotlib.pyplot as plt" << std::endl;
    pythonFile << "import numpy as np" << std::endl;
    pythonFile << "pdf = [";
    // Write PDF
    for (int iy = 0; iy < imageSize; ++iy) {
        pythonFile << "[";
        for (int ix = 0; ix < imageSize; ++ix) {
            pythonFile << result.at(iy * imageSize + ix);
            if (ix != imageSize - 1)
                pythonFile << ",";
            else
                pythonFile << "]" << std::endl;
        }
        if (iy != imageSize - 1)
            pythonFile << ",";
        else
            pythonFile << "]" << std::endl;
    }
    pythonFile << "x = np.linspace(" << -1 << ", " << 1 << ", " << imageSize
               << ")" << std::endl;
    pythonFile << "y = np.linspace(" << -1 << ", " << 1 << ", " << imageSize
               << ")" << std::endl;
    pythonFile << "fig, axs = plt.subplots(1,1, figsize=(5, 5))" << std::endl;
    pythonFile << "pdf = np.array(pdf)" << std::endl;
    pythonFile << "vmaxval = pdf.max()" << std::endl;
    pythonFile << "pdf_plot = axs.imshow(pdf, "
                  "interpolation=\'sinc\', vmin=0, vmax=vmaxval, extent=[-"
               << 1 << ", " << 1 << ", -" << 1 << ", " << 1
               << "], origin=\'lower\')" << std::endl;
    pythonFile << "axs.set_axis_off()" << std::endl;
    // pythonFile << "axs.title.set_text(\'NDF\')" << std::endl;
    // pythonFile << "props = dict(fraction=0.046, pad=0.04)" << std::endl;
    // pythonFile << "fig.colorbar(pdf_plot, ax=axs, **props)" << std::endl;
    pythonFile << "plt.tight_layout()" << std::endl;
    pythonFile << "plt.savefig('" << outFilenameWithoutExtension
               << ".png', bbox_inches='tight')" << std::endl;
    // pythonFile << "plt.show()" << std::endl;

    pbrtCleanup();
    return 0;
}

// ============================================================================
// BEGIN Chi2Test from
// pbrt-v3/src/tests/bsdfs.cpp

/// Regularized lower incomplete gamma function (based on code from Cephes)
double RLGamma(double a, double x) {
    const double epsilon = 0.000000000000001;
    const double big = 4503599627370496.0;
    const double bigInv = 2.22044604925031308085e-16;
    if (a < 0 || x < 0)
        throw std::runtime_error("LLGamma: invalid arguments range!");

    if (x == 0) return 0.0f;

    double ax = (a * std::log(x)) - x - std::lgamma(a);
    if (ax < -709.78271289338399) return a < x ? 1.0 : 0.0;

    if (x <= 1 || x <= a) {
        double r2 = a;
        double c2 = 1;
        double ans2 = 1;

        do {
            r2 = r2 + 1;
            c2 = c2 * x / r2;
            ans2 += c2;
        } while ((c2 / ans2) > epsilon);

        return std::exp(ax) * ans2 / a;
    }

    int c = 0;
    double y = 1 - a;
    double z = x + y + 1;
    double p3 = 1;
    double q3 = x;
    double p2 = x + 1;
    double q2 = z * x;
    double ans = p2 / q2;
    double error;

    do {
        c++;
        y += 1;
        z += 2;
        double yc = y * c;
        double p = (p2 * z) - (p3 * yc);
        double q = (q2 * z) - (q3 * yc);

        if (q != 0) {
            double nextans = p / q;
            error = std::abs((ans - nextans) / nextans);
            ans = nextans;
        } else {
            // zero div, skip
            error = 1;
        }

        // shift
        p3 = p2;
        p2 = p;
        q3 = q2;
        q2 = q;

        // normalize fraction when the numerator becomes large
        if (std::abs(p) > big) {
            p3 *= bigInv;
            p2 *= bigInv;
            q3 *= bigInv;
            q2 *= bigInv;
        }
    } while (error > epsilon);

    return 1.0 - (std::exp(ax) * ans);
}

/// Chi^2 distribution cumulative distribution function
double Chi2CDF(double x, int dof) {
    if (dof < 1 || x < 0) {
        return 0.0;
    } else if (dof == 2) {
        return 1.0 - std::exp(-0.5 * x);
    } else {
        return (Float)RLGamma(0.5 * dof, 0.5 * x);
    }
}

/// Run A Chi^2 test based on the given frequency tables
std::pair<bool, std::string> Chi2Test(const Float *frequencies,
                                      const Float *expFrequencies, int res,
                                      int sampleCount, Float minExpFrequency,
                                      Float significanceLevel, int numTests) {
    struct Cell {
        Float expFrequency;
        size_t index;
    };

    /* Sort all cells by their expected frequencies */
    std::vector<Cell> cells(res);
    for (size_t i = 0; i < cells.size(); ++i) {
        cells[i].expFrequency = expFrequencies[i];
        cells[i].index = i;
    }
    std::sort(cells.begin(), cells.end(), [](const Cell &a, const Cell &b) {
        return a.expFrequency < b.expFrequency;
    });

    /* Compute the Chi^2 statistic and pool cells as necessary */
    Float pooledFrequencies = 0, pooledExpFrequencies = 0, chsq = 0;
    int pooledCells = 0, dof = 0;

    for (const Cell &c : cells) {
        if (expFrequencies[c.index] == 0) {
            if (frequencies[c.index] > sampleCount * 1e-5f) {
                /* Uh oh: samples in a c that should be completely empty
                   according to the probability density function.
                   Ordinarily, even a single sample requires immediate
                   rejection of the null hypothesis. But due to
                   finite-precision computations and rounding errors, this
                   can occasionally happen without there being an actual
                   bug. Therefore, the criterion here is a bit more
                   lenient. */

                std::string result = StringPrintf(
                    "Encountered %f samples in a c with expected "
                    "frequency 0. Rejecting the null hypothesis!",
                    frequencies[c.index]);
                return std::make_pair(false, result);
            }
        } else if (expFrequencies[c.index] < minExpFrequency) {
            /* Pool cells with low expected frequencies */
            pooledFrequencies += frequencies[c.index];
            pooledExpFrequencies += expFrequencies[c.index];
            pooledCells++;
        } else if (pooledExpFrequencies > 0 &&
                   pooledExpFrequencies < minExpFrequency) {
            /* Keep on pooling cells until a sufficiently high
               expected frequency is achieved. */
            pooledFrequencies += frequencies[c.index];
            pooledExpFrequencies += expFrequencies[c.index];
            pooledCells++;
        } else {
            Float diff = frequencies[c.index] - expFrequencies[c.index];
            chsq += (diff * diff) / expFrequencies[c.index];
            ++dof;
        }
    }

    if (pooledExpFrequencies > 0 || pooledFrequencies > 0) {
        Float diff = pooledFrequencies - pooledExpFrequencies;
        chsq += (diff * diff) / pooledExpFrequencies;
        ++dof;
    }

    /* All parameters are assumed to be known, so there is no
       additional DF reduction due to model parameters */
    dof -= 1;

    if (dof <= 0) {
        std::string result = StringPrintf(
            "The number of degrees of freedom %d is too low!", dof);
        return std::make_pair(false, result);
    }

    /* Probability of obtaining a test statistic at least
       as extreme as the one observed under the assumption
       that the distributions match */
    Float pval = 1 - (Float)Chi2CDF(chsq, dof);

    /* Apply the Sidak correction term, since we'll be conducting multiple
       independent
       hypothesis tests. This accounts for the fact that the probability of
       a failure
       increases quickly when several hypothesis tests are run in sequence.
     */
    Float alpha = 1.0f - std::pow(1.0f - significanceLevel, 1.0f / numTests);

    if (pval < alpha || !std::isfinite(pval)) {
        std::string result = StringPrintf(
            "Rejected the null hypothesis (p-value = %f, "
            "significance level = %f",
            pval, alpha);
        return std::make_pair(false, result);
    } else {
        return std::make_pair(true, std::string(""));
    }
}

// END Chi2Test from
// pbrt-v3/src/tests/bsdfs.cpp
// ============================================================================

int chisquaretestglitteringvndf(int argc, char *argv[]) {
    MemoryArena arena;
    Options opt;
    pbrtInit(opt);

    // Total number of samples to be generated.
    int nSampleHisto = 1000000;
    // Resolution of the histogram.
    int res = 512;

    Float thetaO = 1.5;
    Float phiO = 0.;

    // Number of levels in the dictionary
    int nLevels = 8;
    // Number of 1D marginal distributions in the dictionary
    int N = 768;
    // Alpha roughness used by the dictionary generator
    float alphaDict = 0.5;
    // Logarithm microfacet density
    float logMicrofacetDensity = 20.;
    // Alpha roughness of the material
    Float alphaX = 0.3;
    Float alphaY = 0.3;
    // Correlation factor
    Float rho = 0.;
    // Center of the pixel footprint
    Point2f st(0, 0);
    // Partial derivatives defining the extent of the pixel footprint
    Vector2f dstdx(0.001, 0.);
    Vector2f dstdy(0., 0.001);

    Float alphaXBaseMaterial = 0.1;
    Float alphaYBaseMaterial = 0.1;
    Float rhoBaseMaterial = 0.;
    Float microfacetRelativeArea = 1.;
    Float densityRandomization = 0.01;
    bool sampleVisibleArea = true;
    bool sampleApproximation = false;

    int i;
    auto parseArg = [&]() -> std::pair<std::string, double> {
        const char *ptr = argv[i];
        // Skip over a leading dash or two.
        CHECK_EQ(*ptr, '-');
        ++ptr;
        if (*ptr == '-') ++ptr;

        // Copy the flag name to the string.
        std::string flag;
        while (*ptr && *ptr != '=') flag += *ptr++;
        if (!*ptr && i + 1 == argc)
            usage("missing value after %s flag", argv[i]);
        const char *value = (*ptr == '=') ? (ptr + 1) : argv[++i];
        return {flag, atof(value)};
    };

    std::pair<std::string, double> arg;
    for (i = 0; i < argc; ++i) {
        if (argv[i][0] != '-')
            break;

        else {
            std::pair<std::string, double> arg = parseArg();
            if (std::get<0>(arg) == "nsamplehisto") {
                nSampleHisto = int(std::get<1>(arg));
            } else if (std::get<0>(arg) == "res") {
                res = std::get<1>(arg);
            } else if (std::get<0>(arg) == "alphax") {
                alphaX = std::get<1>(arg);
            } else if (std::get<0>(arg) == "alphay") {
                alphaY = std::get<1>(arg);
            } else if (std::get<0>(arg) == "rho") {
                rho = std::get<1>(arg);
            } else if (std::get<0>(arg) == "stx") {
                st.x = std::get<1>(arg);
            } else if (std::get<0>(arg) == "dstdxx") {
                dstdx.x = std::get<1>(arg);
            } else if (std::get<0>(arg) == "dstdyy") {
                dstdy.y = std::get<1>(arg);
            } else if (std::get<0>(arg) == "thetao") {
                thetaO = std::get<1>(arg);
            } else if (std::get<0>(arg) == "phio") {
                phiO = std::get<1>(arg);
            } else
                usage();
        }
    }

    if (i + 1 >= argc)
        usage("missing second filename for \"chisquaretestglitteringvndf\"");
    else if (i >= argc)
        usage("missing filenames for \"chisquaretestglitteringvndf\"");

    const char *dictFilename = argv[i], *outFilename = argv[i + 1];
    std::string dataPath(dictFilename);
    Fresnel *fresnel = new FresnelNoOp();

    RNG rngSample;

    // Contain the piecewise linear distributions
    // Load the texture containing the dictionary
    std::unique_ptr<TextureMapping2D> map;
    Float su = 1.;
    Float sv = 1.;
    Float du = 0.;
    Float dv = 0.;
    map.reset(new UVMapping2D(su, sv, du, dv));
    ImageTexture<RGBSpectrum, Spectrum> distributionsDict(
        std::move(map), dictFilename, false, 8.f, ImageWrap::Clamp, 1.f, false);

    // Will contain the N distributions
    std::vector<PiecewiseLinearDistribution1D> distributions;
    int distResolution = distributionsDict.Width() / (N / 3);

    // Populate the 1D distributions by using the dictionary stored in a
    // texture
    std::vector<Float> dist1D(distResolution);
    for (int lvl = 0; lvl < nLevels; lvl++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < distResolution; j++) {
                dist1D[j] = distributionsDict.Texel(
                    0, i / 3 * distResolution + j, lvl)[i % 3];
            }
            distributions.push_back(
                PiecewiseLinearDistribution1D(dist1D.data(), distResolution));
        }
    }
    BxDF *bxdf = createGlitteringBRDF(
        arena, alphaX, alphaY, rho, alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial, st, dstdx, dstdy, nLevels, N, logMicrofacetDensity,
        0.5, &distributions, microfacetRelativeArea, densityRandomization,
        sampleVisibleArea, sampleApproximation);
    GlitteringConductorReflection *sparklingReflection =
        dynamic_cast<GlitteringConductorReflection *>(bxdf);

    std::cout << "Number of microfacets within pixel footprint: "
              << sparklingReflection->GetNMicrofacetWithinP(st, dstdx, dstdy)
              << std::endl;

    Float xYMax = 1.;

    Vector3f wo = SphericalDirection(std::sin(thetaO), std::cos(thetaO), phiO);

    std::vector<Float> pdf(res * res);
    Float integralPDF = 0;
    for (int ixy = 0; ixy < res * res; ++ixy) {
        Float x = ((Mod(ixy, res) + 0.5) / Float(res) * 2. - 1.) * xYMax;
        Float y = (((ixy / res) + 0.5) / Float(res) * 2. - 1.) * xYMax;

        Vector2f projectedNormal(x, y);
        float normalZSqr = 1. - projectedNormal.x * projectedNormal.x -
                           projectedNormal.y * projectedNormal.y;
        if (normalZSqr <= 0.0001) continue;
        Vector3f microNormal(projectedNormal.x, projectedNormal.y,
                             std::sqrt(normalZSqr));

        Float density = sparklingReflection->Visible_D_P(microNormal, wo);
        // 4. * xYMax * xYMax / (res * res): sub interval area
        // 1. / microNormal.z: Jacobian of the projected normal domain to
        // normal domain
        pdf.at(ixy) = density * 4. * xYMax * xYMax / (res * res) *
                      nSampleHisto / microNormal.z;

        // integral
        // integralPDF += pdf.at(ixy) * 4. * xYMax * xYMax / (res * res);
        integralPDF += pdf.at(ixy) / nSampleHisto;
    }
    std::cout << "Integral PDF: " << integralPDF << std::endl;

    // Extract the output file name without the extension
    std::string outputfilename(outFilename);
    std::size_t found = outputfilename.rfind(".");
    CHECK(found != std::string::npos);
    std::string outFilenameWithoutExtension =
        outputfilename.substr(std::size_t(0), found);
    std::string extension = outputfilename.substr(found, outputfilename.size());

    RNG rngSampling;

    for (int k = 0; k < CHI2_RUNS; ++k) {
        // 1D histogram
        std::vector<Float> histogram(res * res);

        // Accumulate samples in the histogram
        for (int i = 0; i < nSampleHisto; ++i) {
            Vector3f sampleValue = sparklingReflection->Sample_Visible_D_P(wo);

            // [-1, 1]
            Vector2f sampleValueNormalized(sampleValue.x, sampleValue.y);
            // [0, 1]
            sampleValueNormalized =
                (sampleValueNormalized + Vector2f(1., 1.)) / 2.;
            Vector2i histogramIndexSample(sampleValueNormalized.x * res,
                                          sampleValueNormalized.y * res);
            if (histogramIndexSample.x < 0 || histogramIndexSample.y < 0 ||
                histogramIndexSample.x >= res || histogramIndexSample.y >= res)
                continue;

            histogram.at(histogramIndexSample.y * res +
                         histogramIndexSample.x) += 1.;
        }

        Float integralHistogram = 0.;
        for (int ixy = 0; ixy < res * res; ++ixy) {
            Float density = histogram.at(ixy);
            integralHistogram += density * 4. * xYMax * xYMax / (res * res);
        }
        // std::cout << "Integral histogram: " << integralHistogram <<
        // std::endl;

        auto result =
            Chi2Test(histogram.data(), pdf.data(), res * res, nSampleHisto,
                     CHI2_MINFREQ, CHI2_SLEVEL, CHI2_RUNS);

        // Create python file containing the plot
        std::ofstream pythonFile(outFilenameWithoutExtension + "_test" +
                                 std::to_string(k + 1) + extension);
        pythonFile << "import matplotlib.pyplot as plt" << std::endl;
        pythonFile << "import numpy as np" << std::endl;
        pythonFile << "histogram = [";
        // Write histogram
        for (int iy = 0; iy < res; ++iy) {
            pythonFile << "[";
            for (int ix = 0; ix < res; ++ix) {
                pythonFile << histogram.at(iy * res + ix);
                if (ix != res - 1)
                    pythonFile << ",";
                else
                    pythonFile << "]" << std::endl;
            }
            if (iy != res - 1)
                pythonFile << ",";
            else
                pythonFile << "]" << std::endl;
        }
        pythonFile << "pdf = [";
        // Write PDF
        for (int iy = 0; iy < res; ++iy) {
            pythonFile << "[";
            for (int ix = 0; ix < res; ++ix) {
                pythonFile << pdf.at(iy * res + ix);
                if (ix != res - 1)
                    pythonFile << ",";
                else
                    pythonFile << "]" << std::endl;
            }
            if (iy != res - 1)
                pythonFile << ",";
            else
                pythonFile << "]" << std::endl;
        }
        pythonFile << "x = np.linspace(" << -xYMax << ", " << xYMax << ", "
                   << res << ")" << std::endl;
        pythonFile << "y = np.linspace(" << -xYMax << ", " << xYMax << ", "
                   << res << ")" << std::endl;
        pythonFile << "fig, axs = plt.subplots(1,3, figsize=(15, 5))"
                   << std::endl;
        pythonFile << "pdf = np.array(pdf)" << std::endl;
        pythonFile << "histogram = np.array(histogram)" << std::endl;
        pythonFile << "diff=histogram - pdf" << std::endl;
        pythonFile << "absdiff=np.abs(diff)" << std::endl;
        pythonFile << "a = pdf.shape[1] / pdf.shape[0]" << std::endl;
        pythonFile << "vmaxval = pdf.max()" << std::endl;
        pythonFile << "vmaxval = max(vmaxval, histogram.max())" << std::endl;
        pythonFile
            << "pdf_plot = axs[0].imshow(pdf, aspect=a, "
               "interpolation=\'nearest\', vmin=0, vmax=vmaxval, extent=[-"
            << xYMax << ", " << xYMax << ", -" << xYMax << ", " << xYMax
            << "], origin=\'lower\')" << std::endl;
        pythonFile
            << "hist_plot = axs[1].imshow(histogram, aspect=a, "
               "interpolation=\'nearest\', vmin=0, vmax=vmaxval, extent=[-"
            << xYMax << ", " << xYMax << ", -" << xYMax << ", " << xYMax
            << "], origin=\'lower\')" << std::endl;
        pythonFile << "diff_plot = axs[2].imshow(absdiff, aspect=a, vmin=0, "
                      "vmax=vmaxval, interpolation=\'nearest\', extent=[-"
                   << xYMax << ", " << xYMax << ", -" << xYMax << ", " << xYMax
                   << "], origin=\'lower\')" << std::endl;
        pythonFile << "axs[0].title.set_text(\'PDF\')" << std::endl;
        pythonFile << "axs[1].title.set_text(\'Histogram\')" << std::endl;
        pythonFile << "axs[2].title.set_text(\'Absolute difference\')"
                   << std::endl;
        pythonFile << "props = dict(fraction=0.046, pad=0.04)" << std::endl;
        pythonFile << "fig.colorbar(pdf_plot, ax=axs[0], **props)" << std::endl;
        pythonFile << "fig.colorbar(hist_plot, ax=axs[1], **props)"
                   << std::endl;
        pythonFile << "fig.colorbar(diff_plot, ax=axs[2], **props)"
                   << std::endl;
        pythonFile << "plt.tight_layout()" << std::endl;
        pythonFile << "plt.savefig(\'"
                   << outFilenameWithoutExtension + "_test" +
                          std::to_string(k + 1)
                   << ".png\')" << std::endl;
        pythonFile << "plt.show()" << std::endl;

        if (result.first)
            std::cout << "Test " << k + 1 << ": success" << std::endl;
        else {
            std::cout << "Test " << k + 1 << ": failure. " << result.first
                      << ". " << result.second << std::endl;
            ParallelCleanup();
            return 1;
        }
    }

    pbrtCleanup();
    return 0;
}

int chisquaretestglitteringbrdf(int argc, char *argv[]) {
    MemoryArena arena;
    Options opt;
    pbrtInit(opt);

    // Total number of samples to be generated.
    int nSampleHisto = 1000000;
    // Resolution of the histogram.
    int res = 512;

    Float thetaO = 0.2;
    Float phiO = 0.;

    // Number of levels in the dictionary
    int nLevels = 8;
    // Number of 1D marginal distributions in the dictionary
    int N = 768;
    // Alpha roughness used by the dictionary generator
    float alphaDict = 0.5;
    // Logarithm microfacet density
    float logMicrofacetDensity = 20.;
    // Microfacet relative area
    float microfacetRelativeArea = 1.;
    // Alpha roughness of the material
    Float alphaX = 0.3;
    Float alphaY = 0.3;
    // Correlation factor
    Float rho = 0.;
    // Center of the pixel footprint
    Point2f st(0, 0);
    // Partial derivatives defining the extent of the pixel footprint
    Vector2f dstdx(0.001, 0.);
    Vector2f dstdy(0., 0.001);
    // Vector2f dstdx(0.1, 0.);
    // Vector2f dstdy(0., 0.1);

    Float alphaXBaseMaterial = 0.1;
    Float alphaYBaseMaterial = 0.1;
    Float rhoBaseMaterial = 0.;
    Float densityRandomization = 0.01;
    bool sampleVisibleArea = true;
    bool sampleApproximation = false;

    int i;
    auto parseArg = [&]() -> std::pair<std::string, double> {
        const char *ptr = argv[i];
        // Skip over a leading dash or two.
        CHECK_EQ(*ptr, '-');
        ++ptr;
        if (*ptr == '-') ++ptr;

        // Copy the flag name to the string.
        std::string flag;
        while (*ptr && *ptr != '=') flag += *ptr++;
        if (!*ptr && i + 1 == argc)
            usage("missing value after %s flag", argv[i]);
        const char *value = (*ptr == '=') ? (ptr + 1) : argv[++i];
        return {flag, atof(value)};
    };

    std::pair<std::string, double> arg;
    for (i = 0; i < argc; ++i) {
        if (argv[i][0] != '-')
            break;

        else {
            std::pair<std::string, double> arg = parseArg();
            if (std::get<0>(arg) == "nsamplehisto") {
                nSampleHisto = int(std::get<1>(arg));
            } else if (std::get<0>(arg) == "res") {
                res = std::get<1>(arg);
            } else if (std::get<0>(arg) == "alphax") {
                alphaX = std::get<1>(arg);
            } else if (std::get<0>(arg) == "alphay") {
                alphaY = std::get<1>(arg);
            } else if (std::get<0>(arg) == "rho") {
                rho = std::get<1>(arg);
            } else if (std::get<0>(arg) == "stx") {
                st.x = std::get<1>(arg);
            } else if (std::get<0>(arg) == "dstdxx") {
                dstdx.x = std::get<1>(arg);
            } else if (std::get<0>(arg) == "dstdyy") {
                dstdy.y = std::get<1>(arg);
            } else if (std::get<0>(arg) == "thetao") {
                thetaO = std::get<1>(arg);
            } else if (std::get<0>(arg) == "phio") {
                phiO = std::get<1>(arg);
            } else if (std::get<0>(arg) == "mra") {
                microfacetRelativeArea = std::get<1>(arg);
            } else
                usage();
        }
    }

    if (i + 1 >= argc)
        usage("missing second filename for \"chisquaretestglitteringbrdf\"");
    else if (i >= argc)
        usage("missing filenames for \"chisquaretestglitteringbrdf\"");

    const char *dictFilename = argv[i], *outFilename = argv[i + 1];
    std::string dataPath(dictFilename);
    Fresnel *fresnel = new FresnelNoOp();

    // Contain the piecewise linear distribution
    // Load the texture containing the dictionary
    std::unique_ptr<TextureMapping2D> map;
    Float su = 1.;
    Float sv = 1.;
    Float du = 0.;
    Float dv = 0.;
    map.reset(new UVMapping2D(su, sv, du, dv));
    ImageTexture<RGBSpectrum, Spectrum> distributionsDict(
        std::move(map), dictFilename, false, 8.f, ImageWrap::Clamp, 1.f, false);

    // Will contain the N distributions
    std::vector<PiecewiseLinearDistribution1D> distributions;
    int distResolution = distributionsDict.Width() / (N / 3);

    // Populate the 1D distributions by using the dictionary stored in a
    // texture
    std::vector<Float> dist1D(distResolution);
    for (int lvl = 0; lvl < nLevels; lvl++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < distResolution; j++) {
                dist1D[j] = distributionsDict.Texel(
                    0, i / 3 * distResolution + j, lvl)[i % 3];
            }
            distributions.push_back(
                PiecewiseLinearDistribution1D(dist1D.data(), distResolution));
        }
    }
    BxDF *bxdf = createGlitteringBRDF(
        arena, alphaX, alphaY, rho, alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial, st, dstdx, dstdy, nLevels, N, logMicrofacetDensity,
        0.5, &distributions, microfacetRelativeArea, densityRandomization,
        sampleVisibleArea, sampleApproximation);

    GlitteringConductorReflection *sparklingReflection =
        dynamic_cast<GlitteringConductorReflection *>(bxdf);

    std::cout << "Number of microfacets within the pixel footprint: "
              << sparklingReflection->GetNMicrofacetWithinP(st, dstdx, dstdy)
              << std::endl;

    Float xYMax = 1.;

    Vector3f wo = SphericalDirection(std::sin(thetaO), std::cos(thetaO), phiO);

    std::vector<Float> pdf(res * res);
    Float integralPDF = 0;
    for (int ixy = 0; ixy < res * res; ++ixy) {
        Float x = ((Mod(ixy, res) + 0.5) / Float(res) * 2. - 1.) * xYMax;
        Float y = (((ixy / res) + 0.5) / Float(res) * 2. - 1.) * xYMax;

        Vector2f projectedDirection(x, y);
        float normalZSqr = 1. - projectedDirection.x * projectedDirection.x -
                           projectedDirection.y * projectedDirection.y;
        if (normalZSqr <= 0.) continue;
        Vector3f wi(projectedDirection.x, projectedDirection.y,
                    std::sqrt(normalZSqr));

        Float density = sparklingReflection->Pdf(wo, wi);

        CHECK_GE(density, 0.);
        // 4. * xYMax * xYMax / (res * res): sub interval area
        // 1. / microNormal.z: Jacobian of the projected normal domain to
        // normal domain
        pdf.at(ixy) =
            density * 4. * xYMax * xYMax / (res * res) * nSampleHisto / wi.z;

        // integral
        // integralPDF += pdf.at(ixy) * 4. * xYMax * xYMax / (res * res);
        integralPDF += pdf.at(ixy) / nSampleHisto;
    }
    std::cout << "Integral PDF: " << integralPDF << std::endl;

    // Extract the output file name without the extension
    std::string outputfilename(outFilename);
    std::size_t found = outputfilename.rfind(".");
    CHECK(found != std::string::npos);
    std::string outFilenameWithoutExtension =
        outputfilename.substr(std::size_t(0), found);
    std::string extension = outputfilename.substr(found, outputfilename.size());

    RNG rngSampling;

    for (int k = 0; k < CHI2_RUNS; ++k) {
        // 1D histogram
        std::vector<Float> histogram(res * res);

        // Accumulate samples in the histogram
        for (int i = 0; i < nSampleHisto; ++i) {
            Vector3f wi;
            Point2f u(rngSampling.UniformFloat(), rngSampling.UniformFloat());
            Float pdf = 0.;
            BxDFType flags;
            Spectrum BRDFValue =
                sparklingReflection->Sample_f(wo, &wi, u, &pdf, &flags);
            if (wi.z < 0.001) continue;

            // [-1, 1]
            Vector2f sampleValueNormalized(wi.x, wi.y);
            // [0, 1]
            sampleValueNormalized =
                (sampleValueNormalized + Vector2f(1., 1.)) / 2.;
            Vector2i histogramIndexSample(sampleValueNormalized.x * res,
                                          sampleValueNormalized.y * res);
            if (histogramIndexSample.x < 0 || histogramIndexSample.y < 0 ||
                histogramIndexSample.x >= res || histogramIndexSample.y >= res)
                continue;

            histogram.at(histogramIndexSample.y * res +
                         histogramIndexSample.x) += 1.;
        }

        Float integralHistogram = 0.;
        for (int ixy = 0; ixy < res * res; ++ixy) {
            Float density = histogram.at(ixy);
            integralHistogram += density * 4. * xYMax * xYMax / (res * res);
        }
        // std::cout << "Integral histogram: " << integralHistogram <<
        // std::endl;

        auto result =
            Chi2Test(histogram.data(), pdf.data(), res * res, nSampleHisto,
                     CHI2_MINFREQ, CHI2_SLEVEL, CHI2_RUNS);

        // Create python file containing the plot
        std::ofstream pythonFile(outFilenameWithoutExtension + "_test" +
                                 std::to_string(k + 1) + extension);
        pythonFile << "import matplotlib.pyplot as plt" << std::endl;
        pythonFile << "import numpy as np" << std::endl;
        pythonFile << "histogram = [";
        // Write histogram
        for (int iy = 0; iy < res; ++iy) {
            pythonFile << "[";
            for (int ix = 0; ix < res; ++ix) {
                pythonFile << histogram.at(iy * res + ix);
                if (ix != res - 1)
                    pythonFile << ",";
                else
                    pythonFile << "]" << std::endl;
            }
            if (iy != res - 1)
                pythonFile << ",";
            else
                pythonFile << "]" << std::endl;
        }
        pythonFile << "pdf = [";
        // Write PDF
        for (int iy = 0; iy < res; ++iy) {
            pythonFile << "[";
            for (int ix = 0; ix < res; ++ix) {
                pythonFile << pdf.at(iy * res + ix);
                if (ix != res - 1)
                    pythonFile << ",";
                else
                    pythonFile << "]" << std::endl;
            }
            if (iy != res - 1)
                pythonFile << ",";
            else
                pythonFile << "]" << std::endl;
        }
        pythonFile << "x = np.linspace(" << -xYMax << ", " << xYMax << ", "
                   << res << ")" << std::endl;
        pythonFile << "y = np.linspace(" << -xYMax << ", " << xYMax << ", "
                   << res << ")" << std::endl;
        pythonFile << "fig, axs = plt.subplots(1,3, figsize=(15, 5))"
                   << std::endl;
        pythonFile << "pdf = np.array(pdf)" << std::endl;
        pythonFile << "histogram = np.array(histogram)" << std::endl;
        pythonFile << "diff=histogram - pdf" << std::endl;
        pythonFile << "absdiff=np.abs(diff)" << std::endl;
        pythonFile << "a = pdf.shape[1] / pdf.shape[0]" << std::endl;
        pythonFile << "vmaxval = pdf.max()" << std::endl;
        pythonFile << "vmaxval = max(vmaxval, histogram.max())" << std::endl;
        pythonFile
            << "pdf_plot = axs[0].imshow(pdf, aspect=a, "
               "interpolation=\'nearest\', vmin=0, vmax=vmaxval, extent=[-"
            << xYMax << ", " << xYMax << ", -" << xYMax << ", " << xYMax
            << "], origin=\'lower\')" << std::endl;
        pythonFile
            << "hist_plot = axs[1].imshow(histogram, aspect=a, "
               "interpolation=\'nearest\', vmin=0, vmax=vmaxval, extent=[-"
            << xYMax << ", " << xYMax << ", -" << xYMax << ", " << xYMax
            << "], origin=\'lower\')" << std::endl;
        pythonFile << "diff_plot = axs[2].imshow(absdiff, aspect=a, vmin=0, "
                      "vmax=vmaxval, interpolation=\'nearest\', extent=[-"
                   << xYMax << ", " << xYMax << ", -" << xYMax << ", " << xYMax
                   << "], origin=\'lower\')" << std::endl;
        pythonFile << "axs[0].title.set_text(\'PDF\')" << std::endl;
        pythonFile << "axs[1].title.set_text(\'Histogram\')" << std::endl;
        pythonFile << "axs[2].title.set_text(\'Absolute difference\')"
                   << std::endl;
        pythonFile << "props = dict(fraction=0.046, pad=0.04)" << std::endl;
        pythonFile << "fig.colorbar(pdf_plot, ax=axs[0], **props)" << std::endl;
        pythonFile << "fig.colorbar(hist_plot, ax=axs[1], **props)"
                   << std::endl;
        pythonFile << "fig.colorbar(diff_plot, ax=axs[2], **props)"
                   << std::endl;
        pythonFile << "plt.tight_layout()" << std::endl;
        pythonFile << "plt.savefig(\'"
                   << outFilenameWithoutExtension + "_test" +
                          std::to_string(k + 1)
                   << ".png\')" << std::endl;
        pythonFile << "plt.show()" << std::endl;

        if (result.first)
            std::cout << "Test " << k + 1 << ": success" << std::endl;
        else {
            std::cout << "Test " << k + 1 << ": failure. " << result.first
                      << ". " << result.second << std::endl;
            ParallelCleanup();
            return 1;
        }
    }

    pbrtCleanup();
    return 0;
}

int chisquaretestglitteringbsdf(int argc, char *argv[]) {
    MemoryArena arena;
    Options opt;
    pbrtInit(opt);

    // Total number of samples to be generated.
    int nSampleHisto = 256000000;
    // Resolution of the histogram.
    int res = 4096;

    Float thetaO = 0.2;
    Float phiO = 0.;

    // Number of levels in the dictionary
    int nLevels = 8;
    // Number of 1D marginal distributions in the dictionary
    int N = 768;
    // Alpha roughness used by the dictionary generator
    float alphaDict = 0.5;
    // Logarithm microfacet density
    float logMicrofacetDensity = 20.;
    // Microfacet relative area
    float microfacetRelativeArea = 1.;
    // Alpha roughness of the material
    Float alphaX = 0.25;
    Float alphaY = 0.25;
    // Correlation factor
    Float rho = 0.;
    // Center of the pixel footprint
    Point2f st(0, 0);
    // Partial derivatives defining the extent of the pixel footprint
    Vector2f dstdx(0.001, 0.);
    Vector2f dstdy(0., 0.001);

    Float alphaXBaseMaterial = 0.1;
    Float alphaYBaseMaterial = 0.1;
    Float rhoBaseMaterial = 0.;
    Float densityRandomization = 0.01;
    bool sampleVisibleArea = true;
    bool sampleApproximation = false;

    int i;
    auto parseArg = [&]() -> std::pair<std::string, double> {
        const char *ptr = argv[i];
        // Skip over a leading dash or two.
        CHECK_EQ(*ptr, '-');
        ++ptr;
        if (*ptr == '-') ++ptr;

        // Copy the flag name to the string.
        std::string flag;
        while (*ptr && *ptr != '=') flag += *ptr++;
        if (!*ptr && i + 1 == argc)
            usage("missing value after %s flag", argv[i]);
        const char *value = (*ptr == '=') ? (ptr + 1) : argv[++i];
        return {flag, atof(value)};
    };

    std::pair<std::string, double> arg;
    for (i = 0; i < argc; ++i) {
        if (argv[i][0] != '-')
            break;

        else {
            std::pair<std::string, double> arg = parseArg();
            if (std::get<0>(arg) == "nsamplehisto") {
                nSampleHisto = int(std::get<1>(arg));
            } else if (std::get<0>(arg) == "res") {
                res = std::get<1>(arg);
            } else if (std::get<0>(arg) == "alphax") {
                alphaX = std::get<1>(arg);
            } else if (std::get<0>(arg) == "alphay") {
                alphaY = std::get<1>(arg);
            } else if (std::get<0>(arg) == "rho") {
                rho = std::get<1>(arg);
            } else if (std::get<0>(arg) == "thetao") {
                thetaO = std::get<1>(arg);
            } else if (std::get<0>(arg) == "phio") {
                phiO = std::get<1>(arg);
            } else if (std::get<0>(arg) == "stx") {
                st.x = std::get<1>(arg);
            } else if (std::get<0>(arg) == "dstdxx") {
                dstdx.x = std::get<1>(arg);
            } else if (std::get<0>(arg) == "dstdyy") {
                dstdy.y = std::get<1>(arg);
            } else if (std::get<0>(arg) == "mra") {
                microfacetRelativeArea = std::get<1>(arg);
            } else
                usage();
        }
    }

    if (i + 1 >= argc)
        usage("missing second filename");
    else if (i >= argc)
        usage("missing filenames");

    const char *dictFilename = argv[i], *outFilename = argv[i + 1];
    std::string dataPath(dictFilename);

    // Contain the piecewise linear distribution
    // Load the texture containing the dictionary
    std::unique_ptr<TextureMapping2D> map;
    Float su = 1.;
    Float sv = 1.;
    Float du = 0.;
    Float dv = 0.;
    map.reset(new UVMapping2D(su, sv, du, dv));
    ImageTexture<RGBSpectrum, Spectrum> distributionsDict(
        std::move(map), dictFilename, false, 8.f, ImageWrap::Clamp, 1.f, false);

    // Will contain the N distributions
    std::vector<PiecewiseLinearDistribution1D> distributions;
    int distResolution = distributionsDict.Width() / (N / 3);

    // Populate the 1D distributions by using the dictionary stored in a
    // texture
    std::vector<Float> dist1D(distResolution);
    for (int lvl = 0; lvl < nLevels; lvl++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < distResolution; j++) {
                dist1D[j] = distributionsDict.Texel(
                    0, i / 3 * distResolution + j, lvl)[i % 3];
            }
            distributions.push_back(
                PiecewiseLinearDistribution1D(dist1D.data(), distResolution));
        }
    }

    BxDF *brdf = createGlitteringBRDF(
        arena, alphaX, alphaY, rho, alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial, st, dstdx, dstdy, nLevels, N, logMicrofacetDensity,
        0.5, &distributions, microfacetRelativeArea, densityRandomization,
        sampleVisibleArea, sampleApproximation);

    BxDF *bsdf = createGlitteringBSDF(
        arena, alphaX, alphaY, rho, alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial, st, dstdx, dstdy, nLevels, N, logMicrofacetDensity,
        0.5, &distributions, microfacetRelativeArea, densityRandomization,
        sampleVisibleArea, sampleApproximation);

    RNG rngSample;

    GlitteringDielectricScattering *sparklingScattering =
        dynamic_cast<GlitteringDielectricScattering *>(bsdf);
    GlitteringConductorReflection *sparklingReflection =
        dynamic_cast<GlitteringConductorReflection *>(brdf);

    std::cout << "Number of microfacets within the pixel footprint: "
              << sparklingReflection->GetNMicrofacetWithinP(st, dstdx, dstdy)
              << std::endl;

    Float xYMax = 1.;

    Vector3f wo = SphericalDirection(std::sin(thetaO), std::cos(thetaO), phiO);

    std::vector<Float> pdfZPos(res * res);
    std::vector<Float> pdfZNeg(res * res);
    Float integralPDFZPos = 0;
    Float integralPDFZNeg = 0;
    for (int ixy = 0; ixy < res * res; ++ixy) {
        Float x = ((Mod(ixy, res) + 0.5) / Float(res) * 2. - 1.) * xYMax;
        Float y = (((ixy / res) + 0.5) / Float(res) * 2. - 1.) * xYMax;

        Vector2f projectedDirection(x, y);
        float normalZSqr = 1. - projectedDirection.x * projectedDirection.x -
                           projectedDirection.y * projectedDirection.y;
        if (normalZSqr <= 0.) continue;
        Vector3f wi(projectedDirection.x, projectedDirection.y,
                    std::sqrt(normalZSqr));

        Float density = sparklingScattering->Pdf(wo, wi);

        CHECK_GE(density, 0.);

        pdfZPos.at(ixy) =
            density * 4. * xYMax * xYMax / (res * res) * nSampleHisto / wi.z;

        // integralPDFZPos += pdfZPos.at(ixy) * 4. * xYMax * xYMax / (res *
        // res);
        integralPDFZPos += pdfZPos.at(ixy) / nSampleHisto;
    }
    std::cout << "Integral PDF (upper hemisphere): " << integralPDFZPos
              << std::endl;
    for (int ixy = 0; ixy < res * res; ++ixy) {
        Float x = ((Mod(ixy, res) + 0.5) / Float(res) * 2. - 1.) * xYMax;
        Float y = (((ixy / res) + 0.5) / Float(res) * 2. - 1.) * xYMax;

        Vector2f projectedDirection(x, y);
        float normalZSqr = 1. - projectedDirection.x * projectedDirection.x -
                           projectedDirection.y * projectedDirection.y;
        if (normalZSqr <= 0.) continue;
        Vector3f wi(projectedDirection.x, projectedDirection.y,
                    -std::sqrt(normalZSqr));

        Float density = sparklingScattering->Pdf(wo, wi);

        CHECK_GE(density, 0.);

        pdfZNeg.at(ixy) = density * 4. * xYMax * xYMax / (res * res) *
                          nSampleHisto / std::abs(wi.z);

        // integralPDFZNeg += pdfZNeg.at(ixy) * 4. * xYMax * xYMax / (res *
        // res);
        integralPDFZNeg += pdfZNeg.at(ixy) / nSampleHisto;
    }
    std::cout << "Integral PDF (lower hemisphere): " << integralPDFZNeg
              << std::endl;

    // Extract the output file name without the extension
    std::string outputfilename(outFilename);
    std::size_t found = outputfilename.rfind(".");
    CHECK(found != std::string::npos);
    std::string outFilenameWithoutExtension =
        outputfilename.substr(std::size_t(0), found);
    std::string extension = outputfilename.substr(found, outputfilename.size());

    RNG rngSampling;

    for (int k = 0; k < CHI2_RUNS; ++k) {
        std::vector<Float> histogramZPos(res * res);
        std::vector<Float> histogramZNeg(res * res);
        int nSampleHistoZPos = 0;
        int nSampleHistoZNeg = 0;

        // Accumulate samples in the histograms
        for (int i = 0; i < nSampleHisto; ++i) {
            Vector3f wi;
            Point2f u(rngSampling.UniformFloat(), rngSampling.UniformFloat());
            Float pdf = 0.;
            BxDFType flags;
            Spectrum BRDFValue =
                sparklingScattering->Sample_f(wo, &wi, u, &pdf, &flags);
            if (BRDFValue.IsBlack() || wi.z == 0.) continue;

            // [-1, 1]
            Vector2f sampleValueNormalized(wi.x, wi.y);
            // [0, 1]
            sampleValueNormalized =
                (sampleValueNormalized + Vector2f(1., 1.)) / 2.;
            Vector2i histogramIndexSample(sampleValueNormalized.x * res,
                                          sampleValueNormalized.y * res);
            if (histogramIndexSample.x < 0 || histogramIndexSample.y < 0 ||
                histogramIndexSample.x >= res || histogramIndexSample.y >= res)
                continue;
            if (wi.z > 0) {
                histogramZPos.at(histogramIndexSample.y * res +
                                 histogramIndexSample.x) += 1.;
                nSampleHistoZPos++;
            } else {
                histogramZNeg.at(histogramIndexSample.y * res +
                                 histogramIndexSample.x) += 1.;
                nSampleHistoZNeg++;
            }
        }

        Float integralHistogramZPos = 0.;
        Float integralHistogramZNeg = 0.;
        for (int ixy = 0; ixy < res * res; ++ixy) {
            Float density = histogramZPos.at(ixy);
            integralHistogramZPos += density * 4. * xYMax * xYMax / (res * res);
        }
        for (int ixy = 0; ixy < res * res; ++ixy) {
            Float density = histogramZNeg.at(ixy);
            integralHistogramZNeg += density * 4. * xYMax * xYMax / (res * res);
        }
        /*std::cout << "Integral histogram positive Z: " <<
        integralHistogramZPos
                  << std::endl;
        std::cout << "Integral histogram negative Z: " << integralHistogramZNeg
                  << std::endl;*/

        auto resultZPos =
            Chi2Test(histogramZPos.data(), pdfZPos.data(), res * res,
                     nSampleHistoZPos, CHI2_MINFREQ, CHI2_SLEVEL, CHI2_RUNS);

        auto resultZNeg =
            Chi2Test(histogramZNeg.data(), pdfZNeg.data(), res * res,
                     nSampleHistoZNeg, CHI2_MINFREQ, CHI2_SLEVEL, CHI2_RUNS);

        // Create python file containing the plot
        std::ofstream pythonFile(outFilenameWithoutExtension + "_test" +
                                 std::to_string(k + 1) + extension);
        pythonFile << "import matplotlib.pyplot as plt" << std::endl;
        pythonFile << "import numpy as np" << std::endl;
        pythonFile << "histogramZPos = [";
        // Write histogram positive Z
        for (int iy = 0; iy < res; ++iy) {
            pythonFile << "[";
            for (int ix = 0; ix < res; ++ix) {
                pythonFile << histogramZPos.at(iy * res + ix);
                if (ix != res - 1)
                    pythonFile << ",";
                else
                    pythonFile << "]" << std::endl;
            }
            if (iy != res - 1)
                pythonFile << ",";
            else
                pythonFile << "]" << std::endl;
        }
        pythonFile << "pdfZPos = [";
        // Write PDF positive Z
        for (int iy = 0; iy < res; ++iy) {
            pythonFile << "[";
            for (int ix = 0; ix < res; ++ix) {
                pythonFile << pdfZPos.at(iy * res + ix);
                if (ix != res - 1)
                    pythonFile << ",";
                else
                    pythonFile << "]" << std::endl;
            }
            if (iy != res - 1)
                pythonFile << ",";
            else
                pythonFile << "]" << std::endl;
        }
        pythonFile << "histogramZNeg = [";
        // Write histogram positive Z
        for (int iy = 0; iy < res; ++iy) {
            pythonFile << "[";
            for (int ix = 0; ix < res; ++ix) {
                pythonFile << histogramZNeg.at(iy * res + ix);
                if (ix != res - 1)
                    pythonFile << ",";
                else
                    pythonFile << "]" << std::endl;
            }
            if (iy != res - 1)
                pythonFile << ",";
            else
                pythonFile << "]" << std::endl;
        }
        pythonFile << "pdfZNeg = [";
        // Write PDF positive Z
        for (int iy = 0; iy < res; ++iy) {
            pythonFile << "[";
            for (int ix = 0; ix < res; ++ix) {
                pythonFile << pdfZNeg.at(iy * res + ix);
                if (ix != res - 1)
                    pythonFile << ",";
                else
                    pythonFile << "]" << std::endl;
            }
            if (iy != res - 1)
                pythonFile << ",";
            else
                pythonFile << "]" << std::endl;
        }
        pythonFile << "x = np.linspace(" << -xYMax << ", " << xYMax << ", "
                   << res << ")" << std::endl;
        pythonFile << "y = np.linspace(" << -xYMax << ", " << xYMax << ", "
                   << res << ")" << std::endl;
        /*pythonFile << "fig, axs = plt.subplots(1,3, figsize=(15, 5))"
                   << std::endl;*/
        pythonFile << "fig, axs = plt.subplots(2, 3, sharex=True, sharey=True)"
                   << std::endl;

        pythonFile << "pdfZPos = np.array(pdfZPos)" << std::endl;
        pythonFile << "histogramZPos = np.array(histogramZPos)" << std::endl;
        pythonFile << "diffZPos=histogramZPos - pdfZPos" << std::endl;
        pythonFile << "absdiffZPos=np.abs(diffZPos)" << std::endl;
        pythonFile << "a = pdfZPos.shape[1] / pdfZPos.shape[0]" << std::endl;
        pythonFile << "vmaxvalZPos = pdfZPos.max()" << std::endl;
        pythonFile << "vmaxvalZPos = max(vmaxvalZPos, histogramZPos.max())"
                   << std::endl;
        pythonFile
            << "pdfZPos_plot = axs[0,0].imshow(pdfZPos, aspect=a, "
               "interpolation=\'nearest\', vmin=0, vmax=vmaxvalZPos, extent=[-"
            << xYMax << ", " << xYMax << ", -" << xYMax << ", " << xYMax
            << "], origin=\'lower\')" << std::endl;
        pythonFile
            << "histZPos_plot = axs[0,1].imshow(histogramZPos, aspect=a, "
               "interpolation=\'nearest\', vmin=0, vmax=vmaxvalZPos, extent=[-"
            << xYMax << ", " << xYMax << ", -" << xYMax << ", " << xYMax
            << "], origin=\'lower\')" << std::endl;
        pythonFile
            << "diffZPos_plot = axs[0,2].imshow(absdiffZPos, aspect=a, vmin=0, "
               "vmax=vmaxvalZPos, interpolation=\'nearest\', extent=[-"
            << xYMax << ", " << xYMax << ", -" << xYMax << ", " << xYMax
            << "], origin=\'lower\')" << std::endl;
        pythonFile << "axs[0,0].title.set_text(\'PDF\')" << std::endl;
        pythonFile << "axs[0,0].set_ylabel(\'Upper hemisphere\')" << std::endl;
        pythonFile << "axs[0,1].title.set_text(\'Histogram\')" << std::endl;
        pythonFile << "axs[0,2].title.set_text(\'Absolute difference\')"
                   << std::endl;
        pythonFile << "props = dict(fraction=0.046, pad=0.04)" << std::endl;
        pythonFile << "fig.colorbar(pdfZPos_plot, ax=axs[0,0], **props)"
                   << std::endl;
        pythonFile << "fig.colorbar(histZPos_plot, ax=axs[0,1], **props)"
                   << std::endl;
        pythonFile << "fig.colorbar(diffZPos_plot, ax=axs[0,2], **props)"
                   << std::endl;

        pythonFile << "pdfZNeg = np.array(pdfZNeg)" << std::endl;
        pythonFile << "histogramZNeg = np.array(histogramZNeg)" << std::endl;
        pythonFile << "diffZNeg=histogramZNeg - pdfZNeg" << std::endl;
        pythonFile << "absdiffZNeg=np.abs(diffZNeg)" << std::endl;
        pythonFile << "vmaxvalZNeg = pdfZNeg.max()" << std::endl;
        pythonFile << "vmaxvalZNeg = max(vmaxvalZNeg, histogramZNeg.max())"
                   << std::endl;
        pythonFile
            << "pdfZNeg_plot = axs[1,0].imshow(pdfZNeg, aspect=a, "
               "interpolation=\'nearest\', vmin=0, vmax=vmaxvalZNeg, extent=[-"
            << xYMax << ", " << xYMax << ", -" << xYMax << ", " << xYMax
            << "], origin=\'lower\')" << std::endl;
        pythonFile
            << "histZNeg_plot = axs[1,1].imshow(histogramZNeg, aspect=a, "
               "interpolation=\'nearest\', vmin=0, vmax=vmaxvalZNeg, extent=[-"
            << xYMax << ", " << xYMax << ", -" << xYMax << ", " << xYMax
            << "], origin=\'lower\')" << std::endl;
        pythonFile
            << "diffZNeg_plot = axs[1,2].imshow(absdiffZNeg, aspect=a, vmin=0, "
               "vmax=vmaxvalZNeg, interpolation=\'nearest\', extent=[-"
            << xYMax << ", " << xYMax << ", -" << xYMax << ", " << xYMax
            << "], origin=\'lower\')" << std::endl;
        pythonFile << "axs[1,0].title.set_text(\'PDF\')" << std::endl;
        pythonFile << "axs[1,0].set_ylabel(\'Lower hemisphere\')" << std::endl;
        pythonFile << "axs[1,1].title.set_text(\'Histogram\')" << std::endl;
        pythonFile << "axs[1,2].title.set_text(\'Absolute difference\')"
                   << std::endl;
        pythonFile << "fig.colorbar(pdfZNeg_plot, ax=axs[1,0], **props)"
                   << std::endl;
        pythonFile << "fig.colorbar(histZNeg_plot, ax=axs[1,1], **props)"
                   << std::endl;
        pythonFile << "fig.colorbar(diffZNeg_plot, ax=axs[1,2], **props)"
                   << std::endl;

        pythonFile << "plt.tight_layout()" << std::endl;
        pythonFile << "plt.savefig(\'"
                   << outFilenameWithoutExtension + "_test" +
                          std::to_string(k + 1)
                   << ".png\')" << std::endl;
        pythonFile << "plt.show()" << std::endl;

        if (resultZPos.first)
            std::cout << "Test Upper Hemisphere " << k + 1 << ": success"
                      << std::endl;
        else {
            std::cout << "Test Upper Hemisphere " << k + 1 << ": failure. "
                      << resultZPos.first << ". " << resultZPos.second
                      << std::endl;
            /*ParallelCleanup();
            return 1;*/
        }

        if (resultZNeg.first)
            std::cout << "Test Lower Hemisphere " << k + 1 << ": success"
                      << std::endl;
        else {
            std::cout << "Test Lower Hemisphere " << k + 1 << ": failure. "
                      << resultZNeg.first << ". " << resultZNeg.second
                      << std::endl;
            /*ParallelCleanup();
            return 1;*/
        }
    }

    pbrtCleanup();
    return 0;
}

int convergencecomparisons(int argc, char *argv[]) {
    MemoryArena arena;
    Options opt;
    pbrtInit(opt);

    // Total number of samples
    int nSamples = 10000;
    // Number of runs
    int nRuns = 8;

    Float thetaO = 0.;
    Float phiO = 0.;

    // Number of levels in the dictionary
    int nLevels = 8;
    // Number of 1D marginal distributions in the dictionary
    int N = 768;
    // Alpha roughness used by the dictionary generator
    float alphaDict = 0.5;
    // Logarithm microfacet density
    float logMicrofacetDensity = 21.;
    // Microfacet relative area
    float microfacetRelativeArea = 1.;
    // Alpha roughness of the material
    Float alphaX = 0.6;
    Float alphaY = 0.6;
    // Correlation factor
    Float rho = 0.;
    // Center of the pixel footprint
    Point2f st(0, 0);
    // Partial derivatives defining the extent of the pixel footprint
    Vector2f dstdx(0.01, 0.);
    Vector2f dstdy(0., 0.01);

    Float alphaXBaseMaterial = 0.1;
    Float alphaYBaseMaterial = 0.1;
    Float rhoBaseMaterial = 0.;
    Float densityRandomization = 2.;
    bool sampleVisibleArea = true;
    bool pgf = false;

    int i;
    auto parseArg = [&]() -> std::pair<std::string, double> {
        const char *ptr = argv[i];
        // Skip over a leading dash or two.
        CHECK_EQ(*ptr, '-');
        ++ptr;
        if (*ptr == '-') ++ptr;

        // Copy the flag name to the string.
        std::string flag;
        while (*ptr && *ptr != '=') flag += *ptr++;
        if (!*ptr && i + 1 == argc)
            usage("missing value after %s flag", argv[i]);
        const char *value = (*ptr == '=') ? (ptr + 1) : argv[++i];
        return {flag, atof(value)};
    };

    std::pair<std::string, double> arg;
    for (i = 0; i < argc; ++i) {
        if (argv[i][0] != '-')
            break;

        else {
            std::pair<std::string, double> arg = parseArg();
            if (std::get<0>(arg) == "nsamples") {
                nSamples = int(std::get<1>(arg));
            } else if (std::get<0>(arg) == "nruns") {
                nRuns = std::get<1>(arg);
            } else if (std::get<0>(arg) == "pgf") {
                pgf = std::get<1>(arg);
            } else if (std::get<0>(arg) == "alphax") {
                alphaX = std::get<1>(arg);
            } else if (std::get<0>(arg) == "alphay") {
                alphaY = std::get<1>(arg);
            } else if (std::get<0>(arg) == "rho") {
                rho = std::get<1>(arg);
            } else if (std::get<0>(arg) == "stx") {
                st.x = std::get<1>(arg);
            } else if (std::get<0>(arg) == "dstdxx") {
                dstdx.x = std::get<1>(arg);
            } else if (std::get<0>(arg) == "dstdyy") {
                dstdy.y = std::get<1>(arg);
            } else if (std::get<0>(arg) == "thetao") {
                thetaO = std::get<1>(arg);
            } else if (std::get<0>(arg) == "phio") {
                phiO = std::get<1>(arg);
            } else if (std::get<0>(arg) == "mra") {
                microfacetRelativeArea = std::get<1>(arg);
            } else
                usage();
        }
    }

    if (i + 1 >= argc)
        usage("missing second filename for \"convergencecomparisons\"");
    else if (i >= argc)
        usage("missing filenames for \"convergencecomparisons\"");

    const char *dictFilename = argv[i], *outFilename = argv[i + 1];
    std::string dataPath(dictFilename);
    Fresnel *fresnel = new FresnelNoOp();

    nRuns = Clamp(nRuns, 8, 9999999);

    // Load the texture containing the dictionary
    std::unique_ptr<TextureMapping2D> map;
    Float su = 1.;
    Float sv = 1.;
    Float du = 0.;
    Float dv = 0.;
    map.reset(new UVMapping2D(su, sv, du, dv));
    ImageTexture<RGBSpectrum, Spectrum> distributionsDict(
        std::move(map), dictFilename, false, 8.f, ImageWrap::Clamp, 1.f, false);

    // Will contain the N distributions
    std::vector<PiecewiseLinearDistribution1D> distributions;
    int distResolution = distributionsDict.Width() / (N / 3);

    // Populate the 1D distributions by using the dictionary stored in a
    // texture
    std::vector<Float> dist1D(distResolution);
    for (int lvl = 0; lvl < nLevels; lvl++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < distResolution; j++) {
                dist1D[j] = distributionsDict.Texel(
                    0, i / 3 * distResolution + j, lvl)[i % 3];
            }
            distributions.push_back(
                PiecewiseLinearDistribution1D(dist1D.data(), distResolution));
        }
    }

    BxDF *bxdf = createGlitteringBRDF(
        arena, alphaX, alphaY, rho, alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial, st, dstdx, dstdy, nLevels, N, logMicrofacetDensity,
        0.5, &distributions, microfacetRelativeArea, densityRandomization,
        sampleVisibleArea, false);

    GlitteringConductorReflection *glitteringConductorReflection =
        dynamic_cast<GlitteringConductorReflection *>(bxdf);

    BxDF *bsdfApprox = createGlitteringBSDF(
        arena, alphaX, alphaY, rho, alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial, st, dstdx, dstdy, nLevels, N, logMicrofacetDensity,
        0.5, &distributions, microfacetRelativeArea, densityRandomization,
        sampleVisibleArea, true);

    GlitteringDielectricScattering *glitteringDielectricScatteringApprox =
        dynamic_cast<GlitteringDielectricScattering *>(bsdfApprox);

    BxDF *bsdf = createGlitteringBSDF(
        arena, alphaX, alphaY, rho, alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial, st, dstdx, dstdy, nLevels, N, logMicrofacetDensity,
        0.5, &distributions, microfacetRelativeArea, densityRandomization,
        sampleVisibleArea, false);

    GlitteringDielectricScattering *glitteringDielectricScattering =
        dynamic_cast<GlitteringDielectricScattering *>(bsdf);

    int nMicrofacetWithinP =
        glitteringConductorReflection->GetNMicrofacetWithinP(st, dstdx, dstdy);

    std::cout << "Number of microfacets within the pixel footprint: "
              << nMicrofacetWithinP << std::endl;

    Float xYMax = 1.;

    Vector3f wo = SphericalDirection(std::sin(thetaO), std::cos(thetaO), phiO);

    // Will contain the average value of the shadowing term
    Float avgShadowing = 0.;
    int res = 512;
    for (int ixy = 0; ixy < res * res; ++ixy) {
        Float x = ((Mod(ixy, res) + 0.5) / Float(res) * 2. - 1.) * xYMax;
        Float y = (((ixy / res) + 0.5) / Float(res) * 2. - 1.) * xYMax;

        Vector2f projectedDirection(x, y);
        float normalZSqr = 1. - projectedDirection.x * projectedDirection.x -
                           projectedDirection.y * projectedDirection.y;
        if (normalZSqr <= 0.) continue;
        Vector3f wi(projectedDirection.x, projectedDirection.y,
                    std::sqrt(normalZSqr));

        Float BRDFValue = glitteringDielectricScattering->f(wo, wi)[0] * wi.z;
        Float BRDFValueZNeg =
            glitteringDielectricScattering->f(wo, -wi)[0] * wi.z;

        CHECK_GE(BRDFValue, 0.);
        CHECK_GE(BRDFValueZNeg, 0.);
        // 4. * xYMax * xYMax / (res * res): sub interval area
        // 1. / wi.z: Jacobian of the projected normal domain to
        // normal domain
        avgShadowing += (BRDFValue + BRDFValueZNeg) * 4. * xYMax * xYMax /
                        (res * res) / wi.z;
    }
    std::cout << "Average value of the shadowing term: " << avgShadowing
              << std::endl;

    // Extract the output file name without the extension
    std::string outputfilename(outFilename);
    std::size_t found = outputfilename.rfind(".");
    CHECK(found != std::string::npos);
    std::string outFilenameWithoutExtension =
        outputfilename.substr(std::size_t(0), found);
    std::string extension = outputfilename.substr(found, outputfilename.size());

    Float xMinTime = -MaxFloat;
    Float xMaxTime = MaxFloat;

    for (int k = 0; k < nRuns; ++k) {
        {
            // std::cout << "Run number " << k + 1 << std::endl;
            // Table of weights
            std::vector<Float> weights(nSamples);
            weights[0] = 0.;

            // Table of averaged weights
            std::vector<Float> avgWeights(nSamples + 1);

            // Accumulate samples
            for (int i = 0; i < nSamples; ++i) {
                Vector3f wi;
                // Not used by Sample_f
                Point2f u(0., 0.);
                Float pdf = 0.;
                BxDFType flags;
                Spectrum BRDFValue = glitteringDielectricScattering->Sample_f(
                    wo, &wi, u, &pdf, &flags);
                Float weight = 0.;
                if (pdf > 0.) weight = BRDFValue[0] * std::abs(wi.z) / pdf;
                weights[i] = weight;

                avgWeights[i + 1] =
                    weight / (i + 1) + avgWeights[i] * i / (i + 1);
            }

            // std::cout << "Average of the weights: " << avgWeights[nSamples]
            //           << std::endl;

            // Write the average of the weights to a file
            std::ofstream avgWeightsFile(outFilenameWithoutExtension +
                                         "_avgw_run" + std::to_string(k + 1) +
                                         ".dat");
            for (int j = 1; j <= nSamples; ++j)
                avgWeightsFile << avgWeights[j] << std::endl;
            avgWeightsFile.close();

            // Write the weights to a file
            std::ofstream weightsFile(outFilenameWithoutExtension + "_w_run" +
                                      std::to_string(k + 1) + ".dat");
            for (int j = 0; j < nSamples; ++j)
                weightsFile << weights[j] << std::endl;
            weightsFile.close();
        }
        {
            // Table of weights
            std::vector<Float> weights(nSamples);
            weights[0] = 0.;

            // Table of averaged weights
            std::vector<Float> avgWeights(nSamples + 1);

            // Accumulate samples
            for (int i = 0; i < nSamples; ++i) {
                Vector3f wi;
                // Not used by Sample_f
                Point2f u(0., 0.);
                Float pdf = 0.;
                BxDFType flags;
                Spectrum BRDFValue =
                    glitteringDielectricScatteringApprox->Sample_f(
                        wo, &wi, u, &pdf, &flags);
                Float weight = 0.;
                if (pdf > 0.) weight = BRDFValue[0] * std::abs(wi.z) / pdf;
                weights[i] = weight;

                avgWeights[i + 1] =
                    weight / (i + 1) + avgWeights[i] * i / (i + 1);
            }

            // std::cout
            //     << "Average weight using Gaussian approximation sampling: "
            //     << avgWeights[nSamples] << std::endl;

            // Write the average of the weights to a file
            std::ofstream avgWeightsFile(outFilenameWithoutExtension +
                                         "_approx_avgw_run" +
                                         std::to_string(k + 1) + ".dat");
            for (int j = 1; j <= nSamples; ++j)
                avgWeightsFile << avgWeights[j] << std::endl;
            avgWeightsFile.close();

            // Write the weights to a file
            std::ofstream weightsFile(outFilenameWithoutExtension +
                                      "_approx_w_run" + std::to_string(k + 1) +
                                      ".dat");
            for (int j = 0; j < nSamples; ++j)
                weightsFile << weights[j] << std::endl;
            weightsFile.close();
        }
    }
    // Radiances (weights) against number of samples (only 8 runs)
    {
        std::ofstream matplotlibfile(outFilenameWithoutExtension +
                                     "_radiance_n.py");

        matplotlibfile << "import matplotlib.pyplot as plt" << std::endl;
        matplotlibfile << "import matplotlib" << std::endl;
        matplotlibfile << "import numpy as np" << std::endl;
        matplotlibfile << std::endl;

        if (pgf) {
            matplotlibfile << "matplotlib.use(\"pgf\")" << std::endl;
            matplotlibfile << "matplotlib.rcParams.update({" << std::endl;
            matplotlibfile << "    \"pgf.texsystem\" : \"pdflatex\","
                           << std::endl;
            matplotlibfile << "    'font.family' : 'serif'," << std::endl;
            matplotlibfile << "    'text.usetex' : True," << std::endl;
            matplotlibfile << "    'pgf.rcfonts' : False," << std::endl;
            matplotlibfile << "})" << std::endl;
        }

        matplotlibfile << "def filterCurve(curve):" << std::endl;
        matplotlibfile << "    c_0_100 = curve[0:100]" << std::endl;
        matplotlibfile << "    c_100_1000 = curve[100:1000:10]" << std::endl;
        matplotlibfile << "    c_1000_10000 = curve[1000:10000:100]"
                       << std::endl;
        matplotlibfile << "    c_10000 = curve[9999]" << std::endl;
        matplotlibfile << "    filtered_curve = np.concatenate((c_0_100, "
                          "c_100_1000, c_1000_10000))"
                       << std::endl;
        matplotlibfile
            << "    filtered_curve = np.append(filtered_curve, c_10000)"
            << std::endl;
        matplotlibfile << "    return filtered_curve" << std::endl;

        matplotlibfile << "convergence_va_list = []" << std::endl;
        matplotlibfile << "for index in range(8):" << std::endl;
        matplotlibfile << "    convergence_va_list.append(np.loadtxt('"
                       << outFilenameWithoutExtension
                       << "_approx_avgw_run'+str(index+1)+'.dat'))"
                       << std::endl;

        matplotlibfile << "convergence_ve_list = []" << std::endl;
        matplotlibfile << "for index in range(8):" << std::endl;
        matplotlibfile << "    convergence_ve_list.append(np.loadtxt('"
                       << outFilenameWithoutExtension
                       << "_avgw_run'+str(index+1)+'.dat'))" << std::endl;

        matplotlibfile << "fig, ax = plt.subplots()" << std::endl;

        matplotlibfile << "x_values = filterCurve(np.arange(1,10001))"
                       << std::endl;

        matplotlibfile
            << "ax.plot(x_values, filterCurve(convergence_va_list[0]), "
               "'#ff7f0e', label='previous', linewidth=1)"
            << std::endl;
        matplotlibfile << "for index in range(1, 8):" << std::endl;
        matplotlibfile
            << "    ax.plot(x_values, filterCurve(convergence_va_list[index]), "
               "'#ff7f0e', linewidth=1)"
            << std::endl;

        matplotlibfile
            << "ax.plot(x_values, filterCurve(convergence_ve_list[0]), "
               "'#1f77b4', label='our', linewidth=1)"
            << std::endl;
        matplotlibfile << "for index in range(1, 8):" << std::endl;
        matplotlibfile
            << "    ax.plot(x_values, filterCurve(convergence_ve_list[index]), "
               "'#1f77b4', linewidth=1)"
            << std::endl;

        matplotlibfile
            << "ax.set(xlabel='Number of samples', title='$\\\\theta_o="
            << thetaO << "$, $\\\\alpha=" << alphaX << "$')" << std::endl;
        matplotlibfile << "ax.set_xlim(1,10000)" << std::endl;
        matplotlibfile << "ax.grid()" << std::endl;
        matplotlibfile << "ax.legend()" << std::endl;
        matplotlibfile << "plt.xscale('log')" << std::endl;
        matplotlibfile << "fig.set_size_inches(w=2.42, h=1.815)" << std::endl;
        matplotlibfile << "plt.tight_layout()" << std::endl;
        if (pgf)
            matplotlibfile << "fig.savefig(\"" << outFilenameWithoutExtension
                           << "_n.pgf\")" << std::endl;
        else
            matplotlibfile << "plt.show()" << std::endl;
    }

    // Pointwise boxplots
    {
        std::ofstream matplotlibfile(outFilenameWithoutExtension +
                                     "_pointwise_boxplots.py");

        matplotlibfile << "import matplotlib.pyplot as plt" << std::endl;
        matplotlibfile << "import matplotlib" << std::endl;
        matplotlibfile << "import numpy as np" << std::endl;
        matplotlibfile << std::endl;

        if (pgf) {
            matplotlibfile << "matplotlib.use(\"pgf\")" << std::endl;
            matplotlibfile << "matplotlib.rcParams.update({" << std::endl;
            matplotlibfile << "    \"pgf.texsystem\" : \"pdflatex\","
                           << std::endl;
            matplotlibfile << "    'font.family' : 'serif'," << std::endl;
            matplotlibfile << "    'text.usetex' : True," << std::endl;
            matplotlibfile << "    'pgf.rcfonts' : False," << std::endl;
            matplotlibfile << "})" << std::endl;
        }

        matplotlibfile << "def filterCurve(curve):" << std::endl;
        matplotlibfile << "    c_0_100 = curve[0:100]" << std::endl;
        matplotlibfile << "    c_100_1000 = curve[100:1000:10]" << std::endl;
        matplotlibfile << "    c_1000_10000 = curve[1000:10000:100]"
                       << std::endl;
        matplotlibfile << "    c_10000 = curve[9999]" << std::endl;
        matplotlibfile << "    filtered_curve = np.concatenate((c_0_100, "
                          "c_100_1000, c_1000_10000))"
                       << std::endl;
        matplotlibfile
            << "    filtered_curve = np.append(filtered_curve, c_10000)"
            << std::endl;
        matplotlibfile << "    return filtered_curve" << std::endl;

        matplotlibfile << "convergence_va_list = []" << std::endl;
        matplotlibfile << "for index in range(" << nRuns << "):" << std::endl;
        matplotlibfile << "    convergence_va_list.append(np.loadtxt('"
                       << outFilenameWithoutExtension
                       << "_approx_avgw_run'+str(index+1)+'.dat'))"
                       << std::endl;

        matplotlibfile << "convergence_ve_list = []" << std::endl;
        matplotlibfile << "for index in range(" << nRuns << "):" << std::endl;
        matplotlibfile << "    convergence_ve_list.append(np.loadtxt('"
                       << outFilenameWithoutExtension
                       << "_avgw_run'+str(index+1)+'.dat'))" << std::endl;

        matplotlibfile << "fig, ax = plt.subplots()" << std::endl;

        matplotlibfile << "F0_approx = np.zeros(10000)" << std::endl;
        matplotlibfile << "F25_approx = np.zeros(10000)" << std::endl;
        matplotlibfile << "F50_approx = np.zeros(10000)" << std::endl;
        matplotlibfile << "F75_approx = np.zeros(10000)" << std::endl;
        matplotlibfile << "F100_approx = np.zeros(10000)" << std::endl;
        matplotlibfile << "F0 = np.zeros(10000)" << std::endl;
        matplotlibfile << "F25 = np.zeros(10000)" << std::endl;
        matplotlibfile << "F50 = np.zeros(10000)" << std::endl;
        matplotlibfile << "F75 = np.zeros(10000)" << std::endl;
        matplotlibfile << "F100 = np.zeros(10000)" << std::endl;
        matplotlibfile << "" << std::endl;
        matplotlibfile << "# Go through the number of samples" << std::endl;
        matplotlibfile << "for N in range(10000):" << std::endl;
        matplotlibfile << "    # Go through each realisations" << std::endl;
        matplotlibfile << "    array_size = " << nRuns << std::endl;
        matplotlibfile << "    array_size_m1 = array_size - 1" << std::endl;
        matplotlibfile << "    F_r_N = np.zeros(array_size)" << std::endl;
        matplotlibfile << "    F_r_N_approx = np.zeros(array_size)"
                       << std::endl;
        matplotlibfile << "    for r in range(" << nRuns << "):" << std::endl;
        matplotlibfile << "        F_r_N[r] = convergence_ve_list[r][N]"
                       << std::endl;
        matplotlibfile << "        F_r_N_approx[r] = convergence_va_list[r][N]"
                       << std::endl;
        matplotlibfile << "    F_r_N_sorted = np.sort(F_r_N)" << std::endl;
        matplotlibfile << "    F_r_N_approx_sorted = np.sort(F_r_N_approx)"
                       << std::endl;
        matplotlibfile << "    F0[N] = F_r_N_sorted[0]" << std::endl;
        matplotlibfile << "    F0_approx[N] = F_r_N_approx_sorted[0]"
                       << std::endl;
        matplotlibfile << "    F100[N] = F_r_N_sorted[array_size_m1]"
                       << std::endl;
        matplotlibfile
            << "    F100_approx[N] = F_r_N_approx_sorted[array_size_m1]"
            << std::endl;
        matplotlibfile << "    F25[N] = F_r_N_sorted[int(0.25*array_size_m1)]"
                       << std::endl;
        matplotlibfile << "    F25_approx[N] = "
                          "F_r_N_approx_sorted[int(0.25*array_size_m1)]        "
                       << std::endl;
        matplotlibfile << "    F50[N] = F_r_N_sorted[int(0.5*array_size_m1)]"
                       << std::endl;
        matplotlibfile
            << "    F50_approx[N] = F_r_N_approx_sorted[int(0.5*array_size_m1)]"
            << std::endl;
        matplotlibfile << "    F75[N] = F_r_N_sorted[int(0.75*array_size_m1)]"
                       << std::endl;
        matplotlibfile << "    F75_approx[N] = "
                          "F_r_N_approx_sorted[int(0.75*array_size_m1)]"
                       << std::endl;

        matplotlibfile << "x_values = filterCurve(np.arange(1,10001))"
                       << std::endl;

        matplotlibfile << "ax.plot(x_values, filterCurve(F0), '#2ca02c', "
                          "linewidth=1, linestyle=(0, (1, 1)))"
                       << std::endl;
        matplotlibfile << "ax.plot(x_values, filterCurve(F100), '#2ca02c', "
                          "linewidth=1, linestyle=(0, (1, 1)))"
                       << std::endl;
        matplotlibfile << "ax.plot(x_values, filterCurve(F25), '#2ca02c', "
                          "linewidth=1, linestyle=(0, (5, 1)))"
                       << std::endl;
        matplotlibfile << "ax.plot(x_values, filterCurve(F75), '#2ca02c', "
                          "linewidth=1, linestyle=(0, (5, 1)))"
                       << std::endl;
        matplotlibfile
            << "ax.plot(x_values, filterCurve(F50), '#2ca02c', linewidth=1)"
            << std::endl;
        matplotlibfile << "" << std::endl;
        matplotlibfile << "ax.plot(x_values, filterCurve(F0_approx), "
                          "'#d62728', linewidth=1, linestyle=(0, (1, 1)))"
                       << std::endl;
        matplotlibfile << "ax.plot(x_values, filterCurve(F100_approx), "
                          "'#d62728', linewidth=1, linestyle=(0, (1, 1)))"
                       << std::endl;
        matplotlibfile << "ax.plot(x_values, filterCurve(F25_approx), "
                          "'#d62728', linewidth=1, linestyle=(0, (5, 1)))"
                       << std::endl;
        matplotlibfile << "ax.plot(x_values, filterCurve(F75_approx), "
                          "'#d62728', linewidth=1, linestyle=(0, (5, 1)))"
                       << std::endl;
        matplotlibfile << "ax.plot(x_values, filterCurve(F50_approx), "
                          "'#d62728', linewidth=1)"
                       << std::endl;

        matplotlibfile << "ax.set(title='$\\\\theta_o=" << thetaO
                       << "$, $\\\\alpha=" << alphaX << "$')" << std::endl;
        matplotlibfile << "ax.set_xlim(1,10000)" << std::endl;
        matplotlibfile << "ax.set_ylim(0,2)" << std::endl;
        matplotlibfile << "ax.grid()" << std::endl;
        matplotlibfile << "plt.xscale('log')" << std::endl;
        matplotlibfile << "fig.set_size_inches(w=2.42, h=1.815)" << std::endl;
        matplotlibfile << "plt.tight_layout()" << std::endl;
        if (pgf)
            matplotlibfile << "fig.savefig(\"" << outFilenameWithoutExtension
                           << "_pointwise_boxplot.pgf\")" << std::endl;
        else
            matplotlibfile << "plt.show()" << std::endl;
    }

    pbrtCleanup();
    return 0;
}

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = 1;  // Warning and above.

    if (argc < 2) usage();

    if (!strcmp(argv[1], "plotglitteringndf"))
        return plotglitteringndf(argc - 2, argv + 2);
    else if (!strcmp(argv[1], "plotglitteringsdf"))
        return plotglitteringsdf(argc - 2, argv + 2);
    else if (!strcmp(argv[1], "chisquaretestglitteringvndf"))
        return chisquaretestglitteringvndf(argc - 2, argv + 2);
    else if (!strcmp(argv[1], "chisquaretestglitteringbrdf"))
        return chisquaretestglitteringbrdf(argc - 2, argv + 2);
    else if (!strcmp(argv[1], "chisquaretestglitteringbsdf"))
        return chisquaretestglitteringbsdf(argc - 2, argv + 2);
    else if (!strcmp(argv[1], "convergencecomparisons"))
        return convergencecomparisons(argc - 2, argv + 2);
    else
        usage("unknown command \"%s\"", argv[1]);

    return 0;
}
