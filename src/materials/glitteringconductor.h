
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_GLITTERINGCONDUCTOR_H
#define PBRT_MATERIALS_GLITTERINGCONDUCTOR_H

// The BRDF is built from a virtual mipmap. The number of levels must be set to
// a high value in order to avoid negative levels when evaluating the BRDF
#define N_LEVELS_VIRTUAL_MIPMAP 32

#include "material.h"
#include "microfacet.h"
#include "pbrt.h"
#include "reflection.h"
#include "rng.h"
#include "sampling.h"
#include "spectrum.h"
#include "textures/imagemap.h"

namespace pbrt {

// Clamp the ray footprint defined by center st and vectors dst0 and dst1
inline void ClampRayFootprint(Point2f st, Vector2f &dst0, Vector2f &dst1) {
    Float maxAnisotropy = 8.;
    if (dst0.LengthSquared() < dst1.LengthSquared()) std::swap(dst0, dst1);
    // Compute ellipse minor and major axes
    Float majorLength = dst0.Length();
    Float minorLength = dst1.Length();

    // Clamp ellipse eccentricity if too large
    if (minorLength * maxAnisotropy < majorLength && minorLength > 0) {
        Float scale = majorLength / (minorLength * maxAnisotropy);
        dst1 *= scale;
    }
}

// Sparkling material
class GlitteringConductorMaterial : public Material {
  public:
    GlitteringConductorMaterial(
        const std::shared_ptr<Texture<Spectrum>> &R,
        const std::shared_ptr<Texture<Spectrum>> &eta,
        const std::shared_ptr<Texture<Spectrum>> &k,
        const std::shared_ptr<Texture<Float>> &alphaXGlitteringMaterial,
        const std::shared_ptr<Texture<Float>> &alphaYGlitteringMaterial,
        const std::shared_ptr<Texture<Float>> &rhoGlitteringMaterial,
        const std::shared_ptr<Texture<Float>> &alphaXBaseMaterial,
        const std::shared_ptr<Texture<Float>> &alphaYBaseMaterial,
        const std::shared_ptr<Texture<Float>> &rhoBaseMaterial,
        const std::shared_ptr<Texture<Spectrum>> &dictionary, int nLevels,
        int N, float alphaDict,
        const std::shared_ptr<Texture<Float>> &logMicrofacetDensity,
        const std::shared_ptr<Texture<Float>> &microfacetRelativeArea,
        float densityRandomisation, Float su, Float sv, bool fresnelNoOp,
        bool remapParameters, bool sampleVisibleArea, bool sampleApproximation);

    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // The coefficient of specular reflection.
    std::shared_ptr<Texture<Spectrum>> R;

    /*
        eta: Index of refraction to use in computing the material's reflectance.
        k: Absorption coefficient to use in computing the material's
        reflectance.
    */
    std::shared_ptr<Texture<Spectrum>> eta, k;

    /*
       alphaXGlitteringMaterial: Roughness in the x direction for the
        glittering material.
       alphaYGlitteringMaterial: Roughness in the y direction for the
        glittering material.
       rhoGlitteringMaterial: Slope correlation factor for the
        glittering material.
    */
    std::shared_ptr<Texture<Float>> alphaXGlitteringMaterial,
        alphaYGlitteringMaterial, rhoGlitteringMaterial;

    /*
        logMicrofacetDensity: Logarithm of the microfacet density, without the
            microfacet relative area applied.
        microfacetRelativeArea: Percentage of the surface without microfacets.

        Note: **Effective** microfacet density = exp(logMicrofacetDensity) *
                                                    microfacetRelativeArea
    */
    std::shared_ptr<Texture<Float>> logMicrofacetDensity;
    std::shared_ptr<Texture<Float>> microfacetRelativeArea;

    /*
       alphaXBaseMaterial: Roughness in the x direction for the base material.
       alphaYBaseMaterial: Roughness in the y direction for the base material.
       rhoBaseMaterial: Slope correlation factor for the base material.

       Note: these parameters are only used when the microfacetRelativeArea <
        1, i.e., when (1 - microfacetRelativeArea) of the surface is covered by
        a base material.
    */
    std::shared_ptr<Texture<Float>> alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial;

    // Dictionary of piecewise linear 1D distributions
    std::shared_ptr<std::vector<PiecewiseLinearDistribution1D>> distributions;
    /*
        nLevels: Number of levels of detail for the 1D distributions
        N: Number of multi-scale 1D distributions in the dictionary
    */
    int nLevels, N;

    // Roughness used to generate the dictionary.
    float alphaDict;

    // Randomly changes the density of microfacets per cell. More precisely,
    // this parameter is the standard deviation of a normal distribution which
    // is sampled to randomize the microfacet density.
    float densityRandomisation;

    // Resolution of the a 1D distribution in the dictionary.
    int distResolution;

    // su: scale the x/u direction of the surface
    // sv: scale the y/v direction of the surface
    Float su, sv;

    // If true, the material is non-absorptive
    bool fresnelNoOp;

    // If true, material parameters are remapped from [0-1] to adequate
    // intervals. Usefull when using textures.
    bool remapParameters;

    // If true, samples the visible area of the normal distribution function.
    bool sampleVisibleArea;

    // If true, samples the Gaussian approximation of the normal distribution
    // function.
    bool sampleApproximation;
};

GlitteringConductorMaterial *CreateGlitteringConductorMaterial(
    const TextureParams &mp);

// Sample a normal distribution with mean mu and standard deviation sigma by
// using a uniform random number U.
inline Float SampleNormalDistribution(Float U, Float mu, Float sigma) {
    Float x = sigma * 1.414213f * ErfInv(2.0f * U - 1.0f) + mu;
    return x;
}

// Evaluate a bivariate normal distribution with position (x,y), standard
// deviations sigmaX and sigmaY, and slope correlation factor rho.
inline Float BivariateNormalDistribution(Float x, Float y, Float sigmaX,
                                         Float sigmaY, Float rho) {
    Float xSqr = x * x;
    Float ySqr = y * y;
    Float sigmaXSqr = sigmaX * sigmaX;
    Float sigmaYSqr = sigmaY * sigmaY;

    Float z = ((xSqr / sigmaXSqr) - ((2. * rho * x * y) / (sigmaX * sigmaY)) +
               (ySqr / sigmaYSqr));
    return std::exp(-z / (2. * (1. - rho * rho))) /
           (2. * Pi * sigmaX * sigmaY * std::sqrt(1. - rho * rho));
}

// Sample a bivariate normal distribution with uniform random position u,
// standard deviations sigmaX and sigmaY, and slope correlation factor rho.
inline Vector2f SampleBivariateNormalDistribution(Point2f u, Float sigmaX,
                                                  Float sigmaY, Float rho) {
    // Sample an original 2D centered unit-variance Gaussian
    Float sampleXO = SampleNormalDistribution(u.x, 0., 1.);
    Float sampleYO = SampleNormalDistribution(u.y, 0., 1.);

    Vector2f sampleO(sampleXO, sampleYO);

    // Matrix M transformed the original unit-variance Gaussian
    // into a non-axis aligned Gaussian with standard deviations sigmaX and
    // sigmaY and correlation factor rho
    // See Chermain et al. 2021 for more information
    Float M00 = sigmaX * std::sqrt(1. - rho * rho);
    Float M01 = rho * sigmaX;
    Float M10 = 0.;
    Float M11 = sigmaY;

    // Transform the sample using the linear transformation M
    Vector2f sample(sampleO.x * M00 + sampleO.y * M01,
                    sampleO.x * M10 + sampleO.y * M11);

    return sample;
}

// Sparkling BRDF class
class GlitteringConductorReflection : public BxDF {
  public:
    GlitteringConductorReflection(
        const Spectrum &R, Fresnel *fresnel, RNG *rngEval, RNG *rngSample,
        Point2f st, Vector2f dstdx, Vector2f dstdy,
        Float alphaXGlitteringMaterial, Float alphaYGlitteringMaterial,
        Float rhoGlitteringMaterial, Float alphaXBaseMaterial,
        Float alphaYBaseMaterial, Float rhoBaseMaterial, int distResolution,
        std::vector<PiecewiseLinearDistribution1D> *distributions, int nLevels,
        int N, Float alphaDict, Float logMicrofacetDensity,
        Float microfacetRelativeArea, Float densityRandomisation,
        bool sampleVisibleArea, bool sampleApproximation)
        : BxDF(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION)),
          R(R),
          fresnel(fresnel),
          rngEval(rngEval),
          rngSample(rngSample),
          st(st),
          dstdx(dstdx),
          dstdy(dstdy),
          alphaXGlitteringMaterial(alphaXGlitteringMaterial),
          alphaYGlitteringMaterial(alphaYGlitteringMaterial),
          rhoGlitteringMaterial(rhoGlitteringMaterial),
          alphaXBaseMaterial(alphaXBaseMaterial),
          alphaYBaseMaterial(alphaYBaseMaterial),
          rhoBaseMaterial(rhoBaseMaterial),
          distributions(distributions),
          nLevels(nLevels),
          N(N),
          alphaDict(alphaDict),
          logMicrofacetDensity(logMicrofacetDensity),
          microfacetRelativeArea(microfacetRelativeArea),
          densityRandomisation(densityRandomisation),
          distResolution(distResolution),
          sampleVisibleArea(sampleVisibleArea),
          sampleApproximation(sampleApproximation) {}

    // It returns the value of the distribution function for a given pair of
    // directions wo and wi.
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    // It samples the distribution function for the given direction wo
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    // It evaluates the PDF of the sampling procedure Sample_f
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;

    // It evaluates the normal distribution function (NDF) within ray
    // footprint P (st, dst0, dst1)
    Float D_P(const Vector3f &wm, Point2f st, Vector2f dst0,
              Vector2f dst1) const;
    Float P22_P(const Vector2f &slope_h, Point2f st, Vector2f dst0,
                Vector2f dst1) const;
    // It evaluates the visible normal distribution function (VNDF) within the
    // ray footprint P
    Float Visible_D_P(const Vector3f &wm, const Vector3f &wo) const;

    // It samples the normal distribution function (NDF) within the ray
    // footprint
    // P
    Vector3f Sample_D_P() const;
    // It samples the visible normal distribution function (NDF) within the ray
    // footprint P
    Vector3f Sample_Visible_D_P(const Vector3f &wo) const;

    // It **evaluates** the 2D slope distribution function with slope_h in the
    // ray footprint P (st, dst0, dst1) at the discrete LOD l in the virtual MIP
    // hierarchy. It uses two 1D distributions at level lDist.
    Float P22_P_discreteLOD(int l, const Vector2f &slope_h, Point2f st,
                            Vector2f dst0, Vector2f dst1, Float lDist) const;
    // It **samples** the 2D slope distribution function in the ray
    // footprint P (st, dst0, dst1) at the discrete LOD l in the virtual MIP
    // hierarchy. It uses two 1D distributions at level lDist.
    Vector2f Sample_P22_P_discreteLOD(int l, Point2f st, Vector2f dst0,
                                      Vector2f dst1, Float lDist) const;

    // It **evaluates** the transformed slope distribution function with
    // slope_h, at the discrete LOD l within the virtual MIP hierarchy. This
    // distribution is associated to cell (s0, t0) in the MIP hierarchy. It uses
    // two 1D distributions at level lDist.
    Float P22_M(int l, uint64_t s0, uint64_t t0, const Vector2f &slopeH,
                Float lDist) const;
    // It **samples** the transformed slope distribution function at the
    // discrete LOD l on the virtual MIP hierarchy. This distribution is
    // associated to cell (s0, t0) in the MIP hierarchy. It uses two 1D
    // distributions at level lDist.
    Vector2f Sample_P22_M(int l, uint64_t s0, uint64_t t0, Float lDist) const;

    // It **evaluates** the ith 1D distribution with x at level l.
    Float P(Float x, int i, int l) const;
    // It **samples** the ith 1D distribution at level l. It uses two random
    // numbers: one for sampling the positive 1D distribution, the other to
    // inverse the sign of the sampled value.
    Float Sample_P(const Point2f &u, int i, int l) const;

    // It gives the ray footprint area
    inline Float GetRayFootprintArea(Point2f st, Vector2f dst0,
                                     Vector2f dst1) const {
        Float minorLength = std::min(dst0.Length(), dst1.Length());
        Float l = N_LEVELS_VIRTUAL_MIPMAP - 1. + std::log2(minorLength);
        CHECK_GE(l, 0.);
        l = std::max((Float)0., l);
        int il = int(floor(l));
        Float w = l - il;

        uint64_t nCellsOnTheSide =
            std::pow(2, (N_LEVELS_VIRTUAL_MIPMAP - 1 - il));
        double cellArea = 1. / (nCellsOnTheSide * nCellsOnTheSide);

        return Lerp(w, cellArea, cellArea * 4.);
    }

    // It gives the number of microfacet within the ray footprint P
    inline int GetNMicrofacetWithinP(Point2f st, Vector2f dst0,
                                     Vector2f dst1) const {
        Float nParticlesWithinRayFootprint =
            GetRayFootprintArea(st, dst0, dst1) * microfacetRelativeArea *
            std::exp(logMicrofacetDensity);
        return nParticlesWithinRayFootprint;
    }
    // It evaluates the SDF at the last LOD
    inline Float SDFLastLOD(const Vector2f &slopeH) const {
        Float slopeDensity = BivariateNormalDistribution(
            slopeH.x, slopeH.y, alphaXGlitteringMaterial / Sqrt2,
            alphaYGlitteringMaterial / Sqrt2, rhoGlitteringMaterial);

        if (microfacetRelativeArea != 1.)
            // The distribution of the last LOD is mixed with the
            // distribution of the base material.
            slopeDensity =
                slopeDensity * microfacetRelativeArea +
                (1. - microfacetRelativeArea) *
                    BivariateNormalDistribution(
                        slopeH.x, slopeH.y, alphaXBaseMaterial / Sqrt2,
                        alphaYBaseMaterial / Sqrt2, rhoBaseMaterial);
        return slopeDensity;
    }
    // It samples the SDF at the last LOD
    inline Vector2f SampleSDFLastLOD() const {
        Float uMaterial = rngSample->UniformFloat();
        Point2f u(rngSample->UniformFloat(), rngSample->UniformFloat());
        Vector2f sampledSlope;
        // Samples the glittering SDF at the last LOD
        if (uMaterial < microfacetRelativeArea) {
            sampledSlope = SampleBivariateNormalDistribution(
                u, alphaXGlitteringMaterial / Sqrt2,
                alphaYGlitteringMaterial / Sqrt2, rhoGlitteringMaterial);
        }
        // Samples the base material SDF
        else {
            sampledSlope = SampleBivariateNormalDistribution(
                u, alphaXBaseMaterial / Sqrt2, alphaYBaseMaterial / Sqrt2,
                rhoBaseMaterial);
        }
        return sampledSlope;
    }

    std::string ToString() const;

    // The coefficient of specular reflection
    const Spectrum R;

    // Fresnel reflection
    const Fresnel *fresnel;

    // Center of the ray footprint
    const Point2f st;
    // Derivatives of the intersection point on surface with respect to
    // pixel coordinates. Used to construct the ray footprint.
    const Vector2f dstdx, dstdy;

    // Random number generator used during evaluation
    RNG *rngEval;
    // Random number generator used during sampling
    RNG *rngSample;

    // Roughness in the x direction
    const Float alphaXGlitteringMaterial;
    // Roughness in the y direction
    const Float alphaYGlitteringMaterial;
    // Slope correlation factor
    const Float rhoGlitteringMaterial;

    /*
        logMicrofacetDensity: Logarithm of the microfacet density, without the
            microfacet relative area applied.
        microfacetRelativeArea: Percentage of the surface without microfacets.

        Note: **Effective** microfacet density = exp(logMicrofacetDensity) *
                                                    microfacetRelativeArea
    */
    Float logMicrofacetDensity;
    Float microfacetRelativeArea;

    /*
       alphaXBaseMaterial: Roughness in the x direction for the base material.
       alphaYBaseMaterial: Roughness in the y direction for the base material.
       rhoBaseMaterial: Slope correlation factor for the base material.

       Note: these parameters are only used when the microfacetRelativeArea <
        1, i.e., when (1 - microfacetRelativeArea) of the surface is covered by
        a base material.
    */
    const Float alphaXBaseMaterial, alphaYBaseMaterial, rhoBaseMaterial;

    // Randomly changes the density of microfacets per cell. More precisely,
    // this parameter is the standard deviation of a normal distribution which
    // is sampled to randomize the microfacet density.
    Float densityRandomisation;

    // Roughness used to generate the dictionary
    const Float alphaDict;

    // Dictionary of piecewise linear 1D distributions
    std::vector<PiecewiseLinearDistribution1D> *distributions;

    /*
        nLevels: Number of level of detail for the 1D distributions
        N: Number of multi-scale 1D distributions in the dictionary
    */
    int nLevels, N, distResolution;

    // If true, samples the visible area of the normal distribution function.
    bool sampleVisibleArea;

    // If true, samples the Gaussian approximation of the normal distribution
    // function.
    bool sampleApproximation;
};

}  // namespace pbrt

#endif