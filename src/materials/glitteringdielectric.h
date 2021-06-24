
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

#ifndef PBRT_MATERIALS_GLITTERINGDIELECTRIC_H
#define PBRT_MATERIALS_GLITTERINGDIELECTRIC_H

#include "material.h"
#include "pbrt.h"
#include "reflection.h"
#include "spectrum.h"

#include "sampling.h"
#include "textures/imagemap.h"

namespace pbrt {

class GlitteringDielectricMaterial : public Material {
  public:
    GlitteringDielectricMaterial(
        const std::shared_ptr<Texture<Spectrum>> &Kr,
        const std::shared_ptr<Texture<Spectrum>> &Kt,
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
        float densityRandomisation, Float su, Float sv, bool remapParameters,
        bool sampleVisibleArea, bool sampleApproximation,
        const std::shared_ptr<Texture<Float>> &index)
        : Kr(Kr),
          Kt(Kt),
          alphaXGlitteringMaterial(alphaXGlitteringMaterial),
          alphaYGlitteringMaterial(alphaYGlitteringMaterial),
          rhoGlitteringMaterial(rhoGlitteringMaterial),
          alphaXBaseMaterial(alphaXBaseMaterial),
          alphaYBaseMaterial(alphaYBaseMaterial),
          rhoBaseMaterial(rhoBaseMaterial),
          nLevels(nLevels),
          N(N),
          alphaDict(alphaDict),
          logMicrofacetDensity(logMicrofacetDensity),
          microfacetRelativeArea(microfacetRelativeArea),
          densityRandomisation(densityRandomisation),
          su(su),
          sv(sv),
          remapParameters(remapParameters),
          sampleVisibleArea(sampleVisibleArea),
          sampleApproximation(sampleApproximation),
          index(index) {
        distributions =
            std::make_shared<std::vector<PiecewiseLinearDistribution1D>>();
        auto *dict = dynamic_cast<ImageTexture<RGBSpectrum, Spectrum> *>(
            dictionary.get());
        distResolution = dict->Width() / (N / 3);
        std::vector<Float> dist1D(distResolution);

        // Convert the 1D distributions stored into the dictionary into
        // a vector of PiecewiseLinearDistribution1D. This structure enables
        // sampling of the piecewise linear 1D distribution.
        for (int lvl = 0; lvl < nLevels; lvl++) {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < distResolution; j++) {
                    dist1D[j] =
                        dict->Texel(0, i / 3 * distResolution + j, lvl)[i % 3];
                }
                (*distributions)
                    .push_back(PiecewiseLinearDistribution1D(&dist1D[0],
                                                             distResolution));
            }
        }
    }
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    std::shared_ptr<Texture<Spectrum>> Kr, Kt;
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

    // If true, material parameters are remapped from [0-1] to adequate
    // intervals. Usefull when using textures.
    bool remapParameters;

    // If true, samples the visible area of the normal distribution function.
    bool sampleVisibleArea;

    // If true, samples the Gaussian approximation of the normal distribution
    // function.
    bool sampleApproximation;

    std::shared_ptr<Texture<Float>> index;
};

GlitteringDielectricMaterial *CreateGlitteringDielectricMaterial(
    const TextureParams &mp);

class GlitteringDielectricScattering : public BxDF {
  public:
    GlitteringDielectricScattering(
        const Spectrum &R, const Spectrum &T, Float etaA, Float etaB,
        TransportMode mode, RNG *rngEval, RNG *rngSample, Point2f st,
        Vector2f dstdx, Vector2f dstdy, Float alphaXGlitteringMaterial,
        Float alphaYGlitteringMaterial, Float rhoGlitteringMaterial,
        Float alphaXBaseMaterial, Float alphaYBaseMaterial,
        Float rhoBaseMaterial, int distResolution,
        std::vector<PiecewiseLinearDistribution1D> *distributions, int nLevels,
        int N, Float alphaDict, Float logMicrofacetDensity,
        Float microfacetRelativeArea, Float densityRandomisation,
        bool sampleVisibleArea, bool sampleApproximation)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_GLOSSY)),
          R(R),
          T(T),
          etaA(etaA),
          etaB(etaB),
          mode(mode),
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
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;

    std::string ToString() const;

  private:
    const Spectrum R, T;
    const Float etaA, etaB;
    const TransportMode mode;

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
