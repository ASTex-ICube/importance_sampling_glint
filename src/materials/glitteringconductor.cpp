
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

#include "materials/glitteringconductor.h"

#include <memory>
#include <utility>

#include "interaction.h"
#include "paramset.h"
#include "reflection.h"
#include "texture.h"

namespace pbrt {

GlitteringConductorMaterial::GlitteringConductorMaterial(
    const std::shared_ptr<Texture<Spectrum>> &R,
    const std::shared_ptr<Texture<Spectrum>> &eta,
    const std::shared_ptr<Texture<Spectrum>> &k,
    const std::shared_ptr<Texture<Float>> &alphaXGlitteringMaterial,
    const std::shared_ptr<Texture<Float>> &alphaYGlitteringMaterial,
    const std::shared_ptr<Texture<Float>> &rhoGlitteringMaterial,
    const std::shared_ptr<Texture<Float>> &alphaXBaseMaterial,
    const std::shared_ptr<Texture<Float>> &alphaYBaseMaterial,
    const std::shared_ptr<Texture<Float>> &rhoBaseMaterial,
    const std::shared_ptr<Texture<Spectrum>> &dictionary, int nLevels, int N,
    float alphaDict,
    const std::shared_ptr<Texture<Float>> &logMicrofacetDensity,
    const std::shared_ptr<Texture<Float>> &microfacetRelativeArea,
    float densityRandomisation, Float su, Float sv, bool fresnelNoOp,
    bool remapParameters, bool sampleVisibleArea, bool sampleApproximation)
    : R(R),
      eta(eta),
      k(k),
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
      fresnelNoOp(fresnelNoOp),
      remapParameters(remapParameters),
      sampleVisibleArea(sampleVisibleArea),
      sampleApproximation(sampleApproximation) {
    distributions =
        std::make_shared<std::vector<PiecewiseLinearDistribution1D>>();
    auto *dict =
        dynamic_cast<ImageTexture<RGBSpectrum, Spectrum> *>(dictionary.get());
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
                .push_back(
                    PiecewiseLinearDistribution1D(&dist1D[0], distResolution));
        }
    }
}

void GlitteringConductorMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
    bool allowMultipleLobes) const {
    Vector2f dstdx(su * si->dudx, sv * si->dvdx);
    Vector2f dstdy(su * si->dudy, sv * si->dvdy);
    Point2f st(su * si->uv[0], sv * si->uv[1]);

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);

    Fresnel *frMf;
    if (fresnelNoOp)
        frMf = ARENA_ALLOC(arena, FresnelNoOp)();
    else
        frMf = ARENA_ALLOC(arena, FresnelConductor)(1., eta->Evaluate(*si),
                                                    k->Evaluate(*si));

    RNG *rngEval = ARENA_ALLOC(arena, RNG)();
    RNG *rngSample = ARENA_ALLOC(arena, RNG)();

    Spectrum R = this->R->Evaluate(*si).Clamp();

    // Retrieve the parameters of the BRDF
    Float alphaXP = alphaXGlitteringMaterial->Evaluate(*si);
    Float alphaYP = alphaYGlitteringMaterial->Evaluate(*si);
    Float rhoP = rhoGlitteringMaterial->Evaluate(*si);
    Float logMicrofacetDensityP = logMicrofacetDensity->Evaluate(*si);
    Float microfacetRelativeAreaP = microfacetRelativeArea->Evaluate(*si);
    Float alphaXBaseMaterialP = alphaXBaseMaterial->Evaluate(*si);
    Float alphaYBaseMaterialP = alphaYBaseMaterial->Evaluate(*si);
    Float rhoBaseMaterialP = rhoBaseMaterial->Evaluate(*si);

    if (remapParameters) {
        // Alpha roughness in [0.1, 0.6]
        alphaXP = Clamp(alphaXP, 0., 1.) * 0.5 + 0.1;
        alphaYP = Clamp(alphaYP, 0., 1.) * 0.5 + 0.1;
        // Slope correlation factor in [-0.99, 0.99]
        rhoP = Clamp(rhoP, 0., 1.) * 1.98 - 0.99;
        // Logarithm microfacet density in [17., 25.]
        logMicrofacetDensityP = Clamp(logMicrofacetDensityP, 0., 1.) * 8. + 17.;
        // Microfacet relative area in [0.01, 1.]
        microfacetRelativeAreaP =
            Clamp(microfacetRelativeAreaP, 0., 1.) * 0.99 + 0.01;
    } else {
        // Alpha roughness in [0.02, 0.6]
        alphaXP = Clamp(alphaXP, 0.02, 0.6);
        alphaYP = Clamp(alphaYP, 0.02, 0.6);
        // Slope correlation factor in [-0.99, 0.99]
        rhoP = Clamp(rhoP, -0.99, 0.99);
        // Logarithm microfacet density in [10., 50.]
        logMicrofacetDensityP = Clamp(logMicrofacetDensityP, 10., 50.);
        // Microfacet relative area in [0.01, 1.]
        microfacetRelativeAreaP = Clamp(microfacetRelativeAreaP, 0.01, 1.);
    }

    // Limits the ray footprint to maximum anisotropy
    ClampRayFootprint(st, dstdx, dstdy);

    Float u = st.x * 457. + 269. * st.y + 599. * dstdx.x + 397. * dstdx.y +
              29. * dstdy.x + 283. * dstdy.y;
    Float du = std::abs(u) - std::floor(std::abs(u));
    rngSample->SetSequence(du * std::numeric_limits<int>::max());

    si->bsdf->Add(ARENA_ALLOC(arena, GlitteringConductorReflection)(
        R, frMf, rngEval, rngSample, st, dstdx, dstdy, alphaXP, alphaYP, rhoP,
        alphaXBaseMaterialP, alphaYBaseMaterialP, rhoBaseMaterialP,
        distResolution, distributions.get(), nLevels, N, alphaDict,
        logMicrofacetDensityP, microfacetRelativeAreaP, densityRandomisation,
        sampleVisibleArea, sampleApproximation));
}

const int CopperSamples = 56;
const Float CopperWavelengths[CopperSamples] = {
    298.7570554, 302.4004341, 306.1337728, 309.960445,  313.8839949,
    317.9081487, 322.036826,  326.2741526, 330.6244747, 335.092373,
    339.6826795, 344.4004944, 349.2512056, 354.2405086, 359.374429,
    364.6593471, 370.1020239, 375.7096303, 381.4897785, 387.4505563,
    393.6005651, 399.9489613, 406.5055016, 413.2805933, 420.2853492,
    427.5316483, 435.0322035, 442.8006357, 450.8515564, 459.2006593,
    467.8648226, 476.8622231, 486.2124627, 495.936712,  506.0578694,
    516.6007417, 527.5922468, 539.0616435, 551.0407911, 563.5644455,
    576.6705953, 590.4008476, 604.8008683, 619.92089,   635.8162974,
    652.5483053, 670.1847459, 688.8009889, 708.4810171, 729.3186941,
    751.4192606, 774.9011125, 799.8979226, 826.5611867, 855.0632966,
    885.6012714};

const Float CopperN[CopperSamples] = {
    1.400313, 1.38,  1.358438, 1.34,  1.329063, 1.325, 1.3325,   1.34,
    1.334375, 1.325, 1.317812, 1.31,  1.300313, 1.29,  1.281563, 1.27,
    1.249062, 1.225, 1.2,      1.18,  1.174375, 1.175, 1.1775,   1.18,
    1.178125, 1.175, 1.172812, 1.17,  1.165312, 1.16,  1.155312, 1.15,
    1.142812, 1.135, 1.131562, 1.12,  1.092437, 1.04,  0.950375, 0.826,
    0.645875, 0.468, 0.35125,  0.272, 0.230813, 0.214, 0.20925,  0.213,
    0.21625,  0.223, 0.2365,   0.25,  0.254188, 0.26,  0.28,     0.3};

const Float CopperK[CopperSamples] = {
    1.662125, 1.687, 1.703313, 1.72,  1.744563, 1.77,  1.791625, 1.81,
    1.822125, 1.834, 1.85175,  1.872, 1.89425,  1.916, 1.931688, 1.95,
    1.972438, 2.015, 2.121562, 2.21,  2.177188, 2.13,  2.160063, 2.21,
    2.249938, 2.289, 2.326,    2.362, 2.397625, 2.433, 2.469187, 2.504,
    2.535875, 2.564, 2.589625, 2.605, 2.595562, 2.583, 2.5765,   2.599,
    2.678062, 2.809, 3.01075,  3.24,  3.458187, 3.67,  3.863125, 4.05,
    4.239563, 4.43,  4.619563, 4.817, 5.034125, 5.26,  5.485625, 5.717};

GlitteringConductorMaterial *CreateGlitteringConductorMaterial(
    const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> R =
        mp.GetSpectrumTexture("R", Spectrum(1.f));

    static Spectrum copperN =
        Spectrum::FromSampled(CopperWavelengths, CopperN, CopperSamples);
    std::shared_ptr<Texture<Spectrum>> eta =
        mp.GetSpectrumTexture("eta", copperN);
    static Spectrum copperK =
        Spectrum::FromSampled(CopperWavelengths, CopperK, CopperSamples);
    std::shared_ptr<Texture<Spectrum>> k = mp.GetSpectrumTexture("k", copperK);

    std::shared_ptr<Texture<Float>> alphaX = mp.GetFloatTexture("alphax", 0.5);
    std::shared_ptr<Texture<Float>> alphaY = mp.GetFloatTexture("alphay", 0.5);
    std::shared_ptr<Texture<Float>> rho = mp.GetFloatTexture("rho", 0.);
    std::shared_ptr<Texture<Float>> logMicrofacetDensity =
        mp.GetFloatTexture("logmicrofacetdensity", 20.f);
    std::shared_ptr<Texture<Float>> microfacetRelativeArea =
        mp.GetFloatTexture("microfacetrelativearea", 1.);

    std::shared_ptr<Texture<Float>> alphaXBaseMaterial =
        mp.GetFloatTexture("alphaxbasematerial", 0.01);
    std::shared_ptr<Texture<Float>> alphaYBaseMaterial =
        mp.GetFloatTexture("alphaybasematerial", 0.01);
    std::shared_ptr<Texture<Float>> rhoBaseMaterial =
        mp.GetFloatTexture("rhobasematerial", 0.);

    std::shared_ptr<Texture<Spectrum>> dictionary =
        mp.GetSpectrumTextureOrNull("dictionary");

    int nLevels = mp.FindInt("nlevels", 8);
    int N = mp.FindInt("N", 768);
    Float alphaDict = mp.FindFloat("alpha_dict", 0.5f);

    Float densityRandomisation = mp.FindFloat("densityrandomisation", 2.f);

    Float su = mp.FindFloat("su", 1.f);
    Float sv = mp.FindFloat("sv", 1.f);

    bool fresnelNoOp = mp.FindBool("fresnelnoop", false);
    bool remapParameters = mp.FindBool("remapparameters", false);
    bool sampleVisibleArea = mp.FindBool("samplevisiblearea", true);
    bool sampleApproximation = mp.FindBool("sampleapproximation", false);

    return new GlitteringConductorMaterial(
        R, eta, k, alphaX, alphaY, rho, alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial, dictionary, nLevels, N, alphaDict,
        logMicrofacetDensity, microfacetRelativeArea, densityRandomisation, su,
        sv, fresnelNoOp, remapParameters, sampleVisibleArea,
        sampleApproximation);
}

// It evaluates the normal distribution function (NDF) into the ray
// footprint P (st, dst0, dst1)
Float GlitteringConductorReflection::D_P(const Vector3f &wm, Point2f st,
                                         Vector2f dst0, Vector2f dst1) const {
    Vector2f slopeH(-wm.x / wm.z, -wm.y / wm.z);

    // Without footprint, we evaluate the SDF at the last LOD
    Float minorLength = std::min(dst0.Length(), dst1.Length());
    if (minorLength == 0.) {
        return SDFLastLOD(slopeH) / (wm.z * wm.z * wm.z * wm.z);
    } else {
        // LOD of the MIP hierarchy
        Float l = N_LEVELS_VIRTUAL_MIPMAP - 1. + std::log2(minorLength);
        // Negative LODs are not allowed. Set N_LEVELS_VIRTUAL_MIPMAP to greater
        // value to avoid this case
        CHECK_GE(l, 0.);

        // Discrete MIP LOD
        int il = int(std::floor(l));

        // Weight of the two adjacent LODs
        Float w = l - il;

        // Number of cells on the side at current level
        uint64_t nCellsOnTheSide =
            std::pow(2, (N_LEVELS_VIRTUAL_MIPMAP - 1 - il));
        double cellArea = 1. / (nCellsOnTheSide * nCellsOnTheSide);
        // Number of microfacets in the current LOD
        double n_il = cellArea * std::exp(logMicrofacetDensity);

        // Corresponding continuous distribution LOD
        Float LOD_dist_il = std::log(n_il) / 1.38629;  // 2. * log(2) = 1.38629

        // The model is not defined for negative distribution LOD
        // At this scale, we begin to see each cell individually and each
        // distribution of 4 lobes within these cells.
        if (LOD_dist_il < 0.) {
            // We will increase the LOD of the distribution by
            // std::abs(std::floor(LOD_dist_il)), so we need to increase the LOD
            // of the MIP hierarchy by std::abs(std::floor(LOD_dist_il)) in
            // order to keep consistency between the two hierarchies
            int floor_LOD_dist_il = std::abs(std::floor(LOD_dist_il));
            il += floor_LOD_dist_il;
            w = 0.;

            LOD_dist_il = 0.;
        }

        // Corresponding continuous distribution LOD + 1
        Float LOD_dist_ilp1 = LOD_dist_il + 1.;

        // Default value at distribution LOD >= nLevels
        Float P22_P_il = SDFLastLOD(slopeH);
        Float P22_P_ilp1 = P22_P_il;

        // When the distribution level is greater than nLevels, we don't need to
        // go through the pixel footprint's cells.
        if (LOD_dist_il < nLevels && w != 1.)
            P22_P_il =
                P22_P_discreteLOD(il, slopeH, st, dst0, dst1, LOD_dist_il);
        if (LOD_dist_ilp1 < nLevels && w != 0.)
            P22_P_ilp1 = P22_P_discreteLOD(il + 1, slopeH, st, dst0, dst1,
                                           LOD_dist_ilp1);

        // The value of the P-SDF
        Float P22_P = Lerp(w, P22_P_il, P22_P_ilp1);

        // The value of the P-NDF (slope to normal transformation)
        return P22_P / (wm.z * wm.z * wm.z * wm.z);
    }
}
Float GlitteringConductorReflection::P22_P(const Vector2f &slope_h, Point2f st,
                                           Vector2f dst0, Vector2f dst1) const {
    Vector2f slopeH = slope_h;

    // Without footprint, we evaluate the SDF at the last LOD
    Float minorLength = std::min(dst0.Length(), dst1.Length());
    if (minorLength == 0.) {
        return SDFLastLOD(slopeH);
    } else {
        // LOD of the MIP hierarchy
        Float l = N_LEVELS_VIRTUAL_MIPMAP - 1. + std::log2(minorLength);
        // Negative LODs are not allowed. Set N_LEVELS_VIRTUAL_MIPMAP to greater
        // value to avoid this case
        CHECK_GE(l, 0.);

        // Discrete MIP LOD
        int il = int(std::floor(l));

        // Weight of the two adjacent LODs
        Float w = l - il;

        // Number of cells on the side at current level
        uint64_t nCellsOnTheSide =
            std::pow(2, (N_LEVELS_VIRTUAL_MIPMAP - 1 - il));
        double cellArea = 1. / (nCellsOnTheSide * nCellsOnTheSide);
        // Number of microfacets in the current LOD
        double n_il = cellArea * std::exp(logMicrofacetDensity);

        // Corresponding continuous distribution LOD
        Float LOD_dist_il = std::log(n_il) / 1.38629;  // 2. * log(2) = 1.38629

        // The model is not defined for negative distribution LOD
        // At this scale, we begin to see each cell individually and each
        // distribution of 4 lobes within these cells.
        if (LOD_dist_il < 0.) {
            // We will increase the LOD of the distribution by
            // std::abs(std::floor(LOD_dist_il)), so we need to increase the LOD
            // of the MIP hierarchy by std::abs(std::floor(LOD_dist_il)) in
            // order to keep consistency between the two hierarchies
            int floor_LOD_dist_il = std::abs(std::floor(LOD_dist_il));
            il += floor_LOD_dist_il;
            w = 0.;

            LOD_dist_il = 0.;
        }

        // Corresponding continuous distribution LOD + 1
        Float LOD_dist_ilp1 = LOD_dist_il + 1.;

        // Default value at distribution LOD >= nLevels
        Float P22_P_il = SDFLastLOD(slopeH);
        Float P22_P_ilp1 = P22_P_il;

        // When the distribution level is greater than nLevels, we don't need to
        // go through the pixel footprint's cells.
        if (LOD_dist_il < nLevels && w != 1.)
            P22_P_il =
                P22_P_discreteLOD(il, slopeH, st, dst0, dst1, LOD_dist_il);
        if (LOD_dist_ilp1 < nLevels && w != 0.)
            P22_P_ilp1 = P22_P_discreteLOD(il + 1, slopeH, st, dst0, dst1,
                                           LOD_dist_ilp1);

        // The value of the P-SDF
        Float P22_P = Lerp(w, P22_P_il, P22_P_ilp1);

        return P22_P;
    }
}

// It evaluates the visible normal distribution function (VNDF) into the ray
// footprint P
Float GlitteringConductorReflection::Visible_D_P(const Vector3f &wm,
                                                 const Vector3f &wo) const {
    Float dotWOWM = Dot(wo, wm);
    // Back facing surfaces return 0
    if (dotWOWM <= 0. || wo.z <= 0. || wm.z <= 0.) return 0.;
    // P-NDF value
    Float D_PValue = D_P(wm, st, dstdx, dstdy);
    // P-VNDF value
    Float density = G1VCavity(wm, wo) * dotWOWM * D_PValue / wo.z;
    return density;
}

// It samples the normal distribution function (NDF) within the ray footprint
// P
Vector3f GlitteringConductorReflection::Sample_D_P() const {
    Float minorLength = std::min(dstdx.Length(), dstdy.Length());

    // If the intersection has not partial derivatives or if we sample the
    // approximated SDF, we sample the SDF at the last LOD
    if (minorLength == 0. || sampleApproximation) {
        Vector2f sampledSlope = SampleSDFLastLOD();

        return SlopeToNormal(sampledSlope);
    } else {
        // LOD of the MIP hierarchy
        Float l = N_LEVELS_VIRTUAL_MIPMAP - 1. + std::log2(minorLength);
        // Negative LODs are not allowed. Set N_LEVELS_VIRTUAL_MIPMAP to greater
        // value to avoid this case
        CHECK_GE(l, 0.);

        // Discrete MIP LOD
        int il = int(floor(l));

        // Weight of the two adjacent LODs
        float w = l - float(il);

        // Number of cells on the side at current level
        uint64_t nCellsOnTheSide =
            std::pow(2, (N_LEVELS_VIRTUAL_MIPMAP - 1 - il));
        double cellArea = 1. / (nCellsOnTheSide * nCellsOnTheSide);
        // Number of microfacets in the current LOD
        double n_il = cellArea * std::exp(logMicrofacetDensity);
        // Corresponding continuous distribution LOD
        Float LOD_dist_il = std::log(n_il) / 1.38629;  // 2. * log(2) = 1.38629

        // The model is not defined for negative distribution LOD
        // At this scale, we begin to see each cell individually and each
        // distribution of 4 lobes within these cells.
        if (LOD_dist_il < 0.) {
            // We will increase the LOD of the distribution by
            // std::abs(std::floor(LOD_dist_il)), so we need to increase the LOD
            // of the MIP hierarchy by std::abs(std::floor(LOD_dist_il)) in
            // order to keep consistency between the two hierarchies
            int floor_LOD_dist_il = std::abs(std::floor(LOD_dist_il));
            il += floor_LOD_dist_il;
            w = 0.;

            LOD_dist_il = 0.;
        }

        // Sample an adjacent LOD based on w value
        Float uLOD = rngSample->UniformFloat();
        if (uLOD < w) {
            il += 1;
            LOD_dist_il += 1.;
        }

        Vector2f sampledSlope;
        // Samples a cell within the pixel footprint then its associated SDF
        if (LOD_dist_il < nLevels) {
            sampledSlope =
                Sample_P22_P_discreteLOD(il, st, dstdx, dstdy, LOD_dist_il);
        }
        // When the distribution level is greater than nLevels, we don't need to
        // sample the pixel footprint In this case, we sample the SDF of the
        // last LOD
        else {
            sampledSlope = SampleSDFLastLOD();
        }
        return SlopeToNormal(sampledSlope);
    }
}

// It samples the visible normal distribution function (NDF) within the ray
// footprint P
Vector3f GlitteringConductorReflection::Sample_Visible_D_P(
    const Vector3f &wo) const {
    // See Heitz and d'Eon 2014 for more details

    // Samples the NDF
    Vector3f wm = Sample_D_P();

    // Computes the dual of the sampled micronormal (the other facet of the
    // V-cavity)
    Vector3f wmP(-wm.x, -wm.y, wm.z);
    // Just in case
    CHECK_GT(wm.z, 0.);

    Float dotWOWMP = std::max(Dot(wo, wmP), Float(0.));
    Float dotWOWM = std::max(Dot(wo, wm), Float(0.));
    Float deno = dotWOWM + dotWOWMP;
    // Just in case
    CHECK_GT(deno, 0.);
    // Projected area of the other facet of the V-cavity
    Float projectedAreaWMP = dotWOWMP / deno;
    // Samples a facet of the V-cavity
    if (rngSample->UniformFloat() < projectedAreaWMP)
        return wmP;
    else
        return wm;
}

// It **evaluates** the ith 1D distribution with x at level l.
Float GlitteringConductorReflection::P(Float x, int i, int l) const {
    // The distribution equals 0 after 4 standard deviations.
    // Slope standard deviation = alpha / sqrt(2)
    // 0.707106 \approx 1 / sqrt(2)
    Float alpha_dict_isqrt2_4 = alphaDict * 0.707106 * 4.;
    if (x >= alpha_dict_isqrt2_4) {
        return 0.;
    }

    Float tableIndex_f = (x / alpha_dict_isqrt2_4) * (distResolution - 1);
    int tableIndex_i = std::floor(tableIndex_f);

    Float d_index = tableIndex_f - tableIndex_i;

    // Linear interpolation between to discrete values of the 1D distribution
    Float density =
        (1 - d_index) * (*distributions)[l * N + i].func[tableIndex_i] +
        (d_index) * (*distributions)[l * N + i].func[tableIndex_i + 1];

    // We re-normalize just in case the dictionary is not perfectly normalized
    Float normalization =
        (*distributions)[l * N + i].funcInt * 2. * alpha_dict_isqrt2_4;
    return density / normalization;
}

// It **samples** the ith 1D distribution at level l. It uses two random
// numbers: one for sampling the positive 1D distribution, the other to
// inverse the sign of the sampled value.
Float GlitteringConductorReflection::Sample_P(const Point2f &u, int i,
                                              int l) const {
    // The distribution equals 0 after 4 standard deviations.
    // Slope standard deviation = alpha / sqrt(2)
    // 0.707106 \approx 1 / sqrt(2)
    Float alpha_dist_isqrt2_4 = alphaDict * 0.707106 * 4.;

    // Samples the positive part of the 1D distribution
    // Sampled value is in [0, 1)
    Float sample = (*distributions)[l * N + i].SampleContinuous(u.x, nullptr);
    // The 1D distribution is an even function
    // [-1, 1]
    if (u.y < 0.5) sample *= -1.;

    // [-alpha_dist_isqrt2_4, alpha_dist_isqrt2_4]
    return sample * alpha_dist_isqrt2_4;
}

// It **evaluates** the transformed slope distribution function with
// slope_h, at the discrete LOD l within the virtual MIP hierarchy. This
// distribution is associated to cell (s0, t0) in the MIP hierarchy. It uses
// two 1D distributions at level lDist.
Float GlitteringConductorReflection::P22_M(int l, uint64_t s0, uint64_t t0,
                                           const Vector2f &slopeH,
                                           Float lDist) const {
    // Coherent index between MIPMAP LODs
    uint64_t s064 = s0 * 1 << l;
    uint64_t t064 = t0 * 1 << l;

    // Cantor pairing function
    uint64_t seed = ((s064 + t064) * (s064 + t064 + 1)) / 2 + t064;

    // Seeds pseudo random generator with the coherent index
    rngEval->SetSequence(seed);

    // Recovers RMS roughness
    Float sigmaX = alphaXGlitteringMaterial / Sqrt2;
    Float sigmaY = alphaYGlitteringMaterial / Sqrt2;

    // If the current cell has no microfacets, a default material is evaluated
    // to avoid "black holes" on the surface.
    Float uMicrofacetRelativeArea = rngEval->UniformFloat();
    if (uMicrofacetRelativeArea > microfacetRelativeArea)
        return BivariateNormalDistribution(slopeH.x, slopeH.y,
                                           alphaXBaseMaterial,
                                           alphaYBaseMaterial, rhoBaseMaterial);

    // Sample a Gaussian to randomize the distribution LOD around the
    // distribution level lDist.
    Float uDensityRandomisation = rngEval->UniformFloat();
    lDist = SampleNormalDistribution(uDensityRandomisation, lDist,
                                     densityRandomisation);
    int iLDist = Clamp(int(std::round(lDist)), 0, nLevels);

    // If the current distribution level is greater than or equal to nLevels,
    // then the SDF is the SDF of the last LOD
    if (iLDist == nLevels) {
        return SDFLastLOD(slopeH);
    }

    // Inverse linear transformation represented by matrix invM
    Float SIGMA_DICT = alphaDict / Sqrt2;
    Float invM00 =
        SIGMA_DICT / (sigmaX * std::sqrt(1. - rhoGlitteringMaterial *
                                                  rhoGlitteringMaterial));
    Float invM01 = -SIGMA_DICT * rhoGlitteringMaterial /
                   (sigmaY * std::sqrt(1. - rhoGlitteringMaterial *
                                                rhoGlitteringMaterial));
    Float invM10 = 0.;
    Float invM11 = SIGMA_DICT / sigmaY;

    // Samples a random rotation to disalign the glints
    Float uTheta = rngEval->UniformFloat();
    Float theta = 2. * Pi * uTheta;
    Float cosTheta = std::cos(theta);
    Float sinTheta = std::sin(theta);

    // Inverse linear transformation represented by matrix invRM
    // (incorporate a random rotation)
    Float invRM00 = cosTheta * invM00 + sinTheta * invM10;
    Float invRM01 = cosTheta * invM01 + sinTheta * invM11;
    Float invRM10 = -sinTheta * invM00 + cosTheta * invM10;
    Float invRM11 = -sinTheta * invM01 + cosTheta * invM11;

    // Get back to original (dictionary) slope space
    Vector2f slopeHO(slopeH.x * invRM00 + slopeH.y * invRM01,
                     slopeH.x * invRM10 + slopeH.y * invRM11);
    // 1D distribution are even functions
    Vector2f absSlopeHO(std::abs(slopeHO.x), std::abs(slopeHO.y));

    // The distribution equals 0 after 4 standard deviations.
    // Slope standard deviation = alpha / sqrt(2)
    // 0.707106 \approx 1 / sqrt(2)
    Float alpha_dist_isqrt2_4 = alphaDict * 0.707106 * 4.;
    if (absSlopeHO.x > alpha_dist_isqrt2_4 ||
        absSlopeHO.y > alpha_dist_isqrt2_4)
        return 0.;

    // Sample the two 1D distributions
    Float u1 = rngEval->UniformFloat();
    Float u2 = rngEval->UniformFloat();
    int i = int(u1 * Float(N));
    int j = int(u2 * Float(N));

    // Evaluate the two 1D distributions in the original space
    Float P_i = P(absSlopeHO.x, i, iLDist);
    Float P_j = P(absSlopeHO.y, j, iLDist);

    // The value of the transformed SDF
    return SIGMA_DICT * SIGMA_DICT * P_i * P_j /
           (sigmaX * sigmaY *
            std::sqrt(1. - rhoGlitteringMaterial * rhoGlitteringMaterial));
}

// It **samples** the transformed slope distribution function at the
// discrete LOD l on the virtual MIP hierarchy. This distribution is
// associated to cell (s0, t0) in the MIP hierarchy. It uses two 1D
// distributions at level lDist.
Vector2f GlitteringConductorReflection::Sample_P22_M(int l, uint64_t s0,
                                                     uint64_t t0,
                                                     Float lDist) const {
    // Coherent index
    uint64_t s064 = s0 * 1 << l;
    uint64_t t064 = t0 * 1 << l;

    // https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
    uint64_t seed = ((s064 + t064) * (s064 + t064 + 1)) / 2 + t064;

    // Seed pseudo random generator to retrieve cell material properties
    // (microfacet relative area, density randomisation, marginal distributions)
    rngEval->SetSequence(seed);

    // Recover RMS roughness
    Float sigmaX = alphaXGlitteringMaterial / Sqrt2;
    Float sigmaY = alphaYGlitteringMaterial / Sqrt2;

    // If the current cell has no microfacets, a default material is sampled
    Float uMicrofacetRelativeArea = rngEval->UniformFloat();
    if (uMicrofacetRelativeArea > microfacetRelativeArea) {
        Point2f uSampling(rngSample->UniformFloat(), rngSample->UniformFloat());
        return SampleBivariateNormalDistribution(
            uSampling, alphaXBaseMaterial, alphaYBaseMaterial, rhoBaseMaterial);
    }

    // Sample a Gaussian to randomise the distribution LOD around the
    // distribution level lDist
    Float uDensityRandomisation = rngEval->UniformFloat();
    lDist = SampleNormalDistribution(uDensityRandomisation, lDist,
                                     densityRandomisation);
    int iLDist = Clamp(int(std::round(lDist)), 0, nLevels);

    // If the current distribution level is greater than or equal to nLevels,
    // then the SDF of the last LOD is sampled
    if (iLDist == nLevels) {
        return SampleSDFLastLOD();
    }

    // Linear transformation represented by matrix M
    Float SIGMA_DICT = alphaDict / Sqrt2;
    Float M00 = 1. / SIGMA_DICT *
                (sigmaX *
                 std::sqrt(1. - rhoGlitteringMaterial * rhoGlitteringMaterial));
    Float M01 = 1. / SIGMA_DICT * rhoGlitteringMaterial * sigmaX;
    Float M10 = 0.;
    Float M11 = 1. / SIGMA_DICT * sigmaY;

    // Samples a random rotation to disalign the glints
    Float uTheta = rngEval->UniformFloat();
    Float theta = 2. * Pi * uTheta;
    Float cosTheta = std::cos(theta);
    Float sinTheta = std::sin(theta);

    // Linear transformation represented by matrix MR
    // (incorporate a random rotation)
    Float MR00 = cosTheta * M00 + sinTheta * M01;
    Float MR01 = -sinTheta * M00 + cosTheta * M01;
    Float MR10 = cosTheta * M10 + sinTheta * M11;
    Float MR11 = -sinTheta * M10 + cosTheta * M11;

    // Compute the indices of the two marginal distributions
    Float u1 = rngEval->UniformFloat();
    Float u2 = rngEval->UniformFloat();
    int i = int(u1 * Float(N));
    int j = int(u2 * Float(N));

    // Sample original (dictionary) slope space
    Point2f uSamplingI(rngSample->UniformFloat(), rngSample->UniformFloat());
    Float slopeHXO = Sample_P(uSamplingI, i, iLDist);
    Point2f uSamplingJ(rngSample->UniformFloat(), rngSample->UniformFloat());
    Float slopeHYO = Sample_P(uSamplingJ, j, iLDist);

    // Transforms the sample
    Vector2f slopeH(slopeHXO * MR00 + slopeHYO * MR01,
                    slopeHXO * MR10 + slopeHYO * MR11);

    return slopeH;
}

// It **evaluates** the 2D slope distribution function with slope_h in the
// ray footprint P (st, dst0, dst1) at the discrete LOD l in the virtual MIP
// hierarchy. It uses two 1D distributions at level lDist.
Float GlitteringConductorReflection::P22_P_discreteLOD(
    int l, const Vector2f &slope_h, Point2f st, Vector2f dst0, Vector2f dst1,
    Float lDist) const {
    int pyrSize = std::pow(2, N_LEVELS_VIRTUAL_MIPMAP - 1 - l);

    st[0] = st[0] * pyrSize - 0.5f;
    st[1] = st[1] * pyrSize - 0.5f;

    dst0[0] *= pyrSize;
    dst0[1] *= pyrSize;
    dst1[0] *= pyrSize;
    dst1[1] *= pyrSize;

    // Compute ellipse coefficients to bound filter region
    Float A = dst0[1] * dst0[1] + dst1[1] * dst1[1] + 1;
    Float B = -2 * (dst0[0] * dst0[1] + dst1[0] * dst1[1]);
    Float C = dst0[0] * dst0[0] + dst1[0] * dst1[0] + 1;
    Float invF = 1 / (A * C - B * B * 0.25f);
    A *= invF;
    B *= invF;
    C *= invF;

    // Compute the ellipse's bounding box in texture space
    Float det = -B * B + 4 * A * C;
    Float invDet = 1 / det;
    Float uSqrt = std::sqrt(det * C), vSqrt = std::sqrt(A * det);
    int s0 = std::ceil(st[0] - 2 * invDet * uSqrt);
    int s1 = std::floor(st[0] + 2 * invDet * uSqrt);
    int t0 = std::ceil(st[1] - 2 * invDet * vSqrt);
    int t1 = std::floor(st[1] + 2 * invDet * vSqrt);

    // Scan over ellipse bound and compute quadratic equation
    Float sum(0.f);
    Float sumWts = 0;
    for (int it = t0; it <= t1; ++it) {
        Float tt = it - st[1];
        int itIndex = it - t0;
        for (int is = s0; is <= s1; ++is) {
            Float ss = is - st[0];
            int isIndex = is - s0;
            // Compute squared radius and evaluate the SDF if inside ellipse
            Float r2 = A * ss * ss + B * ss * tt + C * tt * tt;
            if (r2 < 1) {
                Float alpha = 2;
                // Weighting function used in pbrt-v3 EWA function
                Float W_P = std::exp(-alpha * r2) - std::exp(-alpha);

                sum += P22_M(l, is, it, slope_h, lDist) * W_P;
                sumWts += W_P;
            }
        }
    }
    return sum / sumWts;
}

// It **samples** the 2D slope distribution function in the ray
// footprint P (st, dst0, dst1) at the discrete LOD l in the virtual MIP
// hierarchy. It uses two 1D distributions at level lDist.
Vector2f GlitteringConductorReflection::Sample_P22_P_discreteLOD(
    int l, Point2f st, Vector2f dst0, Vector2f dst1, Float lDist) const {
    int pyrSize = std::pow(2, N_LEVELS_VIRTUAL_MIPMAP - 1 - l);
    st[0] = st[0] * pyrSize - 0.5f;
    st[1] = st[1] * pyrSize - 0.5f;
    dst0[0] *= pyrSize;
    dst0[1] *= pyrSize;
    dst1[0] *= pyrSize;
    dst1[1] *= pyrSize;

    // Compute ellipse coefficients to bound filter region
    Float A = dst0[1] * dst0[1] + dst1[1] * dst1[1] + 1;
    Float B = -2 * (dst0[0] * dst0[1] + dst1[0] * dst1[1]);
    Float C = dst0[0] * dst0[0] + dst1[0] * dst1[0] + 1;
    Float invF = 1 / (A * C - B * B * 0.25f);
    A *= invF;
    B *= invF;
    C *= invF;

    // Compute the ellipse's bounding box in texture space
    Float det = -B * B + 4 * A * C;
    Float invDet = 1 / det;
    Float uSqrt = std::sqrt(det * C), vSqrt = std::sqrt(A * det);
    int s0 = std::ceil(st[0] - 2 * invDet * uSqrt);
    int s1 = std::floor(st[0] + 2 * invDet * uSqrt);
    int t0 = std::ceil(st[1] - 2 * invDet * vSqrt);
    int t1 = std::floor(st[1] + 2 * invDet * vSqrt);

    // Width of the bounding box of the ray footprint
    int resX = s1 - s0;
    // Just in case
    CHECK_GT(resX, 0);
    // Height of the bounding box of the ray footprint
    int resY = t1 - t0;
    // Just in case
    CHECK_GT(resY, 0);

    // Will contain the weights associated to the cells within the ray
    // footprint
    std::vector<Float> WPValues((resX + 1) * (resY + 1));

    // Scan over ellipse bound and compute quadratic equation
    Float sum(0.f);
    Float sumWts = 0;
    for (int it = t0; it <= t1; ++it) {
        Float tt = it - st[1];
        int itIndex = it - t0;
        for (int is = s0; is <= s1; ++is) {
            Float ss = is - st[0];
            int isIndex = is - s0;
            // Compute squared radius and compute the weight if inside ellipse
            Float r2 = A * ss * ss + B * ss * tt + C * tt * tt;
            if (r2 < 1) {
                Float alpha = 2;
                // Weighting function used in pbrt-v3 EWA function
                Float W_P = std::exp(-alpha * r2) - std::exp(-alpha);

                WPValues[itIndex * (resX + 1) + isIndex] = W_P;

                sumWts += W_P;
            } else {
                WPValues[itIndex * (resX + 1) + isIndex] = 0.;
            }
        }
    }
    // Normalizes the weights
    for (int i = 0; i < WPValues.size(); ++i) WPValues[i] /= sumWts;
    // Create a piecewise constant 2D distribution with the weights of the cells
    // within the ray footprint
    Distribution2D WPLOD(WPValues.data(), resX + 1, resY + 1);

    // Sample a cell within the ray footprint
    Point2f u(rngSample->UniformFloat(), rngSample->UniformFloat());
    Point2i sampledCell = WPLOD.SampleDiscrete(u);
    sampledCell += Vector2i(s0, t0);

    // Samples the transformed SDF within the sampled cell
    return Sample_P22_M(l, sampledCell.x, sampledCell.y, lDist);
}

// It returns the value of the distribution function for a given pair of
// directions wo and wi.
Spectrum GlitteringConductorReflection::f(const Vector3f &wo,
                                          const Vector3f &wi) const {
    VLOG(1) << "wo.z: " << wo.z
            << " -> number of microfacets within ray footprint: "
            << GetNMicrofacetWithinP(st, dstdx, dstdy);

    Float cosThetaO = CosTheta(wo), cosThetaI = CosTheta(wi);
    Vector3f wh = wi + wo;
    // Handle degenerate cases for microfacet reflection
    if (cosThetaI <= 0 || cosThetaO <= 0) return Spectrum(0.);
    if (wh.x == 0 && wh.y == 0 && wh.z <= 0) return Spectrum(0.);

    wh = Normalize(wh);

    // Local masking shadowing
    if (Dot(wo, wh) <= 0. || Dot(wi, wh) <= 0.) return Spectrum(0.);

    float D_P_value = D_P(wh, st, dstdx, dstdy);

    Spectrum F = fresnel->Evaluate(Dot(wi, Faceforward(wh, Vector3f(0, 0, 1))));

    Float G1wowh = G1VCavity(wh, wo);
    Float G1wiwh = G1VCavity(wh, wi);
    float G = G1wowh * G1wiwh;

    return R * (F * G * D_P_value) / (4. * wo.z * wi.z);
}

// It samples the distribution function for the given direction wo
Spectrum GlitteringConductorReflection::Sample_f(const Vector3f &wo,
                                                 Vector3f *wi, const Point2f &u,
                                                 Float *pdf,
                                                 BxDFType *sampledType) const {
    // Sample microfacet orientation $\wh$ and reflected direction $\wi$
    if (wo.z <= 0) return 0.;
    Vector3f wh;
    if (sampleVisibleArea)
        wh = Sample_Visible_D_P(wo);
    else
        wh = Sample_D_P();
    Float dotWOWM = Dot(wo, wh);
    // We cannot reflect with a back facing normal
    if (dotWOWM <= 0.) return 0.;

    *wi = Reflect(wo, wh);
    if (wi->z <= 0.) return 0.;

    // Compute PDF of _wi_ for microfacet reflection
    Spectrum F =
        fresnel->Evaluate(Dot(*wi, Faceforward(wh, Vector3f(0, 0, 1))));
    if (sampleVisibleArea) {
        //---------------------------------------------------------------------
        // PDF value
        Float visible_D_P = 0.;
        if (!sampleApproximation) {
            visible_D_P = Visible_D_P(wh, wo);
            *pdf = visible_D_P / (4 * Dot(wo, wh));
        } else {
            Vector2f slopeH(-wh.x / wh.z, -wh.y / wh.z);
            Float D_PValue = SDFLastLOD(slopeH) / (wh.z * wh.z * wh.z * wh.z);
            Float density = G1VCavity(wh, wo) * dotWOWM * D_PValue / wo.z;
            Float visible_D_PApprox = density;
            *pdf = visible_D_PApprox / (4 * Dot(wo, wh));
        }
        //---------------------------------------------------------------------

        //---------------------------------------------------------------------
        // BRDF value
        Spectrum BRDFValue;
        // We have already compute D_P
        if (!sampleApproximation)
            BRDFValue = R * visible_D_P * G1VCavity(wh, *wi) * F /
                        (4. * Dot(wo, wh) * wi->z);
        // We do not have already compute D_P
        else
            BRDFValue = R * D_P(wh, st, dstdx, dstdy) * G1VCavity(wh, wo) *
                        G1VCavity(wh, *wi) * F / (4. * wo.z * wi->z);
        return BRDFValue;
        //---------------------------------------------------------------------
    } else {
        //---------------------------------------------------------------------
        // PDF value
        Float D_PValue;
        if (!sampleApproximation) {
            D_PValue = D_P(wh, st, dstdx, dstdy);
            *pdf = D_PValue * wh.z / (4 * Dot(wo, wh));
        } else {
            Vector2f slopeH(-wh.x / wh.z, -wh.y / wh.z);
            Float D_PValueApprox =
                SDFLastLOD(slopeH) / (wh.z * wh.z * wh.z * wh.z);
            *pdf = D_PValueApprox * wh.z / (4 * Dot(wo, wh));
        }
        //---------------------------------------------------------------------

        //---------------------------------------------------------------------
        // BRDF value
        Float G1wowh = G1VCavity(wh, wo);
        Float G1wiwh = G1VCavity(wh, *wi);
        float G = G1wowh * G1wiwh;

        Spectrum BRDFValue;
        // We have already compute D_P
        if (!sampleApproximation)
            BRDFValue = R * D_PValue * G * F / (4. * wo.z * wi->z);
        // We do not have already compute D_P
        else
            BRDFValue =
                R * D_P(wh, st, dstdx, dstdy) * G * F / (4. * wo.z * wi->z);
        return BRDFValue;
        //---------------------------------------------------------------------
    }
}

// It evaluates the PDF of the sampling procedure Sample_f
Float GlitteringConductorReflection::Pdf(const Vector3f &wo,
                                         const Vector3f &wi) const {
    if (wo.z <= 0. || wi.z <= 0.) return 0.;
    Vector3f wh = Normalize(wo + wi);
    Float pdf;
    Float dotWOWM = Dot(wo, wh);
    // Should be very rare
    if (dotWOWM <= 0.) return 0.;
    if (sampleVisibleArea) {
        if (!sampleApproximation)
            pdf = Visible_D_P(wh, wo) / (4 * Dot(wo, wh));
        else {
            Vector2f slopeH(-wh.x / wh.z, -wh.y / wh.z);
            Float D_PValue = SDFLastLOD(slopeH) / (wh.z * wh.z * wh.z * wh.z);
            Float density = G1VCavity(wh, wo) * dotWOWM * D_PValue / wo.z;
            Float visible_D_PApprox = density;
            pdf = visible_D_PApprox / (4 * Dot(wo, wh));
        }
    } else {
        if (!sampleApproximation)
            pdf = D_P(wh, st, dstdx, dstdy) * wh.z / (4 * Dot(wo, wh));
        else {
            Vector2f slopeH(-wh.x / wh.z, -wh.y / wh.z);
            Float D_PValueApprox =
                SDFLastLOD(slopeH) / (wh.z * wh.z * wh.z * wh.z);
            pdf = D_PValueApprox * wh.z / (4 * Dot(wo, wh));
        }
    }
    return pdf;
}

// ToString utility function
std::string GlitteringConductorReflection::ToString() const {
    // TODO
    return std::string("[ SparklingReflection R: ") + R.ToString() +
           std::string(" fresnel: ") + fresnel->ToString() + std::string(" ]");
}

}  // namespace pbrt