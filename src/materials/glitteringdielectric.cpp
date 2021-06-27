
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

#include "materials/glitteringdielectric.h"

#include "interaction.h"
#include "materials/glitteringconductor.h"
#include "paramset.h"
#include "reflection.h"
#include "spectrum.h"
#include "texture.h"

namespace pbrt {

void GlitteringDielectricMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
    bool allowMultipleLobes) const {
    Vector2f dstdx(su * si->dudx, sv * si->dvdx);
    Vector2f dstdy(su * si->dudy, sv * si->dvdy);
    Point2f st(su * si->uv[0], sv * si->uv[1]);

    Float eta = index->Evaluate(*si);

    // Retrieve the parameters of the BSDF
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

    Spectrum R = Kr->Evaluate(*si).Clamp();
    Spectrum T = Kt->Evaluate(*si).Clamp();
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, eta);

    RNG *rngEval = ARENA_ALLOC(arena, RNG)();
    RNG *rngSample = ARENA_ALLOC(arena, RNG)();

    Float u = st.x * 457. + 269. * st.y + 599. * dstdx.x + 397. * dstdx.y +
              29. * dstdy.x + 283. * dstdy.y;
    Float du = std::abs(u) - std::floor(std::abs(u));
    rngSample->SetSequence(du * std::numeric_limits<int>::max());

    if (R.IsBlack() && T.IsBlack()) return;

    bool isSpecular = alphaXP == 0 && alphaYP == 0;
    if (isSpecular && allowMultipleLobes) {
        si->bsdf->Add(
            ARENA_ALLOC(arena, FresnelSpecular)(R, T, 1.f, eta, mode));
    } else {
        si->bsdf->Add(ARENA_ALLOC(arena, GlitteringDielectricScattering)(
            R, T, 1.f, eta, mode, rngEval, rngSample, st, dstdx, dstdy, alphaXP,
            alphaYP, rhoP, alphaXBaseMaterialP, alphaYBaseMaterialP,
            rhoBaseMaterialP, distResolution, distributions.get(), nLevels, N,
            alphaDict, logMicrofacetDensityP, microfacetRelativeAreaP,
            densityRandomisation, sampleVisibleArea, sampleApproximation));
    }
}

GlitteringDielectricMaterial *CreateGlitteringDielectricMaterial(
    const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> Kr =
        mp.GetSpectrumTexture("Kr", Spectrum(1.f));
    std::shared_ptr<Texture<Spectrum>> Kt =
        mp.GetSpectrumTexture("Kt", Spectrum(1.f));
    std::shared_ptr<Texture<Float>> eta = mp.GetFloatTextureOrNull("eta");
    if (!eta) eta = mp.GetFloatTexture("index", 1.5f);

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

    bool remapParameters = mp.FindBool("remapparameters", false);
    bool sampleVisibleArea = mp.FindBool("samplevisiblearea", true);
    bool sampleApproximation = mp.FindBool("sampleapproximation", false);

    return new GlitteringDielectricMaterial(
        Kr, Kt, alphaX, alphaY, rho, alphaXBaseMaterial, alphaYBaseMaterial,
        rhoBaseMaterial, dictionary, nLevels, N, alphaDict,
        logMicrofacetDensity, microfacetRelativeArea, densityRandomisation, su,
        sv, remapParameters, sampleVisibleArea, sampleApproximation, eta);
}

Spectrum GlitteringDielectricScattering::f(const Vector3f &wo,
                                           const Vector3f &wi) const {
    if (SameHemisphere(wo, wi)) {
        // Reflection
        Vector3f woZPos = wo;
        Vector3f wiZPos = wi;
        // Flip the directions
        if (wo.z < 0.) {
            woZPos = -wo;
            wiZPos = -wi;
        }

        Float cosThetaO = CosTheta(woZPos), cosThetaI = CosTheta(wiZPos);
        Vector3f wh = wiZPos + woZPos;
        // Handle degenerate cases for microfacet reflection
        if (cosThetaI <= 0 || cosThetaO <= 0) return Spectrum(0.);
        wh = Normalize(wh);
        if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);
        if (Dot(woZPos, wh) <= 0.f || Dot(wiZPos, wh) <= 0.f) return 0.f;

        Float F = FrDielectric(wo.z, etaA, etaB);

        Fresnel *fr = new FresnelNoOp();
        GlitteringConductorReflection glitteringConductorReflection(
            R, fr, rngEval, rngSample, st, dstdx, dstdy,
            alphaXGlitteringMaterial, alphaYGlitteringMaterial,
            rhoGlitteringMaterial, alphaXBaseMaterial, alphaYBaseMaterial,
            rhoBaseMaterial, distResolution, distributions, nLevels, N,
            alphaDict, logMicrofacetDensity, microfacetRelativeArea,
            densityRandomisation, sampleVisibleArea, sampleApproximation);
        Float DValue = glitteringConductorReflection.D_P(wh, st, dstdx, dstdy);
        delete fr;
        CHECK_GE(DValue, 0.);

        Float G1wowh = G1VCavity(wh, woZPos);
        Float G1wiwh = G1VCavity(wh, wiZPos);
        Float G = G1wowh * G1wiwh;
        CHECK_GE(G, 0.);
        CHECK_LE(G, 1.);

        return R * DValue * G * F / (4 * cosThetaI * cosThetaO);
    } else {
        // Transmission

        Float cosThetaO = CosTheta(wo);
        Float cosThetaI = CosTheta(wi);
        if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0);

        Float eta = CosTheta(wo) > 0 ? (etaB / etaA) : (etaA / etaB);
        Vector3f wh = Normalize(wo + wi * eta);

        if (wh.z == 0.) return 0.;
        if (wh.z < 0) wh = -wh;

        // Same side?
        if (Dot(wo, wh) * Dot(wi, wh) > 0) return Spectrum(0);

        Float F = FrDielectric(wo.z, etaA, etaB);

        Fresnel *fr = new FresnelNoOp();
        GlitteringConductorReflection glitteringConductorReflection(
            R, fr, rngEval, rngSample, st, dstdx, dstdy,
            alphaXGlitteringMaterial, alphaYGlitteringMaterial,
            rhoGlitteringMaterial, alphaXBaseMaterial, alphaYBaseMaterial,
            rhoBaseMaterial, distResolution, distributions, nLevels, N,
            alphaDict, logMicrofacetDensity, microfacetRelativeArea,
            densityRandomisation, sampleVisibleArea, sampleApproximation);
        Float DValue = glitteringConductorReflection.D_P(wh, st, dstdx, dstdy);
        delete fr;
        CHECK_GE(DValue, 0.);

        Float G1wowh = G1VCavity(wh, wo.z < 0. ? -wo : wo);
        Float G1wiwh = G1VCavity(wh, wi.z < 0. ? -wi : wi);
        Float G = G1wowh * G1wiwh;
        CHECK_GE(G, 0.);
        CHECK_LE(G, 1.);

        Float sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
        // mode == TransportMode::Importance
        Float factor = 1.f;

        return (Spectrum(1.f) - F) * T *
               std::abs(DValue * G * eta * eta * AbsDot(wi, wh) *
                        AbsDot(wo, wh) * factor * factor /
                        (cosThetaI * cosThetaO * sqrtDenom * sqrtDenom));
    }
}

Spectrum GlitteringDielectricScattering::Sample_f(const Vector3f &wo,
                                                  Vector3f *wi,
                                                  const Point2f &u, Float *pdf,
                                                  BxDFType *sampledType) const {
    if (wo.z == 0.) return 0.;
    //-------------------------------------------------------------------------
    Fresnel *fr = new FresnelNoOp();
    GlitteringConductorReflection glitteringConductorReflection(
        R, fr, rngEval, rngSample, st, dstdx, dstdy, alphaXGlitteringMaterial,
        alphaYGlitteringMaterial, rhoGlitteringMaterial, alphaXBaseMaterial,
        alphaYBaseMaterial, rhoBaseMaterial, distResolution, distributions,
        nLevels, N, alphaDict, logMicrofacetDensity, microfacetRelativeArea,
        densityRandomisation, sampleVisibleArea, sampleApproximation);
    Vector3f wm = glitteringConductorReflection.Sample_D_P();
    delete fr;

    Vector3f wmP(-wm.x, -wm.y, wm.z);

    bool flip = wo.z < 0.;

    Float dotWOWMP = std::max(Dot(flip ? -wo : wo, wmP), Float(0.));
    Float dotWOWM = std::max(Dot(flip ? -wo : wo, wm), Float(0.));
    Float deno = dotWOWM + dotWOWMP;
    Float projectedAreaWMP = dotWOWMP / deno;
    Vector3f wh;
    Float uFacet = rngSample->UniformFloat();
    if (uFacet < projectedAreaWMP && sampleVisibleArea)
        wh = wmP;
    else
        wh = wm;

    if (Dot(flip ? -wo : wo, wh) < 0) return 0.;  // Should be rare
    //-------------------------------------------------------------------------

    Float F = FrDielectric(wo.z, etaA, etaB);
    if (rngSample->UniformFloat() < F) {
        // Reflection
        Vector3f wiTmp = Reflect(flip ? -wo : wo, wh);
        *wi = flip ? -wiTmp : wiTmp;
        if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);
        if (Dot(wiTmp, wh) <= 0.f) return 0.f;
        if (sampledType) *sampledType = BxDFType(BSDF_GLOSSY | BSDF_REFLECTION);

        *pdf = Pdf(wo, *wi);
        return f(wo, *wi);
    } else {
        bool entering = CosTheta(wo) > 0;
        Float etaI = entering ? etaA : etaB;
        Float etaT = entering ? etaB : etaA;

        if (!Refract(wo, Faceforward(Normal3f(wh.x, wh.y, wh.z), wo),
                     etaI / etaT, wi))
            return 0;

        if (Dot(wo, wh) * Dot(*wi, wh) >= 0.) return 0;

        if (sampledType)
            *sampledType = BxDFType(BSDF_GLOSSY | BSDF_TRANSMISSION);

        *pdf = Pdf(wo, *wi);
        return f(wo, *wi);
    }
}

Float GlitteringDielectricScattering::Pdf(const Vector3f &wo,
                                          const Vector3f &wi) const {
    if (SameHemisphere(wo, wi)) {
        // Reflection
        Vector3f woZPos = wo;
        Vector3f wiZPos = wi;
        // Flip the directions
        if (wo.z < 0.) {
            woZPos = -wo;
            wiZPos = -wi;
        }

        Float cosThetaO = CosTheta(woZPos), cosThetaI = CosTheta(wiZPos);
        Vector3f wh = wiZPos + woZPos;
        // Handle degenerate cases for microfacet reflection
        if (cosThetaI <= 0 || cosThetaO <= 0) return 0.;
        wh = Normalize(wh);
        if (wh.x == 0 && wh.y == 0 && wh.z == 0) return 0.;
        Float dotWOWH = Dot(woZPos, wh);
        if (dotWOWH <= 0.f || Dot(wiZPos, wh) <= 0.f) return 0.f;

        Float F = FrDielectric(wo.z, etaA, etaB);

        Fresnel *fr = new FresnelNoOp();
        GlitteringConductorReflection glitteringConductorReflection(
            R, fr, rngEval, rngSample, st, dstdx, dstdy,
            alphaXGlitteringMaterial, alphaYGlitteringMaterial,
            rhoGlitteringMaterial, alphaXBaseMaterial, alphaYBaseMaterial,
            rhoBaseMaterial, distResolution, distributions, nLevels, N,
            alphaDict, logMicrofacetDensity, microfacetRelativeArea,
            densityRandomisation, sampleVisibleArea, sampleApproximation);
        Float PDF_D = glitteringConductorReflection.Pdf(woZPos, wiZPos);
        delete fr;
        CHECK_GE(PDF_D, 0.);

        return PDF_D * F;
    } else {
        // Transmission

        Float cosThetaO = CosTheta(wo);
        Float cosThetaI = CosTheta(wi);
        if (cosThetaI == 0 || cosThetaO == 0) return 0.;

        Float eta = CosTheta(wo) > 0 ? (etaB / etaA) : (etaA / etaB);
        Vector3f wh = Normalize(wo + wi * eta);

        if (wh.z == 0.) return 0.;
        if (wh.z < 0) wh = -wh;

        // Same side?
        if (Dot(wo, wh) * Dot(wi, wh) > 0) return 0.;

        Float F = FrDielectric(CosTheta(wo), etaA, etaB);

        Fresnel *fr = new FresnelNoOp();
        GlitteringConductorReflection glitteringConductorReflection(
            R, fr, rngEval, rngSample, st, dstdx, dstdy,
            alphaXGlitteringMaterial, alphaYGlitteringMaterial,
            rhoGlitteringMaterial, alphaXBaseMaterial, alphaYBaseMaterial,
            rhoBaseMaterial, distResolution, distributions, nLevels, N,
            alphaDict, logMicrofacetDensity, microfacetRelativeArea,
            densityRandomisation, sampleVisibleArea, sampleApproximation);

        Float pdfD;
        Float dotWOWM = Dot(wo.z < 0 ? -wo : wo, wh);
        if (sampleVisibleArea) {
            if (!sampleApproximation)
                pdfD = glitteringConductorReflection.Visible_D_P(
                    wh, wo.z < 0 ? -wo : wo);
            else {
                Vector2f slopeH(-wh.x / wh.z, -wh.y / wh.z);
                Float D_PValue =
                    glitteringConductorReflection.SDFLastLOD(slopeH) /
                    (wh.z * wh.z * wh.z * wh.z);
                Float density = G1VCavity(wh, wo.z < 0 ? -wo : wo) * dotWOWM *
                                D_PValue / std::abs(wo.z);
                Float visible_D_PApprox = density;
                pdfD = visible_D_PApprox;
            }
        } else {
            if (!sampleApproximation)
                pdfD = glitteringConductorReflection.D_P(wh, st, dstdx, dstdy) *
                       wh.z;
            else {
                Vector2f slopeH(-wh.x / wh.z, -wh.y / wh.z);
                Float D_PValueApprox =
                    glitteringConductorReflection.SDFLastLOD(slopeH) /
                    (wh.z * wh.z * wh.z * wh.z);
                pdfD = D_PValueApprox * wh.z;
            }
        }
        delete fr;

        Float sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
        // mode == TransportMode::Importance
        Float dwh_dwi =
            std::abs((eta * eta * Dot(wi, wh)) / (sqrtDenom * sqrtDenom));

        Float pdf = pdfD * dwh_dwi * (1. - F);
        CHECK_GE(pdf, 0.);

        return pdf;
    }
}

std::string GlitteringDielectricScattering::ToString() const {
    return std::string("[ FresnelSpecular R: ") + R.ToString() +
           std::string(" T: ") + T.ToString() +
           StringPrintf(" etaA: %f etaB: %f ", etaA, etaB) +
           std::string(" mode : ") +
           (mode == TransportMode::Radiance ? std::string("RADIANCE")
                                            : std::string("IMPORTANCE")) +
           std::string(" ]");
}

}  // namespace pbrt
