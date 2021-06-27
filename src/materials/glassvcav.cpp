
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

#include "materials/glassvcav.h"

#include "interaction.h"
#include "materials/glitteringconductor.h"
#include "paramset.h"
#include "reflection.h"
#include "spectrum.h"
#include "texture.h"

namespace pbrt {

void GlassVCavMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
    bool allowMultipleLobes) const {
    Float eta = index->Evaluate(*si);
    Float alphaXP = alphaX->Evaluate(*si);
    Float alphaYP = alphaY->Evaluate(*si);
    Float rhoP = rho->Evaluate(*si);
    Spectrum R = Kr->Evaluate(*si).Clamp();
    Spectrum T = Kt->Evaluate(*si).Clamp();
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, eta);

    RNG *rngSample = ARENA_ALLOC(arena, RNG)();
    Vector2f dstdx(si->dudx, si->dvdx);
    Vector2f dstdy(si->dudy, si->dvdy);
    Point2f st(si->uv[0], si->uv[1]);
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
        si->bsdf->Add(ARENA_ALLOC(arena, FresnelMicrofacetScattering)(
            R, T, 1.f, eta, mode, rngSample, alphaXP, alphaYP, rhoP));
    }
}

GlassVCavMaterial *CreateGlassVCavMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> Kr =
        mp.GetSpectrumTexture("Kr", Spectrum(1.f));
    std::shared_ptr<Texture<Spectrum>> Kt =
        mp.GetSpectrumTexture("Kt", Spectrum(1.f));
    std::shared_ptr<Texture<Float>> eta = mp.GetFloatTextureOrNull("eta");
    if (!eta) eta = mp.GetFloatTexture("index", 1.5f);
    std::shared_ptr<Texture<Float>> alphaX = mp.GetFloatTexture("alphax", 0.1f);
    std::shared_ptr<Texture<Float>> alphaY = mp.GetFloatTexture("alphay", 0.1f);
    std::shared_ptr<Texture<Float>> rho = mp.GetFloatTexture("rho", 0.f);
    return new GlassVCavMaterial(Kr, Kt, alphaX, alphaY, rho, eta);
}

Spectrum FresnelMicrofacetScattering::f(const Vector3f &wo,
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

        Vector2f slope_h(-wh.x / wh.z, -wh.y / wh.z);
        Float slopeDensity = BivariateNormalDistribution(
            slope_h.x, slope_h.y, alphaX / Sqrt2, alphaY / Sqrt2, rho);
        Float DValue = slopeDensity / (wh.z * wh.z * wh.z * wh.z);
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

        Vector2f slope_h(-wh.x / wh.z, -wh.y / wh.z);
        Float slopeDensity = BivariateNormalDistribution(
            slope_h.x, slope_h.y, alphaX / Sqrt2, alphaY / Sqrt2, rho);
        Float DValue = slopeDensity / (wh.z * wh.z * wh.z * wh.z);
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

Spectrum FresnelMicrofacetScattering::Sample_f(const Vector3f &wo, Vector3f *wi,
                                               const Point2f &u, Float *pdf,
                                               BxDFType *sampledType) const {
    if (wo.z == 0.) return 0.;
    //-------------------------------------------------------------------------
    // Sample the visible normal distribution function
    Vector2f slopem = SampleBivariateNormalDistribution(
        Point2f(rngSample->UniformFloat(), rngSample->UniformFloat()),
        alphaX / Sqrt2, alphaY / Sqrt2, rho);
    Vector3f wm = SlopeToNormal(slopem);

    Vector3f wmP(-wm.x, -wm.y, wm.z);

    bool flip = wo.z < 0.;

    Float dotWOWMP = std::max(Dot(flip ? -wo : wo, wmP), Float(0.));
    Float dotWOWM = std::max(Dot(flip ? -wo : wo, wm), Float(0.));
    Float deno = dotWOWM + dotWOWMP;
    Float projectedAreaWMP = dotWOWMP / deno;
    Vector3f wh;
    Float uFacet = rngSample->UniformFloat();
    if (uFacet < projectedAreaWMP)
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

Float FresnelMicrofacetScattering::Pdf(const Vector3f &wo,
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

        Vector2f slope_h(-wh.x / wh.z, -wh.y / wh.z);
        Float slopeDensity = BivariateNormalDistribution(
            slope_h.x, slope_h.y, alphaX / Sqrt2, alphaY / Sqrt2, rho);
        Float DValue = slopeDensity / (wh.z * wh.z * wh.z * wh.z);
        CHECK_GE(DValue, 0.);

        Float visibleD = G1VCavity(wh, woZPos) * dotWOWH * DValue / woZPos.z;

        return visibleD * F / (4 * dotWOWH);
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

        Vector2f slope_h(-wh.x / wh.z, -wh.y / wh.z);
        Float slopeDensity = BivariateNormalDistribution(
            slope_h.x, slope_h.y, alphaX / Sqrt2, alphaY / Sqrt2, rho);
        Float DValue = slopeDensity / (wh.z * wh.z * wh.z * wh.z);
        CHECK_GE(DValue, 0.);

        Float visibleD = G1VCavity(wh, wo.z < 0 ? -wo : wo) *
                         Dot(wo.z < 0 ? -wo : wo, wh) * DValue / std::abs(wo.z);

        Float sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
        // mode == TransportMode::Importance
        Float dwh_dwi =
            std::abs((eta * eta * Dot(wi, wh)) / (sqrtDenom * sqrtDenom));

        return visibleD * dwh_dwi * (1. - F);
    }
}

std::string FresnelMicrofacetScattering::ToString() const {
    return std::string("[ FresnelSpecular R: ") + R.ToString() +
           std::string(" T: ") + T.ToString() +
           StringPrintf(" etaA: %f etaB: %f ", etaA, etaB) +
           std::string(" mode : ") +
           (mode == TransportMode::Radiance ? std::string("RADIANCE")
                                            : std::string("IMPORTANCE")) +
           std::string(" ]");
}

}  // namespace pbrt
