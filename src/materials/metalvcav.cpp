
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

#include "materials/metalvcav.h"

#include "interaction.h"
#include "materials/glitteringconductor.h"
#include "paramset.h"
#include "reflection.h"
#include "texture.h"

namespace pbrt {

// MetalVCavMaterial Method Definitions
MetalVCavMaterial::MetalVCavMaterial(
    const std::shared_ptr<Texture<Spectrum>> &R,
    const std::shared_ptr<Texture<Spectrum>> &eta,
    const std::shared_ptr<Texture<Spectrum>> &k,
    const std::shared_ptr<Texture<Float>> &alphaX,
    const std::shared_ptr<Texture<Float>> &alphaY,
    const std::shared_ptr<Texture<Float>> &rho, bool fresnelNoOp)
    : R(R),
      eta(eta),
      k(k),
      alphaX(alphaX),
      alphaY(alphaY),
      rho(rho),
      fresnelNoOp(fresnelNoOp) {}

void MetalVCavMaterial::ComputeScatteringFunctions(
    SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
    bool allowMultipleLobes) const {
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);

    Spectrum RValue = R->Evaluate(*si).Clamp();
    Float alphaXP = alphaX->Evaluate(*si);
    Float alphaYP = alphaY->Evaluate(*si);
    Float rhoP = rho->Evaluate(*si);

    CHECK_GT(rhoP, -1.) << "Slope correlation factor must be between (-1, 1)";
    CHECK_LT(rhoP, 1.) << "Slope correlation factor must be between (-1, 1)";

    Fresnel *frMf;
    if (fresnelNoOp)
        frMf = ARENA_ALLOC(arena, FresnelNoOp)();
    else
        frMf = ARENA_ALLOC(arena, FresnelConductor)(1., eta->Evaluate(*si),
                                                    k->Evaluate(*si));

    RNG *rngSample = ARENA_ALLOC(arena, RNG)();
    Float u = si->uv.x * 457. + 101. * si->uv.y + 857. * si->dudx + 263. * si->dudy +
              311. * si->dvdx + 149. * si->dvdx;
    Float du = std::abs(u) - std::floor(std::abs(u));
    rngSample->SetSequence(du * std::numeric_limits<int>::max());

    si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetVCavReflection)(
        RValue, frMf, alphaXP, alphaYP, rhoP, rngSample));
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

MetalVCavMaterial *CreateMetalVCavMaterial(const TextureParams &mp) {
    static Spectrum copperN =
        Spectrum::FromSampled(CopperWavelengths, CopperN, CopperSamples);
    std::shared_ptr<Texture<Spectrum>> eta =
        mp.GetSpectrumTexture("eta", copperN);
    static Spectrum copperK =
        Spectrum::FromSampled(CopperWavelengths, CopperK, CopperSamples);
    std::shared_ptr<Texture<Spectrum>> k = mp.GetSpectrumTexture("k", copperK);
    std::shared_ptr<Texture<Spectrum>> R = mp.GetSpectrumTexture("Rs", 1.);
    std::shared_ptr<Texture<Float>> alphaX = mp.GetFloatTexture("alphax", 0.5);
    std::shared_ptr<Texture<Float>> alphaY = mp.GetFloatTexture("alphay", 0.5);
    std::shared_ptr<Texture<Float>> rho = mp.GetFloatTexture("rho", 0.);
    bool fresnelNoOp = mp.FindBool("fresnelnoop", false);

    return new MetalVCavMaterial(R, eta, k, alphaX, alphaY, rho, fresnelNoOp);
}

Spectrum MicrofacetVCavReflection::f(const Vector3f &wo,
                                     const Vector3f &wi) const {
    Float cosThetaO = CosTheta(wo), cosThetaI = CosTheta(wi);
    Vector3f wh = wi + wo;
    // Handle degenerate cases for microfacet reflection
    if (cosThetaI <= 0 || cosThetaO <= 0) return Spectrum(0.);
    if (wh.x == 0 && wh.y == 0 && wh.z <= 0) return Spectrum(0.);
    wh = Normalize(wh);
    // For the Fresnel call, make sure that wh is in the same hemisphere
    // as the surface normal, so that TIR is handled correctly.
    Spectrum F = fresnel->Evaluate(Dot(wi, Faceforward(wh, Vector3f(0, 0, 1))));

    // Local masking shadowing
    if (Dot(wo, wh) <= 0. || Dot(wi, wh) <= 0.) return Spectrum(0.);
    Float G1wowh = G1VCavity(wh, wo);
    Float G1wiwh = G1VCavity(wh, wi);

    Vector2f slope_h(-wh.x / wh.z, -wh.y / wh.z);
    Float slopeDensity = BivariateNormalDistribution(
        slope_h.x, slope_h.y, alphaX / Sqrt2, alphaY / Sqrt2, rho);
    Float DValue = slopeDensity / (wh.z * wh.z * wh.z * wh.z);

    return R * DValue * G1wowh * G1wiwh * F / (4 * cosThetaI * cosThetaO);
}

Spectrum MicrofacetVCavReflection::Sample_f(const Vector3f &wo, Vector3f *wi,
                                            const Point2f &u, Float *pdf,
                                            BxDFType *sampledType) const {
    // Sample microfacet orientation $\wh$ and reflected direction $\wi$
    if (wo.z <= 0) return 0.;
    Vector2f slopem = SampleBivariateNormalDistribution(
        Point2f(rngSample->UniformFloat(), rngSample->UniformFloat()),
        alphaX / Sqrt2, alphaY / Sqrt2, rho);
    Vector3f wm = SlopeToNormal(slopem);

    Vector3f wmP(-wm.x, -wm.y, wm.z);
    Float dotWOWMP = std::max(Dot(wo, wmP), Float(0.));
    Float dotWOWM = std::max(Dot(wo, wm), Float(0.));
    Float deno = dotWOWM + dotWOWMP;
    Float projectedAreaWMP = dotWOWMP / deno;
    Vector3f wh;
    if (rngSample->UniformFloat() < projectedAreaWMP)
        wh = wmP;
    else
        wh = wm;

    if (Dot(wo, wh) < 0) return 0.;  // Should be rare
    *wi = Reflect(wo, wh);
    if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);

    Float dotWOWH = Dot(wo, wh);
    if (dotWOWH <= 0. || wo.z <= 0. || wh.z <= 0.) return 0.;
    Float D_Value =
        BivariateNormalDistribution(-wh.x / wh.z, -wh.y / wh.z, alphaX / Sqrt2,
                                    alphaY / Sqrt2, rho) /
        (wh.z * wh.z * wh.z * wh.z);

    Float visibleD = G1VCavity(wh, wo) * dotWOWH * D_Value / wo.z;
    *pdf = visibleD / (4 * Dot(wo, wh));

    Spectrum F =
        fresnel->Evaluate(Dot(*wi, Faceforward(wh, Vector3f(0, 0, 1))));
    Spectrum BRDFValue =
        R * visibleD * G1VCavity(wh, *wi) * F / (4. * Dot(wo, wh) * wi->z);
    return BRDFValue;
}

Float MicrofacetVCavReflection::Pdf(const Vector3f &wo,
                                    const Vector3f &wi) const {
    if (!SameHemisphere(wo, wi)) return 0;
    Vector3f wh = Normalize(wo + wi);
    Float dotWOWH = Dot(wo, wh);
    if (dotWOWH <= 0. || wo.z <= 0. || wh.z <= 0.) return 0.;
    Float D_Value =
        BivariateNormalDistribution(-wh.x / wh.z, -wh.y / wh.z, alphaX / Sqrt2,
                                    alphaY / Sqrt2, rho) /
        (wh.z * wh.z * wh.z * wh.z);
    Float visibleD = G1VCavity(wh, wo) * dotWOWH * D_Value / wo.z;
    return visibleD / (4 * Dot(wo, wh));
}

std::string MicrofacetVCavReflection::ToString() const {
    return std::string("[ MicrofacetVCavReflection R: ") + R.ToString() +
           std::string(" fresnel: ") + fresnel->ToString() + std::string(" ]");
}

}  // namespace pbrt
