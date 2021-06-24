
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

#ifndef PBRT_MATERIALS_METALVCAV_H
#define PBRT_MATERIALS_METALVCAV_H

#include "material.h"
#include "pbrt.h"
#include "reflection.h"
#include "spectrum.h"

namespace pbrt {

class MetalVCavMaterial : public Material {
  public:
    MetalVCavMaterial(const std::shared_ptr<Texture<Spectrum>> &R,
                      const std::shared_ptr<Texture<Spectrum>> &eta,
                      const std::shared_ptr<Texture<Spectrum>> &k,
                      const std::shared_ptr<Texture<Float>> &alphaX,
                      const std::shared_ptr<Texture<Float>> &alphaY,
                      const std::shared_ptr<Texture<Float>> &rho,
                      bool fresnelNoOp);
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    std::shared_ptr<Texture<Spectrum>> eta, k, R;
    std::shared_ptr<Texture<Float>> alphaX, alphaY, rho;
    bool fresnelNoOp;
};

MetalVCavMaterial *CreateMetalVCavMaterial(const TextureParams &mp);

class MicrofacetVCavReflection : public BxDF {
  public:
    MicrofacetVCavReflection(const Spectrum &R, Fresnel *fresnel, Float alphaX,
                             Float alphaY, Float rho, RNG *rngSample)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
          R(R),
          alphaX(alphaX),
          alphaY(alphaY),
          rho(rho),
          fresnel(fresnel),
          rngSample(rngSample) {}

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // MicrofacetReflection Private Data
    const Spectrum R;
    const Fresnel *fresnel;
    const Float alphaX, alphaY, rho;
    RNG *rngSample;
};

}  // namespace pbrt

#endif  // PBRT_MATERIALS_METALVCAV_H
