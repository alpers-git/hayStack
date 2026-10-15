// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hayStack/HayStack.h"
#include "hayStack/TransferFunction.h"
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace hs {

  /*! the default colour map for a signed field: a diverging ramp that is
      *transparent at zero*, so vacuum stays out of the way and the two phases
      of an orbital lobe read as opposite colours. the usual one-sided default
      makes the whole domain opaque, which for a field that is near zero almost
      everywhere just fills the screen. */
  inline std::vector<vec4f> makeSignedColorMap(int N = 128, float alphaGamma = 2.f,
                                              float peakAlpha = 1.f)
  {
    const vec3f neg(0.23f,0.30f,0.75f);
    const vec3f mid(0.87f,0.87f,0.87f);
    const vec3f pos(0.71f,0.02f,0.15f);
    std::vector<vec4f> cm(N);
    for (int i=0;i<N;i++) {
      const float t = (N == 1) ? 0.f : i/float(N-1);
      const float v = 2.f*t - 1.f;
      const float a = std::fabs(v);
      const vec3f c = (v < 0.f)
        ? (neg + (mid-neg)*(1.f+v))
        : (mid + (pos-mid)*v);
      cm[i] = vec4f(c.x,c.y,c.z,peakAlpha*std::pow(a,alphaGamma));
    }
    return cm;
  }

  /*! the transfer function an orbital field should start from: diverging,
      transparent at zero, symmetric domain, and scaled so that crossing the
      whole box at full classification is about one optical depth rather than
      the opaque wall a one-sided default produces. */
  inline TransferFunction makeSignedTransferFunction(const range1f &valueRange,
                                                     const box3f &bounds,
                                                     float alphaGamma = 2.f,
                                                     float opacityScale = 1.f,
                                                     float domainScale = 0.15f)
  {
    TransferFunction xf;
    const float hi = std::fabs(valueRange.upper);
    const float lo = std::fabs(valueRange.lower);
    /* a two-sided field needs the domain pulled well inside its peak, because
       the peak sits on a nucleus and would otherwise take all the contrast; a
       one-sided field varies across its whole range, so it keeps all of it */
    const bool oneSided = std::min(lo,hi) <= 0.05f*std::max(lo,hi);
    const float a = (oneSided ? 1.f : domainScale) * std::max(lo,hi);
    const float A = (a > 0.f) ? a : 1.f;
    xf.relDomain = range1f(0.f,100.f);
    xf.baseDensity = 100.f;
    (void)bounds;

    /* a field that never changes sign - an electron density, say - gets a
       one-sided ramp: a diverging map would put all of it on one half and read
       as a single flat colour */
    if (!oneSided) {
      xf.absDomain = range1f(-A,+A);
      xf.colorMap = makeSignedColorMap(128,alphaGamma,opacityScale);
    } else {
      const bool negative = lo > hi;
      xf.absDomain = negative ? range1f(-A,0.f) : range1f(0.f,A);
      const vec3f end = negative ? vec3f(0.23f,0.30f,0.75f)
                                 : vec3f(0.71f,0.02f,0.15f);
      const vec3f start(0.95f,0.90f,0.55f);
      const int N = 128;
      xf.colorMap.resize(N);
      for (int i=0;i<N;i++) {
        const float t = i/float(N-1);
        /* reversed for a negative field, so the far end of the domain is the
           strong colour either way */
        const float u = negative ? (1.f-t) : t;
        const vec3f c = start + (end-start)*u;
        xf.colorMap[i]
          = vec4f(c.x,c.y,c.z,opacityScale*std::pow(u,alphaGamma));
      }
    }
    return xf;
  }

  /*! what the viewer pushes down when the user changes the orbital-channel
      selection. Barney sums the enabled channels per voxel and rebuilds the
      macro-cell ranges and majorants from that, on the set. */
  struct ChannelFilter {
    /*! bit c enables orbital channel c */
    uint32_t channelMask = 0xffffffffu;
    /*! set when channelMask changed and needs pushing */
    bool     channelDirty = true;
  };

  /*! the orbital-current field on a dense regular grid.

    One current vector per channel per voxel; the enabled channels are summed
    and the magnitude is taken at the sample point. Vectors rather than
    magnitudes because |sum_c j_c| != sum_c |j_c| wherever channels oppose. */
  struct OrbGridVolume {
    typedef std::shared_ptr<OrbGridVolume> SP;
    static SP create() { return std::make_shared<OrbGridVolume>(); }

    /*! numChannels * dims.x*dims.y*dims.z vec3, channel-major, x fastest */
    std::vector<vec3f>    channelVec;
    vec3i                 dims { 0,0,0 };
    vec3f                 gridOrigin  { 0,0,0 };
    vec3f                 gridSpacing { 1,1,1 };

    int                   numChannels = 1;
    uint32_t              channelMask = 0xffffffffu;
    std::vector<std::string> channelNames;
    /*! peak |j| of each channel on its own. These span a factor of ten, so a
        viewer that keeps one fixed transfer-function domain across a channel
        toggle draws the weak ones nearly invisible - which reads as the filter
        being broken rather than the channel being weak. */
    std::vector<float>    channelPeak;
    /*! exact peak of |sum_{c in S} j_c| per channel subset, indexed by the
        mask. Summing channelPeak instead overshoots wherever channels cancel. */
    std::vector<float>    subsetPeak;

    box3f                 bounds;
    range1f               valueRange { 0.f, 0.f };
    float alphaGamma   = 1.f;
    float opacityScale = 1.f;
    float domainScale  = 1.f;

    size_t numVoxels() const
    { return (size_t)dims.x*(size_t)dims.y*(size_t)dims.z; }
  };

}
