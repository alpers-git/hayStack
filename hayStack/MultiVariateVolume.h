// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hayStack/HayStack.h"
#include "hayStack/TransferFunction.h"
#include <cmath>
#include <limits>
#include <string>
#include <vector>

/*! \file MultiVariateVolume.h multi-variate volumes

  A multi-variate volume is ONE field on a dense regular grid whose voxels
  carry several co-located variables ("channels"), each a scalar or a
  3-vector. The renderer mixes the channels with per-channel weights and
  reduces the mix to a single scalar - its magnitude, one component, or its
  dot product with a direction - before one transfer function classifies it.

  This is deliberately a different mode from loading several volume files
  into the same space. Those stay independent volumes, each with its own
  transfer function, composited after classification; here the variables
  combine before classification, so they can reinforce or cancel.

  Rendered by barney's "multiVariateGrid" spatial field. Loaders fill this
  struct from their own file formats (see OrbGridFromFile for .orbg).
*/

namespace hs {

  enum class Reduction : int { Magnitude = 0, Component = 1, Dot = 2 };

  inline const char *reductionName(Reduction r)
  {
    switch (r) {
    case Reduction::Component: return "component";
    case Reduction::Dot:       return "dot";
    default:                   return "magnitude";
    }
  }

  inline bool parseReduction(const std::string &s, Reduction &r)
  {
    for (int i=0;i<3;i++)
      if (s == reductionName(Reduction(i))) { r = Reduction(i); return true; }
    return false;
  }

  /*! how the channels are combined; the viewer pushes it whenever it changes */
  struct ChannelMix {
    std::vector<float> weights;
    Reduction reduction = Reduction::Magnitude;
    int       component = 0;
    vec3f     direction { 1.f,0.f,0.f };
    bool      dirty     = true;

    /*! whether the reduced value can take both signs */
    bool isSigned() const { return reduction != Reduction::Magnitude; }
  };

  struct MultiVariateVolume {
    typedef std::shared_ptr<MultiVariateVolume> SP;
    static SP create() { return std::make_shared<MultiVariateVolume>(); }

    /*! numChannels*numVoxels*numComponents floats, channel-major, x fastest */
    std::vector<float>       channelData;
    vec3i                    dims        { 0,0,0 };
    vec3f                    gridOrigin  { 0,0,0 };
    vec3f                    gridSpacing { 1,1,1 };
    int                      numChannels   = 1;
    int                      numComponents = 1;
    std::vector<std::string> channelNames;

    /*! optional: peak |channel| per channel, and the exact peak magnitude of
        every 0/1 channel subset, indexed by its bit mask */
    std::vector<float>       channelPeak;
    std::vector<float>       subsetPeak;
    std::string              unit;

    /*! the mix the volume starts with */
    ChannelMix               mix;

    box3f   bounds;
    range1f valueRange { 0.f, 0.f };
    float   alphaGamma   = 1.f;
    float   opacityScale = 1.f;
    float   domainScale  = 1.f;

    size_t numVoxels() const
    { return (size_t)dims.x*(size_t)dims.y*(size_t)dims.z; }

    /*! upper bound on |reduced value| under a mix; exact when all enabled
        channels share one weight and subsetPeak covers them, 0 if unknown */
    float peakFor(const ChannelMix &m) const
    {
      float scale = 1.f;
      if (m.reduction == Reduction::Dot)
        scale = std::sqrt(dot(m.direction,m.direction));
      uint64_t mask = 0;
      float shared = 0.f;
      bool uniform = numChannels <= 32 && (int)m.weights.size() >= numChannels;
      for (int c=0;uniform && c<numChannels;c++) {
        const float w = m.weights[c];
        if (w == 0.f) continue;
        if (shared == 0.f) shared = w;
        uniform = (w == shared);
        mask |= (1ull<<c);
      }
      if (uniform && mask < subsetPeak.size())
        return scale*std::fabs(shared)*subsetPeak[mask];
      float sum = 0.f;
      for (int c=0;c<numChannels && c<(int)channelPeak.size();c++)
        sum += std::fabs(c < (int)m.weights.size() ? m.weights[c] : 1.f)
          * channelPeak[c];
      return scale*sum;
    }
  };

  /*! diverging colour map that is transparent at zero; N should be odd so
      that zero lands exactly on the transparent centre entry */
  inline std::vector<vec4f> makeSignedColorMap(int N = 129, float alphaGamma = 2.f,
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

  /*! starting transfer function: symmetric and diverging for a two-sided
      field, a one-sided ramp otherwise */
  inline TransferFunction makeSignedTransferFunction(const range1f &valueRange,
                                                     const box3f &bounds,
                                                     float alphaGamma = 2.f,
                                                     float opacityScale = 1.f,
                                                     float domainScale = 0.15f)
  {
    TransferFunction xf;
    const float hi = std::fabs(valueRange.upper);
    const float lo = std::fabs(valueRange.lower);
    const bool oneSided = std::min(lo,hi) <= 0.05f*std::max(lo,hi);
    const float a = (oneSided ? 1.f : domainScale) * std::max(lo,hi);
    const float A = (a > 0.f) ? a : 1.f;
    xf.relDomain = range1f(0.f,100.f);
    xf.baseDensity = 100.f;
    (void)bounds;

    if (!oneSided) {
      xf.absDomain = range1f(-A,+A);
      xf.colorMap = makeSignedColorMap(129,alphaGamma,opacityScale);
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
        const float u = negative ? (1.f-t) : t;
        const vec3f c = start + (end-start)*u;
        xf.colorMap[i]
          = vec4f(c.x,c.y,c.z,opacityScale*std::pow(u,alphaGamma));
      }
    }
    return xf;
  }

}
