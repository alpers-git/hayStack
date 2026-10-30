// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hayStack/MultiVariateVolume.h"
#include <string>
#include <vector>
#include <algorithm>

namespace hs {

  /*! side panel for a multi-variate volume: the transfer-function domain in
      the field's own units, and the channel mix applied before classification */
  struct MultiVariatePanel {
    MultiVariatePanel(const MultiVariateVolume::SP &volume);

    /*! returns true if the transfer function needs updating */
    bool run_ui();

    /*! true when a mix control was released this frame; pushing only on
        release keeps a drag from re-mixing every frame */
    bool mixCommitted() const { return committed; }

    const ChannelMix &getMix() const { return mix; }

    range1f domain() const
    {
      if (signedDomain()) return { -absMax, +absMax };
      if (negativeField)  return { -absMax, 0.f };
      return { std::min(absMin,absMax*0.999f), +absMax };
    }
    range1f relDomain() const { return { relLo, relHi }; }
    float   baseDensity() const { return density; }

    /*! labelled value axis to draw under the colour-map gizmo */
    void drawRuler(float width = 0.f) const;

  private:
    bool  signedDomain() const { return twoSided || mix.isSigned(); }
    float mixDomain() const;
    void  fitDomain(float peak);
    const char *valueLabel() const;

    MultiVariateVolume::SP volume;
    ChannelMix mix;
    std::vector<float> restoreWeight;
    bool committed = false;

    /*! the slider edits log10(absMax), bounded around the field's own peak */
    float logAbsMax = 0.f;
    float absMax    = 1.f;
    float absMin    = 0.f;
    float relLo     = 0.f;
    float relHi     = 100.f;
    float density   = 100.f;
    float logDomainLo = -6.f, logDomainHi = 12.f;

    /*! rescale the domain to the mix's peak whenever the mix changes */
    bool autoFit = true;
    bool twoSided = false;
    bool negativeField = false;
  };

} // ::hs
