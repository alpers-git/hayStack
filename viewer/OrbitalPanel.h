// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hayStack/OrbGridVolume.h"
#include <string>
#include <vector>
#include <algorithm>

namespace hs {

  /*! the side panel for an orbital-current grid field. it owns the
      transfer-function domain - in the field's own units, with a value axis
      under the colour map - and the per-channel selection that the field
      filters on.

      it replaces the generic volume and iso ui rather than sitting beside it:
      the generic controls all assume a domain that is the world scalar range,
      which is the wrong framing here. */
  struct OrbitalPanel {
    OrbitalPanel(const OrbGridVolume::SP &volume);

    /*! draw the panel. returns true if anything the renderer needs changed. */
    bool run_ui();

    /*! true when a control was released this frame - the caller uses this to
        hold off the majorant rebuild until a drag is over */
    bool filterCommitted() const { return committed; }

    const ChannelFilter &getFilter() const { return filter; }

    /*! the transfer-function domain the user picked. only a field that
        actually takes both signs gets a symmetric domain - forcing one on a
        non-negative field (a magnitude, say) puts every value in
        the upper half of the colour map, so the lower half is dead and the
        data reads as skewed to the top with no low end visible at all. */
    range1f domain() const
    {
      if (twoSided)      return { -absMax, +absMax };
      if (negativeField) return { -absMax, 0.f };
      /* absMin lets the bottom of the domain be lifted off zero, which is how
         you push the weak tail out of the picture without touching opacity */
      return { std::min(absMin,absMax*0.999f), +absMax };
    }

    /*! the secondary 0..100 window, same control the generic volume ui has */
    range1f relDomain() const { return { relLo, relHi }; }
    float   baseDensity() const { return density; }

    /*! physical unit of the sampled scalar: a current DENSITY (amperes per
      cubic Angstrom), not a current. That is a different quantity from the
      .cube's microamps, which is an interpolated per-atom current and never
      divided by a volume. */
    static const char *valueUnit() { return "A/A^3"; }

    /*! draw a labelled value axis under the colour-map gizmo, so the ramp can
        be read against real numbers instead of a 0..1 bar */
    void drawRuler(float width = 0.f) const;

  private:
    OrbGridVolume::SP volume;
    ChannelFilter     filter;
    bool              committed = false;

    // ---- transfer function ----
    /*! the slider state is the LOG, not absMax: absMax is derived from it.
        round-tripping through pow/log every frame loses the value the user
        picked, and an orbital field's peak (~1e-9) lives nowhere near any
        range that could be hardcoded. */
    float logAbsMax = 0.f;
    float absMax    = 1.f;
    float absMin    = 0.f;
    float relLo     = 0.f;
    float relHi     = 100.f;
    float density   = 100.f;

    /*! slider bounds, derived from the field's own peak. a hardcoded range is
        the bug this replaces: a value outside it gets clamped the instant the
        slider is touched, which silently rescaled the domain by ~250x and
        blanked the volume. */
    float logDomainLo = -6.f, logDomainHi = 12.f;

    // ---- orbital channels ----
    /*! one checkbox per (l_A,l_B) channel the file carries */
    std::vector<char> chanOn;
    bool hasChannels = false;
    /*! rescale the transfer-function domain to the enabled channels' own peak.
        On by default: the per-channel peaks span a factor of ten, so holding
        one domain across a toggle makes a weak channel look like an empty
        filter instead of a dim one. */
    bool autoFitChannels = true;
    /*! peak of the enabled channels' combined field, scaled like the startup
        domain */
    float channelDomain() const;

    /*! whether the field really straddles zero, by the same 5% test the
        default transfer function uses */
    bool twoSided = false;
    bool negativeField = false;
  };

} // ::hs
