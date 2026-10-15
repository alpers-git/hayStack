// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#include "OrbitalPanel.h"
#include "imgui.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace hs {

  namespace {
    /*! the domain slider has to span many decades, so it lives in log10 */
    float safeLog10(float v) { return std::log10(std::max(v,1e-30f)); }
  }

  OrbitalPanel::OrbitalPanel(const OrbGridVolume::SP &volume)
    : volume(volume)
  {
    if (!volume) return;

    /* orbital channels come straight from the file, and start from whatever
       the url asked for - so a build launched with channels=s-d opens on that
       selection rather than resetting it on the first frame */
    hasChannels = !volume->channelNames.empty();
    if (hasChannels) {
      chanOn.assign(volume->channelNames.size(),0);
      for (size_t c=0;c<chanOn.size();c++)
        chanOn[c] = (volume->channelMask >> c) & 1u;
      filter.channelMask = volume->channelMask;
    }

    /* the converter's sampled range gives the scale */
    absMax = std::max(std::fabs(volume->valueRange.lower),
                      std::fabs(volume->valueRange.upper));
    if (absMax <= 0.f) absMax = 1.f;
    logAbsMax = safeLog10(absMax);
    {
      const float lo = std::fabs(volume->valueRange.lower);
      const float hi = std::fabs(volume->valueRange.upper);
      twoSided = std::min(lo,hi) > 0.05f*std::max(lo,hi);
      negativeField = !twoSided && lo > hi;
    }

    /* bound the slider around the data, not around a guess. an orbital current
       field peaks near 1e-9 and an electron density near 1e0, so no fixed
       range can hold both: whichever one falls outside gets clamped the moment
       the slider is touched. 5 decades down is fully transparent and 3 up is
       fully saturated, which is all the travel that means anything. */
    logDomainLo = logAbsMax - 5.f;
    logDomainHi = logAbsMax + 3.f;
  }

  /*! the exact peak of the enabled channels' combined field when the file
      carries it per subset, else the sum of their individual peaks. |sum| can
      never exceed that sum (triangle inequality) and usually falls well short,
      because these channels are anti-aligned at better than half the sample
      points - so it is an upper bound that still lands in the right decade. */
  float OrbitalPanel::channelDomain() const
  {
    if (!hasChannels || volume->channelPeak.empty()) return 0.f;
    uint32_t mask = 0;
    for (size_t c=0;c<chanOn.size();c++) if (chanOn[c]) mask |= (1u<<c);
    /* the exact peak of this selection, when the field knows it. Summing the
       per-channel peaks assumes the channels reinforce; they partially cancel,
       so for the orbital-current field the sum overshoots the real maximum by
       over 3x on the full set - enough to push everything to the bottom of the
       colour ramp and make a toggle look like the data disappeared. */
    /* domainScale is what the startup transfer function is built with, so the
       auto-fit has to apply it too - otherwise the first toggle jumps to a
       different exposure than the one the scene opened on. */
    const float scale = volume->domainScale > 0.f ? volume->domainScale : 1.f;
    if (mask < volume->subsetPeak.size())
      return volume->subsetPeak[mask] * scale;
    float sum = 0.f;
    for (size_t c=0;c<chanOn.size() && c<volume->channelPeak.size();c++)
      if (chanOn[c]) sum += volume->channelPeak[c];
    return sum * scale;
  }


  /*! a tick axis under the colour ramp.

    The gizmo draws the map over a unit bar with no numbers on it, so there is
    no way to tell what value a colour corresponds to. This lays five labelled
    ticks across the current domain, in the field's own units. */
  void OrbitalPanel::drawRuler(float width) const
  {
    const range1f d = domain();
    if (!(std::isfinite(d.lower) && std::isfinite(d.upper))) return;
    if (width <= 0.f) width = ImGui::GetContentRegionAvail().x;
    if (width < 32.f) return;

    const float h = 18.f;
    ImDrawList *dl = ImGui::GetWindowDrawList();
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const ImU32 col = ImGui::GetColorU32(ImGuiCol_Text);

    dl->AddLine(ImVec2(p0.x, p0.y + 2.f),
                ImVec2(p0.x + width, p0.y + 2.f), col, 1.f);
    const int nTicks = 5;
    for (int i = 0; i < nTicks; i++) {
      const float f = i / float(nTicks - 1);
      const float x = p0.x + f * width;
      dl->AddLine(ImVec2(x, p0.y + 2.f), ImVec2(x, p0.y + 7.f), col, 1.f);
      char lab[32];
      snprintf(lab, sizeof(lab), "%.2e", d.lower + f * (d.upper - d.lower));
      const ImVec2 ts = ImGui::CalcTextSize(lab);
      float tx = x - ts.x * 0.5f;                 // keep the ends inside
      tx = std::max(p0.x, std::min(tx, p0.x + width - ts.x));
      dl->AddText(ImVec2(tx, p0.y + 8.f), col, lab);
    }
    ImGui::Dummy(ImVec2(width, h));
    ImGui::TextDisabled("|j|  (%s)", valueUnit());
  }

  bool OrbitalPanel::run_ui()
  {
    committed = false;
    if (!volume) return false;

    bool changed = false;
    bool release = false;

    ImGui::Begin("Orbital field");

    ImGui::Text("bounds  %.1f %.1f %.1f  ..  %.1f %.1f %.1f",
                volume->bounds.lower.x,volume->bounds.lower.y,volume->bounds.lower.z,
                volume->bounds.upper.x,volume->bounds.upper.y,volume->bounds.upper.z);
    ImGui::Text("field   %.3g .. %.3g",
                volume->valueRange.lower,volume->valueRange.upper);

    // ================================================================
    // transfer function domain, in the field's own units
    // ================================================================
    if (ImGui::CollapsingHeader("Volume",ImGuiTreeNodeFlags_DefaultOpen)) {
      if (ImGui::SliderFloat("|domain| (1e^)",&logAbsMax,
                             logDomainLo,logDomainHi,"%.2f")) {
        absMax = std::pow(10.f,logAbsMax);
        changed = true;
      }
      release |= ImGui::IsItemDeactivatedAfterEdit();
      ImGui::SameLine();
      ImGui::Text("= %.3g",absMax);
      ImGui::SameLine();

      /* the slider alone cannot express an exact value, and at 1e-10 the
         default %.3f formatting of every imgui float widget prints 0.000 -
         which is why the domain read as empty. Type the number instead. */
      {
        float typed = absMax;
        ImGui::SetNextItemWidth(160.f);
        if (ImGui::InputFloat(("max (" + std::string(valueUnit()) + ")").c_str(),
                              &typed, 0.f, 0.f, "%.4e",
                              ImGuiInputTextFlags_EnterReturnsTrue)) {
          if (typed > 0.f && std::isfinite(typed)) {
            absMax = typed;
            logAbsMax = safeLog10(absMax);
            changed = true; release = true;
          }
        }
        const range1f d = domain();
        ImGui::Text("domain  %.4e .. %.4e  %s", d.lower, d.upper, valueUnit());
      }
      ImGui::SameLine();
      /* a way back to the field's own scale, so a domain dragged somewhere
         useless is one click to recover rather than a guess */
      if (ImGui::SmallButton("fit")) {
        /* fit to what is actually shown: snapping to the file's whole-field
           range would leave a single weak channel 8x under-scaled, which is
           the same trap the auto-rescale below exists to avoid */
        const float ch = channelDomain();
        logAbsMax = safeLog10(ch > 0.f ? ch
                              : std::max(std::fabs(volume->valueRange.lower),
                                         std::fabs(volume->valueRange.upper)));
        absMax = std::pow(10.f,logAbsMax);
        changed = release = true;
      }
      /* the same domain controls the generic volume ui offers, in the units
         this field actually uses rather than a 0..1 bar */
      {
        float lo = absMin;
        ImGui::SetNextItemWidth(160.f);
        if (ImGui::InputFloat(("min (" + std::string(valueUnit()) + ")").c_str(),
                              &lo, 0.f, 0.f, "%.4e",
                              ImGuiInputTextFlags_EnterReturnsTrue)) {
          if (std::isfinite(lo) && lo >= 0.f) {
            absMin = lo; changed = true; release = true;
          }
        }
        ImGui::SameLine();
        if (ImGui::SmallButton("0##absmin")) {
          absMin = 0.f; changed = release = true;
        }
      }
      float rl = relLo, rh = relHi;
      if (ImGui::DragFloatRange2("rel domain",&rl,&rh,0.5f,0.f,100.f,"%.1f")) {
        relLo = std::min(rl,rh); relHi = std::max(rl,rh);
        changed = true;
      }
      release |= ImGui::IsItemDeactivatedAfterEdit();
      if (ImGui::SliderFloat("density",&density,0.f,200.f))
        changed = true;
      release |= ImGui::IsItemDeactivatedAfterEdit();
      if (twoSided)
        ImGui::TextDisabled("domain is symmetric: the field takes both signs, "
                            "so the colour map's midpoint is zero");
      else
        ImGui::TextDisabled("domain is one-sided %s: this field never changes "
                            "sign, so a symmetric domain would leave half the "
                            "colour map unused",
                            negativeField ? "(negative)" : "(positive)");
    }

    // ================================================================
    // orbital channels
    // ================================================================
    if (hasChannels
        && ImGui::CollapsingHeader("Orbital channels",
                                   ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::TextDisabled("which (l_A,l_B) overlap carries the current");
      if (ImGui::SmallButton("all##ch")) {
        for (auto &c : chanOn) c = 1;
        filter.channelDirty = changed = release = true;
      }
      ImGui::SameLine();
      if (ImGui::SmallButton("none##ch")) {
        for (auto &c : chanOn) c = 0;
        filter.channelDirty = changed = release = true;
      }
      for (size_t c=0;c<chanOn.size();c++) {
        if (c % 3) ImGui::SameLine();
        bool on = chanOn[c] != 0;
        const char *nm = c < volume->channelNames.size()
          ? volume->channelNames[c].c_str() : "?";
        if (ImGui::Checkbox(nm,&on)) {
          chanOn[c] = on ? 1 : 0;
          filter.channelDirty = changed = release = true;
        }
      }
      uint32_t mask = 0;
      for (size_t c=0;c<chanOn.size();c++) if (chanOn[c]) mask |= (1u<<c);
      filter.channelMask = mask;

      ImGui::Checkbox("rescale colours to selection",&autoFitChannels);
      if (autoFitChannels && filter.channelDirty) {
        const float d = channelDomain();
        if (d > 0.f) {
          logAbsMax = safeLog10(d);
          absMax = std::pow(10.f,logAbsMax);
        }
      }
      if (!volume->channelPeak.empty()) {
        std::string s;
        for (size_t c=0;c<chanOn.size() && c<volume->channelPeak.size();c++)
          if (chanOn[c]) {
            char b[64];
            snprintf(b,sizeof(b),"%s %.2g  ",
                     volume->channelNames[c].c_str(),volume->channelPeak[c]);
            s += b;
          }
        ImGui::TextDisabled("peak |j| per enabled channel: %s",s.c_str());
        ImGui::TextDisabled("these are NOT additive: |sum| != sum of |parts|, "
                            "and these channels oppose each other at over half "
                            "the sample points");
      }
    }

    ImGui::End();

    /* the majorant grid is rebuilt on every push, so only push on release */
    committed = release;
    return changed;
  }

} // ::hs
