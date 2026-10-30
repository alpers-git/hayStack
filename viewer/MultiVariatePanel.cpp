// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#include "MultiVariatePanel.h"
#include "imgui.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace hs {

  namespace {
    float safeLog10(float v) { return std::log10(std::max(v,1e-30f)); }
  }

  MultiVariatePanel::MultiVariatePanel(const MultiVariateVolume::SP &volume)
    : volume(volume)
  {
    if (!volume) return;
    mix = volume->mix;
    mix.weights.resize(volume->numChannels,1.f);
    restoreWeight.resize(volume->numChannels,1.f);
    for (int c=0;c<volume->numChannels;c++)
      if (mix.weights[c] != 0.f) restoreWeight[c] = mix.weights[c];

    absMax = std::max(std::fabs(volume->valueRange.lower),
                      std::fabs(volume->valueRange.upper));
    if (absMax <= 0.f) absMax = 1.f;
    logAbsMax = safeLog10(absMax);
    if (!mix.isSigned()) {
      const float lo = std::fabs(volume->valueRange.lower);
      const float hi = std::fabs(volume->valueRange.upper);
      twoSided = std::min(lo,hi) > 0.05f*std::max(lo,hi);
      negativeField = !twoSided && lo > hi;
    }
    /* 5 decades below the peak is transparent, 3 above saturated */
    logDomainLo = logAbsMax - 5.f;
    logDomainHi = logAbsMax + 3.f;
  }

  float MultiVariatePanel::mixDomain() const
  {
    const float scale = volume->domainScale > 0.f ? volume->domainScale : 1.f;
    return volume->peakFor(mix) * scale;
  }

  void MultiVariatePanel::fitDomain(float peak)
  {
    if (!(peak > 0.f)) return;
    logAbsMax = safeLog10(peak);
    absMax = std::pow(10.f,logAbsMax);
  }

  const char *MultiVariatePanel::valueLabel() const
  {
    const bool scalar = volume->numComponents == 1;
    switch (mix.reduction) {
    case Reduction::Component:
      if (scalar) return "s";
      return mix.component == 0 ? "v.x" : (mix.component == 1 ? "v.y" : "v.z");
    case Reduction::Dot:
      return "v . d";
    default:
      return scalar ? "|s|" : "|v|";
    }
  }

  void MultiVariatePanel::drawRuler(float width) const
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
      float tx = x - ts.x * 0.5f;
      tx = std::max(p0.x, std::min(tx, p0.x + width - ts.x));
      dl->AddText(ImVec2(tx, p0.y + 8.f), col, lab);
    }
    ImGui::Dummy(ImVec2(width, h));
    ImGui::TextDisabled("%s  (%s)", valueLabel(), volume->unit.c_str());
  }

  bool MultiVariatePanel::run_ui()
  {
    committed = false;
    if (!volume) return false;

    bool changed = false;
    bool mixChanged = false;
    const std::string unit = volume->unit.empty() ? "" : " ("+volume->unit+")";

    ImGui::Begin("Multi-variate field");
    ImGui::Text("%d channels x %d components on %d x %d x %d voxels",
                volume->numChannels,volume->numComponents,
                volume->dims.x,volume->dims.y,volume->dims.z);
    ImGui::Text("bounds  %.1f %.1f %.1f  ..  %.1f %.1f %.1f",
                volume->bounds.lower.x,volume->bounds.lower.y,volume->bounds.lower.z,
                volume->bounds.upper.x,volume->bounds.upper.y,volume->bounds.upper.z);

    if (ImGui::CollapsingHeader("Volume",ImGuiTreeNodeFlags_DefaultOpen)) {
      if (ImGui::SliderFloat("|domain| (1e^)",&logAbsMax,
                             logDomainLo,logDomainHi,"%.2f")) {
        absMax = std::pow(10.f,logAbsMax);
        changed = true;
      }
      ImGui::SameLine();
      ImGui::Text("= %.3g",absMax);

      float typed = absMax;
      ImGui::SetNextItemWidth(160.f);
      if (ImGui::InputFloat(("max"+unit).c_str(),&typed,0.f,0.f,"%.4e",
                            ImGuiInputTextFlags_EnterReturnsTrue)
          && typed > 0.f && std::isfinite(typed)) {
        absMax = typed;
        logAbsMax = safeLog10(absMax);
        changed = true;
      }
      ImGui::SameLine();
      if (ImGui::SmallButton("fit")) {
        const float m = mixDomain();
        fitDomain(m > 0.f ? m : std::max(std::fabs(volume->valueRange.lower),
                                         std::fabs(volume->valueRange.upper)));
        changed = true;
      }

      float lo = absMin;
      ImGui::SetNextItemWidth(160.f);
      if (ImGui::InputFloat(("min"+unit).c_str(),&lo,0.f,0.f,"%.4e",
                            ImGuiInputTextFlags_EnterReturnsTrue)
          && std::isfinite(lo) && lo >= 0.f) {
        absMin = lo;
        changed = true;
      }
      ImGui::SameLine();
      if (ImGui::SmallButton("0##absmin")) {
        absMin = 0.f;
        changed = true;
      }
      const range1f d = domain();
      ImGui::Text("domain  %.4e .. %.4e%s",d.lower,d.upper,unit.c_str());

      float rl = relLo, rh = relHi;
      if (ImGui::DragFloatRange2("rel domain",&rl,&rh,0.5f,0.f,100.f,"%.1f")) {
        relLo = std::min(rl,rh); relHi = std::max(rl,rh);
        changed = true;
      }
      if (ImGui::SliderFloat("density",&density,0.f,200.f))
        changed = true;
      ImGui::TextDisabled(signedDomain()
                          ? "signed values: the colour map is centred on zero"
                          : "one-sided values: the colour map starts at zero");
    }

    if (ImGui::CollapsingHeader("Channel mix",ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::TextDisabled("weighted sum of the channels, reduced to one scalar");
      if (ImGui::SmallButton("all##ch")) {
        for (int c=0;c<volume->numChannels;c++) mix.weights[c] = restoreWeight[c];
        mixChanged = true;
      }
      ImGui::SameLine();
      if (ImGui::SmallButton("none##ch")) {
        std::fill(mix.weights.begin(),mix.weights.end(),0.f);
        mixChanged = true;
      }
      for (int c=0;c<volume->numChannels;c++) {
        ImGui::PushID(c);
        const char *name = c < (int)volume->channelNames.size()
          ? volume->channelNames[c].c_str() : "?";
        bool on = mix.weights[c] != 0.f;
        if (ImGui::Checkbox("##on",&on)) {
          mix.weights[c] = on ? restoreWeight[c] : 0.f;
          mixChanged = true;
        }
        ImGui::SameLine();
        ImGui::SetNextItemWidth(160.f);
        if (ImGui::SliderFloat(name,&mix.weights[c],-2.f,2.f,"%.2f")
            && mix.weights[c] != 0.f)
          restoreWeight[c] = mix.weights[c];
        mixChanged |= ImGui::IsItemDeactivatedAfterEdit();
        ImGui::PopID();
      }

      const bool scalar = volume->numComponents == 1;
      static const char *vectorItems[] = { "magnitude", "component", "dot" };
      static const char *scalarItems[] = { "|value|", "value" };
      int r = (int)mix.reduction;
      if (scalar && r > 1) r = 1;
      if (ImGui::Combo("reduction",&r,scalar ? scalarItems : vectorItems,
                       scalar ? 2 : 3)) {
        mix.reduction = Reduction(r);
        mix.component = scalar ? 0 : mix.component;
        mixChanged = true;
      }
      if (mix.reduction == Reduction::Component && !scalar) {
        static const char *axes[] = { "x", "y", "z" };
        for (int k=0;k<3;k++) {
          if (k) ImGui::SameLine();
          if (ImGui::RadioButton(axes[k],mix.component == k)) {
            mix.component = k;
            mixChanged = true;
          }
        }
      }
      if (mix.reduction == Reduction::Dot) {
        float d[3] = { mix.direction.x, mix.direction.y, mix.direction.z };
        if (ImGui::InputFloat3("direction",d,"%.3f",
                               ImGuiInputTextFlags_EnterReturnsTrue)) {
          mix.direction = vec3f(d[0],d[1],d[2]);
          mixChanged = true;
        }
      }
      ImGui::Checkbox("rescale colours to mix",&autoFit);

      if (!volume->channelPeak.empty()) {
        std::string s;
        for (int c=0;c<volume->numChannels && c<(int)volume->channelPeak.size();c++)
          if (mix.weights[c] != 0.f) {
            char b[64];
            snprintf(b,sizeof(b),"%s %.2g  ",
                     c < (int)volume->channelNames.size()
                     ? volume->channelNames[c].c_str() : "?",
                     volume->channelPeak[c]);
            s += b;
          }
        ImGui::TextDisabled("peak per enabled channel: %s",s.c_str());
        ImGui::TextDisabled("peaks are not additive: channels can cancel");
      }
    }
    ImGui::End();

    if (mixChanged) {
      mix.dirty = true;
      committed = true;
      if (autoFit) fitDomain(mixDomain());
      changed = true;
    }
    return changed;
  }

} // ::hs
