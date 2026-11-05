// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#include "ClipPlanePanel.h"
#include "imgui.h"
#include <algorithm>
#include <cmath>

namespace hs {

  namespace {
    const float degToRad = float(M_PI)/180.f;
  }

  ClipPlanePanel::ClipPlanePanel(const box3f &bounds)
    : bounds(bounds),
      pivot(bounds.center())
  {}

  vec3f ClipPlanePanel::normal() const
  {
    const float az = azimuth*degToRad, el = elevation*degToRad;
    return vec3f(cosf(el)*cosf(az), cosf(el)*sinf(az), sinf(el));
  }

  void ClipPlanePanel::setNormal(const vec3f &n)
  {
    const vec3f N = normalize(n);
    azimuth   = atan2f(N.y,N.x)/degToRad;
    elevation = asinf(std::clamp(N.z,-1.f,1.f))/degToRad;
  }

  float ClipPlanePanel::halfExtent(const vec3f &N) const
  {
    const vec3f s = bounds.size();
    return .5f*(fabsf(N.x)*s.x + fabsf(N.y)*s.y + fabsf(N.z)*s.z);
  }

  void ClipPlanePanel::set(const vec3f &n, float d)
  {
    const float len = length(n);
    if (len == 0.f) return;
    setNormal(n);
    const vec3f N = n/len;
    const vec3f c = bounds.center();
    pivot   = c - (dot(N,c) + d/len)*N;
    enabled = true;
  }

  hm::ClipPlane ClipPlanePanel::plane() const
  {
    hm::ClipPlane p;
    p.enabled = enabled;
    p.N = normal();
    p.d = -dot(p.N,pivot);
    return p;
  }

  void ClipPlanePanel::drag(float dx, float dy, bool tilt,
                            const Camera &camera, int windowHeight)
  {
    const vec3f dir   = normalize(camera.vi-camera.vp);
    const vec3f right = normalize(cross(dir,camera.vu));
    const vec3f up    = cross(right,dir);
    const float h     = float(std::max(windowHeight,1));
    const vec3f N     = normal();
    const vec3f c     = bounds.center();

    if (tilt) {
      /* trackball: the side facing the camera follows the mouse, half a
         turn per window height */
      const vec3f axis = cross(-dir,dx*right-dy*up);
      if (length(axis) == 0.f) return;
      const vec3f a     = normalize(axis);
      const float angle = float(M_PI)*sqrtf(dx*dx+dy*dy)/h;
      setNormal(N*cosf(angle) + cross(a,N)*sinf(angle)
                + a*dot(a,N)*(1.f-cosf(angle)));
    } else {
      /* world units per pixel at the plane's depth, so the plane keeps up
         with the mouse along its normal's on-screen direction */
      const float depth
        = std::max(dot(pivot-camera.vp,dir),1e-3f*length(bounds.size()));
      const float perPixel = 2.f*depth*tanf(.5f*camera.fovy*degToRad)/h;
      const float sx = dot(N,right), sy = -dot(N,up);
      const float len = sqrtf(sx*sx+sy*sy);
      /* a normal facing the camera has no screen direction: dragging up then
         pushes the plane away */
      const float shift
        = len > .2f
        ? (dx*sx+dy*sy)/len * perPixel/len
        : -dy*(dot(N,dir) < 0.f ? -1.f : 1.f) * perPixel;
      pivot = pivot + shift*N;
    }

    /* keep the plane within the scene's extent along its (new) normal */
    const vec3f M   = normal();
    const float r   = halfExtent(M);
    const float pos = dot(M,pivot-c);
    pivot = pivot + (std::clamp(pos,-r,r)-pos)*M;
  }

  bool ClipPlanePanel::run_ui(const vec3f &viewDir)
  {
    /* a section of the caller's window, closed until the user opens it */
    if (!ImGui::CollapsingHeader("Clip plane")) return false;
    bool changed = false;
    changed |= ImGui::Checkbox("clip the scene",&enabled);

    ImGui::BeginDisabled(!enabled);
    changed |= ImGui::SliderFloat("azimuth",&azimuth,-180.f,180.f,"%.1f deg");
    changed |= ImGui::SliderFloat("elevation",&elevation,-90.f,90.f,"%.1f deg");

    /* signed distance of the plane from the scene centre */
    const vec3f N = normal();
    const vec3f c = bounds.center();
    const float r = halfExtent(N);
    float pos = dot(N,pivot-c);
    if (ImGui::SliderFloat("position",&pos,-r,r,"%.3f")) {
      pivot = pivot + (pos - dot(N,pivot-c))*N;
      changed = true;
    }

    if (ImGui::Button("+x")) { setNormal(vec3f(1.f,0.f,0.f)); changed = true; }
    ImGui::SameLine();
    if (ImGui::Button("+y")) { setNormal(vec3f(0.f,1.f,0.f)); changed = true; }
    ImGui::SameLine();
    if (ImGui::Button("+z")) { setNormal(vec3f(0.f,0.f,1.f)); changed = true; }
    ImGui::SameLine();
    if (ImGui::Button("flip")) { setNormal(-N); changed = true; }
    ImGui::SameLine();
    if (ImGui::Button("face camera")) { setNormal(viewDir); changed = true; }
    ImGui::SameLine();
    if (ImGui::Button("centre")) { pivot = c; changed = true; }

    ImGui::Text("normal (%.3f %.3f %.3f) points to the side that stays",
                N.x,N.y,N.z);
    ImGui::TextDisabled("in the render: ctrl+drag tilts, shift+drag slides");
    ImGui::EndDisabled();
    return changed;
  }

} // ::hs
