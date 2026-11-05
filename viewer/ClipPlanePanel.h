// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hayMaker/RenderEngineInterface.h"

namespace hs {

  /*! one clip plane the user can switch on, tilt, and slide through the scene;
      the side its normal points to stays visible */
  struct ClipPlanePanel {
    ClipPlanePanel(const box3f &bounds);

    /*! draws a collapsed "Clip plane" section into the current window;
        returns true if the plane changed. viewDir is for "face camera" */
    bool run_ui(const vec3f &viewDir);

    /*! switch the plane on as dot(n,p)+d >= 0, e.g. from the command line */
    void set(const vec3f &n, float d);

    hm::ClipPlane plane() const;
    bool isEnabled() const { return enabled; }

    /*! a mouse drag of (dx,dy) pixels in the render window: tilt turns the
        normal like a trackball, otherwise the plane slides along its normal */
    void drag(float dx, float dy, bool tilt,
              const Camera &camera, int windowHeight);

  private:
    vec3f normal() const;
    void  setNormal(const vec3f &n);
    /*! half the bounds' extent along N, the range the plane can slide over */
    float halfExtent(const vec3f &N) const;

    box3f bounds;
    bool  enabled   = false;
    /*! normal direction in degrees; 0/0 is +x */
    float azimuth   = 0.f;
    float elevation = 0.f;
    /*! the plane tilts about this point and slides it along the normal */
    vec3f pivot;
  };

} // ::hs
