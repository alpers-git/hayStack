// SPDX-FileCopyrightText: Copyright (c) 2023-2026 Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hayStack/TransferFunction.h"
#include "hayStack/NanoVDBVolume.h"
#include "hayStack/MultiVariateVolume.h"

/* parallel renderer abstraction */
namespace hm {
  
  /*! base abstraction for any renderer - no matter whether it's a
      single node or multiple workers on the back */
  struct RenderEngineInterface {

    virtual void setTransferFunction(const hs::TransferFunction &xf) {}
    /*! re-set the iso-values of all 'gpuIso=' iso-surfaces, in
        the order they were created */
    virtual void setIsoValues(const std::vector<float> &values) {}
    /*! re-set the channel mix of every multi-variate volume; only pushed
        when the caller marks it dirty */
    virtual void setChannelMix(const hs::ChannelMix &mix) {}
    /*! swap which file channels the multi-variate volumes hold; mix.weights
        has to match the new set */
    virtual void setLoadedChannels(const std::vector<int> &channels,
                                   const hs::ChannelMix &mix) {}
    virtual void setVolumeScatterSettings(const hs::VolumeScatterSettings &settings) {}
    virtual hs::VolumeScatterSettings getVolumeScatterSettings() const { return {}; }
    
    virtual void renderFrame() {}
    virtual void resize(const vec2i &fbSize, uint32_t *hostRgba) {}
    virtual void resetAccumulation() {}
    virtual void setCamera(const hs::Camera &camera) {}
    // virtual void setXF(const range1f &domain,
    //                    const std::vector<vec4f> &colors) {}
    virtual void screenShot() {}
    virtual void terminate() {}
    virtual void setLights(float ambient,
                           const std::vector<hs::PointLight> &pointLights,
                           const std::vector<hs::DirLight> &dirLights) {}
  };

}
