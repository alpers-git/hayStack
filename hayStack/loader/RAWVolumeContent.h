// SPDX-FileCopyrightText: Copyright (c) 2023-2026 Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hayStack/loader/DataLoader.h"
#include "hayStack/StructuredVolume.h"

namespace hs {
  namespace loader {

    /*! parse the comma-separated 'gpuIso=' value list off a resource
        specifier. Shared with the viewer, which needs to know up
        front how many iso-surfaces there will be (and at what values)
        so it can put a slider on each of them. */
    std::vector<float> parseGpuIsoValues(const ResourceSpecifier &dataURL);

    /*! a file of 'raw' spheres */
    struct RAWVolumeContent : public LoadableContent {
      // using ScalarType = StructuredVolume::ScalarType;
    
      RAWVolumeContent(const std::string &fileName,
                       int thisPartID,
                       const box3i &cellRange,
                       vec3i fullVolumeDims,
                       const std::string &texelFormat,
                       int numChannels,
                       /*! if not NaN, we'll actually not store the
                         volume, but run iso-value extraction and use
                         the resulting surface(s) */
                       const float isoValue,
                       /*! if non-empty, we store the volume as usual,
                         but tag it to be rendered as GPU-side
                         ray-traced iso-surfaces at these values
                         instead of being volume-rendered */
                       const std::vector<float> &gpuIsoValues,
                       const vec3f &isoColor);
    
      static void create(DataLoader *loader,
                         const ResourceSpecifier &dataURL);
      size_t projectedSize() override;
      void   executeLoad(OnePartition &dataGroup) override;

      std::string toString() override;

      const std::string   fileName;
      const int           thisPartID;
      const vec3i         fullVolumeDims;
      const box3i         cellRange;
      const int           numChannels;
      const std::string   texelFormat;
      const float         isoValue;
      const std::vector<float> gpuIsoValues;
      const vec3f         isoColor;

      /*! world placement of the grid. by default the volume lives in VOXEL
          space (origin at the cell range, spacing 1), which forces every other
          object in the scene to be converted to voxel coordinates to line up.
          'gridOrigin='/'gridSpacing=' place it in world units instead - for a
          Gaussian cube, the header's origin and voxel vectors in Angstrom - so
          spheres, cylinders and orbital fields can all stay in Angstrom. */
      bool  havePlacement = false;
      vec3f gridOriginOpt { 0.f };
      vec3f gridSpacingOpt { 1.f };
    };
  
  }
}
