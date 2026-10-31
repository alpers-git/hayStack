// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hayStack/loader/DataLoader.h"
#include "hayStack/MultiVariateVolume.h"
#include <fstream>

namespace hs {
  namespace loader {

    /*! loads a .orbg file (per-channel orbital-current vectors on a dense grid)
        as a 3-component multi-variate volume; ':load=' limits what is read */
    struct OrbGridFromFile : public LoadableContent {
      OrbGridFromFile(const ResourceSpecifier &data, int thisPartID);
      static void create(DataLoader *loader, const ResourceSpecifier &dataURL);
      size_t projectedSize() override;
      void   executeLoad(OnePartition &dataGroup) override;
      std::string toString() override;

      const ResourceSpecifier data;
      const size_t fileSize;
      const int thisPartID;
    };

  }
}
