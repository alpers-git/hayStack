// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hayStack/loader/DataLoader.h"
#include "hayStack/MultiVariateVolume.h"
#include <functional>

namespace hs {
  namespace loader {

    /*! extra channel-name spellings a loader accepts in 'channels=' */
    typedef std::function<bool(const std::string &name,
                               const std::string &token)> NameMatch;

    /*! the channel-mix url options shared by every multi-variate loader:
        channels=a,b | weights=w0,w1,.. , reduction=, component=, direction= */
    ChannelMix parseChannelMix(const ResourceSpecifier &url,
                               const std::vector<std::string> &channelNames,
                               const NameMatch &matches = {});

    /*! 'multivariate://a.raw,b.raw:dims=X,Y,Z[:components=1|3:names=..]' loads
        each float32 raw file as one channel of ONE volume (not one volume each) */
    struct MultiVariateFromRaw : public LoadableContent {
      MultiVariateFromRaw(const ResourceSpecifier &data);
      static void create(DataLoader *loader, const ResourceSpecifier &dataURL);
      size_t projectedSize() override;
      void   executeLoad(OnePartition &dataGroup) override;
      std::string toString() override;

      const ResourceSpecifier  data;
      std::vector<std::string> files;
    };

  }
}
