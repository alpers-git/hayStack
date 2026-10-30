// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#include "hayStack/loader/MultiVariateContent.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace hs {
  namespace loader {

    static std::vector<std::string> splitCommas(const std::string &s)
    {
      std::vector<std::string> out;
      std::stringstream ss(s);
      for (std::string tok; std::getline(ss,tok,','); )
        if (!tok.empty()) out.push_back(tok);
      return out;
    }

    ChannelMix parseChannelMix(const ResourceSpecifier &url,
                               const std::vector<std::string> &names,
                               const NameMatch &matches)
    {
      ChannelMix mix;
      mix.weights.assign(names.size(),1.f);

      const std::string sel = url.get("channels");
      if (!sel.empty() && sel != "all") {
        std::fill(mix.weights.begin(),mix.weights.end(),0.f);
        for (auto &tok : splitCommas(sel)) {
          bool found = false;
          for (size_t c=0;c<names.size();c++)
            if (names[c] == tok || (matches && matches(names[c],tok))) {
              mix.weights[c] = 1.f;
              found = true;
            }
          if (!found) {
            std::string all; for (auto &n : names) all += n+" ";
            throw std::runtime_error("#hs.mv: unknown channel '"+tok
                                     +"'; this field has: "+all);
          }
        }
      }

      const std::string w = url.get("weights");
      if (!w.empty()) {
        const auto toks = splitCommas(w);
        if (toks.size() != names.size())
          throw std::runtime_error("#hs.mv: 'weights' needs one value per channel ("
                                   +std::to_string(names.size())+")");
        for (size_t c=0;c<toks.size();c++)
          mix.weights[c] = std::stof(toks[c]);
      }

      const std::string r = url.get("reduction");
      if (!r.empty() && !parseReduction(r,mix.reduction))
        throw std::runtime_error("#hs.mv: unknown reduction '"+r
                                 +"'; use magnitude, component or dot");
      mix.component = url.get_int("component",0);
      mix.direction = url.get_vec3f("direction",mix.direction);
      return mix;
    }

    MultiVariateFromRaw::MultiVariateFromRaw(const ResourceSpecifier &data)
      : data(data), files(splitCommas(data.where))
    {
      if (files.empty())
        throw std::runtime_error("#hs.mv: no files in '"+data.where+"'");
    }

    void MultiVariateFromRaw::create(DataLoader *loader,
                                     const ResourceSpecifier &dataURL)
    { loader->addContent(new MultiVariateFromRaw(dataURL)); }

    std::string MultiVariateFromRaw::toString()
    { return "MultiVariateFromRaw{"+std::to_string(files.size())+" channels, proj size "
        +prettyNumber(projectedSize())+"B}"; }

    size_t MultiVariateFromRaw::projectedSize()
    {
      size_t sum = 0;
      for (auto &f : files) sum += getFileSize(f);
      return 2*sum;
    }

    void MultiVariateFromRaw::executeLoad(OnePartition &dataGroup)
    {
      const std::string format = data.get("format","float");
      if (format != "float" && format != "float32")
        throw std::runtime_error("#hs.mv: only float32 raw channels are supported");

      MultiVariateVolume::SP vol = MultiVariateVolume::create();
      vol->numChannels   = (int)files.size();
      vol->numComponents = data.get_int("components",1);
      if (vol->numComponents != 1 && vol->numComponents != 3)
        throw std::runtime_error("#hs.mv: 'components' must be 1 or 3");
      const std::string dimsString = data.get("dims");
      if (sscanf(dimsString.c_str(),"%i,%i,%i",
                 &vol->dims.x,&vol->dims.y,&vol->dims.z) != 3)
        throw std::runtime_error("#hs.mv: needs dims=X,Y,Z");
      vol->gridOrigin  = data.get_vec3f("gridOrigin",vec3f(0.f));
      vol->gridSpacing = data.get_vec3f("gridSpacing",vec3f(1.f));
      vol->bounds.lower = vol->gridOrigin;
      vol->bounds.upper = vol->gridOrigin
        + vol->gridSpacing*vec3f(vol->dims-vec3i(1));

      vol->channelNames = splitCommas(data.get("names"));
      if (vol->channelNames.empty())
        for (auto &f : files) {
          const size_t slash = f.find_last_of('/');
          vol->channelNames.push_back(slash == std::string::npos ? f : f.substr(slash+1));
        }
      if ((int)vol->channelNames.size() != vol->numChannels)
        throw std::runtime_error("#hs.mv: 'names' needs one name per file");

      const size_t perChannel = vol->numVoxels()*vol->numComponents;
      vol->channelData.resize(perChannel*vol->numChannels);
      vol->channelPeak.assign(vol->numChannels,0.f);
      for (int c=0;c<vol->numChannels;c++) {
        std::ifstream in(files[c],std::ios::binary);
        float *dst = vol->channelData.data() + c*perChannel;
        in.read((char*)dst,perChannel*sizeof(float));
        if (!in.good())
          throw std::runtime_error("#hs.mv: could not read "+std::to_string(perChannel)
                                   +" floats from "+files[c]);
        for (size_t v=0;v<vol->numVoxels();v++) {
          const float *x = dst + v*vol->numComponents;
          const float m = vol->numComponents == 1 ? std::fabs(x[0])
            : std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
          vol->channelPeak[c] = std::max(vol->channelPeak[c],m);
        }
      }

      vol->mix = parseChannelMix(data,vol->channelNames);
      const float p = vol->peakFor(vol->mix);
      vol->valueRange = vol->mix.isSigned() ? range1f(-p,p) : range1f(0.f,p);
      vol->alphaGamma   = data.get_float("alphaGamma",1.f);
      vol->opacityScale = data.get_float("opacity",1.f);
      vol->domainScale  = data.get_float("domainScale",1.f);

      std::cout << "#hs.mv: " << vol->numChannels << " channels x "
                << vol->numComponents << " components on " << vol->dims.x << "x"
                << vol->dims.y << "x" << vol->dims.z << ", reduction "
                << reductionName(vol->mix.reduction) << ", range ["
                << vol->valueRange.lower << " " << vol->valueRange.upper << "]"
                << std::endl;
      dataGroup.multiVariates.push_back(vol);
    }

  }
}
