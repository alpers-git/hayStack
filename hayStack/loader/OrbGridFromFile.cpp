// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#include "hayStack/loader/OrbGridFromFile.h"
#include <iostream>
#include <stdexcept>

namespace hs {
  namespace loader {

    OrbGridFromFile::OrbGridFromFile(const ResourceSpecifier &data, int thisPartID)
      : data(data), fileSize(getFileSize(data.where)), thisPartID(thisPartID)
    {}

    void OrbGridFromFile::create(DataLoader *loader, const ResourceSpecifier &dataURL)
    { loader->addContent(new OrbGridFromFile(dataURL,0)); }

    std::string OrbGridFromFile::toString()
    { return "OrbGrid{fileName="+data.where+", proj size "
        +prettyNumber(projectedSize())+"B}"; }

    size_t OrbGridFromFile::projectedSize()
    { return 2 * fileSize; }

    template<typename T>
    static void readArray(std::ifstream &in, std::vector<T> &v, const char *what)
    {
      uint64_t n = 0;
      in.read((char *)&n,8);
      if (!in.good())
        throw std::runtime_error(std::string("#hs.orbg: truncated reading ")+what);
      v.resize(n);
      if (n) in.read((char *)v.data(),n*sizeof(T));
      if (!in.good())
        throw std::runtime_error(std::string("#hs.orbg: truncated reading ")+what);
    }

    static uint32_t parseChannelMask(const std::string &spec,
                                     const std::vector<std::string> &names)
    {
      if (spec.empty() || spec == "all") return 0xffffffffu;
      uint32_t mask = 0;
      for (size_t pos=0; pos<spec.size(); ) {
        size_t comma = spec.find(',',pos);
        if (comma == std::string::npos) comma = spec.size();
        const std::string tok = spec.substr(pos,comma-pos);
        pos = comma+1;
        if (tok.empty()) continue;
        bool found = false;
        for (size_t c=0;c<names.size();c++)
          if (names[c] == tok ||
              names[c] == std::string(1,tok[tok.size()-1])+"-"+tok.substr(0,1)) {
            mask |= (1u<<c); found = true;
          }
        if (!found) {
          std::string all; for (auto &n : names) all += n+" ";
          throw std::runtime_error("#hs.orbg: unknown channel '"+tok
                                   +"'; this field has: "+all);
        }
      }
      return mask;
    }

    void OrbGridFromFile::executeLoad(OnePartition &dataGroup)
    {
      std::ifstream in(data.where,std::ios::binary);
      if (!in.good())
        throw std::runtime_error("#hs.orbg: could not open "+data.where);
      char magic[8]; in.read(magic,8);
      if (std::string(magic,7) != "ORBGRID")
        throw std::runtime_error("#hs.orbg: "+data.where+" is not a .orbg file");

      uint32_t version = 0, numChannels = 0;
      in.read((char *)&version,4);
      in.read((char *)&numChannels,4);
      int d[3]; float o[3], sp[3], vr[2];
      in.read((char *)d,12); in.read((char *)o,12);
      in.read((char *)sp,12); in.read((char *)vr,8);

      OrbGridVolume::SP vol = OrbGridVolume::create();
      vol->numChannels = (int)numChannels;
      vol->dims = vec3i(d[0],d[1],d[2]);
      vol->gridOrigin = vec3f(o[0],o[1],o[2]);
      vol->gridSpacing = vec3f(sp[0],sp[1],sp[2]);
      vol->valueRange = range1f(vr[0],vr[1]);
      vol->bounds.lower = vol->gridOrigin;
      vol->bounds.upper = vol->gridOrigin
        + vol->gridSpacing*vec3f(float(d[0]-1),float(d[1]-1),float(d[2]-1));

      readArray(in,vol->channelVec,"channel vectors");
      readArray(in,vol->channelPeak,"channel peaks");
      readArray(in,vol->subsetPeak,"subset peaks");
      {
        uint64_t n = 0; in.read((char *)&n,8);
        std::string blob(n,'\0'); if (n) in.read(&blob[0],n);
        size_t pos = 0;
        while (pos <= blob.size() && !blob.empty()) {
          size_t nl = blob.find('\n',pos);
          if (nl == std::string::npos) { vol->channelNames.push_back(blob.substr(pos)); break; }
          vol->channelNames.push_back(blob.substr(pos,nl-pos));
          pos = nl+1;
        }
      }
      if (vol->channelVec.size() != (size_t)numChannels*vol->numVoxels())
        throw std::runtime_error("#hs.orbg: channel array is "
                                 +std::to_string(vol->channelVec.size())
                                 +", expected numChannels*numVoxels");

      vol->channelMask = parseChannelMask(data.get("channels"),vol->channelNames);
      if (vol->channelMask < vol->subsetPeak.size()) {
        const float p = vol->subsetPeak[vol->channelMask];
        if (p > 0.f) vol->valueRange.upper = p;
      }
      vol->alphaGamma   = data.get_float("alphaGamma",1.f);
      vol->opacityScale = data.get_float("opacity",1.f);
      vol->domainScale  = data.get_float("domainScale",1.f);

      std::string chans;
      for (int c=0;c<vol->numChannels;c++)
        if (vol->channelMask & (1u<<c)) chans += vol->channelNames[c]+" ";
      std::cout << "#hs.orbg: " << vol->dims.x << "x" << vol->dims.y << "x"
                << vol->dims.z << " = " << prettyNumber(vol->numVoxels())
                << " voxels, " << numChannels << " channels [" << chans << "]"
                << ", bounds " << vol->bounds
                << ", range [0 " << vol->valueRange.upper << "]" << std::endl;

      dataGroup.orbGrids.push_back(vol);
    }

  }
}
