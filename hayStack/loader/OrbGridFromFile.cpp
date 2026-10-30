// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#include "hayStack/loader/OrbGridFromFile.h"
#include "hayStack/loader/MultiVariateContent.h"
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

    /*! a uint64 item count followed by the items; floatsPerItem > 1 reads
        e.g. an array of vec3f into a flat float vector */
    static void readArray(std::ifstream &in, std::vector<float> &v,
                          const char *what, int floatsPerItem = 1)
    {
      uint64_t n = 0;
      in.read((char *)&n,8);
      if (!in.good())
        throw std::runtime_error(std::string("#hs.orbg: truncated reading ")+what);
      v.resize(n*floatsPerItem);
      if (n) in.read((char *)v.data(),v.size()*sizeof(float));
      if (!in.good())
        throw std::runtime_error(std::string("#hs.orbg: truncated reading ")+what);
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

      MultiVariateVolume::SP vol = MultiVariateVolume::create();
      vol->numChannels   = (int)numChannels;
      vol->numComponents = 3;
      vol->unit          = "A/A^3";
      vol->dims          = vec3i(d[0],d[1],d[2]);
      vol->gridOrigin    = vec3f(o[0],o[1],o[2]);
      vol->gridSpacing   = vec3f(sp[0],sp[1],sp[2]);
      vol->valueRange    = range1f(vr[0],vr[1]);
      vol->bounds.lower  = vol->gridOrigin;
      vol->bounds.upper  = vol->gridOrigin
        + vol->gridSpacing*vec3f(float(d[0]-1),float(d[1]-1),float(d[2]-1));

      readArray(in,vol->channelData,"channel vectors",3);
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
      if (vol->channelData.size() != (size_t)numChannels*vol->numVoxels()*3)
        throw std::runtime_error("#hs.orbg: channel array is "
                                 +std::to_string(vol->channelData.size()/3)
                                 +" vectors, expected numChannels*numVoxels");

      vol->mix = parseChannelMix
        (data,vol->channelNames,
         [](const std::string &name, const std::string &tok) {
           return name == std::string(1,tok.back())+"-"+tok.substr(0,1);
         });
      const float p = vol->peakFor(vol->mix);
      if (p > 0.f)
        vol->valueRange = vol->mix.isSigned() ? range1f(-p,p) : range1f(0.f,p);
      vol->alphaGamma   = data.get_float("alphaGamma",1.f);
      vol->opacityScale = data.get_float("opacity",1.f);
      vol->domainScale  = data.get_float("domainScale",1.f);

      std::string chans;
      for (int c=0;c<vol->numChannels;c++)
        if (vol->mix.weights[c] != 0.f) chans += vol->channelNames[c]+" ";
      std::cout << "#hs.orbg: " << vol->dims.x << "x" << vol->dims.y << "x"
                << vol->dims.z << " = " << prettyNumber(vol->numVoxels())
                << " voxels, " << numChannels << " channels [" << chans << "]"
                << ", reduction " << reductionName(vol->mix.reduction)
                << ", bounds " << vol->bounds
                << ", range [" << vol->valueRange.lower << " "
                << vol->valueRange.upper << "]" << std::endl;

      dataGroup.multiVariates.push_back(vol);
    }

  }
}
