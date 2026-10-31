// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#include "hayStack/loader/OrbGridFromFile.h"
#include "hayStack/loader/MultiVariateContent.h"
#include <cctype>
#include <cstring>
#include <iostream>
#include <map>
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

    /*! "s0-d1" -> "s0","d1"; false if the name is not a pair of orbitals */
    static bool splitOrbitalPair(const std::string &name,
                                 std::string &a, std::string &b)
    {
      const size_t dash = name.find('-');
      if (dash == std::string::npos) return false;
      a = name.substr(0,dash);
      b = name.substr(dash+1);
      for (const std::string &o : { a, b }) {
        if (o.empty() || !std::strchr("spdfg",o[0])) return false;
        for (size_t i=1;i<o.size();i++)
          if (!std::isdigit((unsigned char)o[i])) return false;
      }
      return true;
    }

    /*! 'load=' keeps a pair named in either order, and every pair made only
        of the listed orbitals ("s1") or orbital types ("s") */
    static std::vector<int> pickChannels(const std::string &load,
                                         const std::vector<std::string> &names)
    {
      std::vector<int> pick;
      if (load.empty() || load == "all") {
        for (size_t c=0;c<names.size();c++) pick.push_back((int)c);
        return pick;
      }
      auto names_orbital = [](const std::string &tok, const std::string &orbital) {
        return tok == orbital || (tok.size() == 1 && tok[0] == orbital[0]);
      };
      const std::vector<std::string> tokens = splitCommas(load);
      std::vector<bool> used(tokens.size(),false);
      for (size_t c=0;c<names.size();c++) {
        std::string a, b;
        const bool isPair = splitOrbitalPair(names[c],a,b);
        bool keep = false, hasA = false, hasB = false;
        for (size_t t=0;t<tokens.size();t++) {
          if (tokens[t].find('-') != std::string::npos) {
            if (tokens[t] == names[c] || (isPair && tokens[t] == b+"-"+a))
              keep = used[t] = true;
          } else if (isPair) {
            if (names_orbital(tokens[t],a)) hasA = used[t] = true;
            if (names_orbital(tokens[t],b)) hasB = used[t] = true;
          }
        }
        if (keep || (hasA && hasB)) pick.push_back((int)c);
      }
      for (size_t t=0;t<tokens.size();t++)
        if (!used[t]) {
          std::string all; for (auto &n : names) all += n+" ";
          throw std::runtime_error("#hs.orbg: 'load' entry '"+tokens[t]
                                   +"' matches no channel; this file has: "+all);
        }
      return pick;
    }

    /*! where a .orbg keeps its channels, so any set can be read at any time */
    struct OrbGridIndex {
      std::string              path;
      std::streamoff           channelBase = 0;
      size_t                   floatsPerChannel = 0;
      std::vector<std::string> names;
      std::vector<float>       peak, subsetPeak;
    };

    static void readChannels(const OrbGridIndex &file, MultiVariateVolume &vol,
                             const std::vector<int> &pick)
    {
      std::ifstream in(file.path,std::ios::binary);
      if (!in.good())
        throw std::runtime_error("#hs.orbg: could not open "+file.path);
      std::vector<float> data(pick.size()*file.floatsPerChannel);
      std::vector<std::string> names;
      std::vector<float> peaks;
      for (size_t i=0;i<pick.size();i++) {
        names.push_back(file.names[pick[i]]);
        if ((size_t)pick[i] < file.peak.size()) peaks.push_back(file.peak[pick[i]]);
        in.seekg(file.channelBase
                 +(std::streamoff)(pick[i]*file.floatsPerChannel*sizeof(float)));
        in.read((char *)(data.data()+i*file.floatsPerChannel),
                file.floatsPerChannel*sizeof(float));
        if (!in.good())
          throw std::runtime_error("#hs.orbg: truncated reading channel "
                                   +file.names[pick[i]]);
      }
      vol.channelData.swap(data);
      vol.channelNames   = names;
      vol.channelPeak    = peaks;
      vol.loadedChannels = pick;
      vol.numChannels    = (int)pick.size();
      /* subset masks count file channels, so they only apply to the full set */
      vol.subsetPeak = pick.size() == file.names.size()
        ? file.subsetPeak : std::vector<float>{};
    }

    /*! tags each file channel with its two orbitals and groups them by type, s
        before p before d, so viewers can switch orbitals instead of pairs */
    static void tagOrbitals(MultiVariateVolume &vol)
    {
      static const std::string types = "spdfg";
      std::map<std::pair<int,int>,std::string> labels;
      for (auto &name : vol.fileChannelNames) {
        std::string a, b;
        if (!splitOrbitalPair(name,a,b)) return;
        for (const std::string &o : { a, b })
          labels[{ (int)types.find(o[0]),
                   o.size() > 1 ? std::stoi(o.substr(1)) : -1 }] = o;
      }
      std::map<std::string,int> index;
      int group = -1, prevType = -1;
      for (auto &l : labels) {
        if (l.first.first != prevType) {
          prevType = l.first.first;
          group++;
          vol.groupNames.push_back(std::string(1,(char)std::toupper(types[prevType]))
                                   +"-orbital");
        }
        index[l.second] = (int)vol.labelNames.size();
        vol.labelNames.push_back(l.second);
        vol.labelGroup.push_back(group);
      }
      for (auto &name : vol.fileChannelNames) {
        std::string a, b;
        splitOrbitalPair(name,a,b);
        vol.channelLabels.push_back(a == b ? std::vector<int>{ index[a] }
                                    : std::vector<int>{ index[a], index[b] });
      }
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
      vol->numComponents = 3;
      vol->unit          = "A/A^3";
      vol->dims          = vec3i(d[0],d[1],d[2]);
      vol->gridOrigin    = vec3f(o[0],o[1],o[2]);
      vol->gridSpacing   = vec3f(sp[0],sp[1],sp[2]);
      vol->valueRange    = range1f(vr[0],vr[1]);
      vol->bounds.lower  = vol->gridOrigin;
      vol->bounds.upper  = vol->gridOrigin
        + vol->gridSpacing*vec3f(float(d[0]-1),float(d[1]-1),float(d[2]-1));

      /* the channel block sits in front of the names, so skip over it first */
      uint64_t numVectors = 0;
      in.read((char *)&numVectors,8);
      const std::streamoff channelBase = in.tellg();
      const size_t numVoxels = vol->numVoxels();
      if (!in.good() || numVectors != (uint64_t)numChannels*numVoxels)
        throw std::runtime_error("#hs.orbg: channel array is "
                                 +std::to_string(numVectors)
                                 +" vectors, expected numChannels*numVoxels");
      in.seekg(numVectors*3*sizeof(float),std::ios::cur);

      OrbGridIndex file;
      file.path             = data.where;
      file.channelBase      = channelBase;
      file.floatsPerChannel = numVoxels*3;
      std::vector<std::string> &fileNames = file.names;
      readArray(in,file.peak,"channel peaks");
      readArray(in,file.subsetPeak,"subset peaks");
      {
        uint64_t n = 0; in.read((char *)&n,8);
        std::string blob(n,'\0'); if (n) in.read(&blob[0],n);
        size_t pos = 0;
        while (pos <= blob.size() && !blob.empty()) {
          size_t nl = blob.find('\n',pos);
          if (nl == std::string::npos) { fileNames.push_back(blob.substr(pos)); break; }
          fileNames.push_back(blob.substr(pos,nl-pos));
          pos = nl+1;
        }
      }
      if (fileNames.size() != numChannels)
        throw std::runtime_error("#hs.orbg: "+std::to_string(fileNames.size())
                                 +" channel names for "
                                 +std::to_string(numChannels)+" channels");

      const std::vector<int> pick = pickChannels(data.get("load"),fileNames);
      if (pick.empty())
        throw std::runtime_error("#hs.orbg: 'load' selects no channels");
      vol->fileChannels     = (int)numChannels;
      vol->fileChannelNames = fileNames;
      readChannels(file,*vol,pick);
      vol->loadChannels = [file](MultiVariateVolume &v,const std::vector<int> &channels)
      { readChannels(file,v,channels); };
      tagOrbitals(*vol);

      vol->mix = parseChannelMix
        (data,vol->channelNames,
         [](const std::string &name, const std::string &tok) {
           const size_t dash = tok.find('-');
           return dash == std::string::npos
             ? name == std::string(1,tok.back())+"-"+tok.substr(0,1)
             : name == tok.substr(dash+1)+"-"+tok.substr(0,dash);
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
                << " voxels, " << vol->numChannels << " of " << numChannels
                << " channels loaded, enabled [" << chans << "]"
                << ", reduction " << reductionName(vol->mix.reduction)
                << ", bounds " << vol->bounds
                << ", range [" << vol->valueRange.lower << " "
                << vol->valueRange.upper << "]" << std::endl;

      dataGroup.multiVariates.push_back(vol);
    }

  }
}
