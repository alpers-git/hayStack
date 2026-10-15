// SPDX-FileCopyrightText: Copyright (c) 2023-2026 Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hayStack/OnePartition.h"
#include "hayMaker/TextureLibrary.h"
#include "hayMaker/MaterialLibrary.h"
#include "hayStack/TransferFunction.h"

namespace hs {
  struct NanoVDBVolume;
};

namespace hm {
  using hs::OnePartition;
  using hs::TransferFunction;
  using hs::VolumeScatterParams;
  using hs::VolumeScatterSettings;
  
  struct HayMaker;
  
  /*! implements rendering operations for one logical partition, on
    one single device */
  struct AnariDeviceRenderer {
      
    AnariDeviceRenderer(int gpuID,
                        int tetherIndex,
                        int tetherCount,
                        HayMaker     *hayMaker,
                        OnePartition *myPartition);
    void renderInitialAnariWorld();    
    void renderFrame();

    struct {
      std::vector<affine3f>     xfms;
      std::vector<anari::Group> groups;
    } rootInstances;
    anari::Group rootGroup;

    std::vector<anari::Volume>  rootVolumes;
    std::vector<anari::Surface> rootGeoms;
    std::map<anari::Volume, VolumeScatterParams> principledScatterByVolume;
      
    anari::Light envLight;
    std::vector<anari::Light> lights;

    anari::Group volumeGroup = 0;
    HayMaker     *const hayMaker;
    OnePartition *const myPartition;
    
    void setCamera(const hs::Camera &camera);
    
    void setTransferFunction(const TransferFunction &xf);

    /*! re-set the iso-value of every iso-surface created from a
        'gpuIso=' volume, in the order they were created. Values
        beyond the number of surfaces we have are ignored. */
    void setIsoValues(const std::vector<float> &values);
    void setChannelFilter(const hs::ChannelFilter &filter);
      
    // void renderAll();
    void renderMiniScene(mini::Scene::SP miniScene);
    anari::Group render(const mini::Object::SP &miniObject);

    void applyTransferFunction(const TransferFunction &xf);
    void applyDefaultPrincipledTransferFunction(anari::Volume volume);
    void applyPrincipledScatterParams(anari::Volume volume,
                                      const VolumeScatterParams &scatter);
    void applyVolumeScatterSettings(const VolumeScatterSettings &settings);
    VolumeScatterSettings volumeScatterSettings;
    void setVolumeScatterSettings(const VolumeScatterSettings &settings);

    /*! creates *default* color mapper that maps from global scalar
        min/max (across all per-vertex arrays and scalarfield types),
        to the chosen default color map */
    void createDefaultColorMapper(/*! input range across all objects,
                                    across all rankgs */
                                  const range1f &inRange,
                                  const std::vector<vec4f> &outColors);
    
    /*! default scalar to color mapper - this spans the total input
        range of all scalar fields in the model (globally reduced
        across all ranks), for the selected haystack's default color
        map. if input range is empty (ie, no scalar field or
        per-vertex scalar attributes this becomes 0 */
    anari::Sampler defaultColorMapper = 0;
    
    // ==================================================================
    // create() functions: create anari object for given piece of
    // haystack input //
    // ==================================================================
    
    anari::Group createGroup(const std::vector<anari::Surface> &geoms,
                             const std::vector<anari::Volume>  &volumes);
    /*! returns either a volume (the regular case) or, if the volume
        asked for iso-surfaces, a list of iso-surface surfaces */
    std::pair<anari::Volume,std::vector<anari::Surface>>
    create(const hs::StructuredVolume &vol);
    anari::Volume create(const hs::NanoVDBVolume &vol);
    std::pair<anari::Volume,std::vector<anari::Surface>>
    create(const hs::TAMRVolume &input);
    
    anari::Volume
    create(const std::pair<umesh::UMesh::SP,box3f> &meshAndDomain);
    /*! the orbital-current field on a dense grid: per-channel voxel vectors
      plus the enabled-channel mask */
    anari::Volume create(const hs::OrbGridVolume::SP &orb);

    /*! the orbGrid spatial fields we created, so a channel toggle can be
        pushed without rebuilding the scene */
    std::vector<anari::SpatialField> orbGridFields;

    using Surfaces = std::vector<anari::Surface>;
    Surfaces create(const hs::SphereSet &content);
    Surfaces create(const hs::TriangleMesh &content);
    Surfaces create(const hs::Capsules &caps);
    Surfaces create(const hs::Cylinders &content);
    
    anari::Surface create(const mini::Mesh::SP &mesh);

    void createAndAdd(const mini::QuadLight &ml);
    void createAndAdd(const mini::DirLight &ml);
    void createAndAdd(const mini::EnvMapLight &ml);

    // ==================================================================
    // materials
    // ==================================================================
    std::pair<anari::Material,std::string> create(mini::Metal::SP);
    std::pair<anari::Material,std::string> create(mini::MetallicPaint::SP);
    std::pair<anari::Material,std::string> create(mini::Matte::SP);
    std::pair<anari::Material,std::string> create(mini::ANARIMaterial::SP);
    std::pair<anari::Material,std::string> create(mini::DisneyMaterial::SP);
    std::pair<anari::Material,std::string> create(mini::Dielectric::SP);
    std::pair<anari::Material,std::string> create(mini::Plastic::SP);
    std::pair<anari::Material,std::string> create(mini::Material::SP);

    // anari::Sampler createScalarMapper(const std::vector<float> &scalars);
    

    
    TextureLibrary  textureLibrary;
    MaterialLibrary materialLibrary;

    TransferFunction currentXF;
#if HS_USE_MULTI_SCATTERING
    VolumeScatterSettings volumeScatterSettings;
    void setVolumeScatterSettings(const VolumeScatterSettings &settings)
    {
      volumeScatterSettings = settings;
      for (auto &entry : principledScatterByVolume)
        entry.second = settings.medium;
      this->applyVolumeScatterSettings(settings);
    }
#endif


    void setInstances(const std::vector<anari::Group> &groups,
                      const std::vector<affine3f> &xfms);
    void setLights(anari::Group rootGroup,
                   const std::vector<anari::Light> &lights);

    bool dirty = true;
    vec2i fbSize { -1,-1 };

    /*! every iso-surface geometry we created, in creation order, so
        the UI can re-set their iso-values at runtime */
    std::vector<anari::Geometry> isoGeoms;

    struct {
      anari::Device device;
      anari::World  world;
      anari::Renderer renderer;
      anari::Frame    frame;
      anari::Camera   camera;
    } anari;
  };
    

}
