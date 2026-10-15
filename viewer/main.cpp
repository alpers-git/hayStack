// SPDX-FileCopyrightText: Copyright (c) 2023++ Ingo Wald
// SPDX-License-Identifier: Apache-2.0

#include "hayMaker/HayMaker.h"
#include "hayStack/loader/DataLoader.h"
#include "hayStack/loader/RAWVolumeContent.h"
#if HS_IMGUI
# include "owlViewer/OWLViewer.h"
# include "owlViewer/XFEditor.h"
# include "OrbitalPanel.h"
#else
# define STB_IMAGE_WRITE_IMPLEMENTATION 1
# define STB_IMAGE_IMPLEMENTATION 1
#endif
#include "stb/stb_image.h"
#include "stb/stb_image_write.h"
#if HS_MPI
#include <unistd.h>
#endif
namespace hm {

  double t_last_render;
  
  typedef enum { DPMODE_NOT_SPECIFIED,
                 DPMODE_DATA_PARALLEL,
                 DPMODE_DATA_REPLICATED } DPMode;

  struct CmdLineCamera {
    vec3f vp   = vec3f(0.f);
    vec3f vi   = vec3f(0.f);
    vec3f vu   = vec3f(0.f,1.f,0.f);
    float fovy = 60.f;
  };
  
  struct FromCL {
    /*! data groups per rank. '0' means 'auto - use few as we can, as
      many as we have to fit for given number of ranks */
    int dpr = 0;
    /*! num data groups. '0' will become '1' by default, but allows us
        for specifying 'hasn't been specified', which in turn allows
        the mpi mode to set it to either '1' or '1 per rank' depending
        on choesn dpmode */
    int ndg = 0;

    bool forceSingleGPU = false;
    
    DPMode dpMode = DPMODE_NOT_SPECIFIED;
    
    /*! which color map to use for color mapping (if applicable) */
    int cmID = 0;
    
    bool mergeUnstructuredMeshes = false;
    vec4f bgColor { NAN, NAN, NAN, NAN };
    float ambientRadiance = .6f;
    std::string xfFileName = "";
    std::string outFileName = "hayStack.png";
    vec2i fbSize = { 800,600 };
    bool createHeadNode = false;
    int  numExtraDisplayRanks = 0;
    int  numFramesAccum = 1;
    bool denoise = true;
    int  spp            = 1;
    bool verbose = true;
    CmdLineCamera camera;
    std::vector<CmdLineCamera> cameraPath;
    // struct {
    //   vec3f vp0, vp1;
    //   vec3f vi0, vi1;
    //   int numSteps = 0;
    // } cameraPath;
    bool measure = 0;
    /*! the positional data arguments, verbatim, so the '/' key can print a
        command line that reproduces the scene */
    std::vector<std::string> contentArgs;
    /*! toggle this many channels off and back on before rendering. The filter
        state ends where it started, so the image must be identical to one
        rendered without it - a regression test for state that survives a
        filter change when it should not. */
    int  toggleCycles = 0;
    std::string envMapFileName;
    /*! initial iso-values, harvested from the 'gpuIso=' key of the
        data URLs; one UI slider gets created per entry */
    std::vector<float> isoValues;
  };
  FromCL fromCL;
  
  // std::string FromCL::outFileName = "hayStack.png";
  // bool FromCL::measure = 0;
  // bool FromCL::verbose = true;

  inline bool verbose() { return fromCL.verbose; }
  
  void usage(const std::string &error="")
  {
    std::cout << "./hs{Offline,Viewer,ViewerQT} ... <args>" << std::endl;
    std::cout << "w/ args:" << std::endl;
    std::cout << "-xf file.xf   ; specify transfer function" << std::endl;
    if (!error.empty())
      throw std::runtime_error("fatal error: " +error);
    exit(0);
  }

#if HS_IMGUI
  struct Viewer : public owl::viewer::OWLViewer
  {
    Viewer(BoundsData worldBounds,
           RenderEngineInterface *const renderer,
           hs::mpi::Comm *world,
           const hs::OrbGridVolume::SP &orbGrid = {})
      : renderer(renderer),
        world(world),
        worldBounds(worldBounds)
    {
      xf.absDomain = worldBounds.scalars;
      isoValues = fromCL.isoValues;
      if (orbGrid) {
        this->orbGrid = orbGrid;
        orbPanel = new hs::OrbitalPanel(orbGrid);
        /* seed only the DOMAIN and density from the field: the world scalar
           range is the wrong framing for an orbital field, but the colour map
           is the XFEditor's business. pushing one in here overrode the
           editor's own maps with an invented one. */
        const hs::TransferFunction seed
          = hs::makeSignedTransferFunction(orbGrid->valueRange,
                                           orbGrid->bounds,
                                           orbGrid->alphaGamma,
                                           orbGrid->opacityScale,
                                           orbGrid->domainScale);
        xf.absDomain   = seed.absDomain;
        xf.relDomain   = seed.relDomain;
        xf.baseDensity = seed.baseDensity;
        /* the map the offline renders use, taken from the editor's own library
           rather than rebuilt by hand */
        useBuiltinColorMap("cool_warm_extended");
        filterDirty = true;
      }
    }

    const BoundsData worldBounds;
    hs::OrbGridVolume::SP orbGrid;
    /*! how many scenes '/' has dumped this session, so each gets its own .xf */
    int sceneDumpCount = 0;
    
    bool uiIsOn = true;
    
    void runImgui() override
    {
      if (!uiIsOn) return;
      
      if (orbPanel) {
        /* this mode owns the volume, so the generic scalar-range and iso
           controls below would only duplicate - and mis-frame - what the
           panel already does */
        /* the gizmo draws no window of its own, and without one it lands in
           imgui's implicit Debug window - nowhere near the panel that owns the
           domain. Give it a window and put the value axis directly under the
           ramp, where the colours it labels actually are. run_ui leaves the
           cursor just below its canvas and the canvas spans the full content
           width, so the ticks line up with the ramp. */
        ImGui::Begin("Transfer function");
        if (xfEditor.run_ui())
          xfDirty = true;
        orbPanel->drawRuler();
        ImGui::End();
        if (orbPanel->run_ui()) {
          xf.absDomain = orbPanel->domain();
          xf.relDomain = orbPanel->relDomain();
          xf.baseDensity = orbPanel->baseDensity();
          xfDirty = true;
        }
        if (orbPanel->filterCommitted())
          filterDirty = true;
        return;
      }

      ImGui::Begin("HayStack/Imgui");
      
      if (worldBounds.scalars.empty()) {
        ImGui::Text("(no volume data)");
      } else {
        ImGui::Text("scalar range [%f .. %f]",
                    worldBounds.scalars.lower,
                    worldBounds.scalars.upper);
        if (xfEditor.run_ui())
          xfDirty = true;
        
        ImGui::Text("scalar range [%f .. %f]",
                    worldBounds.scalars.lower,
                    worldBounds.scalars.upper);
        if (ImGui::SliderFloat("xfDomain(abs).x", &xf.absDomain.lower,
                               worldBounds.scalars.lower,
                               worldBounds.scalars.upper))
          xfDirty = true;
        if (ImGui::SliderFloat("xfDomain(abs).y", &xf.absDomain.upper, 
                               worldBounds.scalars.lower,
                               worldBounds.scalars.upper))
          xfDirty = true;
        if (ImGui::SliderFloat("xfDomain(rel).x", &xf.relDomain.lower,
                               0.f, 100.f))
          xfDirty = true;
        if (ImGui::SliderFloat("xfDomain(rel).y", &xf.relDomain.upper,
                               0.f, 100.f))
          xfDirty = true;
        if (ImGui::SliderFloat("baseDensity", &xf.baseDensity, 
                               0.f,200.f))
          xfDirty = true;
      }
      if (!isoValues.empty()) {
        ImGui::Separator();
        ImGui::Text("iso-surfaces");
        /* iso-values live in scalar-field units, so the field's own
           range is the only sensible slider range */
        for (size_t i=0;i<isoValues.size();i++) {
          char label[64];
          snprintf(label,sizeof(label),"isoValue[%zu]",i);
          if (ImGui::SliderFloat(label,&isoValues[i],
                                 worldBounds.scalars.lower,
                                 worldBounds.scalars.upper))
            isoDirty = true;
        }
      }
      ImGui::End();
    }

    // void colorMapChanged();
    // void rangeChanged();
    // void opacityScaleChanged();
#if HS_USE_MULTI_SCATTERING
    void scatterSettingsChanged(const VolumeScatterSettings &settings);
#endif

    void screenShot()
    {
      std::string fileName = fromCL.outFileName;
      std::cout << "#ht: saving screen shot in " << fileName << std::endl;
      stbi_flip_vertically_on_write(true);
      stbi_write_png(fileName.c_str(),fbSize.x,fbSize.y,4,
                     fbPointer,fbSize.x*sizeof(uint32_t));
    }
    
    /*! this gets called when the user presses a key on the keyboard ... */
    /*! Print an hsOffline command line that reproduces exactly what is on
        screen: the same data, the same camera, and the same look.

        The transfer function is written out and referenced with -xf rather than
        re-derived from url options: the colour map, the domain and the density
        are all live editor state, and there is no set of url options that
        reconstructs an edited map. The orbital channel selection *is* a url
        option, so that is appended to the field's url. */
    void dumpSceneConfig()
    {
      hs::Camera cam;
      OWLViewer::getCameraOrientation((owl::common::vec3f&)cam.vp,
                                      (owl::common::vec3f&)cam.vi,
                                      (owl::common::vec3f&)cam.vu,
                                      cam.fovy);

      char xfName[256];
      snprintf(xfName,sizeof(xfName),"scene_%03d.xf",sceneDumpCount);
      xf.save(xfName);

      /* the channel selection, as ':channels=' understands it */
      std::string chanOpt;
      if (orbGrid && !orbGrid->channelNames.empty() && orbPanel) {
        const hs::ChannelFilter &f = orbPanel->getFilter();
        std::string names;
        size_t on = 0;
        for (size_t c=0;c<orbGrid->channelNames.size();c++)
          if (f.channelMask & (1u<<c)) {
            if (!names.empty()) names += ",";
            names += orbGrid->channelNames[c];
            on++;
          }
        /* the .orbg loader has no spelling for an empty selection */
        if (on > 0 && on != orbGrid->channelNames.size())
          chanOpt = ":channels=" + names;
      }

      const owl::common::vec2i sz = getWindowSize();
      std::cout << "\n#hs: ---- scene " << sceneDumpCount
                << " (transfer function in " << xfName << ") ----\n";
      std::cout << "./hsOffline";
      for (auto &a : fromCL.contentArgs) {
        const bool isField
          = a.find(".orbg") != std::string::npos
          || a.find("orbgrid://") != std::string::npos;
        std::cout << " \\\n    \"" << a
                  << (isField ? chanOpt : std::string()) << "\"";
      }
      std::cout << " \\\n    --camera "
                << cam.vp.x << " " << cam.vp.y << " " << cam.vp.z << "  "
                << cam.vi.x << " " << cam.vi.y << " " << cam.vi.z << "  "
                << cam.vu.x << " " << cam.vu.y << " " << cam.vu.z
                << " -fovy " << cam.fovy;
      std::cout << " \\\n    -res " << sz.x << " " << sz.y
                << " -spp " << fromCL.spp
                << " -xf " << xfName;
      if (!std::isnan(fromCL.bgColor.x))
        std::cout << " --bg-color " << fromCL.bgColor.x << " " << fromCL.bgColor.y
                  << " " << fromCL.bgColor.z << " " << fromCL.bgColor.w;
      std::cout << " \\\n    -o scene_" << sceneDumpCount << ".png"
                << std::endl << std::endl;
      sceneDumpCount++;
    }

    void key(char key, const owl::common::vec2i &where) override
    {
      switch(key) {
      case 'u':
        uiIsOn = !uiIsOn;
        break;
      case '!':
        screenShot();
        break;
      case '*': {
        hs::Camera camera;
        OWLViewer::getCameraOrientation
          ((owl::common::vec3f&)camera.vp,
           (owl::common::vec3f&)camera.vi,
           (owl::common::vec3f&)camera.vu,
           camera.fovy);
      } break;
      case 'P': {
          char *fl = getenv("BARNEY_FOCAL_LENGTH");
          std::cout << "export BARNEY_FOCAL_LENGTH=" << (fl != nullptr ? fl : "0") << std::endl;
          char *lr = getenv("BARNEY_LENS_RADIUS");
          std::cout << "export BARNEY_LENS_RADIUS=" << (lr != nullptr ? lr : "0") << std::endl;
          break;
        }
      case '9': case '0':
      case '(': case ')': {
        char* flc = getenv("BARNEY_FOCAL_LENGTH");
        float ratio = (key == '(' || key == '9') ? (0.9f) : (1.1f);
        if (key == '9' || key == '0')
          ratio = (ratio - 1.f) * 0.1f + 1.f; 
        float fl = (flc == nullptr) ? 1.f : std::stof(flc) * ratio;
#ifdef _WIN32
        _putenv_s("BARNEY_FOCAL_LENGTH", std::to_string(fl).c_str());
#else        
        setenv("BARNEY_FOCAL_LENGTH", std::to_string(fl).c_str(), 1);
#endif        
        std::cout << "Focal length is set to: " << fl << '\n';
        renderer->resetAccumulation();
        break;
      }
      case '{': case '}':
      case '[': case ']': {
        char* lrc = getenv("BARNEY_LENS_RADIUS");
        float ratio = (key == '[' || key == '{') ? (0.9f) : (1.1f);
        if (key == '[' || key == ']')
          ratio = (ratio - 1.f) * 0.1f + 1.f; 
        float lr = (lrc == nullptr) ? 1.f : std::stof(lrc) * ratio;
#ifdef _WIN32
        _putenv_s("BARNEY_LENS_RADIUS", std::to_string(lr).c_str());
#else        
        setenv("BARNEY_LENS_RADIUS", std::to_string(lr).c_str(), 1);
#endif        
        std::cout << "Lens radius is set to: " << lr << '\n';
        renderer->resetAccumulation();
        break;
      }

      case '/':
        dumpSceneConfig();
        break;

      case 'T':
        std::cout << "(T) : dumping transfer function" << std::endl;
        xf.save("hayThere.xf");
        break;

      default: OWLViewer::key(key,where);
      };
    }

    /*! window notifies us that we got resized. We HAVE to override
      this to know our actual render dimensions, and get pointer
      to the device frame buffer that the viewer cated for us */
    void resize(const owl::common::vec2i &newSize) override
    {
      OWLViewer::resize(newSize);
      renderer->resize((const mini::common::vec2i&)newSize,fbPointer);
    }

    /*! gets called whenever the viewer needs us to re-render out widget */
    void render() override
    {
      double _t0 = mini::common::getCurrentTime();
      
      if (xfDirty) {
        auto cm = xfEditor.getColorMap();
        xf.set((const mini::vec4f*)cm.data(),
                        cm.size());
        
        
        renderer->setTransferFunction(xf);
        xfDirty = false;
        accumDirty = true;
      }

      if (isoDirty) {
        renderer->setIsoValues(isoValues);
        isoDirty = false;
        accumDirty = true;
      }

      if (filterDirty && orbPanel) {
        renderer->setChannelFilter(orbPanel->getFilter());
        filterDirty = false;
        accumDirty = true;
      }

#if HS_USE_MULTI_SCATTERING
      if (scatterDirty) {
        renderer->setVolumeScatterSettings(scatterSettings);
        scatterDirty = false;
        accumDirty = true;
      }
#endif

      if (accumDirty) {
        renderer->resetAccumulation();
        accumDirty = false;
      }

      static int numFramesRendered = 0;
      const int measure_warmup_frames = 2;
      const int measure_max_frames = 100;
      const float measure_max_seconds = 60.f;
      
      static double measure_t0 = 0.;
      if (numFramesRendered == measure_warmup_frames)
        measure_t0 = mini::common::getCurrentTime();

      static double t0 = mini::common::getCurrentTime();
      renderer->renderFrame();
      ++numFramesRendered;
      double t1 = mini::common::getCurrentTime();

      if (fromCL.measure) {
        int numFramesMeasured = numFramesRendered - measure_warmup_frames;
        float numSecondsMeasured
          = (numFramesMeasured < 1)
          ? 0.f
          : float(t1 - measure_t0);

        if (numFramesMeasured >= measure_max_frames ||
            numSecondsMeasured >= measure_max_seconds) {
          std::cout << "measure: rendered " << numFramesMeasured << " frames in " << numSecondsMeasured << ", that is:" << std::endl;
          std::cout << "FPS " << double(numFramesMeasured/numSecondsMeasured) << std::endl;
          screenShot();
          renderer->terminate();
          
          world->barrier();
          hs::mpi::finalize();
          exit(0);
        }
      }

      static double sum_t = 0.f;
      static double sum_w = 0.f;
      sum_t = 0.8f*sum_t + (t1-t0);
      sum_w = 0.8f*sum_w + 1.f;
      float timePerFrame = float(sum_t / sum_w);
      float fps = 1.f/timePerFrame;
      std::string title = "HayThere ("+mini::common::prettyDouble(fps)+"fps)";
      setTitle(title.c_str());
      t0 = t1;
      double _t1 = mini::common::getCurrentTime();
      t_last_render = _t1-_t0;
    }
    
    void cameraChanged() override
    {  
      hs::Camera camera;
      OWLViewer::getCameraOrientation
        ((owl::common::vec3f&)camera.vp,
         (owl::common::vec3f&)camera.vi,
         (owl::common::vec3f&)camera.vu,
         camera.fovy);
      renderer->setCamera(camera);
      accumDirty = true;
    }

    //    TransferFunction xf;
    bool xfDirty = true;
    /*! live iso-values, one per 'gpuIso=' iso-surface */
    std::vector<float> isoValues;
    bool isoDirty = false;
#if HS_USE_MULTI_SCATTERING
    VolumeScatterSettings scatterSettings;
    bool scatterDirty = true;
    VolumeScatterPanel *scatterPanel = 0;
#endif
    bool accumDirty = true;
    /*! non-null only when the scene has an orbital-current grid field; when
        it is set the panel replaces the generic volume/iso ui entirely */
    hs::OrbitalPanel *orbPanel = nullptr;
    bool filterDirty = false;
    RenderEngineInterface *const renderer;
    hs::mpi::Comm *world;

    // the ui editor, that edits the color and alpha map (and nothing else)
    owl::viewer::XFEditor xfEditor;
    // the ui elements to edit abs/rel domains
    // range1f xfDomainAbs{0.f,1.f};
    // range1f xfDomainRel{0.f,100.f};
    float xfDensity = 1.f;
    // the objec that stores the _haystack_ idea of what a xf is (with
    // domain, color, and alpha)
    hs::TransferFunction xf;
    
    /*! pick one of the XFEditor's own maps by name. The editor ships
        cool_warm_extended among others; seeding it with a ramp written
        elsewhere replaces its library with an invented map that then cannot be
        got back from the combo box. */
    void useBuiltinColorMap(const char *want)
    {
      const std::vector<std::string> names
        = owl::viewer::XFEditor::colorMapLibrary.getNames();
      for (size_t i=0;i<names.size();i++)
        if (names[i] == want) { xfEditor.changeColorScheme((int)i); return; }
    }

    void loadXF() {
      xf.load(fromCL.xfFileName);
      /* an orbital field takes only the DOMAIN and density from the file and
         keeps the editor's own map (see useBuiltinColorMap): a .xf written by
         make_xf.py carries a hand-built ramp, and pushing it in here made the
         gizmo open on a map that is not one of its own. Every other scene gets
         the file's colour map, as before. */
      if (!orbPanel)
        xfEditor.setColorAndAlpha((const owl::common::vec4f*)xf.colorMap.data(),
                                  xf.colorMap.size());
      xfDirty = true;
    }
    
  };
#endif

#if HS_IMGUI
#if HS_USE_MULTI_SCATTERING
  void Viewer::scatterSettingsChanged(const VolumeScatterSettings &settings)
  {
    scatterSettings = settings;
    scatterDirty = true;
  }
#endif
#endif

  inline mini::common::vec3f get3f(char **av, int &i)
  {
    float x = std::stof(av[++i]);
    float y = std::stof(av[++i]);
    float z = std::stof(av[++i]);
    return mini::common::vec3f(x,y,z);
  }

  size_t computeHashFromString(const char *s)
  {
    size_t hash = 0;
    size_t FNV_PRIME = 0x00000100000001b3ull;
    for (int i=0;s[i] != 0;i++)
      hash = hash * FNV_PRIME ^ s[i];
    return hash;
  }
  
  std::vector<int> parseCommaSeparatedListOfInts(std::string s)
  {
    std::vector<std::string> tokens;
    int last = 0;
    while (true) {
      int pos = s.find(",");
      if (pos == s.npos) {
        tokens.push_back(s);
        break;
      } else {
        tokens.push_back(s.substr(0,pos));
        s = s.substr(pos+1);
      }
    }
    
    std::vector<int> result;
    for (auto s : tokens)
      if (s != "")
        result.push_back(std::stoi(s));
    return result;
  }

#if HS_MPI
  void determineLocalProcessID(mpi::Comm &world, int &localRank, int &localSize)
  {
# if HS_FAKE_MPI
    localRank = 0;
    localSize = 1;
# else
    world.barrier();
    std::vector<char> hostName(10000);
    gethostname(hostName.data(),hostName.size());
    size_t hash = computeHashFromString(hostName.data());
    std::vector<size_t> allHostNames(world.size);
    MPI_Allgather(&hash,sizeof(hash),MPI_BYTE,
                  allHostNames.data(),sizeof(hash),MPI_BYTE,
                  world.comm);
    localRank = 0;
    localSize = 0;
    for (int i=0;i<world.size;i++) {
      if (allHostNames[i] != hash) continue;
      localSize++;
      if (i < world.rank) localRank++;
    }

    for (int i=0;i<world.size;i++) {
      world.barrier();
      if (i != world.rank) continue;
      std::cout << "#hs(" << world.rank << "): determined local rank/size as "
                << localRank << "/" << localSize << std::endl;
    }
#endif
  }
#endif

  int getIntFromEnv(const char *varName, int fallback)
  {
    const char *var = getenv(varName);
    if (!var) return fallback;
    return std::stoi(var);
  }

  std::vector<int> selectGPUs(mpi::Comm &world, int localRank, int localSize)
  {
#ifdef __APPLE__
    return { -1 };
#endif
    std::cout << "#hs(" << world.rank << ":" << localRank << "/" << localSize << "): selecting GPUs ... " << std::endl;
    const char *hcd = getenv("HS_CUDA_DEVICES");
    if (hcd) {
      std::cout << "#hs(" << world.rank << "): found HS_CUDA_DEVICES, using this" << std::endl;
      return parseCommaSeparatedListOfInts(hcd);
    }
    const char *cvd = getenv("CUDA_VISIBLE_DEVICES");
    if (hcd) {
      std::cout << "#hs(" << world.rank << "): found CUDA_VISIBLE_DEVICES being set in env, using GPUs 0,1,... etc" << std::endl;
      auto cudaGPUs = parseCommaSeparatedListOfInts(hcd);
      std::vector<int> result;
      for (int i=0;i<cudaGPUs.size();i++)
        result.push_back(i);
      return result;
    }

    if (world.size == 1) {
      return { -1 };
    }
    
    if (localSize == 1) {
      std::cout << "#hs: single data rank per process, no dedicated gpus specified in env -> selecting gpuIDs={-1} to indicate 'any you can'" << std::endl;
      return { -1 };
    }
    
    return { localRank % localSize };
    // return { 0 };
  }
    
  inline float lerp_l(float f, float a, float b)
  { return (1.f-f)*a + f*b; }
  inline mini::common::vec3f lerp_l(float f,
                                    mini::common::vec3f a,
                                    mini::common::vec3f b)
  { return (1.f-f)*a + f*b; }
  
  void addCameraPath(int numSteps,
                     CmdLineCamera c0,
                     CmdLineCamera c1)
  {
    for (int i=0;i<numSteps;i++) {
      float f = i/(numSteps-1.f);
      CmdLineCamera c;
      c.vp = lerp_l(f,c0.vp,c1.vp);
      c.vi = lerp_l(f,c0.vi,c1.vi);
      c.vu = lerp_l(f,c0.vu,c1.vu);
      c.fovy = lerp_l(f,c0.fovy,c1.fovy);
      fromCL.cameraPath.push_back(c);
    }
  }

  void addCamerasFromFile(const std::string &fileName)
  {
    std::ifstream in(fileName.c_str());
    std::string line;
    while (in.good()) {
      std::getline(in,line);
      CmdLineCamera c;
      int rc = sscanf(line.c_str(),"--camera %f %f %f %f %f %f %f %f %f --fovy %f",
                      &c.vp.x,
                      &c.vp.y,
                      &c.vp.z,
                      &c.vi.x,
                      &c.vi.y,
                      &c.vi.z,
                      &c.vu.x,
                      &c.vu.y,
                      &c.vu.z,
                      &c.fovy);
      if (rc != 10)
        break;
      fromCL.cameraPath.push_back(c);
    }
  }
  
}

using namespace hm;

int main(int ac, char **av)
{
  // /*! init ALL gpus - let's do that right away, so gpus are
  //     initailized before mpi even gets to run */
  // hs::initAllGPUs();

  hs::mpi::init(ac, av);
#if HS_FAKE_MPI
  hs::mpi::Comm world;
#else
  hs::mpi::Comm world(MPI_COMM_WORLD);
#endif
  world.barrier();
  if (world.rank == 0) {
    std::cout << "#hv: hsviewer starting up" << std::endl; fflush(0);
  }
  world.barrier();

  bool hanari = true;
  hs::loader::DynamicDataLoader loader(world);
  for (int i=1;i<ac;i++) {
    const std::string arg = av[i];
    if (arg[0] != '-') {
      loader.addContent(arg);
      fromCL.contentArgs.push_back(arg);
      /* peek at the same URL the loader just consumed, so the UI
         knows how many iso-surfaces exist and where they started */
      try {
        hs::loader::ResourceSpecifier url(arg);
        for (float v : hs::loader::parseGpuIsoValues(url))
          fromCL.isoValues.push_back(v);
      } catch (...) { /* not a resource URL we care about */ }
    } else if (arg == "--no-bg") {
      fromCL.bgColor = ::mini::common::vec4f(0.f);
    } else if (arg == "--bg-color") {
      fromCL.bgColor.x = std::stof(av[++i]);
      fromCL.bgColor.y = std::stof(av[++i]);
      fromCL.bgColor.z = std::stof(av[++i]);
      fromCL.bgColor.w = std::stof(av[++i]);
    } else if (arg == "-dp" || arg == "--data-parallel") {
      fromCL.dpMode = DPMODE_DATA_PARALLEL;
    } else if (arg == "-sg" || arg == "--single-gpu") {
      fromCL.forceSingleGPU = true;
    } else if (arg == "-dp1" || arg == "-dpsg" || arg == "--data-parallel-single-gpu") {
      fromCL.dpMode = DPMODE_DATA_PARALLEL;
      fromCL.forceSingleGPU = true;
    } else if (arg == "-dr" || arg == "--data-replicated") {
      fromCL.dpMode = DPMODE_DATA_REPLICATED;
    } else if (arg == "-cm" || arg == "--color-map") {
      fromCL.cmID = std::stoi(av[++i]);
    } else if (arg == "-env" || arg == "--env-map") {
      fromCL.envMapFileName = av[++i];
      loader.sharedLights.envMap = fromCL.envMapFileName;
    } else if (arg == "--no-denoiser" || arg == "--no-denoise") {
      fromCL.denoise = false;
    } else if (arg == "--num-frames") {
      fromCL.numFramesAccum = std::stoi(av[++i]);
    } else if (arg == "--ambient") {
      fromCL.ambientRadiance = std::stof(av[++i]);
    } else if (arg == "-spp" || arg == "-ppp" || arg == "--paths-per-pixel") {
      fromCL.spp = std::stoi(av[++i]);
    } else if (arg == "-mum" || arg == "--merge-unstructured-meshes" || arg == "--merge-umeshes") {
      fromCL.mergeUnstructuredMeshes = true;
    } else if (arg == "--no-mum") {
      fromCL.mergeUnstructuredMeshes = false;
    } else if (arg == "--default-radius") {
      loader.defaultRadius = std::stof(av[++i]);
    } else if (arg == "--toggle-cycles") {
      fromCL.toggleCycles = std::stoi(av[++i]);
    } else if (arg == "--measure") {
      fromCL.measure = true;
    } else if (arg == "-o") {
      fromCL.outFileName = av[++i];
    } else if (arg == "--dir-light") {
      mini::DirLight light;
      light.direction = get3f(av,i);
      light.radiance  = get3f(av,i);
      loader.sharedLights.directional.push_back(light);
    } else if (arg == "--camera-pdu") {
      fromCL.camera.vp = get3f(av,i);
      fromCL.camera.vi = fromCL.camera.vp + get3f(av,i);
      fromCL.camera.vu = get3f(av,i);
    } else if (arg == "--camera") {
      fromCL.camera.vp = get3f(av,i);
      fromCL.camera.vi = get3f(av,i);
      fromCL.camera.vu = get3f(av,i);
    } else if (arg == "--cameras-from-file") {
      hm::addCamerasFromFile(av[++i]);
    } else if (arg == "--camera-path") {
      CmdLineCamera c0, c1;
      int numSteps = std::stoi(av[++i]);
      // --camera
      ++i;
      c0.vp = get3f(av,i);
      c0.vp = get3f(av,i);
      c0.vu = get3f(av,i);
      ++i;// -fovy
      c0.fovy = std::stof(av[++i]);
      ++i;// --camera
      c1.vp = get3f(av,i);
      c1.vp = get3f(av,i);
      c1.vu = get3f(av,i);
      ++i;// -fovy
      c1.fovy = std::stof(av[++i]);
      hm::addCameraPath(numSteps,c0,c1);
    } else if (arg == "-fovy") {
      fromCL.camera.fovy = std::stof(av[++i]);
    } else if (arg == "-xf") {
      fromCL.xfFileName = av[++i];
    } else if (arg == "-res" || arg == "-os" || arg == "--output-size") {
      fromCL.fbSize.x = std::stoi(av[++i]);
      fromCL.fbSize.y = std::stoi(av[++i]);
    } else if (arg == "-ndg") {
      fromCL.ndg = std::stoi(av[++i]);
      fromCL.dpMode
        = (fromCL.ndg == 1)
        ? DPMODE_DATA_REPLICATED
        : DPMODE_DATA_PARALLEL;
    } else if (arg == "-dpr") {
      fromCL.dpr = std::stoi(av[++i]);
    } else if (arg == "-nhn" || arg == "--no-head-node") {
      fromCL.createHeadNode = false;
    } else if (arg == "-hn" || arg == "-chn" ||
               arg == "--head-node" || arg == "--create-head-node") {
      fromCL.createHeadNode = true;
    } else if (arg == "-h" || arg == "--help") {
      usage();
    } else if (arg == "-anari" || arg == "--hanari") {
      hanari = true;
    } else if (arg == "-native" || arg == "--native") {
      hanari = false;
    } else {
      usage("unknown cmd-line argument '"+arg+"'");
    }    
  }

  int localRank = 0, localSize = 1;
#if HS_MPI
  determineLocalProcessID(world,localRank,localSize);
#endif
  std::vector<int> gpuIDs;
  for (int i=0;i<world.size;i++) {
    world.barrier();
    if (i != world.rank) continue;
    gpuIDs = selectGPUs(world,localRank,localSize);
  }
  assert(!gpuIDs.empty());
  world.barrier();
  
  const bool isHeadNode = fromCL.createHeadNode && (world.rank == 0);
  hs::mpi::Comm workers = world.split(!isHeadNode);

  if (world.size > 1 && fromCL.dpMode == DPMODE_NOT_SPECIFIED)
    throw std::runtime_error("you're running haystack in MPI mode, and with more than one rank, but didn't specify num data groups (-ndg <n>), or whether you want to run data parallel (-dp|--data-parallel) or data replicated (-dr|--data-replicated). Just to make sure we're not actually running the wrong mode I'll hereby bail out ...");
  
  if (fromCL.ndg == 0) {
    if (fromCL.dpMode == DPMODE_DATA_PARALLEL && world.size > 1)
      // we _are_ run in mpi mode with more than one rank, and in data
      // _parallel_mode. if not otherwise specified, use one data
      // group per rank))
      fromCL.ndg = world.size;
    else
      fromCL.ndg = 1;
  }
  
  int numDataGroupsGlobally = fromCL.ndg;
  int dataPerRank   = fromCL.dpr;
  //  LocalModel thisRankData;
  LocalPartitions *localPartitions = 0;
  // localPartitions->colorMapIndex = fromCL.cmID;
  if (!isHeadNode) {
    localPartitions = loader.loadData(numDataGroupsGlobally,dataPerRank);
    // loader.loadData(thisRankData,numDataGroupsGlobally,dataPerRank,verbose());
  }
  if (fromCL.mergeUnstructuredMeshes) {
    std::cout << "merging potentially separate unstructured meshes into single mesh" << std::endl;
    localPartitions->mergeUnstructuredMeshes();
    std::cout << "done mergine umeshes..." << std::endl;
  }
  
  int numPartitionsLocally = localPartitions->numPartitionsOnThisRank();
  if (numPartitionsLocally == 0)
    throw std::runtime_error("no partitions on this rank!?");
  world.barrier();

  // std::vector<int> gpuIDs;
  std::vector<DeviceConfig> deviceConfigs;
  for (int i=0;i<gpuIDs.size();i++) {
    DeviceConfig dc;
    dc.gpuID = gpuIDs[i];
    dc.localPartitionIndex = i % numPartitionsLocally;
    deviceConfigs.push_back(dc);
  }

  GlobalRenderSettings globalRenderSettings;
  globalRenderSettings.samplesPerPixel = fromCL.spp;
  globalRenderSettings.ambientRadiance = fromCL.ambientRadiance;
  globalRenderSettings.bgColor = fromCL.bgColor;
  globalRenderSettings.denoise = fromCL.denoise;
  globalRenderSettings.defaultColorMapIndex = fromCL.cmID;
  
  HayMaker *hayMaker
    = new HayMaker(// mpi/peers:
                   world,workers,
                   globalRenderSettings,
                   // the parition(s) we have loaded locally on this rank
                   localPartitions,
                   
                   deviceConfigs);
                   
    // = HayMaker::createAnariImplementation(world,
    //                                       /* the workers */workers,
    //                                       fromCL.spp,
    //                                        fromCL.ambientRadiance,
    //                                       fromCL.bgColor,
    //                                       thisRankData,
    //                                       gpuIDs,verbose());
  
  world.barrier();
  const BoundsData worldBounds = hayMaker->getWorldBounds();
  bool modelHasVolumeData = !worldBounds.scalars.empty();
  
  if (world.rank == 0)
    std::cout << MINI_TERMINAL_CYAN
              << "#hs: world bounds is " << worldBounds
              << MINI_TERMINAL_DEFAULT << std::endl;

  if (fromCL.camera.vp == fromCL.camera.vi) {
    fromCL.camera.vp
      = worldBounds.spatial.center()
      + mini::common::vec3f(-.3f, .7f, +1.f) * worldBounds.spatial.span();
    fromCL.camera.vi = worldBounds.spatial.center();
  }

  world.barrier();
  if (world.rank == 0)
    std::cout << MINI_TERMINAL_CYAN
              << "#hs: creating context"
              << MINI_TERMINAL_DEFAULT << std::endl;
  world.barrier();
  if (world.rank == 0)
    std::cout << MINI_TERMINAL_CYAN
              << "#hs: building data groups"
              << MINI_TERMINAL_DEFAULT << std::endl;
  if (!isHeadNode)
    hayMaker->renderInitialAnariWorld();
  
  world.barrier();

  RenderEngineInterface *renderer = nullptr;
  if (world.size == 1)
    // no MPI, render direcftly
    renderer = hayMaker;
  else if (world.rank == 0)
    // we're in MPI mode, _and_ the rank that runs the viewer
    renderer = new MPIRenderEngine(world,hayMaker);
  else {
    // we're in MPI mode, but one of the passive workers (ie NOT running the viewer)
    MPIRenderEngine::runWorker(world,hayMaker);
    world.barrier();
    hs::mpi::finalize();

    exit(0);
  }

  /* an orbital-current grid field takes over the volume ui and the default
     transfer function, so find one if the scene has one; the head node has no
     partitions of its own, and in that case both simply stay off */
  hs::OrbGridVolume::SP orbGrid;
  if (localPartitions)
    for (auto part : localPartitions->myPartitions)
      if (part && !part->orbGrids.empty()) {
        orbGrid = part->orbGrids[0];
        break;
      }

#if HS_IMGUI
  Viewer viewer(worldBounds,renderer,&world,orbGrid);

  // viewer.show();
  viewer.enableFlyMode();
  viewer.enableInspectMode();
  viewer.setWorldScale(length(worldBounds.spatial.span()));
  viewer.setCameraOrientation
    (/*origin   */(const owl::common::vec3f &)fromCL.camera.vp,
     /*lookat   */(const owl::common::vec3f &)fromCL.camera.vi,
     /*up-vector*/(const owl::common::vec3f &)fromCL.camera.vu,
     /*fovy(deg)*/fromCL.camera.fovy);

  // QMainWindow secondWindow;
  // if (modelHasVolumeData) {
  //   XFEditor *xfEditor = new XFEditor
  //     ((cutee::common::interval<float>&)worldBounds.scalars);
  //   viewer.xfEditor = xfEditor;

  //   QWidget *sidePanel = new QWidget;
  //   auto *sideLayout = new QVBoxLayout(sidePanel);
  //   sideLayout->addWidget(xfEditor);
  
    // QObject::connect(xfEditor,&cutee::XFEditor::colorMapChanged,
    //                  &viewer, &Viewer::colorMapChanged);
    // QObject::connect(xfEditor,&cutee::XFEditor::rangeChanged,
    //                  &viewer, &Viewer::rangeChanged);
    // QObject::connect(xfEditor,&cutee::XFEditor::opacityScaleChanged,
    //                  &viewer, &Viewer::opacityScaleChanged);
  
  if (!fromCL.xfFileName.empty())
    viewer.loadXF();//fromCL.xfFileName);
    // viewer.xf.load(fromCL.xfFileName);
      // viewer.xfEditor.setColorAndAlpha(loadCM(fromCL.xfFileName));

// #if HS_USE_MULTI_SCATTERING
//     {
//       VolumeScatterPanel *scatterPanel = new VolumeScatterPanel;
//       viewer.scatterPanel = scatterPanel;
//       scatterPanel->setSettings(VolumeScatterSettings{});
//       viewer.scatterSettings = scatterPanel->getSettings();
//       sideLayout->addWidget(scatterPanel);
//       QObject::connect(scatterPanel, &VolumeScatterPanel::settingsChanged,
//                        &viewer, &Viewer::scatterSettingsChanged);
//     }
// #endif

  viewer.showAndRun();
#else

  auto &fbSize = fromCL.fbSize;
  std::vector<uint32_t> pixels(fbSize.x*fbSize.y);
  renderer->resize((const mini::common::vec2i&)fbSize,pixels.data());

  hs::Camera camera;
  camera.vp = fromCL.camera.vp;
  camera.vu = fromCL.camera.vu;
  camera.vi = fromCL.camera.vi;
  camera.fovy = fromCL.camera.fovy;
  renderer->setCamera(camera);
  
  if (fromCL.xfFileName.length() > 0) {
    hs::TransferFunction xf;
    xf.load(fromCL.xfFileName);
    renderer->setTransferFunction(xf);
    renderer->resetAccumulation();
  } else if (orbGrid) {
    /* the stock default is opaque white over the whole domain, which on a
       field that is nonzero nearly everywhere is just a white wall */
    renderer->setTransferFunction
      (hs::makeSignedTransferFunction(orbGrid->valueRange,
                                      orbGrid->bounds,
                                      orbGrid->alphaGamma,
                                      orbGrid->opacityScale,
                                      orbGrid->domainScale));
    renderer->resetAccumulation();
  }

  if (!fromCL.cameraPath.empty()) {
    std::cout << "rendering camera path sequence" << std::endl;
    for (int frameID=0;frameID<fromCL.cameraPath.size();frameID++) {
      hs::Camera camera;
      auto c = fromCL.cameraPath[frameID];
      camera.vi = c.vi;
      camera.vp = c.vp;
      camera.vu = c.vu;
      camera.fovy = c.fovy;
      // float f = frameID / (fromCL.cameraPath.numSteps-1.f);
      // camera.vu = fromCL.camera.vu;
      // camera.vp
      //   = (1.f-f)*fromCL.cameraPath.vp0 + f*fromCL.cameraPath.vp1;
      // camera.vi
      //   = (1.f-f)*fromCL.cameraPath.vi0 + f*fromCL.cameraPath.vi1;
      renderer->setCamera(camera);
      renderer->renderFrame();
      stbi_flip_vertically_on_write(true);
      char suffix[100];
      sprintf(suffix,"_frame%05i.png",frameID);
      std::string fileName = fromCL.outFileName+std::string(suffix);
      std::cout << " ... saving frame " << fileName << std::endl;
      stbi_write_png(fileName.c_str(),fbSize.x,fbSize.y,4,
                     pixels.data(),fbSize.x*sizeof(uint32_t));
    }
    renderer->terminate();
    world.barrier();
    hs::mpi::finalize();
    exit(0);
  }
  
  /* exercise the filter path before rendering: every channel off, then every
     channel back on. The mask ends identical to where it started, so anything
     that differs afterwards is state that a filter change corrupted. */
  if (fromCL.toggleCycles != 0 && orbGrid) {
    hs::ChannelFilter f;
    for (int c=0;c<abs(fromCL.toggleCycles);c++) {
      /* drop ONE channel, the way a user ticking a box does. Clearing every
         channel instead empties every macro-cell, which is a different (and
         more forgiving) path than the one that actually gets exercised.
         A negative count keeps the mask UNCHANGED through the cycle: same
         number of commits and frames, no filter change, so it isolates what
         the extra frames alone do from what a filter change does. */
      f.channelMask = (fromCL.toggleCycles < 0)
        ? orbGrid->channelMask
        : (orbGrid->channelMask & ~1u);
      f.channelDirty = true;
      renderer->setChannelFilter(f);
      /* a frame is what actually drives the commit down into barney */
      renderer->renderFrame();
      f.channelMask = orbGrid->channelMask;
      f.channelDirty = true;
      renderer->setChannelFilter(f);
      renderer->renderFrame();
    }
    /* the frames above were rendered with the filter half-way through the
       cycle; without this they would accumulate into the one being compared
       and the test would fail on its own book-keeping rather than on the
       renderer's */
    renderer->resetAccumulation();
    std::cout << "#hs: completed " << fromCL.toggleCycles
              << " filter toggle cycles" << std::endl;
  }

#if 1
  if (fromCL.measure) {
    int numFramesRendered = 0;
    const int measure_warmup_frames = 2;
    const int measure_max_frames = 100;
    const float measure_max_seconds = 60.f;
    double measure_t0 = 0.;

    while (true) {
      if (numFramesRendered == measure_warmup_frames)
        measure_t0 = getCurrentTime();
      
      double t0 = getCurrentTime();
      renderer->renderFrame();
      ++numFramesRendered;
      double t1 = getCurrentTime();
      
      int numFramesMeasured = numFramesRendered - measure_warmup_frames;
      float numSecondsMeasured
        = (numFramesMeasured < 1)
        ? 0.f
        : float(t1 - measure_t0);
      if (numFramesMeasured >= measure_max_frames ||
          numSecondsMeasured >= measure_max_seconds) {
        std::cout << "measure: rendered " << numFramesMeasured << " frames in " << numSecondsMeasured << ", that is:" << std::endl;
        std::cout << "FPS " << double(numFramesMeasured/numSecondsMeasured) << std::endl;
        stbi_flip_vertically_on_write(true);
        std::cout << "saving in " << fromCL.outFileName.c_str() << std::endl;
        stbi_write_png(fromCL.outFileName.c_str(),fbSize.x,fbSize.y,4,
                       pixels.data(),fbSize.x*sizeof(uint32_t));
        
        renderer->terminate();
        world.barrier();
        hs::mpi::finalize();
        exit(0);
      }
    }
  }
#endif
  for (int i=0;i<fromCL.numFramesAccum;i++) 
    renderer->renderFrame();

  stbi_flip_vertically_on_write(true);
  std::cout << "saving in " << fromCL.outFileName.c_str() << std::endl;
  stbi_write_png(fromCL.outFileName.c_str(),fbSize.x,fbSize.y,4,
                 pixels.data(),fbSize.x*sizeof(uint32_t));

  renderer->terminate();
#endif

  world.barrier();
  hs::mpi::finalize();

  return 0;
}
