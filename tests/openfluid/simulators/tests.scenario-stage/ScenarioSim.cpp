
#include <openfluid/ware/PluggableSimulator.hpp>

#include "LandProcessor.hpp"


// =====================================================================
// =====================================================================


BEGIN_SIMULATOR_SIGNATURE("tests.scenario-stage")

  // Informations
  DECLARE_NAME("")
  DECLARE_DESCRIPTION("")
  DECLARE_VERSION("")
  DECLARE_STATUS(openfluid::ware::EXPERIMENTAL)


END_SIMULATOR_SIGNATURE


// =====================================================================
// =====================================================================


/**

*/
class ScenarioSimulator : public openfluid::ware::PluggableSimulator
{

  public:

  
    ScenarioSimulator(): PluggableSimulator()
    {
  
  
    }
  
  
    // =====================================================================
    // =====================================================================
  
  
    ~ScenarioSimulator()
    {
  
  
    }
  

    // =====================================================================
    // =====================================================================

  
    void initParams(const openfluid::ware::WareParams_t& Params)
    {

    }


    // =====================================================================
    // =====================================================================
  
  
    void prepareData()
    {
      std::string InputDir, OutputDir;

      OPENFLUID_GetRunEnvironment("dir.input",InputDir);
      OPENFLUID_GetRunEnvironment("dir.output",OutputDir);

      LandProcessor LP(InputDir+"/gisdata-input",
                       OutputDir+"/gisdata-output",
                       OutputDir+"/gisdata-release");

      if (!LP.isReady())
      {
        OPENFLUID_RaiseError("LandProcessor is not ready");
      }

      try
      {
        LP.extractPlotsLimits();
        LP.attributeLinearStructures();
        LP.createSU();
        LP.createRS();
        LP.createLI();
        LP.setSUParameters();
        LP.setRSParameters();
        LP.setLIParameters();
        LP.releaseFiles();
      }
      catch (std::exception &E)
      {
        OPENFLUID_RaiseError(E.what());
      }
  
    }
  
  
    // =====================================================================
    // =====================================================================
  
  
    void checkConsistency()
    {
  
  
    }
  
  
    // =====================================================================
    // =====================================================================
  
  
    openfluid::base::SchedulingRequest initializeRun()
    {  
      
      return Never();
    }


    // =====================================================================
    // =====================================================================
  
  
    openfluid::base::SchedulingRequest runStep()
    {

      return Never();
    }


    // =====================================================================
    // =====================================================================
  
  
    void finalizeRun()
    {
  
  
    }

};


// =====================================================================
// =====================================================================


DEFINE_SIMULATOR_CLASS(ScenarioSimulator);


DEFINE_WARE_LINKUID(WARE_LINKUID)


