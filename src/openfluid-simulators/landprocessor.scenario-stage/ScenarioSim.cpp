
#include <openfluid/ware/PluggableSimulator.hpp>

#include <LandProcessor/LandProcessor.hpp>


// =====================================================================
// =====================================================================


BEGIN_SIMULATOR_SIGNATURE("landprocessor.scenario-stage")

  // Informations
  DECLARE_NAME("Scenario stage of the LandProcessor process")
  DECLARE_DESCRIPTION("")
  DECLARE_VERSION("")
  DECLARE_STATUS(openfluid::ware::BETA)


END_SIMULATOR_SIGNATURE


// =====================================================================
// =====================================================================


/**

*/
class ScenarioSimulator : public openfluid::ware::PluggableSimulator
{
  private:

    std::string GISdataInputDir = "gisdata-input";
    std::string GISdataOutputDir = "gisdata-output";
    std::string GISdataReleaseDir = "gisdata-release";


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

      LandProcessor LP(InputDir+"/"+GISdataInputDir,
                       OutputDir+"/"+GISdataOutputDir,
                       OutputDir+"/"+GISdataReleaseDir);


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
      catch (std::exception& E)
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
