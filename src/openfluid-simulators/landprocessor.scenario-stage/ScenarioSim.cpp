
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

  DECLARE_USED_PARAMETER("LandUseFieldName","","-");


END_SIMULATOR_SIGNATURE


// =====================================================================
// =====================================================================


/**

*/
class ScenarioSimulator : public openfluid::ware::PluggableSimulator
{
  private:

	std::string m_LandUseFieldName;
    std::string GISdataInputDir = "gisdata-input";
    std::string GISdataOutputDir = "gisdata-output";
    std::string GISdataReleaseDir = "gisdata-release";


  public:


    ScenarioSimulator(): PluggableSimulator(), m_LandUseFieldName("LandUse")
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
    	OPENFLUID_GetSimulatorParameter(Params, "LandUseFieldName", m_LandUseFieldName);
    }


    // =====================================================================
    // =====================================================================


    void prepareData()
    {
      std::string InputDir, OutputDir;

      OPENFLUID_GetRunEnvironment("dir.input",InputDir);
      OPENFLUID_GetRunEnvironment("dir.output",OutputDir);

      try
      {
        LandProcessor LP(InputDir+"/"+GISdataInputDir,
                         OutputDir+"/"+GISdataOutputDir,
                         OutputDir+"/"+GISdataReleaseDir);

        LP.setLandUseFieldName(m_LandUseFieldName);
        LP.createSUVector();
        LP.createRSVector();
        LP.createLIVector();
        LP.setSUAttributes();
        LP.setRSAttributes();
        LP.setLIAttributes();
        LP.releaseSURSLIVectors();
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
