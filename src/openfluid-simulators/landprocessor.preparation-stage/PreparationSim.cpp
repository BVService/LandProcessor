
#include <openfluid/ware/PluggableSimulator.hpp>

#include <LandProcessor/LandProcessor.hpp>


// =====================================================================
// =====================================================================


BEGIN_SIMULATOR_SIGNATURE("landprocessor.preparation-stage")

  // Informations
  DECLARE_NAME("Preparation stage of the LandProcessor process")
  DECLARE_DESCRIPTION("")
  DECLARE_VERSION("")
  DECLARE_STATUS(openfluid::ware::BETA)

  DECLARE_USED_PARAMETER("LandUseFieldName","","-");

END_SIMULATOR_SIGNATURE


// =====================================================================
// =====================================================================


/**

*/
class PreparationSimulator : public openfluid::ware::PluggableSimulator
{

  private:

	  std::string m_LandUseFieldName;
    std::string GISdataInputDir = "gisdata-input";
    std::string GISdataOutputDir = "gisdata-output";
    std::string GISdataReleaseDir = "gisdata-release";


  public:


    PreparationSimulator(): PluggableSimulator(), m_LandUseFieldName("LandUse")
    {


    }


    // =====================================================================
    // =====================================================================


    ~PreparationSimulator()
    {


    }


    // =====================================================================
    // =====================================================================


    void initParams(const openfluid::ware::WareParams_t& Params)
    {
       std::string TmpStr;

       OPENFLUID_GetSimulatorParameter(Params, "LandUseFieldName", TmpStr);

       // use parameter for land use field name only if not en empty string
       if (!TmpStr.empty())
         m_LandUseFieldName = TmpStr;
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
        LP.preprocessVectorData();
        LP.preprocessRasterData();
        LP.createSRFandLNR();
        LP.setSRFParameters();
        LP.setLNRParameters();
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


DEFINE_SIMULATOR_CLASS(PreparationSimulator);


DEFINE_WARE_LINKUID(WARE_LINKUID)
