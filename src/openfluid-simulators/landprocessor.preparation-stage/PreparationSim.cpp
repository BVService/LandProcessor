
#include <openfluid/ware/PluggableSimulator.hpp>

#include <LandProcessor/LandProcessor.hpp>


// =====================================================================
// =====================================================================


BEGIN_SIMULATOR_SIGNATURE("landprocessor.preparation-stage")

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
class PreparationSimulator : public openfluid::ware::PluggableSimulator
{

  public:


    PreparationSimulator(): PluggableSimulator()
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
        LP.preprocessVectorData();
        LP.preprocessRasterData();
        LP.createSRFandLNR();
        LP.setSRFParameters();
        LP.setLNRParameters();
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


DEFINE_SIMULATOR_CLASS(PreparationSimulator);


DEFINE_WARE_LINKUID(WARE_LINKUID)
