/*
 * LandProcessor.hpp
 *
 *  Created on: 8 sept. 2016
 *      Author: fabrejc
 */

#ifndef __LANDPROCESSOR_HPP__
#define __LANDPROCESSOR_HPP__


#include <string>


class LandProcessor
{
  private:

    std::string m_InputPath;

    std::string m_OutputPath;

    const std::string m_VectorDir = "vector";

    const std::string m_RasterDir = "raster";

    const std::string m_InputPlotsFile = "plots.shp";

    const std::string m_InputDitchesFile = "fosses.shp";

    const std::string m_InputHedgesFile = "haies.shp";

    const std::string m_InputGrassBandFile = "bandesenherbees.shp";

    const std::string m_InputRivesrFile = "coursdeau.shp";

    const std::string m_InputThalwegsFile = "talweg.shp";

    const std::string m_InputBenchesFile = "talus.shp";

    const std::string m_InputDEMFile = "DEM.tif";

    bool m_IsReady = true;

    unsigned int minimEntSize = 250; // minimal entity size in square meters

    double snapDistance = 1.0e-08; // snap distance in meters


  public:

	LandProcessor(const std::string& InputPath, const std::string& OutputPath);

	virtual ~LandProcessor();

  bool isReady() const
  { return m_IsReady; }

  std::string getOutputVectorPath(const std::string& Filename = "") const;

  std::string getInputVectorPath(const std::string& Filename = "") const;

  std::string getOutputRasterPath(const std::string& Filename = "") const;

  std::string getInputRasterPath(const std::string& Filename = "") const;

	void preprocessVectorData();

	void preprocessRasterData();

	void createESF();

	void createELN();

	void setESFParameters();

	void setELNParameters();

};


#endif /* __LANDPROCESSOR_HPP__ */
