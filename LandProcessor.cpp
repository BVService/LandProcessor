/*
 * LandProcessor.cpp
 *
 *  Created on: 8 sept. 2016
 *      Author: fabrejc
 */


#include <string>
#include <iostream>

#include <openfluid/tools/Filesystem.hpp>
#include <openfluid/utils/GrassGISProxy.hpp>
#include <openfluid/base/Environment.hpp>


#include "LandProcessor.hpp"


LandProcessor::LandProcessor(const std::string& InputPath, const std::string& OutputPath):
  m_InputPath(InputPath), m_OutputPath(OutputPath)
{
  // set input folder (check if it exist, r/w perm, etc.)
  // create output folders (/output/vector and /output/raster)
  // set variables (root, snapDistance, file names, etc.)


  openfluid::tools::Filesystem::makeDirectory(OutputPath);
  openfluid::tools::Filesystem::makeDirectory(getOutputVectorPath());
  openfluid::tools::Filesystem::makeDirectory(getOutputRasterPath());

  m_IsReady = m_IsReady && openfluid::tools::Filesystem::isDirectory(getInputVectorPath());
  m_IsReady = m_IsReady && openfluid::tools::Filesystem::isDirectory(getInputRasterPath());

  openfluid::tools::Filesystem::makeDirectory(openfluid::base::Environment::getTempDir());
}


// =====================================================================
// =====================================================================


LandProcessor::~LandProcessor()
{

}


// =====================================================================
// =====================================================================


std::string LandProcessor::getOutputVectorPath(const std::string& Filename) const
{
  std::string Path = m_OutputPath+"/"+m_VectorDir;

  if (!Filename.empty())
    Path += "/"+Filename;

  return Path;
}


// =====================================================================
// =====================================================================


std::string LandProcessor::getInputVectorPath(const std::string& Filename) const
{
  std::string Path = m_InputPath+"/"+m_VectorDir;

  if (!Filename.empty())
    Path += "/"+Filename;

  return Path;

}


// =====================================================================
// =====================================================================


std::string LandProcessor::getOutputRasterPath(const std::string& Filename) const
{
  std::string Path = m_OutputPath+"/"+m_RasterDir;

  if (!Filename.empty())
    Path += "/"+Filename;

  return Path;
}


// =====================================================================
// =====================================================================


std::string LandProcessor::getInputRasterPath(const std::string& Filename) const
{
  std::string Path = m_InputPath+"/"+m_RasterDir;

  if (!Filename.empty())
    Path += "/"+Filename;

  return Path;

}


// =====================================================================
// =====================================================================

void LandProcessor::preprocessVectorData()
{
  // import vector data
  // check data conformity and integrity
  // check geometry validity (polygones and linestrings)
  // check for multipart, overlays, empty geometries, etc.
  // check SRS
  // correction (snap, union, MtS if necessary)
  // export vector data
}

// =====================================================================
// =====================================================================

void LandProcessor::preprocessRasterData()
{
  // import raster data
  // check data conformity and integrity
  // check SRS
  // correction (sinks and flats)
  // vector-to-raster
  // raster creation (create-calculate-update-export GRASS and GDAL: (GRASS: dem --> demclean --> drainage) (GDAL: rasterID, plotrast, outlets, receivers, downslope, catchments) )
  // export raster data

  openfluid::tools::Filesystem::makeDirectory("/tmp/bvservice-grass");
  openfluid::utils::GrassGISProxy GRASS("/tmp/bvservice-grass","temp");

  GRASS.setOutputFile("/tmp/bvservice-grass/processrasterdata.out");
  GRASS.setErrorFile("/tmp/bvservice-grass/processrasterdata.err");
  GRASS.appendTask("v.in.ogr",{{"input",QString::fromStdString(getInputVectorPath(m_InputPlotsFile))},
                               {"output","plots"}},{"--o"});
  GRASS.runJob();
}

// =====================================================================
// =====================================================================

void LandProcessor::createESF()
{
  // import results from LandProcessor::preprocessRasterData
  // raster-to-vector
  // re-labelling (mult)
  // re-groupping (mult)
  // export
}

// =====================================================================
// =====================================================================

void LandProcessor::createELN()
{
  // import results from LandProcessor::preprocessRasterData
  // features-to-boundaries
  // re-labelling
  // export
}

// =====================================================================
// =====================================================================

void LandProcessor::setESFParameters()
{
  // import
  // set attributes from raster
  // set attributes from vector
  // export
}

// =====================================================================
// =====================================================================

void LandProcessor::setELNParameters()
{
  // import (SL, ELN)
  // different operation (GDAL-based)
  // set attributes
  // export
}
