/*
 * LandProcessor.cpp
 *
 *  Created on: 8 sept. 2016
 *      Author: ezadonina
 */


#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <list>
#include <map>
#include <cstring>
#include <stdexcept>
#include <stdio.h>
#include <gdal/ogr_core.h>
#include <gdal/gdal_alg.h>
#include <gdal/gdal_priv.h>

#include <openfluid/tools/Filesystem.hpp>
#include <openfluid/utils/GrassGISProxy.hpp>
#include <openfluid/base/Environment.hpp>

#include "LandProcessor.hpp"
#include "Helpers.hpp"


LandProcessor::LandProcessor(const std::string& InputPath, const std::string& OutputPath):
m_InputPath(InputPath), m_OutputPath(OutputPath)
{
  openfluid::tools::Filesystem::makeDirectory(OutputPath);
  openfluid::tools::Filesystem::makeDirectory(getOutputVectorPath());
  openfluid::tools::Filesystem::makeDirectory(getOutputRasterPath());

  m_IsReady = m_IsReady && openfluid::tools::Filesystem::isDirectory(getInputVectorPath());
  m_IsReady = m_IsReady && openfluid::tools::Filesystem::isDirectory(getInputRasterPath());

  openfluid::tools::Filesystem::makeDirectory(openfluid::base::Environment::getTempDir());

  openfluid::tools::Filesystem::makeDirectory("/tmp/bvservice-grass");
  openfluid::utils::GrassGISProxy GRASS("/tmp/bvservice-grass","temp");

  GRASS.createLocation(m_EPSGCode.c_str());

  OGRRegisterAll();
  GDALAllRegister();
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

  BLOCK_ENTER(__PRETTY_FUNCTION__)

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  std::string FIDAB, FIDBA;

  int FIDA, FIDB;

  std::vector<std::string> FieldNamesTab, ClipTab, ClipedTab;
  OGRDataSource *DataSource, *Plots;
  OGRLayer *SQLLayer, *PlotsLayer;
  OGRFeature *PlotsFeature;
  OGRGeometry *PlotsGeometry, *NewGeometry;
  OGRFeatureDefn *FeatureDefn;

  GDALDataset *DEM;

  OGRSpatialReference *VectorSRS, RasterSRS;
  const char *VectorDataEPSGCode, *RasterDataEPSGCode;

  // =====================================================================
  // Data checking (format, SRS, type, etc.)
  // =====================================================================

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  if (!openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputPlotsFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsFile + ": there is no such file input vector directory");
  }

  Plots = OGRSFDriverRegistrar::Open( getInputVectorPath(m_InputPlotsFile).c_str(), TRUE );

  if (Plots->GetDriver()->GetName() != m_VectorDriverName)
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsFile + ": vector file is not in ESRI Shapefile format");
  }

  PlotsLayer = Plots->GetLayer(0);
  VectorSRS = PlotsLayer->GetSpatialRef();

  if (!VectorSRS)
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsFile + ": no spatial reference information is provided for the vector file");
  }

  VectorDataEPSGCode = VectorSRS->GetAttrValue("AUTHORITY",1);

  if (!VectorDataEPSGCode)
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsFile + ": no EPSG code provided for the vector data");
  }

  if (VectorDataEPSGCode != m_EPSGCode.substr(m_EPSGCode.find(":")+1, m_EPSGCode.length()-1))
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsFile + ":  vector data EPSG code does not correspond to the default EPSG code");
  }

  if (!openfluid::tools::Filesystem::isFile(getInputRasterPath(m_InputDEMFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputDEMFile + ": there is no such input raster file in input raster directory");
  }

  DEM = (GDALDataset *) GDALOpen( getInputRasterPath(m_InputDEMFile).c_str(), GA_ReadOnly );

  if (DEM->GetDriver()->GetDescription() != m_RasterDriverName)
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputDEMFile + ": raster file is not in GTiff format");
  }

  if (!DEM->GetProjectionRef())
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputDEMFile + ": no spatial reference information is provided for the raster file");
  }

  RasterSRS.SetFromUserInput(DEM->GetProjectionRef());
  RasterDataEPSGCode = RasterSRS.GetAttrValue("AUTHORITY",1);

  if (!RasterDataEPSGCode)
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputDEMFile + ": no EPSG code provided for the raster data");
  }

  if (RasterDataEPSGCode != m_EPSGCode.substr(m_EPSGCode.find(":")+1, m_EPSGCode.length()-1))
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputDEMFile + ": raster data EPSG code does not correspond to the default EPSG code");
  }

  m_SQLRequest = "REPACK " + m_InputPlotsFile.substr(0, m_InputPlotsFile.find("."));
  Plots->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
  PlotsLayer = Plots->GetLayer(0);
  PlotsLayer->ResetReading();

  while ((PlotsFeature = PlotsLayer->GetNextFeature()) != 0)
  {
    if (PlotsFeature->GetGeometryRef()->IsValid() == 0)
    {
      throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsFile + ": one or several features have invalid geometry");
    }
  }

  PlotsLayer = Plots->GetLayer(0);
  PlotsLayer->ResetReading();

  while ((PlotsFeature = PlotsLayer->GetNextFeature()) != 0)
  {
    if (PlotsFeature->GetGeometryRef()->getGeometryType() != 3)
    {
      throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsFile + ": one or several features have geometry that are not of type 'POLYGON'");
    }
  }

  // =====================================================================
  // Creation or calculation of the unique ID field
  // =====================================================================

  PlotsLayer = Plots->GetLayer(0);
  PlotsLayer->ResetReading();
  FeatureDefn = PlotsLayer->GetLayerDefn();

  FieldNamesTab.clear();

  for (int i = 0; i < FeatureDefn->GetFieldCount(); i++)
  {
    FieldNamesTab.push_back(PlotsLayer->GetLayerDefn()->GetFieldDefn(i)->GetNameRef());
  }

  if(std::find(FieldNamesTab.begin(), FieldNamesTab.end(), m_IDFieldName) != FieldNamesTab.end())
  {
    for (int i = 0; i < PlotsLayer->GetFeatureCount(); i++)
    {
      PlotsFeature = PlotsLayer->GetFeature(i);
      PlotsFeature->SetField(m_IDFieldName.c_str(), i+1);
      PlotsLayer->SetFeature(PlotsFeature);
    }
  }
  else
  {
    createField(PlotsLayer, m_IDFieldName.c_str(), OFTInteger);
    for (int i = 0; i < PlotsLayer->GetFeatureCount(); i++)
    {
      PlotsFeature = PlotsLayer->GetFeature(i);
      PlotsFeature->SetField(m_IDFieldName.c_str(), i+1);
      PlotsLayer->SetFeature(PlotsFeature);
    }
  }

  FieldNamesTab.clear();
  OGRDataSource::DestroyDataSource(Plots);

  // =====================================================================
  // Correcting overlaps and gaps between polygons, if at all possible
  // =====================================================================

  Plots = OGRSFDriverRegistrar::Open( getInputVectorPath(m_InputPlotsFile).c_str(), TRUE );

  PlotsLayer = Plots->GetLayer(0);
  PlotsLayer->ResetReading();

  DataSource = OGRSFDriverRegistrar::Open( getInputVectorPath().c_str(), TRUE );

  ClipTab.push_back("0");

  while(ClipTab.size() != ClipedTab.size())
  {
    ClipTab.clear();
    PlotsLayer = Plots->GetLayer(0);
    PlotsLayer->ResetReading();
    while ((PlotsFeature = PlotsLayer->GetNextFeature()) != nullptr)
    {
      PlotsGeometry = PlotsFeature->GetGeometryRef();
      DataSource = OGRSFDriverRegistrar::Open( getInputVectorPath().c_str(), TRUE );
      FIDA = PlotsFeature->GetFID();
      m_SQLRequest = "SELECT *, ROWID AS RID FROM " + m_InputPlotsFile.substr(0, m_InputPlotsFile.find(".")) + " WHERE "
          "ST_Overlaps(geometry, (SELECT geometry FROM " + m_InputPlotsFile.substr(0, m_InputPlotsFile.find(".")) + " WHERE ROWID == " + std::to_string(PlotsFeature->GetFID()) + ")) AND ROWID != " + std::to_string(PlotsFeature->GetFID()) + " AND "
          " ST_GeometryType(ST_Intersection(geometry, (SELECT geometry FROM " + m_InputPlotsFile.substr(0, m_InputPlotsFile.find(".")) + " WHERE ROWID == " + std::to_string(PlotsFeature->GetFID()) + "))) != " + m_QMark + "LINESTRING" + m_QMark + " AND "
          "ST_GeometryType(ST_Intersection(geometry, (SELECT geometry FROM " + m_InputPlotsFile.substr(0, m_InputPlotsFile.find(".")) + " WHERE ROWID == " + std::to_string(PlotsFeature->GetFID()) + "))) != " + m_QMark + "MULTILINESTRING" + m_QMark;
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        NewGeometry = nullptr;
        for (int i = 0; i < SQLLayer->GetFeatureCount(); i++)
        {
          FIDB = SQLLayer->GetFeature(i)->GetFieldAsInteger("RID");
          FIDBA = std::to_string(FIDB) + "/" + std::to_string(FIDA);
          if (std::find(ClipTab.begin(), ClipTab.end(), FIDBA) == ClipTab.end())
          {
            if (NewGeometry == nullptr)
            {
              if (PlotsGeometry->Difference(SQLLayer->GetFeature(i)->GetGeometryRef())->getGeometryType() == 3)
              {
                NewGeometry = PlotsGeometry->Difference(SQLLayer->GetFeature(i)->GetGeometryRef());
              }
            }
            else
            {
              if (NewGeometry->Difference(SQLLayer->GetFeature(i)->GetGeometryRef())->getGeometryType() == 3)
              {
                NewGeometry = NewGeometry->Difference(SQLLayer->GetFeature(i)->GetGeometryRef());
              }
            }
            FIDAB = std::to_string(FIDA) + "/" + std::to_string(FIDB);
            ClipTab.push_back(FIDAB);
          }
        }
        if (NewGeometry != nullptr)
        {
          PlotsFeature->SetGeometry(NewGeometry);
          PlotsLayer->SetFeature(PlotsFeature);
        }
      }
      DataSource->ReleaseResultSet(SQLLayer);
      OGRDataSource::DestroyDataSource(DataSource);
    }
    PlotsLayer = Plots->GetLayer(0);
    PlotsLayer->ResetReading();
    ClipedTab.clear();
    while ((PlotsFeature = PlotsLayer->GetNextFeature()) != nullptr)
    {
      PlotsGeometry = PlotsFeature->GetGeometryRef();
      DataSource = OGRSFDriverRegistrar::Open( getInputVectorPath().c_str(), TRUE );
      FIDA = PlotsFeature->GetFID();
      m_SQLRequest = "SELECT *, ROWID AS RID FROM " + m_InputPlotsFile.substr(0, m_InputPlotsFile.find(".")) + " WHERE "
          "ST_Overlaps(geometry, (SELECT geometry FROM " + m_InputPlotsFile.substr(0, m_InputPlotsFile.find(".")) + " WHERE ROWID == " + std::to_string(PlotsFeature->GetFID()) + ")) AND ROWID != " + std::to_string(PlotsFeature->GetFID()) + " AND "
          " ST_GeometryType(ST_Intersection(geometry, (SELECT geometry FROM " + m_InputPlotsFile.substr(0, m_InputPlotsFile.find(".")) + " WHERE ROWID == " + std::to_string(PlotsFeature->GetFID()) + "))) != " + m_QMark + "LINESTRING" + m_QMark + " AND "
          "ST_GeometryType(ST_Intersection(geometry, (SELECT geometry FROM " + m_InputPlotsFile.substr(0, m_InputPlotsFile.find(".")) + " WHERE ROWID == " + std::to_string(PlotsFeature->GetFID()) + "))) != " + m_QMark + "MULTILINESTRING" + m_QMark;
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        for (int i = 0; i < SQLLayer->GetFeatureCount(); i++)
        {
          FIDB = SQLLayer->GetFeature(i)->GetFieldAsInteger("RID");
          FIDBA = std::to_string(FIDB) + "/" + std::to_string(FIDA);
          if (std::find(ClipedTab.begin(), ClipedTab.end(), FIDBA) == ClipedTab.end())
          {
            FIDAB = std::to_string(FIDA) + "/" + std::to_string(FIDB);
            ClipedTab.push_back(FIDAB);
          }
        }
      }
      DataSource->ReleaseResultSet(SQLLayer);
      OGRDataSource::DestroyDataSource(DataSource);
    }
  }

  ClipTab.clear();
  ClipedTab.clear();
  OGRDataSource::DestroyDataSource(Plots);

  openfluid::utils::GrassGISProxy GRASS("/tmp/bvservice-grass","temp");
  GRASS.setOutputFile("/tmp/bvservice-grass/procesvectordata.out");
  GRASS.setErrorFile("/tmp/bvservice-grass/processvectordata.err");

  GRASS.appendTask("v.in.ogr", {{"input", QString::fromStdString(getInputVectorPath(m_InputPlotsFile))}, {"output", "plotssnapped"}, {"snap", QString::fromStdString(m_SnapDistance)}}, {"--o"});

  GRASS.appendTask("r.in.gdal",{{"input",QString::fromStdString(getInputRasterPath(m_InputDEMFile))}, {"output","dem"}}, {"--o"});

  GRASS.appendTask("g.region", {{"raster", "dem"}});

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsVectorFile).c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputPlotsVectorFile).c_str());
  }

  GRASS.appendTask("v.out.ogr", {{"input", "plotssnapped"}, {"type", "area"}, {"output", QString::fromStdString(getOutputVectorPath(m_OutputPlotsVectorFile))}}, {"-s"});

  if (GRASS.runJob() != 0)
  {
    exit(-1);
  }

  // =====================================================================
  // Removing attributes that will not be used in the processing
  // =====================================================================

  Plots = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsVectorFile).c_str(), TRUE );
  PlotsLayer = Plots->GetLayer(0);
  FeatureDefn = PlotsLayer->GetLayerDefn();

  FieldNamesTab.clear();

  for (int i = 0; i < FeatureDefn->GetFieldCount(); i++)
  {
    if (FeatureDefn->GetFieldDefn(i)->GetNameRef() != m_IDFieldName)
    {
      FieldNamesTab.push_back(FeatureDefn->GetFieldDefn(i)->GetNameRef());
    }
  }

  if (FieldNamesTab.size() > 0)
  {
    for (int i = 0; i < FieldNamesTab.size(); i++)
    {
      int FieldIndex = FeatureDefn->GetFieldIndex(FieldNamesTab.at(i).c_str());
      PlotsLayer->DeleteField(FieldIndex);
    }
  }

  FieldNamesTab.clear();

  OGRDataSource::DestroyDataSource(Plots);
  GDALClose( DEM );

  BLOCK_EXIT(__PRETTY_FUNCTION__)
}


// =====================================================================
// =====================================================================


void LandProcessor::preprocessRasterData()
{
  BLOCK_ENTER(__PRETTY_FUNCTION__)

    // =====================================================================
    // Variables that are used in this block
    // =====================================================================

      int *Row, *PlotsUpperRow, *PlotsMiddleRow, *PlotsLowerRow,
      *DrainageUpperRow, *DrainageMiddleRow, *DrainageLowerRow, *DrainageRow,
      *RasterIDUpperRow, *RasterIDMiddleRow, *RasterIDLowerRow, *RasterIDRow,
      *OutletsUpperRow, *OutletsMiddleRow, *OutletsLowerRow,
      *CatchmentsUpperRow, *CatchmentsMiddleRow, *CatchmentsLowerRow;

  unsigned int PlotsCount, CatchmentsCount;

  GDALDataset *DEM, *PlotsR, *RasterID, *Drainage, *Outlets, *Receivers, *Downslope, *Catchments;
  GDALRasterBand *DEMBand, *PlotsRBand, *RasterIDBand, *DrainageBand,
  *OutletsBand, *ReceiversBand, *DownslopeBand, *CatchmentsBand;

  OGRDataSource *PlotsV;
  OGRLayer *PlotsVLayer;
  OGRFeature *PlotsVFeature;
  OGRGeometry *PlotsVGeometry, *UnionGeometry;
  OGREnvelope *Envelope;

  char **Options = nullptr;

  // ======================================================================================
  // Checking if necessary vector file is present in the appropriate folder
  // ======================================================================================

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsVectorFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::preprocessRasterData(): " + m_OutputPlotsVectorFile + ": no such file in the output vector directory");
  }

  // ======================================================================================
  // Get DTM transform coefficients that will be used to set the region in GRASS GIS
  // ======================================================================================

  DEM = (GDALDataset *) GDALOpen( getInputRasterPath(m_InputDEMFile).c_str(), GA_ReadOnly );

  mp_RasterDriver = DEM->GetDriver();

  getGeoTransform(DEM);

  GDALClose( DEM );

  PlotsV = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsVectorFile).c_str(), TRUE );

  PlotsVLayer = PlotsV->GetLayer(0);

  mp_SRS = PlotsVLayer->GetSpatialRef();

  UnionGeometry = OGRGeometryFactory::createGeometry(wkbPolygon);

  PlotsVLayer->ResetReading();

  PlotsVFeature = nullptr;

  while ((PlotsVFeature = PlotsVLayer->GetNextFeature()) != nullptr)
  {
    PlotsVGeometry = PlotsVFeature->GetGeometryRef();
    UnionGeometry = UnionGeometry->Union(PlotsVGeometry);
  }

  Envelope = new OGREnvelope();

  UnionGeometry->Buffer(mp_GeoTransformVal[1])->getEnvelope(Envelope);

  // =====================================================================
  // GRASS GIS data processing
  // =====================================================================

  openfluid::utils::GrassGISProxy GRASS("/tmp/bvservice-grass","temp");

  GRASS.setOutputFile("/tmp/bvservice-grass/procesrasterdata.out");
  GRASS.setErrorFile("/tmp/bvservice-grass/processrasterdata.err");

  GRASS.appendTask("v.to.rast", {{"input", QString::fromStdString("plotssnapped")},
                                 {"output", QString::fromStdString("plotsraster")},
                                 {"type", QString::fromStdString("area")},
                                 {"use", QString::fromStdString("attr")},
                                 {"attribute_column", QString::fromStdString(m_IDFieldName)}},
                   {"--o"});

  GRASS.appendTask("r.in.gdal",{{"input",QString::fromStdString(getInputRasterPath(m_InputDEMFile))},
                                {"output",QString::fromStdString("dem")}},
                   {"--o"});

  GRASS.appendTask("g.region", {{"align", QString::fromStdString("dem")},
                                {"n", QString::fromStdString(std::to_string(Envelope->MaxY))},
                                {"s", QString::fromStdString(std::to_string(Envelope->MinY))},
                                {"e", QString::fromStdString(std::to_string(Envelope->MaxX))},
                                {"w", QString::fromStdString(std::to_string(Envelope->MinX))},
                                {"nsres", QString::fromStdString(std::to_string(mp_GeoTransformVal[1]))},
                                {"ewres", QString::fromStdString(std::to_string(mp_GeoTransformVal[1]))}});

  // =====================================================================
  // Vector-to-raster conversion
  // =====================================================================

  GRASS.appendTask("r.to.vect", {{"input", QString::fromStdString("plotsraster")},
                                 {"output", QString::fromStdString("plotsrastervector")},
                                 {"type", QString::fromStdString("area")},
                                 {"column", QString::fromStdString("IDPlot")}},
                   {"--o"});

  GRASS.appendTask("v.db.dropcolumn", {{"map", QString::fromStdString("plotsrastervector")},
                                       {"column", QString::fromStdString("label")}});

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath("plotsrastervector.shp").c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath("plotsrastervector.shp").c_str());
  }

  GRASS.appendTask("v.out.ogr", {{"input", QString::fromStdString("plotsrastervector")},
                                 {"type", QString::fromStdString("area")},
                                 {"output", QString::fromStdString(getOutputVectorPath("plotsrastervector.shp"))},
                                 {"format", QString::fromStdString("ESRI_Shapefile")}},
                   {"-s"});

  // =====================================================================
  // DTM correction, slope and drainage raster calculation
  // =====================================================================

  GRASS.appendTask("g.remove", {{"type", QString::fromStdString("raster")},
                                {"name", QString::fromStdString("drainageraster")}},
                   {"-f"});

  GRASS.appendTask("r.fill.dir", {{"input",QString::fromStdString("dem")},
                                  {"output",QString::fromStdString("demclean")},
                                  {"direction",QString::fromStdString("direction")}},
                   {"--o"});

  GRASS.appendTask("r.slope.aspect", {{"elevation", QString::fromStdString("demclean")},
                                      {"slope", QString::fromStdString("slope")}},
                   {"--o"});

  GRASS.appendTask("r.watershed", {{"elevation", QString::fromStdString("demclean")},
                                   {"drainage", QString::fromStdString("drainageraster")}},
                   {"-sb"});



  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputDEMFile).c_str()))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputDEMFile).c_str());
  }

  GRASS.appendTask("r.out.gdal", {{"input",QString::fromStdString("demclean")},
                                  {"output", QString::fromStdString(getOutputRasterPath(m_OutputDEMFile))},
                                  {"format", QString::fromStdString("GTiff")}},
                   {"--o"});

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputSlopeRasterFile).c_str()))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputSlopeRasterFile).c_str());
  }

  GRASS.appendTask("r.out.gdal", {{"input",QString::fromStdString("slope")},
                                  {"output", QString::fromStdString(getOutputRasterPath(m_OutputSlopeRasterFile))},
                                  {"format", QString::fromStdString("GTiff")}},
                   {"--o"});

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputDrainageRasterFile).c_str()))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputDrainageRasterFile).c_str());
  }

  GRASS.appendTask("r.out.gdal", {{"input",QString::fromStdString("drainageraster")},
                                  {"output", QString::fromStdString(getOutputRasterPath(m_OutputDrainageRasterFile))},
                                  {"format", QString::fromStdString("GTiff")}},
                   {"--o"});

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputPlotsRasterFile).c_str()))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputPlotsRasterFile).c_str());
  }

  GRASS.appendTask("r.out.gdal", {{"input", QString::fromStdString("plotsraster")},
                                  {"output", QString::fromStdString(getOutputRasterPath(m_OutputPlotsRasterFile))},
                                  {"format", QString::fromStdString("GTiff")},
                                  {"type", QString::fromStdString("Int16")},
                                  {"nodata", QString::fromStdString("-9999")}},
                   {"--o"});

  if (GRASS.runJob() != 0)
  {
    exit(-1);
  }

  OGRDataSource::DestroyDataSource(PlotsV);

  // =========================================================================================
  // ID raster creation of  that will be used to set outlets, receivers and catchments IDs
  // =========================================================================================

  DEM = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputDEMFile).c_str(), GA_ReadOnly );
  DEMBand = DEM->GetRasterBand(1);
  getGeoTransform(DEM);

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputIDRasterFile).c_str()))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputIDRasterFile).c_str());
  }

  RasterID = mp_RasterDriver->Create( getOutputRasterPath(m_OutputIDRasterFile).c_str(), DEM->GetRasterXSize(), DEM->GetRasterYSize(), 1, GDT_Int32, Options );
  RasterIDBand = RasterID->GetRasterBand(1);

  Row = (int*) CPLMalloc(sizeof(int)*RasterID->GetRasterXSize());

  for (int i = 0; i < RasterID->GetRasterYSize(); i++)
  {
    for (int j = 0; j < RasterID->GetRasterXSize(); j++)
    {
      if (i == 0)
      {
        Row[j] = j+1;
      }
      else
      {
        Row[j] = RasterID->GetRasterYSize()*i+(j+1);
      }
    }
    RasterIDBand->RasterIO( GF_Write, 0, i, RasterID->GetRasterXSize(), 1,
                            Row, RasterID->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
  }

  CPLFree(Row);
  RasterIDBand->FlushCache();
  RasterIDBand->SetNoDataValue(-9999);
  RasterID->SetGeoTransform(mp_GeoTransformVal);
  RasterID->SetProjection(DEM->GetProjectionRef());
  GDALClose( RasterID );

  // =====================================================================
  // Outlets raster creation
  // =====================================================================

  PlotsR = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputPlotsRasterFile).c_str(), GA_ReadOnly );
  PlotsRBand = PlotsR->GetRasterBand(1);

  Drainage = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputDrainageRasterFile).c_str(), GA_ReadOnly );
  DrainageBand = Drainage->GetRasterBand(1);

  RasterID = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputIDRasterFile).c_str(), GA_ReadOnly );
  RasterIDBand = RasterID->GetRasterBand(1);

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputOutletsRasterFile).c_str()))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputOutletsRasterFile).c_str());
  }

  Outlets = mp_RasterDriver->Create( getOutputRasterPath(m_OutputOutletsRasterFile).c_str(), DEM->GetRasterXSize(), DEM->GetRasterYSize(), 1, GDT_Int32, Options );
  OutletsBand = Outlets->GetRasterBand(1);

  PlotsUpperRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  PlotsMiddleRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  PlotsLowerRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  DrainageRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  RasterIDRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());

  Row = (int*) CPLMalloc(sizeof(int)*DEM->GetRasterXSize());

  for (int i = 0; i < Outlets->GetRasterYSize(); i++)
  {
    for (int j = 0; j < Outlets->GetRasterXSize(); j++)
    {
      if (i == 0 or j == 0 or i == Outlets->GetRasterYSize()-1 or j == DEM->GetRasterXSize()-1)
      {
        Row[j] = -9999;
      }
      else
      {
        PlotsRBand->RasterIO( GF_Read, 0, i-1, DEM->GetRasterXSize(), 1, PlotsUpperRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        PlotsRBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, PlotsMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        PlotsRBand->RasterIO( GF_Read, 0, i+1, DEM->GetRasterXSize(), 1, PlotsLowerRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        DrainageBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, DrainageRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        RasterIDBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, RasterIDRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        if (PlotsMiddleRow[j] == -9999)
        {
          Row[j] = -9999;
        }
        else
        {
          if ((DrainageRow[j] == 1 and PlotsMiddleRow[j]!= PlotsUpperRow[j+1]) or
              (DrainageRow[j] == 2 and PlotsMiddleRow[j] != PlotsUpperRow[j]) or
              (DrainageRow[j] == 3 and PlotsMiddleRow[j] != PlotsUpperRow[j-1]) or
              (DrainageRow[j] == 4 and PlotsMiddleRow[j] != PlotsMiddleRow[j-1]) or
              (DrainageRow[j] == 5 and PlotsMiddleRow[j] != PlotsLowerRow[j-1]) or
              (DrainageRow[j] == 6 and PlotsMiddleRow[j] != PlotsLowerRow[j]) or
              (DrainageRow[j] == 7 and PlotsMiddleRow[j] != PlotsLowerRow[j+1]) or
              (DrainageRow[j] == 8 and PlotsMiddleRow[j] != PlotsMiddleRow[j+1]))
          {
            Row[j] = RasterIDRow[j];
          }
          else
          {
            Row[j] = -9999;
          }
        }
      }
    }
    OutletsBand->RasterIO( GF_Write, 0, i, Outlets->GetRasterXSize(), 1, Row, Outlets->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
  }

  CPLFree(Row);
  CPLFree(PlotsUpperRow);
  CPLFree(PlotsMiddleRow);
  CPLFree(PlotsLowerRow);
  CPLFree(DrainageRow);
  CPLFree(RasterIDRow);

  OutletsBand->FlushCache();
  OutletsBand->SetNoDataValue(-9999);
  Outlets->SetGeoTransform(mp_GeoTransformVal);
  Outlets->SetProjection(DEM->GetProjectionRef());
  GDALClose( Outlets );

  // =====================================================================
  // Receivers raster creation
  // =====================================================================

  Outlets = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputOutletsRasterFile).c_str(), GA_ReadOnly );
  OutletsBand = Outlets->GetRasterBand(1);

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputReceiversRasterFile).c_str()))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputReceiversRasterFile).c_str());
  }

  Receivers = mp_RasterDriver->Create( getOutputRasterPath(m_OutputReceiversRasterFile).c_str(), DEM->GetRasterXSize(), DEM->GetRasterYSize(), 1, GDT_Int32, Options );
  ReceiversBand = Receivers->GetRasterBand(1);

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputDownslopeRasterFile).c_str()))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputDownslopeRasterFile).c_str());
  }

  Downslope = mp_RasterDriver->Create( getOutputRasterPath(m_OutputDownslopeRasterFile).c_str(), DEM->GetRasterXSize(), DEM->GetRasterYSize(), 1, GDT_Int32, Options );
  DownslopeBand = Downslope->GetRasterBand(1);

  OutletsUpperRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  OutletsMiddleRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  OutletsLowerRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  DrainageUpperRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  DrainageMiddleRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  DrainageLowerRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  RasterIDUpperRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  RasterIDMiddleRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  RasterIDLowerRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());

  Row = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());

  for (int i = 0; i<DEM->GetRasterYSize(); i++)
  {
    for (int j = 0; j<DEM->GetRasterXSize(); j++)
    {
      if (i == 0 or j == 0 or i == DEM->GetRasterYSize()-1 or j == DEM->GetRasterXSize()-1)
      {
        Row[j] = -9999;
      }
      else
      {
        OutletsBand->RasterIO( GF_Read, 0, i-1, DEM->GetRasterXSize(), 1, OutletsUpperRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        OutletsBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, OutletsMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        OutletsBand->RasterIO( GF_Read, 0, i+1, DEM->GetRasterXSize(), 1, OutletsLowerRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        DrainageBand->RasterIO( GF_Read, 0, i-1, DEM->GetRasterXSize(), 1, DrainageUpperRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        DrainageBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, DrainageMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        DrainageBand->RasterIO( GF_Read, 0, i+1, DEM->GetRasterXSize(), 1, DrainageLowerRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        RasterIDBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, RasterIDMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        if ((DrainageUpperRow[j-1] == 7 and OutletsUpperRow[j-1] != -9999) or
            (DrainageUpperRow[j] == 6 and OutletsUpperRow[j] != -9999) or
            (DrainageUpperRow[j+1] == 5 and OutletsUpperRow[j+1] != -9999) or
            (DrainageMiddleRow[j+1] == 4 and OutletsMiddleRow[j+1] != -9999) or
            (DrainageLowerRow[j+1] == 3 and OutletsLowerRow[j+1] != -9999) or
            (DrainageLowerRow[j] == 2 and OutletsLowerRow[j] != -9999) or
            (DrainageLowerRow[j-1] == 1 and OutletsLowerRow[j-1] != -9999) or
            (DrainageMiddleRow[j-1] == 8 and OutletsMiddleRow[j-1] != -9999))
        {
          Row[j] = RasterIDMiddleRow[j];
        }
        else
        {
          Row[j] = -9999;
        }
      }
    }
    ReceiversBand->RasterIO( GF_Write, 0, i, DEM->GetRasterXSize(), 1, Row, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
  }

  CPLFree(OutletsUpperRow);
  CPLFree(OutletsMiddleRow);
  CPLFree(OutletsLowerRow);
  CPLFree(DrainageUpperRow);
  CPLFree(DrainageMiddleRow);
  CPLFree(DrainageLowerRow);
  CPLFree(RasterIDMiddleRow);
  CPLFree(Row);

  ReceiversBand->FlushCache();
  ReceiversBand->SetNoDataValue(-9999);
  Receivers->SetGeoTransform(mp_GeoTransformVal);
  Receivers->SetProjection(DEM->GetProjectionRef());
  GDALClose( Receivers );

  // =============================================================================================
  // Creation of the raster that will contain ID of the cell that the cell in question drains into
  // =============================================================================================

  OutletsMiddleRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  DrainageMiddleRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  RasterIDMiddleRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());

  Row = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());

  for (int i = 0; i<DEM->GetRasterYSize(); i++)
  {
    for (int j = 0; j<DEM->GetRasterXSize(); j++)
    {
      if (i == 0 or j == 0 or i == DEM->GetRasterYSize()-1 or j == DEM->GetRasterXSize()-1)
      {
        Row[j] = -9999;
      }
      else
      {
        OutletsBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, OutletsMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        DrainageBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, DrainageMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        RasterIDBand->RasterIO( GF_Read, 0, i-1, DEM->GetRasterXSize(), 1, RasterIDUpperRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        RasterIDBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, RasterIDMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        RasterIDBand->RasterIO( GF_Read, 0, i+1, DEM->GetRasterXSize(), 1, RasterIDLowerRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        if (DrainageMiddleRow[j] == 1 and OutletsMiddleRow[j] != -9999)
        {
          Row[j] = RasterIDUpperRow[j+1];
        }
        else if (DrainageMiddleRow[j] == 2 and OutletsMiddleRow[j] != -9999)
        {
          Row[j] = RasterIDUpperRow[j];
        }
        else if (DrainageMiddleRow[j] == 3 and OutletsMiddleRow[j] != -9999)
        {
          Row[j] = RasterIDUpperRow[j-1];
        }
        else if (DrainageMiddleRow[j] == 4 and OutletsMiddleRow[j] != -9999)
        {
          Row[j] = RasterIDMiddleRow[j-1];
        }
        else if (DrainageMiddleRow[j] == 5 and OutletsMiddleRow[j] != -9999)
        {
          Row[j] = RasterIDLowerRow[j-1];
        }
        else if (DrainageMiddleRow[j] == 6 and OutletsMiddleRow[j] != -9999)
        {
          Row[j] = RasterIDLowerRow[j];
        }
        else if (DrainageMiddleRow[j] == 7 and OutletsMiddleRow[j] != -9999)
        {
          Row[j] = RasterIDLowerRow[j+1];
        }
        else if (DrainageMiddleRow[j] == 8 and OutletsMiddleRow[j] != -9999)
        {
          Row[j] = RasterIDMiddleRow[j+1];
        }
        else
        {
          Row[j] = -9999;
        }
      }
    }
    DownslopeBand->RasterIO( GF_Write, 0, i, DEM->GetRasterXSize(), 1, Row, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
  }

  CPLFree(OutletsMiddleRow);
  CPLFree(DrainageMiddleRow);
  CPLFree(RasterIDUpperRow);
  CPLFree(RasterIDMiddleRow);
  CPLFree(RasterIDLowerRow);
  CPLFree(Row);

  DownslopeBand->FlushCache();
  DownslopeBand->SetNoDataValue(-9999);
  Downslope->SetGeoTransform(mp_GeoTransformVal);
  Downslope->SetProjection(DEM->GetProjectionRef());
  GDALClose( Downslope );

  // ============================================================================================
  // Catchments raster creation
  // First, we copy the outlet points to the new raster
  // Then we fill the cells that belong to the same parcel working outwards using drainage raster
  // ============================================================================================

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputCatchmentsRasterFile).c_str()))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputCatchmentsRasterFile).c_str());
  }

  Catchments = mp_RasterDriver->Create( getOutputRasterPath(m_OutputCatchmentsRasterFile).c_str(), DEM->GetRasterXSize(), DEM->GetRasterYSize(), 1, GDT_Int32, Options );
  CatchmentsBand = Catchments->GetRasterBand(1);

  OutletsMiddleRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());

  Row = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());

  for (int i = 0; i < DEM->GetRasterYSize(); i++)
  {
    for (int j = 0; j < DEM->GetRasterXSize(); j++)
    {
      if (i == 0 or j == 0 or i == DEM->GetRasterYSize()-1 or j == DEM->GetRasterXSize()-1)
      {
        Row[j] = -9999;
      }
      else
      {
        OutletsBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, OutletsMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
        if (OutletsMiddleRow[j] != -9999)
        {
          Row[j] = OutletsMiddleRow[j];
        }
        else
        {
          Row[j] = -9999;
        }
      }
    }
    CatchmentsBand->RasterIO( GF_Write, 0, i, DEM->GetRasterXSize(), 1, Row, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
  }

  CPLFree(OutletsMiddleRow);
  CPLFree(Row);

  CatchmentsBand->FlushCache();
  CatchmentsBand->SetNoDataValue(-9999);
  Catchments->SetGeoTransform(mp_GeoTransformVal);
  Catchments->SetProjection(DEM->GetProjectionRef());
  GDALClose( Catchments );

  Catchments = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputCatchmentsRasterFile).c_str(), GA_Update );
  CatchmentsBand = Catchments->GetRasterBand(1);

  CatchmentsUpperRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  CatchmentsMiddleRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  CatchmentsLowerRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  PlotsUpperRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  PlotsMiddleRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  PlotsLowerRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  DrainageUpperRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  DrainageMiddleRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());
  DrainageLowerRow = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());

  Row = (int*) CPLMalloc(sizeof(int) *DEM->GetRasterXSize());

  PlotsCount = 0;

  for (int i = 0; i < DEM->GetRasterYSize(); i++)
  {
    PlotsRBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, PlotsMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
    for (int j = 0; j < DEM->GetRasterXSize(); j++)
    {
      if (PlotsMiddleRow[j] != -9999)
      {
        PlotsCount += 1;
      }
    }
  }

  CatchmentsCount = 0;

  for (int i = 0; i < DEM->GetRasterYSize(); i++)
  {
    CatchmentsBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, CatchmentsMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
    for (int j = 0; j < DEM->GetRasterXSize(); j++)
    {
      if (CatchmentsMiddleRow[j] != -9999)
      {
        CatchmentsCount += 1;
      }
    }
  }

  while (CatchmentsCount != PlotsCount)
  {
    for (int i = 0; i < DEM->GetRasterYSize(); i++)
    {
      for (int j = 0; j < DEM->GetRasterXSize(); j++)
      {
        if (i == 0 or j == 0 or i == DEM->GetRasterYSize()-1 or j == DEM->GetRasterXSize()-1)
        {
          Row[j] = -9999;
        }
        else
        {
          CatchmentsBand->RasterIO( GF_Read, 0, i-1, DEM->GetRasterXSize(), 1, CatchmentsUpperRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
          CatchmentsBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, CatchmentsMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
          CatchmentsBand->RasterIO( GF_Read, 0, i+1, DEM->GetRasterXSize(), 1, CatchmentsLowerRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
          PlotsRBand->RasterIO( GF_Read, 0, i-1, DEM->GetRasterXSize(), 1, PlotsUpperRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
          PlotsRBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, PlotsMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
          PlotsRBand->RasterIO( GF_Read, 0, i+1, DEM->GetRasterXSize(), 1, PlotsLowerRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
          DrainageBand->RasterIO( GF_Read, 0, i-1, DEM->GetRasterXSize(), 1, DrainageUpperRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
          DrainageBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, DrainageMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
          DrainageBand->RasterIO( GF_Read, 0, i+1, DEM->GetRasterXSize(), 1, DrainageLowerRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
          if (CatchmentsMiddleRow[j] == -9999)
          {
            if (CatchmentsUpperRow[j+1] != -9999 and DrainageMiddleRow[j] == 1 and PlotsMiddleRow[j] == PlotsUpperRow[j+1])
            {
              Row[j] = CatchmentsUpperRow[j+1];
            }
            else if (CatchmentsUpperRow[j] != -9999 and DrainageMiddleRow[j] == 2 and PlotsMiddleRow[j] == PlotsUpperRow[j])
            {
              Row[j] = CatchmentsUpperRow[j];
            }
            else if (CatchmentsUpperRow[j-1] != -9999 and DrainageMiddleRow[j] == 3 and PlotsMiddleRow[j] == PlotsUpperRow[j-1])
            {
              Row[j] = CatchmentsUpperRow[j-1];
            }
            else if (CatchmentsMiddleRow[j-1] != -9999 and DrainageMiddleRow[j] == 4 and PlotsMiddleRow[j] == PlotsMiddleRow[j-1])
            {
              Row[j] = CatchmentsMiddleRow[j-1];
            }
            else if (CatchmentsLowerRow[j-1] != -9999 and DrainageMiddleRow[j] == 5 and PlotsMiddleRow[j] == PlotsLowerRow[j-1])
            {
              Row[j] = CatchmentsLowerRow[j-1];
            }
            else if (CatchmentsLowerRow[j] != -9999 and DrainageMiddleRow[j] == 6 and PlotsMiddleRow[j] == PlotsLowerRow[j])
            {
              Row[j] = CatchmentsLowerRow[j];
            }
            else if (CatchmentsLowerRow[j+1] != -9999 and DrainageMiddleRow[j] == 7 and PlotsMiddleRow[j] == PlotsLowerRow[j+1])
            {
              Row[j] = CatchmentsLowerRow[j+1];
            }
            else if (CatchmentsMiddleRow[j+1] != -9999 and DrainageMiddleRow[j] == 8 and PlotsMiddleRow[j] == PlotsMiddleRow[j+1])
            {
              Row[j] = CatchmentsMiddleRow[j+1];
            }
            else
            {
              Row[j] = -9999;
            }
          }
          else
          {
            Row[j] = CatchmentsMiddleRow[j];
          }
        }
      }
      CatchmentsBand->RasterIO( GF_Write, 0, i, DEM->GetRasterXSize(), 1, Row, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
    }
    PlotsCount = 0;
    for (int i = 0; i < DEM->GetRasterYSize(); i++)
    {
      PlotsRBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, PlotsMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
      for (int j = 0; j < DEM->GetRasterXSize(); j++)
      {
        if (PlotsMiddleRow[j] != -9999)
        {
          PlotsCount += 1;
        }
      }
    }
    CatchmentsCount = 0;
    for (int i = 0; i < DEM->GetRasterYSize(); i++)
    {
      CatchmentsBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, CatchmentsMiddleRow, DEM->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
      for (int j = 0; j < DEM->GetRasterXSize(); j++)
      {
        if (CatchmentsMiddleRow[j] != -9999)
        {
          CatchmentsCount += 1;
        }
      }
    }
  }

  CPLFree(CatchmentsUpperRow);
  CPLFree(CatchmentsMiddleRow);
  CPLFree(CatchmentsLowerRow);
  CPLFree(PlotsUpperRow);
  CPLFree(PlotsMiddleRow);
  CPLFree(PlotsLowerRow);
  CPLFree(DrainageUpperRow);
  CPLFree(DrainageMiddleRow);
  CPLFree(DrainageLowerRow);
  CPLFree(Row);

  CatchmentsBand->FlushCache();
  CatchmentsBand->SetNoDataValue(-9999);
  Catchments->SetGeoTransform(mp_GeoTransformVal);
  Catchments->SetProjection(DEM->GetProjectionRef());
  GDALClose( Catchments );
  GDALClose( Outlets );
  GDALClose( PlotsR );
  GDALClose( Drainage );
  GDALClose( RasterID );
  GDALClose( DEM );

  BLOCK_EXIT(__PRETTY_FUNCTION__)
}


// =====================================================================
// =====================================================================


void LandProcessor::createSRFandLNR()
{
  BLOCK_ENTER(__PRETTY_FUNCTION__)

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int i, FieldIndex, Value, IDCat, IDCat2, IDCat3, IDN, ID2N, ID3N, IDCatR, IDPRO, IDPR, IDPR2,
      IDPlotO, IDPlot, IDPlot2, IDPlot3, IDPlotNew, IDPlot2New, IDPlot3New,
      FID, IDField, IDMax, N, Beginning, End, pos, IDPlotAbove;

  std::string FieldName, FileName, ID, ID2, ID3, IDOld, IDNew, ID2New,
  ID3New, IDAbove, IDTo, FaceA, FaceB, FaceAB, FaceBA;

  std::vector<std::string> FieldNamesTab, IDNewTab, ID2NewTab, ID3NewTab, IDOldTab, ID2OldTab,
  ID3OldTab, IDConnectionProblemTab, IDTab, IDDoublesTab, GivesToTab, FacesTab;

  std::vector<int> FIDTable, IDCatTab, IndexToRemove, IDPlotTab, IDPlotNewTab, IDPlot2NewTab, IDPlot3NewTab,
  IDPlotOldTab, IDPlot2OldTab, IDPlot3OldTab, IDPlotLargeTab, MissingIDTab, FIDSmallTab, FeaturesToDelete,
  IDPlotAllSmallTab, FIDSliverTab, FIDConnectionProblemTab, IDIndMax;

  std::vector < std::vector <int> > OutletsTable, ReceiversTable;

  std::vector< std::vector<int> > Table(5,std::vector<int >(0)), BigTable(5,std::vector<int >(0)), CatchmentsTab(10,std::vector<int >(0));;

  GDALDataset *Dataset;

  OGRDataSource *Plots, *Catchments, *Outlets, *Receivers, *DataSource, *Entities, *Union, *Final, *LNR, *SRF;
  OGRLayer *PlotsLayer, *CatchmentsLayer, *SQLLayer, *OutletsLayer, *ReceiversLayer, *EntitiesLayer,
  *UnionLayer, *FinalLayer, *LNRLayer, *SRFLayer;
  OGRFeature *CatchmentsFeature, *OutletsFeature, *ReceiversFeature,
  *SQLFeature, *EntitiesFeature, *NewFeature, *AggregatingFeature,
  *UnionFeature, *UnionFeatureToUpdate, *FinalFeature,
  *LNRFeature, *SRFFeature;
  OGRGeometry *CatchmentsGeometry, *OutletsGeometry, *ReceiversGeometry, *SQLGeometry,
  *EntitiesGeometry, *NewGeometry, *AggregatingGeometry, *UnionGeometry;
  OGRFeatureDefn *FeatureDefn;

  CPLErr ERR;

  // ======================================================================================
  // Checking if necessary vector files are present in the appropriate folder
  // ======================================================================================

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  if (!openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputCatchmentsRasterFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::createSRFandLNR(): " + m_OutputCatchmentsRasterFile + ": no such file in the output raster directory");
  }

  if (!openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputOutletsRasterFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::createSRFandLNR(): " + m_OutputOutletsRasterFile + ": no such file in the output raster directory");
  }

  if (!openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputReceiversRasterFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::createSRFandLNR(): " + m_OutputReceiversRasterFile + ": no such file in the output raster directory");
  }

  if (!openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputDownslopeRasterFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::createSRFandLNR(): " + m_OutputDownslopeRasterFile + ": no such file in the output raster directory");
  }

  // ================================================================================================
  // Creating a catchment vector file based on the catchment raster (IDCat is the catchment ID)
  // ================================================================================================

  Plots = OGRSFDriverRegistrar::Open( getInputVectorPath(m_InputPlotsFile).c_str(), TRUE );

  mp_SRS = Plots->GetLayer(0)->GetSpatialRef();

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputCatchmentsVectorFile).c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputCatchmentsVectorFile).c_str());
  }

  Catchments = mp_VectorDriver->CreateDataSource(getOutputVectorPath(m_OutputCatchmentsVectorFile).c_str());
  CatchmentsLayer = Catchments->CreateLayer("catchments", mp_SRS, wkbPolygon);

  FieldName = "IDCat";
  createField(CatchmentsLayer, FieldName.c_str(), OFTInteger);

  Dataset = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputCatchmentsRasterFile).c_str(), GA_ReadOnly );

  ERR = GDALPolygonize(getRasterBand(Dataset,1), nullptr, CatchmentsLayer, 0, nullptr, nullptr, nullptr);

  CatchmentsLayer->SyncToDisk();

  CatchmentsLayer->ResetReading();

  CatchmentsFeature = nullptr;

  while ((CatchmentsFeature = CatchmentsLayer->GetNextFeature()) != nullptr)
  {
    if (CatchmentsFeature->GetFieldAsInteger(0) == -9999)
    {
      CatchmentsFeature->SetField(0, 0);
    }
    CatchmentsLayer->SetFeature(CatchmentsFeature);
  }

  OGRDataSource::DestroyDataSource(Catchments);

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputOutletsVectorFile).c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputOutletsVectorFile).c_str());
  }

  Outlets = mp_VectorDriver->CreateDataSource(getOutputVectorPath(m_OutputOutletsVectorFile).c_str());
  OutletsLayer = Outlets->CreateLayer("outlets", mp_SRS, wkbPolygon);

  FieldNamesTab.clear();
  FieldNamesTab = {"IDCatO", "IDPlotO", "IDPRO"};

  for (int i = 0; i < FieldNamesTab.size(); i++)
  {
    createField(OutletsLayer, FieldNamesTab[i].c_str(), OFTInteger);
  }

  GDALClose( Dataset );

  Dataset = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputOutletsRasterFile).c_str(), GA_ReadOnly );
  ERR = GDALPolygonize(getRasterBand(Dataset,1), nullptr, OutletsLayer, 0, nullptr, nullptr, nullptr);

  OutletsLayer->SyncToDisk();

  OutletsLayer->ResetReading();

  OutletsFeature = nullptr;

  while ((OutletsFeature = OutletsLayer->GetNextFeature()) != nullptr)
  {
    if (OutletsFeature->GetFieldAsInteger(0) == -9999)
    {
      OutletsFeature->SetField(0, 0);
    }
    OutletsLayer->SetFeature(OutletsFeature);
  }

  GDALClose( Dataset );

  Dataset = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputPlotsRasterFile).c_str(), GA_ReadOnly );
  getGeoTransform(Dataset);

  OutletsLayer->ResetReading();

  OutletsFeature = nullptr;

  while ((OutletsFeature = OutletsLayer->GetNextFeature()) != nullptr)
  {
    if (OutletsFeature->GetFieldAsInteger(0) == 0)
    {
      OutletsFeature->SetField(1, 0);
    }
    else
    {
      OutletsGeometry = OutletsFeature->GetGeometryRef();
      getCentroidPoint(OutletsGeometry);
      Value = extractFromRasterToPoint(Dataset,1);
      if (Value == -9999)
      {
        OutletsFeature->SetField(1, 0);
      }
      else
      {
        OutletsFeature->SetField(1, Value);
      }
    }
    OutletsLayer->SetFeature(OutletsFeature);
  }

  GDALClose( Dataset );

  Dataset = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputDownslopeRasterFile).c_str(), GA_ReadOnly );
  getGeoTransform(Dataset);

  OutletsLayer->ResetReading();

  OutletsFeature = nullptr;

  while ((OutletsFeature = OutletsLayer->GetNextFeature()) != nullptr)
  {
    if (OutletsFeature->GetFieldAsInteger(0) == 0)
    {
      OutletsFeature->SetField(2, 0);
    }
    else
    {
      OutletsGeometry = OutletsFeature->GetGeometryRef();
      getCentroidPoint(OutletsGeometry);
      Value = extractFromRasterToPoint(Dataset, 1);
      if (Value == -9999)
      {
        OutletsFeature->SetField(2, 0);
      }
      else
      {
        OutletsFeature->SetField(2, Value);
      }
    }
    OutletsLayer->SetFeature(OutletsFeature);
  }

  OGRDataSource::DestroyDataSource(Outlets);

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputReceiversVectorFile).c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputReceiversVectorFile).c_str());
  }

  Receivers = mp_VectorDriver->CreateDataSource(getOutputVectorPath(m_OutputReceiversVectorFile).c_str());
  ReceiversLayer = Receivers->CreateLayer("receivers", mp_SRS, wkbPolygon);

  FieldNamesTab.clear();
  FieldNamesTab = {"IDPRR", "IDPlotR", "IDCatR"};

  for (int i = 0; i < FieldNamesTab.size(); i++)
  {
    createField(ReceiversLayer, FieldNamesTab.at(i).c_str(), OFTInteger);
  }

  FieldNamesTab.clear();

  GDALClose( Dataset );

  Dataset = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputReceiversRasterFile).c_str(), GA_ReadOnly );
  ERR = GDALPolygonize(getRasterBand(Dataset,1), nullptr, ReceiversLayer, 0, nullptr, nullptr, nullptr);

  ReceiversLayer->SyncToDisk();

  ReceiversLayer->ResetReading();

  ReceiversFeature = nullptr;

  while ((ReceiversFeature = ReceiversLayer->GetNextFeature()) != nullptr)
  {
    if (ReceiversFeature->GetFieldAsInteger(0) == -9999)
    {
      ReceiversFeature->SetField(0, 0);
    }
    ReceiversLayer->SetFeature(ReceiversFeature);
  }

  GDALClose( Dataset );

  Dataset = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputPlotsRasterFile).c_str(), GA_ReadOnly );
  getGeoTransform(Dataset);

  ReceiversLayer->ResetReading();

  ReceiversFeature = nullptr;

  while ((ReceiversFeature = ReceiversLayer->GetNextFeature()) != nullptr)
  {
    if (ReceiversFeature->GetFieldAsInteger(0) == 0)
    {
      ReceiversFeature->SetField(1, 0);
    }
    else
    {
      ReceiversGeometry = ReceiversFeature->GetGeometryRef();
      getCentroidPoint(ReceiversGeometry);
      Value = extractFromRasterToPoint(Dataset, 1);
      if (Value == -9999)
      {
        ReceiversFeature->SetField(1, 0);
      }
      else
      {
        ReceiversFeature->SetField(1, Value);
      }
    }
    ReceiversLayer->SetFeature(ReceiversFeature);
  }

  GDALClose( Dataset );

  Dataset = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputCatchmentsRasterFile).c_str(), GA_ReadOnly );
  getGeoTransform(Dataset);

  ReceiversLayer->ResetReading();

  ReceiversFeature = nullptr;

  while ((ReceiversFeature = ReceiversLayer->GetNextFeature()) != nullptr)
  {
    if (ReceiversFeature->GetFieldAsInteger(0) == 0)
    {
      ReceiversFeature->SetField(2, 0);
    }
    else
    {
      ReceiversGeometry = ReceiversFeature->GetGeometryRef();
      getCentroidPoint(ReceiversGeometry);
      Value = extractFromRasterToPoint(Dataset, 1);
      if (Value == -9999)
      {
        ReceiversFeature->SetField(2, 0);
      }
      else
      {
        ReceiversFeature->SetField(2, Value);
      }
    }
    ReceiversLayer->SetFeature(ReceiversFeature);
  }

  OGRDataSource::DestroyDataSource(Receivers);

  GDALClose( Dataset );
  OGRDataSource::DestroyDataSource(Plots);

  // ======================================================================
  // Regrouping catchments that have the same IDCat
  // ======================================================================

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputCatchmentsGroupedVectorFile).c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputCatchmentsGroupedVectorFile).c_str());
  }

  Catchments = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputCatchmentsVectorFile).c_str(), TRUE );
  m_SQLRequest = "SELECT ST_Union(geometry), IDCat FROM " + m_OutputCatchmentsVectorFile.substr(0, m_OutputCatchmentsVectorFile.find(".")) + " GROUP BY IDCat";
  FileName = m_OutputCatchmentsGroupedVectorFile.substr(0, m_OutputCatchmentsGroupedVectorFile.find("."));
  SQLLayer = Catchments->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
  Catchments->CopyLayer(SQLLayer, FileName.c_str(), nullptr);
  Catchments->ReleaseResultSet(SQLLayer);
  OGRDataSource::DestroyDataSource(Catchments);

  // ===============================================================================
  // Populating catchments vector attribute table
  // using outlets and receivers vectors
  // The attributes obtain during this operation will be used for further regrouping
  // ===============================================================================

  Outlets = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputOutletsVectorFile).c_str(), TRUE );
  OutletsLayer = Outlets->GetLayer(0);

  Receivers = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputReceiversVectorFile).c_str(), TRUE );
  ReceiversLayer = Receivers->GetLayer(0);

  OutletsTable.resize(OutletsLayer->GetFeature(0)->GetFieldCount());

  for(int i = 0 ; i < OutletsLayer->GetFeature(0)->GetFieldCount(); ++i)
  {
    OutletsTable[i].resize(OutletsLayer->GetFeatureCount());
  }

  for (int i = 0; i < OutletsLayer->GetFeature(0)->GetFieldCount(); i++)
  {
    for (int j = 0; j < OutletsLayer->GetFeatureCount(); j++)
    {
      OutletsTable[i][j] = OutletsLayer->GetFeature(j)->GetFieldAsInteger(i);
    }
  }

  ReceiversTable.resize(ReceiversLayer->GetFeature(0)->GetFieldCount());

  for(int i = 0 ; i < ReceiversLayer->GetFeature(0)->GetFieldCount(); ++i)
  {
    ReceiversTable[i].resize(ReceiversLayer->GetFeatureCount());
  }

  for (int i = 0; i < ReceiversLayer->GetFeature(0)->GetFieldCount(); i++)
  {
    for (int j = 0; j < ReceiversLayer->GetFeatureCount(); j++)
    {
      ReceiversTable[i][j] = ReceiversLayer->GetFeature(j)->GetFieldAsInteger(i);
    }
  }

  OGRDataSource::DestroyDataSource(Outlets);
  OGRDataSource::DestroyDataSource(Receivers);

  // ==============================================================================
  // Creating necessary fields in the attribute table and populating them
  // ("IDPlot", "IDPR", "IDCat2", "IDPlot2", "IDPR2", "IDCat3", "IDPlot3", "IDPR3")
  // based on previously created OutletsTable and ReceiversTable.
  // ==============================================================================

  Catchments = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputCatchmentsGroupedVectorFile).c_str(), TRUE );
  CatchmentsLayer = Catchments->GetLayer(0);

  FieldNamesTab.clear();
  FieldNamesTab = {"IDPlot", "IDPR", "IDCat2", "IDPlot2", "IDPR2", "IDCat3", "IDPlot3", "IDPR3"};

  for (int i = 0; i < FieldNamesTab.size(); i++)
  {
    createField(CatchmentsLayer, FieldNamesTab[i].c_str(), OFTInteger);
  }

  FieldNamesTab.clear();
  FieldNamesTab = {"ID","ID2","ID3"};

  for (int i = 0; i < FieldNamesTab.size(); i++)
  {
    createField(CatchmentsLayer, FieldNamesTab[i].c_str(),  OFTString);
  }

  FieldNamesTab.clear();

  CatchmentsLayer->SyncToDisk();
  CatchmentsLayer->ResetReading();
  CatchmentsFeature = nullptr;

  while ((CatchmentsFeature = CatchmentsLayer->GetNextFeature()) != nullptr)
  {
    IDCat = CatchmentsFeature->GetFieldAsInteger("IDCat");
    if (IDCat != 0)
    {
      if (std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat) != OutletsTable[0].end())
      {
        CatchmentsFeature->SetField("IDPlot", OutletsTable[1].at(std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat) - OutletsTable[0].begin()));
        CatchmentsFeature->SetField("IDPR", OutletsTable[2].at(std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat) - OutletsTable[0].begin()));
        CatchmentsLayer->SetFeature(CatchmentsFeature);
      }
    }
    else
    {
      CatchmentsFeature->SetField("IDPlot", 0);
      CatchmentsFeature->SetField("IDPR", 0);
      CatchmentsLayer->SetFeature(CatchmentsFeature);
    }
    IDPR = CatchmentsFeature->GetFieldAsInteger("IDPR");
    if (IDPR != 0)
    {
      if (std::find(ReceiversTable[0].begin(), ReceiversTable[0].end(), IDPR) != ReceiversTable[0].end())
      {
        CatchmentsFeature->SetField("IDCat2", ReceiversTable[2].at(std::find(ReceiversTable[0].begin(), ReceiversTable[0].end(), IDPR) - ReceiversTable[0].begin()));
        CatchmentsLayer->SetFeature(CatchmentsFeature);
      }
      else
      {
        CatchmentsFeature->SetField("IDCat2", 0);
        CatchmentsLayer->SetFeature(CatchmentsFeature);
      }
    }
    else
    {
      CatchmentsFeature->SetField("IDCat2", 0);
      CatchmentsLayer->SetFeature(CatchmentsFeature);
    }
    IDCat2 = CatchmentsFeature->GetFieldAsInteger("IDCat2");
    if (IDCat2 != 0)
    {
      if (std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat2) != OutletsTable[0].end())
      {
        CatchmentsFeature->SetField("IDPlot2", OutletsTable[1].at(std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat2) - OutletsTable[0].begin()));
        CatchmentsFeature->SetField("IDPR2", OutletsTable[2].at(std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat2) - OutletsTable[0].begin()));
        CatchmentsLayer->SetFeature(CatchmentsFeature);
      }
      else
      {
        CatchmentsFeature->SetField("IDPlot2", 0);
        CatchmentsFeature->SetField("IDPR2", 0);
        CatchmentsLayer->SetFeature(CatchmentsFeature);
      }
    }
    else
    {
      CatchmentsFeature->SetField("IDPlot2", 0);
      CatchmentsFeature->SetField("IDPR2", 0);
      CatchmentsLayer->SetFeature(CatchmentsFeature);
    }
    IDPR2 = CatchmentsFeature->GetFieldAsInteger("IDPR2");
    if (IDPR2 != 0)
    {
      if (std::find(ReceiversTable[0].begin(), ReceiversTable[0].end(), IDPR2) != ReceiversTable[0].end())
      {
        CatchmentsFeature->SetField("IDCat3", ReceiversTable[2].at(std::find(ReceiversTable[0].begin(), ReceiversTable[0].end(), IDPR2) - ReceiversTable[0].begin()));
        CatchmentsLayer->SetFeature(CatchmentsFeature);
      }
      else
      {
        CatchmentsFeature->SetField("IDCat3", 0);
        CatchmentsLayer->SetFeature(CatchmentsFeature);
      }
    }
    else
    {
      CatchmentsFeature->SetField("IDCat3", 0);
      CatchmentsLayer->SetFeature(CatchmentsFeature);
    }
    IDCat3 = CatchmentsFeature->GetFieldAsInteger("IDCat3");
    if (IDCat3 != 0)
    {
      if (std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat3) != OutletsTable[0].end())
      {
        CatchmentsFeature->SetField("IDPlot3", OutletsTable[1].at(std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat3) - OutletsTable[0].begin()));
        CatchmentsFeature->SetField("IDPR3", OutletsTable[2].at(std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat3) - OutletsTable[0].begin()));
        CatchmentsLayer->SetFeature(CatchmentsFeature);
      }
      else
      {
        CatchmentsFeature->SetField("IDPlot3", 0);
        CatchmentsFeature->SetField("IDPR3", 0);
        CatchmentsLayer->SetFeature(CatchmentsFeature);
      }
    }
    else
    {
      CatchmentsFeature->SetField("IDPlot3", 0);
      CatchmentsFeature->SetField("IDPR3", 0);
      CatchmentsLayer->SetFeature(CatchmentsFeature);
    }
  }

  CatchmentsLayer->SyncToDisk();
  OGRDataSource::DestroyDataSource(Catchments);

  // ==========================================================================================================
  // Populating catchments vector attributes ID, ID2 and ID3
  // with "IDPlot+N+Number"-like identifiers based on their attributes and their
  // (respective) location.
  // Zero polygon is the polygon created outside the original zone of interest.
  // It is labeled as following: ID - "0N1", ID2 - "0N1", ID3 - "0N1".
  // It does not drain anywhere and all the catchments ultimately suppose to drain into it.
  // Labeling begins with the catchments that drain directly into zero polygon.
  // They are labeled according to the same concept: "IDPlot+N+Number"
  // The labeling procedure is similar to that of the catchment raster creation:
  // - first, establishing the "exit points" and working inwards progressively labeling everything that drains
  // in that direction.
  // - catchments that belong to the same plot and drain into the same catchment get same ID label
  // only if they share a limit.
  // ==========================================================================================================

  Catchments = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputCatchmentsGroupedVectorFile).c_str(), TRUE );
  CatchmentsLayer = Catchments->GetLayer(0);
  CatchmentsLayer->ResetReading();

  while ((CatchmentsFeature = CatchmentsLayer->GetNextFeature()) != nullptr)
  {
    CatchmentsTab[0].push_back(CatchmentsFeature->GetFID());
    CatchmentsTab[1].push_back(CatchmentsFeature->GetFieldAsInteger("IDCat"));
    CatchmentsTab[2].push_back(CatchmentsFeature->GetFieldAsInteger("IDPlot"));
    CatchmentsTab[3].push_back(CatchmentsFeature->GetFieldAsInteger("IDCat2"));
    CatchmentsTab[4].push_back(CatchmentsFeature->GetFieldAsInteger("IDPlot2"));
    CatchmentsTab[5].push_back(CatchmentsFeature->GetFieldAsInteger("IDCat3"));
    CatchmentsTab[6].push_back(CatchmentsFeature->GetFieldAsInteger("IDPlot3"));
    CatchmentsTab[7].push_back(0);
    CatchmentsTab[8].push_back(0);
    CatchmentsTab[9].push_back(0);
  }

  OGRDataSource::DestroyDataSource(Catchments);

  Table.clear();

  for (int i = 0; i < CatchmentsTab[0].size(); i++)
  {
    FID = CatchmentsTab[0][i];
    IDPlot = CatchmentsTab[2][i];
    IDPlot2 = CatchmentsTab[4][i];
    if (IDPlot2 == 0)
    {
      if (std::find(Table[0].begin(), Table[0].end(), IDPlot) == Table[0].end())
      {
        if(IDPlot != 0)
        {
          Table[0].push_back(IDPlot);
          Table[1].push_back(1);
          Table[2].push_back(IDPlot2);
          Table[3].push_back(1);
          Table[4].push_back(FID);
          CatchmentsTab[7][i] = 1;
          CatchmentsTab[8][i] = 1;
          CatchmentsTab[9][i] = 1;
        }
        else
        {
          Table[0].push_back(0);
          Table[1].push_back(1);
          Table[2].push_back(0);
          Table[3].push_back(1);
          Table[4].push_back(FID);
          CatchmentsTab[7][i] = 1;
          CatchmentsTab[8][i] = 1;
          CatchmentsTab[9][i] = 1;;
        }
      }
    }
  }

  // ==============================================================================================================================
  // Populating the ID, ID2 and ID3 of those catchments that drain into zero polygon and can be regrouped with already labeled
  // catchments
  // ==============================================================================================================================

  BigTable[0] = Table[0];
  BigTable[1] = Table[1];
  BigTable[2] = Table[2];
  BigTable[3] = Table[3];
  BigTable[4] = Table[4];
  FIDTable.clear();

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );

  while (Table[0].size() != 0)
  {
    i = 0;
    while (i < Table[0].size())
    {
      if (Table[0][i] != 0)
      {
        pos  = std::find(CatchmentsTab[0].begin(), CatchmentsTab[0].end(), Table[4][i]) - CatchmentsTab[0].begin();
        FID = CatchmentsTab[0].at(pos);
        IDPlot = CatchmentsTab[2].at(pos);
        IDN = CatchmentsTab[7].at(pos);
        ID2N = CatchmentsTab[8].at(pos);
        ID3N = CatchmentsTab[9].at(pos);
        m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputCatchmentsGroupedVectorFile.substr(0, m_OutputCatchmentsGroupedVectorFile.find(".")) + " WHERE IDPlot2 == 0 "
            "AND IDPlot == " + std::to_string(IDPlot) + " AND "
            "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + m_OutputCatchmentsGroupedVectorFile.substr(0, m_OutputCatchmentsGroupedVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FID) +"))) > 0";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          FIDTable.clear();
          for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
          {
            if (CatchmentsTab[7][SQLLayer->GetFeature(j)->GetFieldAsInteger("RID")] == 0 )
            {
              FIDTable.push_back(SQLLayer->GetFeature(j)->GetFieldAsInteger("RID"));
            }
          }
          DataSource->ReleaseResultSet(SQLLayer);
          for (int k = 0; k < FIDTable.size(); k++)
          {
            pos = std::find(CatchmentsTab[0].begin(), CatchmentsTab[0].end(), FIDTable[k]) - CatchmentsTab[0].begin();
            CatchmentsTab[7][pos] = IDN;
            CatchmentsTab[8][pos] = ID2N;
            CatchmentsTab[9][pos] = ID3N;
            Table[0].push_back(IDPlot);
            Table[1].push_back(IDN);
            Table[2].push_back(0);
            Table[3].push_back(1);
            Table[4].push_back(FIDTable[k]);
            BigTable[0].push_back(IDPlot);
            BigTable[1].push_back(IDN);
            BigTable[2].push_back(0);
            BigTable[3].push_back(1);
            BigTable[4].push_back(FIDTable[k]);
          }
          FIDTable.clear();
        }
        else
        {
          DataSource->ReleaseResultSet(SQLLayer);
        }
      }
      i = i+1;
    }
    IDMax = Table[1].back();
    Table[0].clear();
    Table[1].clear();
    Table[2].clear();
    Table[3].clear();
    Table[4].clear();
    for (int j = 0; j < CatchmentsTab[0].size(); j++)
    {
      FID = CatchmentsTab[0][j];
      IDN = CatchmentsTab[7][j];
      IDPlot = CatchmentsTab[2][j];
      IDPlot2 = CatchmentsTab[4][j];
      if (IDPlot2 == 0 and IDN == 0)
      {
        if (std::find(Table[0].begin(), Table[0].end(), IDPlot) == Table[0].end())
        {
          IDN = IDMax+1;
          CatchmentsTab[7][j] = IDN;
          CatchmentsTab[8][j] = 1;
          CatchmentsTab[9][j] = 1;
          Table[0].push_back(IDPlot);
          Table[1].push_back(IDMax+1);
          Table[2].push_back(0);
          Table[3].push_back(1);
          Table[4].push_back(FID);
          BigTable[0].push_back(IDPlot);
          BigTable[1].push_back(IDMax+1);
          BigTable[2].push_back(0);
          BigTable[3].push_back(1);
          BigTable[4].push_back(FID);
        }
      }
    }
  }

  FIDTable.clear();

  for (int i = 0; i < CatchmentsTab[0].size(); i++)
  {
    IDCat2 = CatchmentsTab[3][i];
    ID2N = CatchmentsTab[8][i];
    if (ID2N == 0)
    {
      for (int j = 0; j < CatchmentsTab[0].size(); j++)
      {
        if (CatchmentsTab[1][j] == IDCat2 and CatchmentsTab[7][j] != 0)
        {
          CatchmentsTab[8][i] = CatchmentsTab[7][j];
          CatchmentsTab[9][i] = CatchmentsTab[8][j];
        }
      }
    }
  }

  // =========================================================================================
  // Setting the rest of the catchments with ID, ID2 and ID3 using the same principal as above
  // =========================================================================================

  N = 0;

  for (int i = 0; i < CatchmentsTab[0].size(); i++)
  {
    if (CatchmentsTab[7][i] == 0)
    {
      N += 1;
    }
  }

  while (N > 0)
  {
    Table[0].clear();
    Table[1].clear();
    Table[2].clear();
    Table[3].clear();
    Table[4].clear();
    for (int i = 0; i < CatchmentsTab[0].size(); i++)
    {
      FID = CatchmentsTab[0][i];
      IDPlot = CatchmentsTab[2][i];
      IDN = CatchmentsTab[7][i];
      IDPlot2 = CatchmentsTab[4][i];
      ID2N = CatchmentsTab[8][i];
      if (IDN == 0 and ID2N != 0)
      {
        if (std::find(Table[0].begin(), Table[0].end(), IDPlot) == Table[0].end() and
            std::find(BigTable[0].begin(), BigTable[0].end(), IDPlot) == BigTable[0].end())
        {
          CatchmentsTab[7][i] = 1;
          Table[0].push_back(IDPlot);
          Table[1].push_back(1);
          Table[2].push_back(IDPlot2);
          Table[3].push_back(ID2N);
          Table[4].push_back(FID);
          BigTable[0].push_back(IDPlot);
          BigTable[1].push_back(1);
          BigTable[2].push_back(IDPlot2);
          BigTable[3].push_back(ID2N);
          BigTable[4].push_back(FID);
        }
        else if (std::find(Table[0].begin(), Table[0].end(), IDPlot) == Table[0].end() and
            std::find(BigTable[0].begin(), BigTable[0].end(), IDPlot) != BigTable[0].end())
        {
          IDMax = 0;
          for (int j = 0; j < BigTable[0].size(); j++)
          {
            if (BigTable[0][j] == IDPlot and BigTable[1][j] > IDMax)
            {
              IDMax = BigTable[1][j];
            }
          }
          CatchmentsTab[7][i] = IDMax+1;
          Table[0].push_back(IDPlot);
          Table[1].push_back(IDMax+1);
          Table[2].push_back(IDPlot2);
          Table[3].push_back(ID2N);
          Table[4].push_back(FID);
          BigTable[0].push_back(IDPlot);
          BigTable[1].push_back(IDMax+1);
          BigTable[2].push_back(IDPlot2);
          BigTable[3].push_back(ID2N);
          BigTable[4].push_back(FID);
        }
      }
    }
    while (Table[0].size() != 0)
    {
      i = 0;
      while (i < Table[0].size())
      {
        if (Table[0][i] != 0)
        {
          FID = CatchmentsTab[0][Table[4][i]];
          IDPlot = CatchmentsTab[2][Table[4][i]];
          IDN = CatchmentsTab[7][Table[4][i]];
          IDPlot2 = CatchmentsTab[4][Table[4][i]];
          ID2N = CatchmentsTab[8][Table[4][i]];
          ID3N = CatchmentsTab[9][Table[4][i]];
          m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputCatchmentsGroupedVectorFile.substr(0, m_OutputCatchmentsGroupedVectorFile.find(".")) + " "
              "WHERE "
              "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + m_OutputCatchmentsGroupedVectorFile.substr(0, m_OutputCatchmentsGroupedVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FID) + "))) > 0";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() != 0)
          {
            for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
            {
              if (CatchmentsTab[7][SQLLayer->GetFeature(j)->GetFieldAsInteger("RID")] == 0 and
                  CatchmentsTab[2][SQLLayer->GetFeature(j)->GetFieldAsInteger("RID")] == IDPlot and
                  CatchmentsTab[8][SQLLayer->GetFeature(j)->GetFieldAsInteger("RID")] == ID2N and
                  CatchmentsTab[4][SQLLayer->GetFeature(j)->GetFieldAsInteger("RID")] == IDPlot2)
              {
                FIDTable.push_back(SQLLayer->GetFeature(j)->GetFieldAsInteger("RID"));
              }
            }
            DataSource->ReleaseResultSet(SQLLayer);
            for (int k = 0; k < FIDTable.size(); k++)
            {
              CatchmentsTab[7][FIDTable[k]] = IDN;
              Table[0].push_back(IDPlot);
              Table[1].push_back(IDN);
              Table[2].push_back(IDPlot2);
              Table[3].push_back(ID2N);
              Table[4].push_back(FIDTable[k]);
              BigTable[0].push_back(IDPlot);
              BigTable[1].push_back(IDN);
              BigTable[2].push_back(IDPlot2);
              BigTable[3].push_back(ID2N);
              BigTable[4].push_back(FIDTable[k]);
            }
            FIDTable.clear();
          }
          else
          {
            DataSource->ReleaseResultSet(SQLLayer);
          }
        }
        i = i+1;
      }
      Table[0].clear();
      Table[1].clear();
      Table[2].clear();
      Table[3].clear();
      Table[4].clear();
      for (int j = 0; j < CatchmentsTab[0].size(); j++)
      {
        FID = CatchmentsTab[0][j];
        IDPlot = CatchmentsTab[2][j];
        IDN = CatchmentsTab[7][j];
        IDPlot2 = CatchmentsTab[4][j];
        ID2N = CatchmentsTab[8][j];
        if (IDN == 0 and ID2N != 0)
        {
          if (std::find(Table[0].begin(), Table[0].end(), IDPlot) == Table[0].end())
          {
            IDMax = 0;
            for (int k = 0; k < BigTable[0].size(); k++)
            {
              if (BigTable[0][k] == IDPlot and BigTable[1][k] > IDMax)
              {
                IDMax = BigTable[1][k];
              }
            }
            CatchmentsTab[7][j] = IDMax+1;
            Table[0].push_back(IDPlot);
            Table[1].push_back(IDMax+1);
            Table[2].push_back(IDPlot2);
            Table[3].push_back(ID2N);
            Table[4].push_back(FID);
            BigTable[0].push_back(IDPlot);
            BigTable[1].push_back(IDMax+1);
            BigTable[2].push_back(IDPlot2);
            BigTable[3].push_back(ID2N);
            BigTable[4].push_back(FID);
          }
        }
      }
      for (int i = 0; i < CatchmentsTab[0].size(); i++)
      {
        IDCat2 = CatchmentsTab[3][i];
        ID2N = CatchmentsTab[8][i];
        if (ID2N == 0)
        {
          for (int j = 0; j < CatchmentsTab[0].size(); j++)
          {
            if (CatchmentsTab[1][j] == IDCat2 and CatchmentsTab[7][j] != 0)
            {
              CatchmentsTab[8][i] = CatchmentsTab[7][j];
              CatchmentsTab[9][i] = CatchmentsTab[8][j];
            }
          }
        }
      }
    }
    N = 0;
    for (int i = 0; i < CatchmentsTab[0].size(); i++)
    {
      if (CatchmentsTab[7][i] == 0)
      {
        N += 1;
      }
    }
  }

  Table[0].clear();
  Table[1].clear();
  Table[2].clear();
  Table[3].clear();
  Table[4].clear();
  BigTable[0].clear();
  BigTable[1].clear();
  BigTable[2].clear();
  BigTable[3].clear();
  BigTable[4].clear();

  OGRDataSource::DestroyDataSource(DataSource);

  Catchments = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputCatchmentsGroupedVectorFile).c_str(), TRUE );
  CatchmentsLayer = Catchments->GetLayer(0);

  for (int i = 0; i < CatchmentsTab[0].size(); i++)
  {
    IDPlot = CatchmentsTab[2][i];
    IDPlot2 = CatchmentsTab[4][i];
    IDPlot3 = CatchmentsTab[6][i];
    IDN = CatchmentsTab[7][i];
    ID2N = CatchmentsTab[8][i];
    ID3N = CatchmentsTab[9][i];
    CatchmentsFeature = CatchmentsLayer->GetFeature(CatchmentsTab[0][i]);
    IDNew = std::to_string(IDPlot)+"N"+std::to_string(IDN);
    ID2New = std::to_string(IDPlot2)+"N"+std::to_string(ID2N);
    ID3New = std::to_string(IDPlot3)+"N"+std::to_string(ID3N);
    CatchmentsFeature->SetField("ID", IDNew.c_str());
    CatchmentsFeature->SetField("ID2", ID2New.c_str());
    CatchmentsFeature->SetField("ID3", ID3New.c_str());
    CatchmentsLayer->SetFeature(CatchmentsFeature);
  }

  OGRDataSource::DestroyDataSource(Catchments);

  CatchmentsTab[0].clear();
  CatchmentsTab[1].clear();
  CatchmentsTab[2].clear();
  CatchmentsTab[3].clear();
  CatchmentsTab[4].clear();
  CatchmentsTab[5].clear();
  CatchmentsTab[6].clear();
  CatchmentsTab[7].clear();
  CatchmentsTab[8].clear();
  CatchmentsTab[9].clear();

  // ======================================================================
  //  Regrouping catchments that have the same attributes into new entities
  // ======================================================================

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str());
  }

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
  m_SQLRequest = "SELECT ST_Union(geometry), IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3 FROM " + m_OutputCatchmentsGroupedVectorFile.substr(0, m_OutputCatchmentsGroupedVectorFile.find(".")) + " GROUP BY IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3";
  SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
  FileName = m_OutputEntitiesVectorFile.substr(0, m_OutputEntitiesVectorFile.find("."));
  DataSource->CopyLayer(SQLLayer, FileName.c_str(), nullptr);
  DataSource->ReleaseResultSet(SQLLayer);
  OGRDataSource::DestroyDataSource(DataSource);

  // ============================================================================================
  // Re-aggregating sub-entities if their size is below selected seed (user-defied parameter)
  // ============================================================================================

  // ============================================================================================
  //  Re-aggregating sub-entities that have no entity larger then selected seed (all the entities
  //  representing the parcel are smaller the the seed)
  // ============================================================================================

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
  m_SQLRequest = "SELECT ROWID as RID, IDPlot FROM " + m_OutputEntitiesVectorFile.substr(0, m_OutputEntitiesVectorFile.find(".")) + " WHERE ST_Area(geometry) > " + std::to_string(m_MinEntSize);
  SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());

  for(int i = 0; i < SQLLayer->GetFeatureCount(); i++)
  {
    if (std::find(IDPlotLargeTab.begin(), IDPlotLargeTab.end(), SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot")) == IDPlotLargeTab.end())
    {
      IDPlotLargeTab.push_back(SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot"));
    }
  }

  DataSource->ReleaseResultSet(SQLLayer);
  m_SQLRequest = "SELECT ROWID as RID, IDPlot FROM " + m_OutputEntitiesVectorFile.substr(0, m_OutputEntitiesVectorFile.find(".")) + " WHERE ST_Area(geometry) <= " + std::to_string(m_MinEntSize);
  SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());

  for(int i = 0; i < SQLLayer->GetFeatureCount(); i++)
  {
    if (std::find(IDPlotLargeTab.begin(), IDPlotLargeTab.end(), SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot")) == IDPlotLargeTab.end() and std::find(IDPlotAllSmallTab.begin(), IDPlotAllSmallTab.end(), SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot")) == IDPlotAllSmallTab.end())
    {
      IDPlotAllSmallTab.push_back(SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot"));
    }
  }

  DataSource->ReleaseResultSet(SQLLayer);

  IDNewTab.clear();
  ID2NewTab.clear();
  ID3NewTab.clear();
  IDPlot2NewTab.clear();
  IDPlot3NewTab.clear();

  if (IDPlotAllSmallTab.size() != 0)
  {
    for (int i = 0; i < IDPlotAllSmallTab.size(); i++)
    {
      m_SQLRequest = "SELECT IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3 FROM " + m_OutputEntitiesVectorFile.substr(0, m_OutputEntitiesVectorFile.find(".")) + " WHERE IDPlot == " + std::to_string(IDPlotAllSmallTab[i]) + " "
          "AND IDPlot != IDPlot3 ORDER BY ST_Area(geometry) DESC";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        IDNewTab.push_back(SQLLayer->GetFeature(0)->GetFieldAsString("ID"));
        IDPlot2NewTab.push_back(SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot2"));
        ID2NewTab.push_back(SQLLayer->GetFeature(0)->GetFieldAsString("ID2"));
        IDPlot3NewTab.push_back(SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot3"));
        ID3NewTab.push_back(SQLLayer->GetFeature(0)->GetFieldAsString("ID3"));
      }
      DataSource->ReleaseResultSet(SQLLayer);
    }
    Entities = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
    EntitiesLayer = Entities->GetLayer(0);
    FeatureDefn = EntitiesLayer->GetLayerDefn();
    IndexToRemove.clear();
    for (int i = 0; i < IDPlotAllSmallTab.size(); i++)
    {
      NewGeometry = nullptr;
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + m_OutputEntitiesVectorFile.substr(0, m_OutputEntitiesVectorFile.find(".")) + " WHERE IDPlot == " + std::to_string(IDPlotAllSmallTab[i]) + " ORDER BY ST_Area(geometry) DESC";
      SQLLayer = Entities->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      IndexToRemove.clear();
      for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
      {
        SQLFeature = SQLLayer->GetFeature(j);
        SQLGeometry = SQLFeature->GetGeometryRef();
        IndexToRemove.push_back(SQLFeature->GetFieldAsInteger("RID"));
        if (NewGeometry == nullptr)
        {
          NewGeometry = SQLGeometry;
        }
        else
        {
          NewGeometry = SQLGeometry->Union(NewGeometry);
        }
      }
      Entities->ReleaseResultSet(SQLLayer);
      EntitiesLayer = Entities->GetLayer(0);
      for (int j = 0; j < IndexToRemove.size(); j++)
      {
        EntitiesLayer->DeleteFeature(IndexToRemove[j]);
      }
      m_SQLRequest = "REPACK " + m_OutputEntitiesVectorFile.substr(0, m_OutputEntitiesVectorFile.find("."));
      Entities->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
      EntitiesLayer->SyncToDisk();
      IndexToRemove.clear();
      NewFeature = OGRFeature::CreateFeature(FeatureDefn);
      NewFeature->SetField("IDPlot", IDPlotAllSmallTab[i]);
      NewFeature->SetField("ID", IDNewTab[i].c_str());
      NewFeature->SetField("IDPlot2", IDPlot2NewTab[i]);
      NewFeature->SetField("ID2", ID2NewTab[i].c_str());
      NewFeature->SetField("IDPlot3", IDPlot3NewTab[i]);
      NewFeature->SetField("ID3", ID3NewTab[i].c_str());
      NewFeature->SetGeometry(NewGeometry);
      EntitiesLayer->CreateFeature(NewFeature);
      EntitiesLayer->SyncToDisk();
    }
    OGRDataSource::DestroyDataSource(Entities);
    Entities = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
    EntitiesLayer = Entities->GetLayer(0);
    EntitiesLayer->ResetReading();
    while ((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
    {
      if (std::find(IDPlotAllSmallTab.begin(), IDPlotAllSmallTab.end(), EntitiesFeature->GetFieldAsInteger("IDPlot2")) != IDPlotAllSmallTab.end())
      {
        pos = std::find(IDPlotAllSmallTab.begin(), IDPlotAllSmallTab.end(), EntitiesFeature->GetFieldAsInteger("IDPlot2")) - IDPlotAllSmallTab.begin();
        EntitiesFeature->SetField("ID2", IDNewTab[pos].c_str());
        EntitiesFeature->SetField("IDPlot3", IDPlot2NewTab[pos]);
        EntitiesFeature->SetField("ID3", ID2NewTab[pos].c_str());
        EntitiesLayer->SetFeature(EntitiesFeature);
      }
    }
    EntitiesLayer->SyncToDisk();
    EntitiesLayer = Entities->GetLayer(0);
    EntitiesLayer->ResetReading();
    while ((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
    {
      if (std::find(IDPlotAllSmallTab.begin(), IDPlotAllSmallTab.end(), EntitiesFeature->GetFieldAsInteger("IDPlot3")) != IDPlotAllSmallTab.end())
      {
        pos = std::find(IDPlotAllSmallTab.begin(), IDPlotAllSmallTab.end(), EntitiesFeature->GetFieldAsInteger("IDPlot3")) - IDPlotAllSmallTab.begin();
        EntitiesFeature->SetField("ID3", IDNewTab[pos].c_str());
        EntitiesLayer->SetFeature(EntitiesFeature);
      }
    }
    EntitiesLayer->SyncToDisk();
    OGRDataSource::DestroyDataSource(Entities);
  }

  OGRDataSource::DestroyDataSource(DataSource);

  // ========================================================================
  // Re-aggregating sub-entities that are below the seed but there are larger
  // entity for it to be merged with that belong to the same plot
  // ========================================================================

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
  m_SQLRequest = "SELECT ROWID as RID, IDPlot FROM " + m_OutputEntitiesVectorFile.substr(0, m_OutputEntitiesVectorFile.find(".")) + " WHERE ST_Area(geometry) <= " + std::to_string(m_MinEntSize);
  SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());

  for(int i = 0; i < SQLLayer->GetFeatureCount(); i++)
  {
    if (std::find(IDPlotLargeTab.begin(), IDPlotLargeTab.end(), SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot")) != IDPlotLargeTab.end() and
        std::find(IDPlotAllSmallTab.begin(), IDPlotAllSmallTab.end(), SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot")) == IDPlotAllSmallTab.end())
    {
      FIDSmallTab.push_back(SQLLayer->GetFeature(i)->GetFieldAsInteger("RID"));
    }
  }

  DataSource->ReleaseResultSet(SQLLayer);
  OGRDataSource::DestroyDataSource(DataSource);

  IDOldTab.clear();
  IDNewTab.clear();
  ID2NewTab.clear();
  IDPlotNewTab.clear();
  IDPlot2NewTab.clear();

  FeaturesToDelete = FIDSmallTab;
  Entities = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
  EntitiesLayer = Entities->GetLayer(0);

  while(FIDSmallTab.size() != 0)
  {
    Beginning = FIDSmallTab.size();
    for (int i = 0; i < FIDSmallTab.size(); i++)
    {
      DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
      EntitiesFeature = EntitiesLayer->GetFeature(FIDSmallTab[i]);
      ID = EntitiesFeature->GetFieldAsString("ID");
      EntitiesGeometry = EntitiesFeature->GetGeometryRef();
      IDPlot = EntitiesFeature->GetFieldAsInteger("IDPlot");
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + m_OutputEntitiesVectorFile.substr(0, m_OutputEntitiesVectorFile.find(".")) + " WHERE "
          "IDPlot == " + std::to_string(IDPlot) + " AND ID3 != " + m_QMark + ID + m_QMark + " AND "
          "ST_Area(geometry) > " + std::to_string(m_MinEntSize) + " AND "
          "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + m_OutputEntitiesVectorFile.substr(0, m_OutputEntitiesVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDSmallTab[i])+ "))) > 0 ORDER BY ST_Area(geometry) DESC";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        SQLFeature = SQLLayer->GetFeature(0);
        FID = SQLFeature->GetFieldAsInteger("RID");
        IDOldTab.push_back(EntitiesFeature->GetFieldAsString("ID"));
        IDNewTab.push_back(SQLFeature->GetFieldAsString("ID"));
        ID2NewTab.push_back(SQLFeature->GetFieldAsString("ID2"));
        IDPlotNewTab.push_back(SQLFeature->GetFieldAsInteger("IDPlot"));
        IDPlot2NewTab.push_back(SQLFeature->GetFieldAsInteger("IDPlot2"));
        AggregatingFeature = EntitiesLayer->GetFeature(FID);
        AggregatingGeometry = AggregatingFeature->GetGeometryRef();
        AggregatingGeometry = AggregatingGeometry->Union(EntitiesGeometry);
        AggregatingFeature->SetGeometry(AggregatingGeometry);
        EntitiesLayer->SetFeature(AggregatingFeature);
        EntitiesLayer->SyncToDisk();
        IndexToRemove.push_back(FIDSmallTab[i]);
        EntitiesFeature = nullptr;
      }
      DataSource->ReleaseResultSet(SQLLayer);
      OGRDataSource::DestroyDataSource(DataSource);
    }
    for (int k = 0; k < IndexToRemove.size(); k++)
    {
      if (std::find(FIDSmallTab.begin(), FIDSmallTab.end(), IndexToRemove[k]) != FIDSmallTab.end())
      {
        pos = std::find(FIDSmallTab.begin(), FIDSmallTab.end(), IndexToRemove[k]) - FIDSmallTab.begin();
        FIDSmallTab.erase(FIDSmallTab.begin()+pos);
      }
    }
    IndexToRemove.clear();
    End = FIDSmallTab.size();
    if (Beginning == End)
    {
      break;
    }
  }

  for (int i = 0; i < FeaturesToDelete.size(); i++)
  {
    if (std::find(FIDSmallTab.begin(), FIDSmallTab.end(), FeaturesToDelete[i]) == FIDSmallTab.end())
    {
      EntitiesLayer->DeleteFeature(FeaturesToDelete[i]);
    }
  }

  m_SQLRequest = "REPACK " + m_OutputEntitiesVectorFile.substr(0, m_OutputEntitiesVectorFile.find("."));
  Entities->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
  EntitiesLayer->SyncToDisk();
  EntitiesLayer = Entities->GetLayer(0);
  EntitiesLayer->ResetReading();

  while ((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
  {
    ID2 = EntitiesFeature->GetFieldAsString("ID2");
    if (std::find(IDOldTab.begin(), IDOldTab.end(), ID2.c_str()) != IDOldTab.end())
    {
      pos = std::find(IDOldTab.begin(), IDOldTab.end(), ID2.c_str()) - IDOldTab.begin();
      EntitiesFeature->SetField("ID2", IDNewTab[pos].c_str());
      EntitiesFeature->SetField("IDPlot3", IDPlot2NewTab[pos]);
      EntitiesFeature->SetField("ID3", ID2NewTab[pos].c_str());
      EntitiesLayer->SetFeature(EntitiesFeature);
    }
  }

  EntitiesLayer->SyncToDisk();
  EntitiesLayer = Entities->GetLayer(0);
  EntitiesLayer->ResetReading();

  while ((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
  {
    ID3 = EntitiesFeature->GetFieldAsString("ID3");
    if (std::find(IDOldTab.begin(), IDOldTab.end(), ID3.c_str()) != IDOldTab.end())
    {
      pos = std::find(IDOldTab.begin(), IDOldTab.end(), ID3.c_str()) - IDOldTab.begin();
      EntitiesFeature->SetField("ID3", IDNewTab[pos].c_str());
      EntitiesLayer->SetFeature(EntitiesFeature);
    }
  }

  OGRDataSource::DestroyDataSource(Entities);

  IDNewTab.clear();
  ID2NewTab.clear();
  IDPlotNewTab.clear();
  IDPlot2NewTab.clear();
  IDOldTab.clear();
  ID2OldTab.clear();
  IDPlotOldTab.clear();
  IDPlot2OldTab.clear();

  // =====================================================================
  // Union with the original parcel vector layer
  // =====================================================================

  openfluid::utils::GrassGISProxy GRASS("/tmp/bvservice-grass","temp");
  GRASS.setOutputFile("/tmp/bvservice-grass/procesvectordata.out");
  GRASS.setErrorFile("/tmp/bvservice-grass/processvectordata.err");

  GRASS.appendTask("v.in.ogr", {{"input", QString::fromStdString(getOutputVectorPath(m_OutputPlotsVectorFile))},
                                {"output", "plots"},
                                {"snap", QString::fromStdString(m_SnapDistance)}},
                   {"--o"});

  GRASS.appendTask("v.in.ogr", {{"input", QString::fromStdString(getOutputVectorPath(m_OutputEntitiesVectorFile))},
                                {"output", "entities"},
                                {"snap", QString::fromStdString(m_SnapDistance)}},
                   {"--o"});

  GRASS.appendTask("v.overlay", {{"ainput", "plots"},
                                 {"binput", "entities"},
                                 {"operator", "or"},
                                 {"output", "union"},
                                 {"snap", QString::fromStdString(m_SnapDistance)}},
                   {"--o"});

  FieldNamesTab.clear();
  FIDTable.clear();

  FieldNamesTab = {m_IDFieldName,"IDPlot","ID","IDPlot2","ID2","IDPlot3","ID3"};

  GRASS.appendTask("v.db.renamecolumn", {{"map", "union"},
                                         {"column", QString::fromStdString("a_" + FieldNamesTab[0] + "," + FieldNamesTab[0])}});

  for (int i = 1; i < FieldNamesTab.size(); i++)
  {
    GRASS.appendTask("v.db.renamecolumn", {{"map", "union"},
                                           {"column", QString::fromStdString("b_" + FieldNamesTab[i] + "," + FieldNamesTab[i])}});
  }

  GRASS.appendTask("v.db.dropcolumn", {{"map", "union"},
                                       {"column", QString::fromStdString("a_cat,b_cat")}});

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str());
  }

  GRASS.appendTask("v.out.ogr", {{"input", "union"},
                                 {"type", "area"},
                                 {"output", QString::fromStdString(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile))}},
                   {"-s"});

  if (GRASS.runJob() != 0)
  {
    exit(-1);
  }

  // =======================================================================================================
  // Re-attribution of the sliver entities resulted from the union between original parcel vector and the
  // entities vector (i.e., the entities whose original parcel ID attribute is not the same as their IDPlot)
  // =======================================================================================================

  Union = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionFeature = nullptr;
  UnionGeometry = nullptr;
  UnionLayer->ResetReading();

  NewFeature = nullptr;
  NewGeometry = nullptr;

  while ((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    IDField = UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str());
    if (!IDField)
    {
      UnionGeometry = UnionFeature->GetGeometryRef();
      if(NewGeometry == nullptr)
      {
        NewGeometry = UnionGeometry;
      }
      else
      {
        NewGeometry = UnionGeometry->Union(NewGeometry);
      }
      UnionLayer->DeleteFeature(UnionFeature->GetFID());
    }
  }

  m_SQLRequest = "REPACK " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find("."));
  Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
  OGRDataSource::DestroyDataSource(Union);

  Union = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);

  FeatureDefn = UnionLayer->GetLayerDefn();
  NewFeature = OGRFeature::CreateFeature(FeatureDefn);

  for (int i = 0; i < NewFeature->GetFieldCount(); i++)
  {
    if (NewFeature->GetFieldDefnRef(i)->GetType() == OFTString)
    {
      NewFeature->SetField(i, "0N1");
    }
    else
    {
      NewFeature->SetField(i, 0);
    }
  }

  NewFeature->SetGeometry(NewGeometry);
  UnionLayer->CreateFeature(NewFeature);
  UnionLayer->SyncToDisk();
  OGRDataSource::DestroyDataSource(Union);

  Union = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  UnionFeature = nullptr;
  UnionGeometry = nullptr;

  IDPlotTab.clear();
  MissingIDTab.clear();
  FIDSliverTab.clear();

  N = 0;

  while ((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    if (std::find(IDPlotTab.begin(), IDPlotTab.end(), UnionFeature->GetFieldAsInteger("IDPlot")) == IDPlotTab.end())
    {
      IDPlotTab.push_back(UnionFeature->GetFieldAsInteger("IDPlot"));
    }
  }

  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  UnionFeature = nullptr;

  while ((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    if (std::find(IDPlotTab.begin(), IDPlotTab.end(), UnionFeature->GetFieldAsInteger("IDPlot")) == IDPlotTab.end() and
        std::find(MissingIDTab.begin(), MissingIDTab.end(), UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) == MissingIDTab.end())
    {
      MissingIDTab.push_back(UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str()));
    }
  }

  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  UnionFeature = nullptr;

  while ((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    if (UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str()) != UnionFeature->GetFieldAsInteger("IDPlot") and
        std::find(MissingIDTab.begin(), MissingIDTab.end(), UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) == MissingIDTab.end())
    {
      FIDSliverTab.push_back(UnionFeature->GetFID());
    }
  }

  OGRDataSource::DestroyDataSource(Union);

  // ======================================================================================
  // Regrouping entities that have at least one other feature that belong to the same plot,
  // have the same original parcel ID as the entities in question and whose
  // original parcel ID attribute is the same as IDPlot
  // ======================================================================================

  IndexToRemove.clear();

  while (FIDSliverTab.size() != 0)
  {
    Union = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    for (int i = 0; i < FIDSliverTab.size(); i++)
    {
      UnionFeature = UnionLayer->GetFeature(FIDSliverTab[i]);
      UnionGeometry = UnionFeature->GetGeometryRef();
      DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
      AggregatingFeature = nullptr;
      AggregatingGeometry = nullptr;
      m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE " + m_IDFieldName + " == IDPlot "
          "AND " + m_IDFieldName + " == " + std::to_string(UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) + "  AND "
          "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDSliverTab[i]) + "))) > 0 "
          "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDSliverTab[i]) + "))) DESC";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() > 1)
      {
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE " + m_IDFieldName + " == "
            "IDPlot AND " + m_IDFieldName + " == " + std::to_string(UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) + " AND "
            "ID2 == " + m_QMark + UnionFeature->GetFieldAsString("ID") + m_QMark + " AND "
            "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDSliverTab[i]) + "))) > 0 "
            "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDSliverTab[i]) + "))) DESC";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() == 0)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE " + m_IDFieldName + " == IDPlot "
              "AND " + m_IDFieldName + " == " + std::to_string(UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) + " AND "
              "ID == " + m_QMark + UnionFeature->GetFieldAsString("ID2") + m_QMark + " AND "
              "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDSliverTab[i]) + "))) > 0 "
              "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDSliverTab[i]) + "))) DESC";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() == 0)
          {
            DataSource->ReleaseResultSet(SQLLayer);
            m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE " + m_IDFieldName + " == IDPlot "
                "AND " + m_IDFieldName + " == " + std::to_string(UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) + " AND "
                "ID2 == " + m_QMark + UnionFeature->GetFieldAsString("ID2") + m_QMark + " AND "
                "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDSliverTab[i]) + "))) > 0 "
                "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDSliverTab[i]) + "))) DESC";
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            if (SQLLayer->GetFeatureCount() == 0)
            {
              DataSource->ReleaseResultSet(SQLLayer);
              m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE " + m_IDFieldName + " == IDPlot "
                  "AND " + m_IDFieldName + " == " + std::to_string(UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) + " AND "
                  "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDSliverTab[i]) + "))) > 0 "
                  "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDSliverTab[i]) + "))) DESC";
              SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
              FID = SQLLayer->GetFeature(0)->GetFieldAsInteger("RID");
              DataSource->ReleaseResultSet(SQLLayer);
              AggregatingFeature = UnionLayer->GetFeature(FID);
              AggregatingGeometry = AggregatingFeature->GetGeometryRef();
              AggregatingFeature->SetGeometry(AggregatingGeometry->Union(UnionGeometry));
              UnionLayer->SetFeature(AggregatingFeature);
              UnionLayer->SyncToDisk();
              IndexToRemove.push_back(FIDSliverTab[i]);
            }
            else
            {
              FID = SQLLayer->GetFeature(0)->GetFieldAsInteger("RID");
              DataSource->ReleaseResultSet(SQLLayer);
              AggregatingFeature = UnionLayer->GetFeature(FID);
              AggregatingGeometry = AggregatingFeature->GetGeometryRef();
              AggregatingFeature->SetGeometry(AggregatingGeometry->Union(UnionGeometry));
              UnionLayer->SetFeature(AggregatingFeature);
              UnionLayer->SyncToDisk();
              IndexToRemove.push_back(FIDSliverTab[i]);
            }
          }
          else
          {
            FID = SQLLayer->GetFeature(0)->GetFieldAsInteger("RID");
            DataSource->ReleaseResultSet(SQLLayer);
            AggregatingFeature = UnionLayer->GetFeature(FID);
            AggregatingGeometry = AggregatingFeature->GetGeometryRef();
            AggregatingFeature->SetGeometry(AggregatingGeometry->Union(UnionGeometry));
            UnionLayer->SetFeature(AggregatingFeature);
            UnionLayer->SyncToDisk();
            IndexToRemove.push_back(FIDSliverTab[i]);
          }
        }
        else
        {
          FID = SQLLayer->GetFeature(0)->GetFieldAsInteger("RID");
          DataSource->ReleaseResultSet(SQLLayer);
          AggregatingFeature = UnionLayer->GetFeature(FID);
          AggregatingGeometry = AggregatingFeature->GetGeometryRef();
          AggregatingFeature->SetGeometry(AggregatingGeometry->Union(UnionGeometry));
          UnionLayer->SetFeature(AggregatingFeature);
          UnionLayer->SyncToDisk();
          IndexToRemove.push_back(FIDSliverTab[i]);
        }
      }
      else if (SQLLayer->GetFeatureCount() == 1)
      {
        FID = SQLLayer->GetFeature(0)->GetFieldAsInteger("RID");
        DataSource->ReleaseResultSet(SQLLayer);
        AggregatingFeature = UnionLayer->GetFeature(FID);
        AggregatingGeometry = AggregatingFeature->GetGeometryRef();
        AggregatingFeature->SetGeometry(AggregatingGeometry->Union(UnionGeometry));
        UnionLayer->SetFeature(AggregatingFeature);
        UnionLayer->SyncToDisk();
        IndexToRemove.push_back(FIDSliverTab[i]);
      }
      else
      {
        DataSource->ReleaseResultSet(SQLLayer);
      }
      OGRDataSource::DestroyDataSource(DataSource);
    }
    if (IndexToRemove.size() != 0)
    {
      for (int i = 0; i < IndexToRemove.size(); i++)
      {
        UnionLayer->DeleteFeature(IndexToRemove[i]);
      }
    }
    IndexToRemove.clear();
    m_SQLRequest = "REPACK " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find("."));
    Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
    OGRDataSource::DestroyDataSource(Union);
    Union = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    UnionLayer->ResetReading();
    UnionFeature = nullptr;
    FIDSliverTab.clear();
    while ((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
    {
      if (UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str()) != UnionFeature->GetFieldAsInteger("IDPlot") and std::find(MissingIDTab.begin(), MissingIDTab.end(), UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) == MissingIDTab.end())
      {
        FIDSliverTab.push_back(UnionFeature->GetFID());
      }
    }
    OGRDataSource::DestroyDataSource(Union);
  }

  // ==================================================================================
  // Regrouping entities that do not have any entities that they can be aggregated with
  // ==================================================================================

  Union = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  FeatureDefn = UnionLayer->GetLayerDefn();
  NewFeature = nullptr;

  if (MissingIDTab.size() != 0)
  {
    for (int i = 0; i < MissingIDTab.size(); i++)
    {
      SQLFeature = nullptr;
      SQLGeometry = nullptr;
      NewGeometry = nullptr;
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE " + m_IDFieldName + " == " + std::to_string(MissingIDTab[i]) + " ORDER BY ST_Area(geometry) DESC";
      SQLLayer = Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      IDField = SQLLayer->GetFeature(0)->GetFieldAsInteger(m_IDFieldName.c_str());
      IDPlot2 = SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot2");
      ID2 = SQLLayer->GetFeature(0)->GetFieldAsString("ID2");
      IDPlot3 = SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot3");
      ID3 = SQLLayer->GetFeature(0)->GetFieldAsString("ID3");
      IndexToRemove.clear();
      for(int j = 0; j < SQLLayer->GetFeatureCount(); j++)
      {
        SQLFeature = SQLLayer->GetFeature(j);
        SQLGeometry = SQLFeature->GetGeometryRef();
        IndexToRemove.push_back(SQLFeature->GetFieldAsInteger("RID"));
        if (NewGeometry == nullptr)
        {
          NewGeometry = SQLGeometry;
        }
        else
        {
          NewGeometry = SQLGeometry->Union(NewGeometry);
        }
      }
      Union->ReleaseResultSet(SQLLayer);
      IDNew = std::to_string(IDField) + "N1";
      NewFeature = OGRFeature::CreateFeature(FeatureDefn);
      NewFeature->SetField("OFLD_ID", IDField);
      NewFeature->SetField("IDPlot", IDField);
      NewFeature->SetField("ID", IDNew.c_str());
      NewFeature->SetField("IDPlot2", IDPlot2);
      NewFeature->SetField("ID2", ID2.c_str());
      NewFeature->SetField("IDPlot3", IDPlot3);
      NewFeature->SetField("ID3", ID3.c_str());
      NewFeature->SetGeometry(NewGeometry);
      UnionLayer->CreateFeature(NewFeature);
      UnionLayer->SyncToDisk();
      for (int k = 0; k < IndexToRemove.size(); k++)
      {
        UnionLayer->DeleteFeature(IndexToRemove[k]);
      }
      m_SQLRequest = "REPACK " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find("."));
      Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
      Union->SyncToDisk();
    }
  }

  OGRDataSource::DestroyDataSource(Union);

  IDPlotTab.clear();
  MissingIDTab.clear();
  FIDSliverTab.clear();
  IndexToRemove.clear();

  // ====================================================================================
  // Regrouping of those features that have the same attributes and share a common border
  // ====================================================================================

  Union = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  UnionFeature = nullptr;
  UnionGeometry = nullptr;

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
  SQLLayer = nullptr;
  SQLFeature = nullptr;
  SQLGeometry = nullptr;

  IndexToRemove.clear();

  while ((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    UnionGeometry = UnionFeature->GetGeometryRef();
    NewGeometry = nullptr;
    FID = UnionFeature->GetFID();
    IDField = UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str());
    IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
    ID = UnionFeature->GetFieldAsString("ID");
    if (std::find(IndexToRemove.begin(), IndexToRemove.end(), FID) == IndexToRemove.end())
    {
      m_SQLRequest = "SELECT ROWID as RID, geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE " + m_IDFieldName + " == " + std::to_string(IDField) + " AND "
          "IDPlot == " + std::to_string(IDPlot) + " AND ID == " + m_QMark + ID + m_QMark + " AND "
          "RID != " + std::to_string(FID) + " ORDER BY ST_Area(geometry)";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        for (int i = 0; i < SQLLayer->GetFeatureCount(); i++)
        {
          SQLFeature = SQLLayer->GetFeature(i);
          SQLGeometry = SQLFeature->GetGeometryRef();
          IndexToRemove.push_back(SQLFeature->GetFieldAsInteger("RID"));
          if (NewGeometry == nullptr)
          {
            NewGeometry = UnionGeometry->Union(SQLGeometry);
          }
          else
          {
            NewGeometry = SQLGeometry->Union(NewGeometry);
          }
        }
      }
      DataSource->ReleaseResultSet(SQLLayer);
      if (NewGeometry != nullptr)
      {
        UnionFeature->SetGeometry(NewGeometry);
        UnionLayer->SetFeature(UnionFeature);
      }
    }
  }

  UnionLayer->SyncToDisk();
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  if (IndexToRemove.size() != 0)
  {
    for (int i = 0; i < IndexToRemove.size(); i++)
    {
      UnionLayer->DeleteFeature(IndexToRemove[i]);
    }
  }

  IndexToRemove.clear();
  m_SQLRequest = "REPACK " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find("."));
  Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
  OGRDataSource::DestroyDataSource(Union);

  // ===============================================================================
  // Checking for missing connections between entities
  // ===============================================================================

  Union = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  DataSource  = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );

  FIDConnectionProblemTab.clear();

  while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    ID = UnionFeature->GetFieldAsString("ID");
    ID2 = UnionFeature->GetFieldAsString("ID2");
    if (ID != "0N1")
    {
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID2 + m_QMark + " AND "
          "ST_Touches(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + "))";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() == 0)
      {
        FIDConnectionProblemTab.push_back(UnionFeature->GetFID());
      }
      DataSource->ReleaseResultSet(SQLLayer);
    }
  }

  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  if (FIDConnectionProblemTab.size() != 0)
  {
    while (FIDConnectionProblemTab.size() != 0)
    {
      Beginning = FIDConnectionProblemTab.size();
      IndexToRemove.clear();
      DataSource  = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
      for (int i = 0; i < FIDConnectionProblemTab.size(); i++)
      {
        UnionFeature = UnionLayer->GetFeature(FIDConnectionProblemTab[i]);
        IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
        ID = UnionFeature->GetFieldAsString("ID");
        ID2 = UnionFeature->GetFieldAsString("ID2");
        m_SQLRequest = "SELECT IDPlot, ID, IDPlot2, ID2 FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID2 == " + m_QMark + ID2 + m_QMark + " AND "
            "IDPlot != " + std::to_string(IDPlot) + " AND "
            "ST_Touches(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
            "ST_Touches(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")) "
            "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + "))) DESC,"
            " ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID2 + m_QMark + "))) DESC";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          IDPlot2New = SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot");
          ID2New = SQLLayer->GetFeature(0)->GetFieldAsString("ID");
          IDPlot3New = SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot2");
          ID3New = SQLLayer->GetFeature(0)->GetFieldAsString("ID2");
          UnionFeature->SetField("IDPlot2", IDPlot2New);
          UnionFeature->SetField("ID2", ID2New.c_str());
          UnionFeature->SetField("IDPlot3", IDPlot3New);
          UnionFeature->SetField("ID3", ID3New.c_str());
          UnionLayer->SetFeature(UnionFeature);
          IndexToRemove.push_back(UnionFeature->GetFID());
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID2 == " + m_QMark + ID + m_QMark;
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() != 0)
          {
            for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
            {
              FID = SQLLayer->GetFeature(j)->GetFieldAsInteger("RID");
              UnionFeatureToUpdate = UnionLayer->GetFeature(FID);
              UnionFeatureToUpdate->SetField("IDPlot3", IDPlot2New);
              UnionFeatureToUpdate->SetField("ID3", ID2New.c_str());
              UnionLayer->SetFeature(UnionFeatureToUpdate);
            }
          }
          DataSource->ReleaseResultSet(SQLLayer);
        }
        else
        {
          DataSource->ReleaseResultSet(SQLLayer);
        }
      }
      Union->SyncToDisk();
      OGRDataSource::DestroyDataSource(DataSource);
      for (int i = 0; i < IndexToRemove.size(); i++)
      {
        if (std::find(FIDConnectionProblemTab.begin(), FIDConnectionProblemTab.end(), IndexToRemove[i]) != FIDConnectionProblemTab.end())
        {
          pos = std::find(FIDConnectionProblemTab.begin(), FIDConnectionProblemTab.end(), IndexToRemove[i]) - FIDConnectionProblemTab.begin();
          FIDConnectionProblemTab.erase(FIDConnectionProblemTab.begin() + pos);
        }
      }
      End = FIDConnectionProblemTab.size();
      if (Beginning == End)
      {
        break;
      }
    }
  }

  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  if (FIDConnectionProblemTab.size() != 0)
  {
    while (FIDConnectionProblemTab.size() != 0)
    {
      Beginning = FIDConnectionProblemTab.size();
      IndexToRemove.clear();
      DataSource  = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
      for (int i = 0; i < FIDConnectionProblemTab.size(); i++)
      {
        UnionFeature = UnionLayer->GetFeature(FIDConnectionProblemTab[i]);
        IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
        ID = UnionFeature->GetFieldAsString("ID");
        ID2 = UnionFeature->GetFieldAsString("ID2");
        m_SQLRequest = "SELECT IDPlot, ID, IDPlot2, ID2 FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE IDPlot != " + std::to_string(IDPlot) + " AND "
            "ST_Touches(geometry, (SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
            "ID2 != " + m_QMark + ID + m_QMark + " AND ID3 != " + m_QMark + ID + m_QMark + " ORDER BY "
            "ST_Distance(geometry, (SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")) ASC, "
            "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + "))) DESC";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          IDPlot2New = SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot");
          ID2New = SQLLayer->GetFeature(0)->GetFieldAsString("ID");
          IDPlot3New = SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot2");
          ID3New = SQLLayer->GetFeature(0)->GetFieldAsString("ID2");
          UnionFeature->SetField("IDPlot2", IDPlot2New);
          UnionFeature->SetField("ID2", ID2New.c_str());
          UnionFeature->SetField("IDPlot3", IDPlot3New);
          UnionFeature->SetField("ID3", ID3New.c_str());
          UnionLayer->SetFeature(UnionFeature);
          IndexToRemove.push_back(UnionFeature->GetFID());
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID2 == " + m_QMark + ID + m_QMark;
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() != 0)
          {
            for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
            {
              FID = SQLLayer->GetFeature(j)->GetFieldAsInteger("RID");
              UnionFeatureToUpdate = UnionLayer->GetFeature(FID);
              UnionFeatureToUpdate->SetField("IDPlot3", IDPlot2New);
              UnionFeatureToUpdate->SetField("ID3", ID2New.c_str());
              UnionLayer->SetFeature(UnionFeatureToUpdate);
            }
          }
          DataSource->ReleaseResultSet(SQLLayer);
        }
        else
        {
          DataSource->ReleaseResultSet(SQLLayer);
        }
      }
      Union->SyncToDisk();
      for (int i = 0; i < IndexToRemove.size(); i++)
      {
        if (std::find(FIDConnectionProblemTab.begin(), FIDConnectionProblemTab.end(), IndexToRemove[i]) != FIDConnectionProblemTab.end())
        {
          pos = std::find(FIDConnectionProblemTab.begin(), FIDConnectionProblemTab.end(), IndexToRemove[i]) - FIDConnectionProblemTab.begin();
          FIDConnectionProblemTab.erase(FIDConnectionProblemTab.begin() + pos);
        }
      }
      OGRDataSource::DestroyDataSource(DataSource);
      End = FIDConnectionProblemTab.size();
      if (Beginning == End)
      {
        break;
      }
    }
  }

  FIDConnectionProblemTab.clear();
  OGRDataSource::DestroyDataSource(Union);

  // =========================================================================================
  // Finding and aggregating entities that are "locked" inside a parcel
  // (could happen as a result of sliver re-attribution)
  // =========================================================================================

  Union = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  DataSource  = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );

  FIDConnectionProblemTab.clear();

  while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
    ID = UnionFeature->GetFieldAsString("ID");
    ID2 = UnionFeature->GetFieldAsString("ID2");
    if (ID != "0N1")
    {
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE "
          "ST_Touches(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND IDPlot != " + std::to_string(IDPlot);
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() == 0)
      {
        FIDConnectionProblemTab.push_back(UnionFeature->GetFID());
      }
      DataSource->ReleaseResultSet(SQLLayer);
    }
  }

  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  IndexToRemove.clear();

  if (FIDConnectionProblemTab.size() != 0)
  {
    DataSource  = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
    for (int i = 0; i < FIDConnectionProblemTab.size(); i++)
    {
      UnionFeature = UnionLayer->GetFeature(FIDConnectionProblemTab[i]);
      IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
      ID = UnionFeature->GetFieldAsString("ID");
      ID2 = UnionFeature->GetFieldAsString("ID2");
      m_SQLRequest = "SELECT ST_Union(geometry, (SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDConnectionProblemTab[i]) + ")), ROWID AS RID, IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3 FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE"
          " IDPlot == " + std::to_string(IDPlot) + " AND "
          "ST_Touches(geometry, (SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDConnectionProblemTab[i]) + ")) AND "
          "ID2 != " + m_QMark + ID + m_QMark + " AND ID3 != " + m_QMark + ID + m_QMark + " ORDER BY "
          "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FIDConnectionProblemTab[i]) + "))) DESC";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        IDNew = SQLLayer->GetFeature(0)->GetFieldAsString("ID");
        IDPlot2New = SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot2");
        ID2New = SQLLayer->GetFeature(0)->GetFieldAsString("ID2");
        IDPlot3New = SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot3");
        ID3New = SQLLayer->GetFeature(0)->GetFieldAsString("ID3");
        UnionGeometry = UnionFeature->GetGeometryRef();
        NewFeature = OGRFeature::CreateFeature(UnionLayer->GetLayerDefn());
        NewFeature->SetField(m_IDFieldName.c_str(), IDPlot);
        NewFeature->SetField("IDPlot", IDPlot);
        NewFeature->SetField("ID", IDNew.c_str());
        NewFeature->SetField("IDPlot2", IDPlot2New);
        NewFeature->SetField("ID2", ID2New.c_str());
        NewFeature->SetField("IDPlot3", IDPlot3New);
        NewFeature->SetField("ID3", ID3New.c_str());
        NewFeature->SetGeometry(SQLLayer->GetFeature(0)->GetGeometryRef());
        UnionLayer->CreateFeature(NewFeature);
        IndexToRemove.push_back(SQLLayer->GetFeature(0)->GetFieldAsInteger("RID"));
        IndexToRemove.push_back(FIDConnectionProblemTab[i]);
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID2 == " + m_QMark + ID + m_QMark;
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
          {
            FID = SQLLayer->GetFeature(j)->GetFieldAsInteger("RID");
            UnionFeatureToUpdate = UnionLayer->GetFeature(FID);
            UnionFeatureToUpdate->SetField("ID2", IDNew.c_str());
            UnionFeatureToUpdate->SetField("IDPlot3", IDPlot2New);
            UnionFeatureToUpdate->SetField("ID3", ID2New.c_str());
            UnionLayer->SetFeature(UnionFeatureToUpdate);
          }
        }
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID3 == " + m_QMark + ID + m_QMark;
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
          {
            FID = SQLLayer->GetFeature(j)->GetFieldAsInteger("RID");
            UnionFeatureToUpdate = UnionLayer->GetFeature(FID);
            UnionFeatureToUpdate->SetField("ID3", IDNew.c_str());
            UnionLayer->SetFeature(UnionFeatureToUpdate);
          }
        }
        DataSource->ReleaseResultSet(SQLLayer);
      }
      else
      {
        DataSource->ReleaseResultSet(SQLLayer);
      }
    }
    Union->SyncToDisk();
    OGRDataSource::DestroyDataSource(DataSource);
  }

  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  if (IndexToRemove.size() != 0)
  {
    for (int i = 0; i < IndexToRemove.size(); i++)
    {
      UnionLayer->DeleteFeature(IndexToRemove[i]);
    }
  }

  IndexToRemove.clear();
  m_SQLRequest = "REPACK " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find("."));
  Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
  Union->SyncToDisk();

  OGRDataSource::DestroyDataSource(Union);
  OGRDataSource::DestroyDataSource(DataSource);

  // ==========================================================================================
  // Finding and aggregating entities that belong to the same parcel
  // while draining to the same entity and forming a common uninterrupted border
  // (once merged together) with that entity
  // ==========================================================================================

  N = 1;

  while (N != 0)
  {
    IndexToRemove.clear();
    IDOldTab.clear();
    IDNewTab.clear();
    IDPlotNewTab.clear();
    Union = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    UnionLayer->ResetReading();
    DataSource  = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
    NewGeometry = nullptr;
    while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
    {
      FID = UnionFeature->GetFID();
      IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
      ID = UnionFeature->GetFieldAsString("ID");
      ID2 = UnionFeature->GetFieldAsString("ID2");
      if (ID != "0N1" and std::find(IndexToRemove.begin(), IndexToRemove.end(), FID) == IndexToRemove.end())
      {
        m_SQLRequest = "SELECT geometry, ROWID as RID, * FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE "
            "ST_Touches(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FID) + ")) AND "
            "IDPlot == " + std::to_string(IDPlot) + " AND ID2 == " + m_QMark + ID2 + m_QMark + " AND "
            "ROWID != " + std::to_string(FID) + " AND "
            "ST_Touches("
            "ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")), "
            "ST_Intersection((SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FID) + "),"
            "(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")))";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          NewGeometry = nullptr;
          for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
          {
            if (std::find(IndexToRemove.begin(), IndexToRemove.end(), SQLLayer->GetFeature(j)->GetFieldAsInteger("RID")) == IndexToRemove.end())
            {
              IndexToRemove.push_back(SQLLayer->GetFeature(j)->GetFieldAsInteger("RID"));
              IDOldTab.push_back(SQLLayer->GetFeature(j)->GetFieldAsString("ID"));
              IDNewTab.push_back(ID.c_str());
              IDPlotNewTab.push_back(IDPlot);
              if (NewGeometry == nullptr)
              {
                NewGeometry = SQLLayer->GetFeature(j)->GetGeometryRef();
              }
              else
              {
                NewGeometry = NewGeometry->Union(SQLLayer->GetFeature(j)->GetGeometryRef());
              }
            }
          }
          if (NewGeometry != nullptr)
          {
            UnionGeometry = UnionFeature->GetGeometryRef();
            UnionGeometry = UnionGeometry->Union(NewGeometry);
            UnionFeature->SetGeometry(UnionGeometry);
            UnionLayer->SetFeature(UnionFeature);
          }
        }
        DataSource->ReleaseResultSet(SQLLayer);
      }
    }
    OGRDataSource::DestroyDataSource(DataSource);
    UnionLayer = Union->GetLayer(0);
    UnionLayer->ResetReading();
    if (IndexToRemove.size() != 0)
    {
      for (int i = 0; i < IndexToRemove.size(); i++)
      {
        UnionLayer->DeleteFeature(IndexToRemove[i]);
      }
    }
    IndexToRemove.clear();
    m_SQLRequest = "REPACK " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find("."));
    Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
    Union->SyncToDisk();
    OGRDataSource::DestroyDataSource(Union);
    Union  = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    DataSource  = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
    for (int i = 0; i < IDOldTab.size(); i++)
    {
      IDOld = IDOldTab[i];
      IDNew = IDNewTab[i];
      IDPlotNew = IDPlotNewTab[i];
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE "
          "ID2 == " + m_QMark + IDOldTab[i] + m_QMark;
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
        {
          FID = SQLLayer->GetFeature(j)->GetFieldAsInteger("RID");
          UnionFeature = UnionLayer->GetFeature(FID);
          UnionFeature->SetField("IDPlot2",IDPlotNew);
          UnionFeature->SetField("ID2",IDNew.c_str());
          UnionLayer->SetFeature(UnionFeature);
        }
      }
      DataSource->ReleaseResultSet(SQLLayer);
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE "
          "ID3 == " + m_QMark + IDOldTab[i] + m_QMark;
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
        {
          FID = SQLLayer->GetFeature(j)->GetFieldAsInteger("RID");
          UnionFeature = UnionLayer->GetFeature(FID);
          UnionFeature->SetField("IDPlot3",IDPlotNew);
          UnionFeature->SetField("ID3",IDNew.c_str());
          UnionLayer->SetFeature(UnionFeature);
        }
      }
      DataSource->ReleaseResultSet(SQLLayer);
    }
    OGRDataSource::DestroyDataSource(Union);
    OGRDataSource::DestroyDataSource(DataSource);
    Union = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    UnionLayer->ResetReading();
    DataSource  = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
    N = 0;
    while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
    {
      FID = UnionFeature->GetFID();
      IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
      ID = UnionFeature->GetFieldAsString("ID");
      ID2 = UnionFeature->GetFieldAsString("ID2");
      if (ID != "0N1")
      {
        m_SQLRequest = "SELECT geometry, ROWID as RID, * FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE "
            "ST_Touches(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FID) + ")) AND "
            "IDPlot == " + std::to_string(IDPlot) + " AND ID2 == " + m_QMark + ID2 + m_QMark + " AND "
            "ROWID != " + std::to_string(FID) + " AND "
            "ST_Touches("
            "ST_Intersection(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")), "
            "ST_Intersection((SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FID) + "),"
            "(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")))";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          N += 1;
        }
        DataSource->ReleaseResultSet(SQLLayer);
      }
    }
    OGRDataSource::DestroyDataSource(DataSource);
    OGRDataSource::DestroyDataSource(Union);
  }

  IndexToRemove.clear();
  IDOldTab.clear();
  IDNewTab.clear();
  IDPlotNewTab.clear();

  // ==========================================================================
  // Checking for entities that have only two other entities as neighbors and
  // one of them belong to the same plot. In this case, the entity in question
  // will be, if possible, aggregated with the entity of the same plot
  // that it shears the border with in order to avoid fracturing the future LNR
  // ==========================================================================

  Union = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  DataSource  = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );

  IndexToRemove.clear();

  while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    FID = UnionFeature->GetFID();
    IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
    ID = UnionFeature->GetFieldAsString("ID");
    ID2 = UnionFeature->GetFieldAsString("ID2");
    if (ID != "0N1" and ID2 != "0N1")
    {
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE "
          "ST_Touches(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
          "ROWID != " + std::to_string(FID);
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() <= 2)
      {
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT ROWID as RID, * FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE "
            "ST_Touches(geometry,(SELECT geometry FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
            "IDPlot == " + std::to_string(IDPlot) + " AND ROWID != " + std::to_string(FID);
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          UnionGeometry = UnionFeature->GetGeometryRef();
          AggregatingFeature = UnionLayer->GetFeature(SQLLayer->GetFeature(0)->GetFieldAsInteger("RID"));
          AggregatingGeometry = AggregatingFeature->GetGeometryRef();
          AggregatingGeometry = AggregatingGeometry->Union(UnionGeometry);
          AggregatingFeature->SetGeometry(AggregatingGeometry);
          UnionLayer->SetFeature(AggregatingFeature);
          IndexToRemove.push_back(FID);
        }
        DataSource->ReleaseResultSet(SQLLayer);
      }
      else
      {
        DataSource->ReleaseResultSet(SQLLayer);
      }
    }
  }

  OGRDataSource::DestroyDataSource(DataSource);

  if (IndexToRemove.size() != 0)
  {
    for (int i = 0; i < IndexToRemove.size(); i++)
    {
      UnionLayer->DeleteFeature(IndexToRemove[i]);
    }
  }

  IndexToRemove.clear();
  m_SQLRequest = "REPACK " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find("."));
  Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
  OGRDataSource::DestroyDataSource(Union);

  // ============================================================================
  // Final regrouping
  // ============================================================================

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputEntitiesGroupedVectorFile).c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputEntitiesGroupedVectorFile).c_str());
  }

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
  m_SQLRequest = "SELECT ST_Union(geometry), ID, ID2, ID3 FROM " + m_OutputPlotsAndEntitiesUnionVectorFile.substr(0, m_OutputPlotsAndEntitiesUnionVectorFile.find(".")) + " GROUP BY ID, ID2, ID3 ORDER BY IDPlot ASC, ID ASC";
  SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
  FileName = m_OutputEntitiesGroupedVectorFile.substr(0, m_OutputEntitiesGroupedVectorFile.find("."));
  DataSource->CopyLayer(SQLLayer, FileName.c_str(), nullptr);
  DataSource->ReleaseResultSet(SQLLayer);
  OGRDataSource::DestroyDataSource(DataSource);

  // =========================================================================================
  // Making sure that ID are consecutive (i.e., that if there are only three entities that
  // belong to the same plot number 5, their IDs will be 5N1, 5N2 and 5N3, not 5N2, 5N9, 5N15)
  // =========================================================================================

  Final = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputEntitiesGroupedVectorFile).c_str(), TRUE );
  FinalLayer = Final->GetLayer(0);
  FinalLayer->ResetReading();

  IDOldTab.clear();
  IDNewTab.clear();

  while ((FinalFeature = FinalLayer->GetNextFeature()) != nullptr)
  {
    FID = FinalFeature->GetFID();
    ID = FinalFeature->GetFieldAsString("ID");
    IDPlot = std::stoi(ID.substr(0, ID.find("N")));
    if (FID == 0)
    {
      if (ID.substr(ID.find("N")+1, ID.length()-1) != "1")
      {
        IDOldTab.push_back(ID.c_str());
        IDNew = std::to_string(IDPlot) + "N1";
        IDNewTab.push_back(IDNew);
        FinalFeature->SetField("ID", IDNew.c_str());
        FinalLayer->SetFeature(FinalFeature);
      }
    }
    else
    {
      IDAbove = FinalLayer->GetFeature(FID-1)->GetFieldAsString("ID");
      IDPlotAbove = std::stoi(IDAbove.substr(0, IDAbove.find("N")));
      if (IDPlot == IDPlotAbove)
      {
        if(std::stoi(ID.substr(ID.find("N")+1, ID.length()-1)) != std::stoi(IDAbove.substr(IDAbove.find("N")+1, IDAbove.length()-1))+1)
        {
          IDOldTab.push_back(ID.c_str());
          IDNew = std::to_string(IDPlot) + "N" + std::to_string(std::stoi(IDAbove.substr(IDAbove.find("N")+1, IDAbove.length()-1))+1);
          IDNewTab.push_back(IDNew.c_str());
          FinalFeature->SetField("ID", IDNew.c_str());
          FinalLayer->SetFeature(FinalFeature);
        }
      }
      else
      {
        if (ID.substr(ID.find("N")+1, ID.length()-1) != "1")
        {
          IDOldTab.push_back(ID.c_str());
          IDNew = std::to_string(IDPlot) + "N1";
          IDNewTab.push_back(IDNew);
          FinalFeature->SetField("ID", IDNew.c_str());
          FinalLayer->SetFeature(FinalFeature);
        }
      }
    }
  }

  FinalLayer = Final->GetLayer(0);
  FinalLayer->ResetReading();

  while ((FinalFeature = FinalLayer->GetNextFeature()) != nullptr)
  {
    ID2 = FinalFeature->GetFieldAsString("ID2");
    ID3 = FinalFeature->GetFieldAsString("ID3");
    if (std::find(IDOldTab.begin(), IDOldTab.end(), ID2.c_str()) != IDOldTab.end())
    {
      pos = std::find(IDOldTab.begin(), IDOldTab.end(), ID2.c_str()) - IDOldTab.begin();
      FinalFeature->SetField("ID2", IDNewTab[pos].c_str());
      FinalLayer->SetFeature(FinalFeature);
    }
    if (std::find(IDOldTab.begin(), IDOldTab.end(), ID3.c_str()) != IDOldTab.end())
    {
      pos = std::find(IDOldTab.begin(), IDOldTab.end(), ID3.c_str()) - IDOldTab.begin();
      FinalFeature->SetField("ID3", IDNewTab[pos].c_str());
      FinalLayer->SetFeature(FinalFeature);
    }
  }

  IDOldTab.clear();
  IDNewTab.clear();
  OGRDataSource::DestroyDataSource(Final);

  // ================================================================================
  // Creating SRF.shp and LNR.shp files
  // ================================================================================

  // ================================================================================
  // Creating linear entities vector file - LNR.shp
  // ================================================================================

  Final = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputEntitiesGroupedVectorFile).c_str(), TRUE );
  FinalLayer = Final->GetLayer(0);
  FinalLayer->ResetReading();
  mp_SRS = Final->GetLayer(0)->GetSpatialRef();

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str());
  }

  LNR = mp_VectorDriver->CreateDataSource( getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str(), nullptr);
  LNRLayer = LNR->CreateLayer( "LNR", mp_SRS, wkbLineString, nullptr );

  FieldNamesTab.clear();
  FieldNamesTab = {"FaceA","FaceB","ID","IDTo"};

  for (int i = 0; i < FieldNamesTab.size(); i++)
  {
    createField(LNRLayer, FieldNamesTab[i].c_str(), OFTString);
  }

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputEntitiesGroupedVectorFile).c_str(), TRUE );
  LNRFeature = OGRFeature::CreateFeature(LNRLayer->GetLayerDefn());
  FacesTab.clear();

  while ((FinalFeature = FinalLayer->GetNextFeature()) != nullptr)
  {
    FID = FinalFeature->GetFID();
    FaceA = FinalFeature->GetFieldAsString("ID");
    IDPlot = std::stoi(FaceA.substr(0, FaceA.find("N")));
    if (IDPlot != 0)
    {
      m_SQLRequest = "SELECT ST_CollectionExtract(ST_Intersection(geometry, (SELECT geometry FROM " + m_OutputEntitiesGroupedVectorFile.substr(0, m_OutputEntitiesGroupedVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FID) + ")),2), ID "
          "FROM " + m_OutputEntitiesGroupedVectorFile.substr(0, m_OutputEntitiesGroupedVectorFile.find(".")) + " WHERE "
          "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + m_OutputEntitiesGroupedVectorFile.substr(0, m_OutputEntitiesGroupedVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FID) + "))) > 0";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        for (int i = 0; i < SQLLayer->GetFeatureCount(); i++)
        {
          FaceB = SQLLayer->GetFeature(i)->GetFieldAsString("ID");
          FaceBA = FaceB + "-" + FaceA;
          if (std::find(FacesTab.begin(), FacesTab.end(), FaceBA.c_str()) == FacesTab.end())
          {
            LNRFeature->SetField(FieldNamesTab[0].c_str(), FaceA.c_str());
            LNRFeature->SetField(FieldNamesTab[1].c_str(), FaceB.c_str());
            NewGeometry = FinalFeature->GetGeometryRef()->Intersection(SQLLayer->GetFeature(i)->GetGeometryRef());
            LNRFeature->SetGeometry(NewGeometry);
            LNRLayer->CreateFeature(LNRFeature);
            LNR->SyncToDisk();
            FaceAB = FaceA + "-" + FaceB;
            FacesTab.push_back(FaceAB.c_str());
          }
        }
      }
      DataSource->ReleaseResultSet(SQLLayer);
    }
  }

  FieldNamesTab.clear();
  FacesTab.clear();
  OGRFeature::DestroyFeature(LNRFeature);
  OGRDataSource::DestroyDataSource(LNR);
  OGRDataSource::DestroyDataSource(DataSource);
  OGRDataSource::DestroyDataSource(Final);

  // ================================================================================
  // Creating surface entities vector file - SRF.shp
  // ================================================================================

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str());
  }

  Final = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputEntitiesGroupedVectorFile).c_str(), TRUE );
  m_SQLRequest = "SELECT ST_Union(geometry), ID, ID2 FROM " + m_OutputEntitiesGroupedVectorFile.substr(0, m_OutputEntitiesGroupedVectorFile.find(".")) + " GROUP BY ID, ID2";
  FileName = m_OutputSurfaceEntitiesVectorFile.substr(0, m_OutputSurfaceEntitiesVectorFile.find("."));
  SQLLayer = Final->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
  Final->CopyLayer(SQLLayer, FileName.c_str(), nullptr);
  Final->ReleaseResultSet(SQLLayer);
  OGRDataSource::DestroyDataSource(Final);

  FieldNamesTab.clear();
  FieldNamesTab = {"IDTo"};

  SRF = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );
  SRFLayer = SRF->GetLayer(0);
  SRFLayer->ResetReading();

  for (int i = 0; i < FieldNamesTab.size(); i++)
  {
    createField(SRFLayer, FieldNamesTab[i].c_str(), OFTString);
  }

  FieldNamesTab.clear();
  OGRDataSource::DestroyDataSource(SRF);

  // ================================================================================
  // Filling IDTo attribute for the linear entities
  // ================================================================================

  LNR = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str(), TRUE );
  LNRLayer = LNR->GetLayer(0);
  LNRLayer->ResetReading();

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );

  while ((LNRFeature = LNRLayer->GetNextFeature()) != nullptr)
  {
    FaceA = LNRFeature->GetFieldAsString("FaceA");
    FaceB = LNRFeature->GetFieldAsString("FaceB");
    m_SQLRequest = "SELECT ID, ID2 FROM " + m_OutputSurfaceEntitiesVectorFile.substr(0, m_OutputSurfaceEntitiesVectorFile.find(".")) + " WHERE "
        "(ID == " + m_QMark + FaceA + m_QMark + " AND ID2 == " + m_QMark + FaceB + m_QMark + ") OR (ID == " + m_QMark + FaceB + m_QMark + " AND ID2 == " + m_QMark + FaceA + m_QMark + ")";
    SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
    if (SQLLayer->GetFeatureCount() != 0)
    {
      if (FaceA == SQLLayer->GetFeature(0)->GetFieldAsString("ID") and FaceB == SQLLayer->GetFeature(0)->GetFieldAsString("ID2"))
      {
        ID = FaceA + "-" + FaceB;
        LNRFeature->SetField("ID", ID.c_str());
        LNRFeature->SetField("IDTo", FaceB.c_str());
        LNRLayer->SetFeature(LNRFeature);
      }
      else if (FaceB == SQLLayer->GetFeature(0)->GetFieldAsString("ID") and FaceA == SQLLayer->GetFeature(0)->GetFieldAsString("ID2"))
      {
        ID = FaceB + "-" + FaceA;
        LNRFeature->SetField("ID", ID.c_str());
        LNRFeature->SetField("IDTo", FaceA.c_str());
        LNRLayer->SetFeature(LNRFeature);
      }
    }
    else
    {
      ID = FaceA + "-" + FaceB;
      IDTo = "None";
      LNRFeature->SetField("ID", ID.c_str());
      LNRFeature->SetField("IDTo", IDTo.c_str());
      LNRLayer->SetFeature(LNRFeature);
    }
    DataSource->ReleaseResultSet(SQLLayer);
  }

  LNRLayer = LNR->GetLayer(0);
  LNRLayer->ResetReading();

  FieldNamesTab.clear();
  FieldNamesTab = {"FaceA", "FaceB"};

  FeatureDefn = LNRLayer->GetLayerDefn();

  if (FieldNamesTab.size() > 0)
  {
    for (int i = 0; i < FieldNamesTab.size(); i++)
    {
      FieldIndex = FeatureDefn->GetFieldIndex(FieldNamesTab[i].c_str());
      LNRLayer->DeleteField(FieldIndex);
    }
  }

  LNR->SyncToDisk();
  FieldNamesTab.clear();

  LNRLayer = LNR->GetLayer(0);

  FieldNamesTab = {"Length"};

  for (int i = 0; i < FieldNamesTab.size(); i++)
  {
    createField(LNRLayer, FieldNamesTab[i].c_str(), OFTReal);
  }

  FieldNamesTab.clear();

  OGRDataSource::DestroyDataSource(LNR);
  OGRDataSource::DestroyDataSource(DataSource);

  // ================================================================================
  // Filling IDTo attribute for the surface entities
  // ================================================================================

  SRF = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );
  SRFLayer = SRF->GetLayer(0);
  SRFLayer->ResetReading();

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );

  N = 0;

  while ((SRFFeature = SRFLayer->GetNextFeature()) != nullptr)
  {
    FID = SRFFeature->GetFID();
    ID = SRFFeature->GetFieldAsString("ID");
    ID2 = SRFFeature->GetFieldAsString("ID2");
    m_SQLRequest = "SELECT * FROM " + m_OutputSurfaceEntitiesVectorFile.substr(0, m_OutputSurfaceEntitiesVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID2 + m_QMark + " AND "
        "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + m_OutputSurfaceEntitiesVectorFile.substr(0, m_OutputSurfaceEntitiesVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FID) + "))) > 0";
    SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
    if (SQLLayer->GetFeatureCount() != 0)
    {
      IDTo = ID + "-" + ID2;
      SRFFeature->SetField("IDTo", IDTo.c_str());
      SRFLayer->SetFeature(SRFFeature);
    }
    else
    {
      SRFFeature->SetField("IDTo", ID2.c_str());
      SRFLayer->SetFeature(SRFFeature);
      N += 1;
    }
    DataSource->ReleaseResultSet(SQLLayer);
  }

  SRFLayer = SRF->GetLayer(0);
  SRFLayer->ResetReading();

  FieldNamesTab.clear();
  FieldNamesTab = {"ID2"};

  FeatureDefn = SRFLayer->GetLayerDefn();

  if (FieldNamesTab.size() > 0)
  {
    for (int i = 0; i < FieldNamesTab.size(); i++)
    {
      FieldIndex = FeatureDefn->GetFieldIndex(FieldNamesTab[i].c_str());
      SRFLayer->DeleteField(FieldIndex);
    }
  }

  SRF->SyncToDisk();
  FieldNamesTab.clear();

  SRFLayer = SRF->GetLayer(0);
  SRFLayer->ResetReading();

  FieldNamesTab.clear();
  FieldNamesTab = {"LandUse", "Surface",  "SlopeMin", "SlopeMax", "SlopeMean"};

  for (int i = 0; i < FieldNamesTab.size(); i++)
  {
    if (FieldNamesTab[i] == "LandUse")
    {
      createField(SRFLayer, FieldNamesTab[i].c_str(), OFTString);
    }
    else
    {
      createField(SRFLayer, FieldNamesTab[i].c_str(), OFTReal);
    }
  }

  FieldNamesTab.clear();
  OGRDataSource::DestroyDataSource(SRF);
  OGRDataSource::DestroyDataSource(DataSource);

}


// ====================================================================
// ====================================================================


void LandProcessor::setSRFParameters()
{

  BLOCK_ENTER(__PRETTY_FUNCTION__)

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int IDPlot, XOff, YOff, XCount, YCount, NPoints;
  double Surface, XOrigin, YOrigin, XMin, XMax, YMin, YMax, PixelWidth, PixelHeight, SlopeMin, SlopeMax, SlopeMean;

  std::string ID, IDTo, LandUseValue;

  std::vector <std::string> FieldNamesTab;
  std::vector <double> XPointsTab, YPointsTab, SlopeValuesTab;

  OGRDataSource *DataSource, *SRF, *Plots;
  OGRLayer *SQLLayer, *SRFLayer, *PlotsLayer;
  OGRFeature *SRFFeature, *PlotsFeature;
  OGRGeometry *SRFGeometry;
  OGRFeatureDefn *FeatureDefn;

  OGRPoint SRFPoint;
  OGRPolygon* SRFPolygon;
  OGRMultiPolygon* SRFMultiPolygon;

  GDALDataset *Slope, *Target;
  GDALRasterBand *SlopeBand, *TargetBand;

  char **Options = nullptr;
  char *SRS_WKT = nullptr;

  float *SlopeRow;
  int *TargetRow;

  // ======================================================================================
  // Checking if necessary files are present in the appropriate folder
  // ======================================================================================

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::setSRFParameters(): " + m_OutputSurfaceEntitiesVectorFile + ": no such file in the output vector directory");
  }

  if (!openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputSlopeRasterFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::setSRFParameters(): " + m_OutputSlopeRasterFile + ": no such file in the output raster directory");
  }

  // ==================================================================================
  // Setting surface parameter for SRF (surface entities)
  // ==================================================================================

  SRF = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );

  SRFLayer = SRF->GetLayer(0);
  SRFLayer->ResetReading();

  while ((SRFFeature = SRFLayer->GetNextFeature()) != nullptr)
  {
    SRFGeometry = SRFFeature->GetGeometryRef();
    if (SRFGeometry->getGeometryType() == 3)
    {
      SRFFeature->SetField("Surface", ((OGRPolygon *) SRFFeature->GetGeometryRef())->get_Area());
      SRFLayer->SetFeature(SRFFeature);
    }
    else if (SRFGeometry->getGeometryType() == 6)
    {
      Surface = 0;
      SRFMultiPolygon = (OGRMultiPolygon*) SRFGeometry;
      for (int i = 0; i < SRFMultiPolygon->getNumGeometries(); i++)
      {
        Surface += ((OGRPolygon*) SRFMultiPolygon->getGeometryRef(i))->get_Area();
      }
      SRFFeature->SetField("Surface", Surface);
      SRFLayer->SetFeature(SRFFeature);
    }
  }

  SRF->SyncToDisk();
  OGRDataSource::DestroyDataSource(SRF);

  // ==================================================================================
  // Setting land use parameter for SRF, if possible (data available for the plot)
  // ==================================================================================

  SRF = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );
  SRFLayer = SRF->GetLayer(0);
  SRFLayer->ResetReading();

  Plots = OGRSFDriverRegistrar::Open( getInputVectorPath(m_InputPlotsFile).c_str(), TRUE );
  PlotsLayer = Plots->GetLayer(0);
  FeatureDefn = PlotsLayer->GetLayerDefn();

  FieldNamesTab.clear();

  for (int i = 0; i < FeatureDefn->GetFieldCount(); i++)
  {
    FieldNamesTab.push_back(PlotsLayer->GetLayerDefn()->GetFieldDefn(i)->GetNameRef());
  }

  if(std::find(FieldNamesTab.begin(), FieldNamesTab.end(), m_LandUseFieldName) != FieldNamesTab.end())
  {
    while ((SRFFeature = SRFLayer->GetNextFeature()) != nullptr)
    {
      ID = SRFFeature->GetFieldAsString("ID");
      IDPlot = std::stoi(ID.substr(0, ID.find("N")));
      if (IDPlot != 0)
      {
        m_SQLRequest = "SELECT LandUse FROM " + m_InputPlotsFile.substr(0, m_InputPlotsFile.find(".")) + " WHERE " + m_IDFieldName + " == "  + std::to_string(IDPlot);
        SQLLayer = Plots->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        LandUseValue = SQLLayer->GetFeature(0)->GetFieldAsString(m_LandUseFieldName.c_str());
        SRFFeature->SetField("LandUse", LandUseValue.c_str());
        SRFLayer->SetFeature(SRFFeature);
        Plots->ReleaseResultSet(SQLLayer);
      }
      else
      {
        SRFFeature->SetField("LandUse", "None");
        SRFLayer->SetFeature(SRFFeature);
      }
    }
  }
  else
  {
    while ((SRFFeature = SRFLayer->GetNextFeature()) != nullptr)
    {
      SRFFeature->SetField("LandUse", "None");
      SRFLayer->SetFeature(SRFFeature);
    }
  }

  FieldNamesTab.clear();
  OGRDataSource::DestroyDataSource(Plots);
  OGRDataSource::DestroyDataSource(SRF);

  // ================================================================================
  // Minimum, maximum and neam slope calculations for SRF
  // ================================================================================

  SRF = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );
  SRFLayer = SRF->GetLayer(0);

  Slope = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputSlopeRasterFile).c_str(), GA_ReadOnly );
  SlopeBand = Slope->GetRasterBand(1);

  mp_RasterDriver = GetGDALDriverManager()->GetDriverByName("MEM");

  getGeoTransform(Slope);

  XOrigin = mp_GeoTransformVal[0];
  YOrigin = mp_GeoTransformVal[3];
  PixelWidth = mp_GeoTransformVal[1];
  PixelHeight = mp_GeoTransformVal[5];

  SRFLayer = SRF->GetLayer(0);
  SRFLayer->ResetReading();

  while((SRFFeature = SRFLayer->GetNextFeature()) != nullptr)
  {
    ID = SRFFeature->GetFieldAsString("ID");
    IDPlot = std::stoi(ID.substr(0, ID.find("N")));
    if (IDPlot != 0)
    {
      SRFGeometry = SRFFeature->GetGeometryRef();
      XPointsTab.clear();
      YPointsTab.clear();
      if (SRFGeometry->getGeometryType() == 3)
      {
        SRFPolygon = (OGRPolygon*) SRFGeometry;
        NPoints = SRFPolygon->getExteriorRing()->getNumPoints();
        for (int i = 0; i < NPoints; i++)
        {
          SRFPolygon->getExteriorRing()->getPoint(i,&SRFPoint);
          XPointsTab.push_back(SRFPoint.getX());
          YPointsTab.push_back(SRFPoint.getY());
        }
      }
      else if (SRFGeometry->getGeometryType() == 6)
      {
        SRFMultiPolygon = (OGRMultiPolygon*) SRFGeometry;
        for (int i = 0; i < SRFMultiPolygon->getNumGeometries(); i++)
        {
          SRFPolygon = (OGRPolygon*) SRFMultiPolygon->getGeometryRef(i);
          NPoints = SRFPolygon->getExteriorRing()->getNumPoints();
          for (int j = 0; j < NPoints; j++)
          {
            SRFPolygon->getExteriorRing()->getPoint(j,&SRFPoint);
            XPointsTab.push_back(SRFPoint.getX());
            YPointsTab.push_back(SRFPoint.getY());
          }
        }
      }
      XMin = XPointsTab.at(std::distance(std::begin(XPointsTab), std::min_element(std::begin(XPointsTab), std::end(XPointsTab))));
      XMax = XPointsTab.at(std::distance(std::begin(XPointsTab), std::max_element(std::begin(XPointsTab), std::end(XPointsTab))));
      YMin = YPointsTab.at(std::distance(std::begin(YPointsTab), std::min_element(std::begin(YPointsTab), std::end(YPointsTab))));
      YMax = YPointsTab.at(std::distance(std::begin(YPointsTab), std::max_element(std::begin(YPointsTab), std::end(YPointsTab))));
      XPointsTab.clear();
      YPointsTab.clear();
      XOff = int((XMin - XOrigin)/PixelWidth);
      YOff = int((YOrigin - YMax)/PixelWidth);
      XCount = int((XMax - XMin)/PixelWidth)+1;
      YCount = int((YMax - YMin)/PixelWidth)+1;
      Target = mp_RasterDriver->Create("", XCount+1, YCount+1, 1, GDT_Byte, Options);
      double NewGeoTransformVal[6] = { XOrigin+(XOff*PixelWidth), PixelWidth, 0, YOrigin-(YOff*PixelWidth), 0, PixelHeight };
      Target->SetGeoTransform(NewGeoTransformVal);
      Target->SetProjection(Slope->GetProjectionRef());
      int BandList[1] = {1};
      double GeomBurnValue[1] = {1};
      char **NewOptions = nullptr;
      NewOptions = CSLSetNameValue(NewOptions, "ALL_TOUCHED", "TRUE");
      GDALProgressFunc ProgressFunc = nullptr;
      GDALTransformerFunc   TransformerFunc = nullptr;
      CPLErr ERR;
      std::vector <OGRGeometryH> GeometryList;
      GeometryList.push_back((OGRGeometryH) SRFGeometry->clone());
      ERR = GDALRasterizeGeometries(Target,
                                    1,
                                    BandList,
                                    1,
                                    (OGRGeometryH*)&GeometryList[0],
                                    TransformerFunc,
                                    nullptr,
                                    GeomBurnValue,
                                    NewOptions,
                                    ProgressFunc,
                                    nullptr);
      TargetBand = Target->GetRasterBand(1);
      SlopeRow = (float*) CPLMalloc(sizeof(float) *Target->GetRasterXSize());
      TargetRow = (int*) CPLMalloc(sizeof(int) *Target->GetRasterXSize());
      SlopeValuesTab.clear();
      for (int i = 0; i < Target->GetRasterYSize(); i++)
      {
        for (int j = 0; j < Target->GetRasterXSize(); j++)
        {
          SlopeBand->RasterIO( GF_Read, XOff, YOff+i, Target->GetRasterXSize(), 1, SlopeRow, Target->GetRasterXSize(), 1, GDT_Float32, 0, 0 );
          TargetBand->RasterIO( GF_Read, 0, i, Target->GetRasterXSize(), 1, TargetRow, Target->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
          if (TargetRow[j] != 0)
          {
            SlopeValuesTab.push_back(SlopeRow[j]);
          }
        }
      }
      SlopeMin = SlopeValuesTab.at(std::distance(std::begin(SlopeValuesTab), std::min_element(std::begin(SlopeValuesTab), std::end(SlopeValuesTab))));
      SlopeMax = SlopeValuesTab.at(std::distance(std::begin(SlopeValuesTab), std::max_element(std::begin(SlopeValuesTab), std::end(SlopeValuesTab))));
      SlopeMean = (double) std::accumulate(SlopeValuesTab.begin(), SlopeValuesTab.end(), 0)/SlopeValuesTab.size();
      SlopeValuesTab.clear();
      SRFFeature->SetField("SlopeMin", SlopeMin);
      SRFFeature->SetField("SlopeMax", SlopeMax);
      SRFFeature->SetField("SlopeMean", SlopeMean);
      SRFLayer->SetFeature(SRFFeature);
      CPLFree(SlopeRow);
      CPLFree(TargetRow);
      GDALClose( Target );
      CSLDestroy(NewOptions);
    }
    else
    {
      SRFFeature->SetField("SlopeMin", 0);
      SRFFeature->SetField("SlopeMax", 0);
      SRFFeature->SetField("SlopeMean", 0);
      SRFLayer->SetFeature(SRFFeature);
    }
  }

  GDALClose( Slope );
  OGRDataSource::DestroyDataSource(SRF);

  BLOCK_EXIT(__PRETTY_FUNCTION__)
}


// =====================================================================
// =====================================================================


void LandProcessor::setLNRParameters()
{
  BLOCK_ENTER(__PRETTY_FUNCTION__)

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  double Length;

  OGRDataSource *LNR;
  OGRLayer *LNRLayer;
  OGRFeature *LNRFeature;

  OGRMultiLineString* LNRMultiLineString;

  // ======================================================================================
  // Checking if necessary vector files are present in the appropriate folder
  // ======================================================================================

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::setLNRParameters(): " + m_OutputLinearEntitiesVectorFile + ": no such file in the output vector directory");
  }

  // ================================================================================
  //  Calculate length for LNR (linear entities)
  // ================================================================================

  LNR = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str(), TRUE );
  LNRLayer = LNR->GetLayer(0);
  LNRLayer->ResetReading();

  while ((LNRFeature = LNRLayer->GetNextFeature()) != nullptr)
  {
    if(LNRFeature->GetGeometryRef()->getGeometryType() == 2)
    {
      LNRFeature->SetField("Length", ((OGRLineString *) LNRFeature->GetGeometryRef())->get_Length());
      LNRLayer->SetFeature(LNRFeature);
    }
    else if (LNRFeature->GetGeometryRef()->getGeometryType() == 5)
    {
      Length = 0;
      LNRMultiLineString = (OGRMultiLineString*) LNRFeature->GetGeometryRef();
      for (int i = 0; i < LNRMultiLineString->getNumGeometries(); i++)
      {
        Length += ((OGRLineString*) LNRMultiLineString->getGeometryRef(i))->get_Length();
      }
      LNRFeature->SetField("Length", Length);
      LNRLayer->SetFeature(LNRFeature);
    }
  }

  OGRDataSource::DestroyDataSource(LNR);

  BLOCK_EXIT(__PRETTY_FUNCTION__)
}


// =====================================================================
// =====================================================================


void LandProcessor::createSU()
{
  BLOCK_ENTER(__PRETTY_FUNCTION__)

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  std::string FieldName, FileName;

  OGRDataSource *SRF, *SU;
  OGRLayer *SQLLayer, *SRFLayer, *SULayer;
  OGRFeature *SRFFeature, *SUFeature;
  OGRGeometry *SUGeometry;
  OGRFeatureDefn *FeatureDefn;

  // ======================================================================================
  // Checking if necessary vector files are present in the appropriate folder
  // ======================================================================================

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::createSU(): " + m_OutputSurfaceEntitiesVectorFile + ": no such file in the output vector directory");
  }

  // ====================================================================
  // Creating the SU file
  // ====================================================================

  SRF = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputSUVectorFile).c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputSUVectorFile).c_str());
  }

  m_SQLRequest = "SELECT * FROM " + m_OutputSurfaceEntitiesVectorFile.substr(0, m_OutputSurfaceEntitiesVectorFile.find(".")) + " WHERE ID != " + m_QMark + "0N1" + m_QMark;
  FileName = m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find("."));
  SQLLayer = SRF->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
  SRF->CopyLayer(SQLLayer, FileName.c_str(), nullptr);
  SRF->ReleaseResultSet(SQLLayer);

  OGRDataSource::DestroyDataSource(SRF);

  // ====================================================================
  // Adding FlowDist attribute that will contain flow distance
  // ====================================================================

  SU = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSUVectorFile).c_str(), TRUE );
  SULayer = SU->GetLayer(0);

  FieldName = "FlowDist";
  createField(SULayer, FieldName.c_str(), OFTReal);
  SULayer->SyncToDisk();

  SULayer = SU->GetLayer(0);
  SULayer->ResetReading();

  while((SUFeature = SULayer->GetNextFeature()) != nullptr)
  {
    SUFeature->SetField(FieldName.c_str(), 0);
    SULayer->SetFeature(SUFeature);
  }

  SULayer->SyncToDisk();
  OGRDataSource::DestroyDataSource(SU);

  BLOCK_EXIT(__PRETTY_FUNCTION__)
}


// =====================================================================
// =====================================================================


void LandProcessor::createRS()
{
  BLOCK_ENTER(__PRETTY_FUNCTION__)

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int FID;

  double Value, LSLength, Length;

  std::string ID, IDTo;

  std::vector <std::string> FieldNamesTab, FieldNamesLTab, FileNamesTab = {m_InputDitchesFile, m_InputThalwegsFile, m_InputRivesrFile};;

  OGRDataSource *DataSource, *LNR, *RS;
  OGRLayer *SQLLayer, *LNRLayer, *RSLayer;
  OGRFeature *LNRFeature, *RSFeature;
  OGRGeometry *RSGeometry, *IntersectionGeometry;
  OGRFeatureDefn *FeatureDefn;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  // ======================================================================================
  // Checking if necessary vector files are present in the appropriate folder
  // ======================================================================================

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::createRS(): " + m_OutputLinearEntitiesVectorFile + ": no such file in the output vector directory");
  }

  if (!openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputDitchesFile).c_str()) and
      !openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputThalwegsFile).c_str()) and
      !openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputRivesrFile).c_str()))
  {
    std::cout << "LandProcessor::createRS(): There are no linear structure files in the input vector directory that could be used to create RS vector" << std::endl;
  }
  else
  {

    // ====================================================================
    // Creating the RS file
    // ====================================================================

    LNR = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str(), TRUE );

    if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputRSVectorFile).c_str()))
    {
      mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputRSVectorFile).c_str());
    }

    RS = mp_VectorDriver->CopyDataSource(LNR, getOutputVectorPath(m_OutputRSVectorFile).c_str(), nullptr);

    OGRDataSource::DestroyDataSource(RS);
    OGRDataSource::DestroyDataSource(LNR);

    // ==========================================================================
    // Creating attributes and populating them
    // ==========================================================================

    RS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );

    FieldNamesTab.clear();

    RSLayer = RS->GetLayer(0);
    RSLayer->ResetReading();

    FieldNamesTab = {"FlowDist", "Ditches", "Thalwegs", "WaterCs", "DitchesL", "ThalwegsL", "WaterCsL", "SurfToLen"};

    for (int i = 0; i < FieldNamesTab.size(); i++)
    {
      createField(RSLayer, FieldNamesTab[i].c_str(), OFTReal);
    }

    RSLayer = RS->GetLayer(0);
    RSLayer->ResetReading();

    while ((RSFeature = RSLayer->GetNextFeature()) != nullptr)
    {
      for (int i = 0; i < FieldNamesTab.size(); i++)
      {
        RSFeature->SetField(FieldNamesTab[i].c_str(), 0);
        RSLayer->SetFeature(RSFeature);
      }
    }

    FieldNamesTab.clear();

    OGRDataSource::DestroyDataSource(RS);

    // ================================================================================
    // Deleting all RS that have 'None' in IDTo field (not crossed by the surface flow)
    // ================================================================================

    RS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );
    RSLayer = RS->GetLayer(0);
    RSLayer->ResetReading();

    while ((RSFeature = RSLayer->GetNextFeature()) != nullptr)
    {
      FID = RSFeature->GetFID();
      IDTo = "None";
      if (RSFeature->GetFieldAsString("IDTo") == IDTo)
      {
        RSLayer->DeleteFeature(FID);
      }
    }

    m_SQLRequest = "REPACK " + m_OutputRSVectorFile.substr(0, m_OutputRSVectorFile.find("."));
    RS->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
    OGRDataSource::DestroyDataSource(RS);

    // =========================================================================================
    // Setting IDTo to 'None' since (by convention) RS# segments are exits
    // =========================================================================================

    RS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );
    RSLayer = RS->GetLayer(0);
    RSLayer->ResetReading();

    while ((RSFeature = RSLayer->GetNextFeature()) != nullptr)
    {
      RSFeature->SetField("IDTo", "None");
      RSLayer->SetFeature(RSFeature);
    }

    RSLayer->SyncToDisk();
    OGRDataSource::DestroyDataSource(RS);

    // ================================================================================
    // Populating attributes that suppose to contain information about linear structures
    // ================================================================================

    FieldNamesTab.clear();
    FieldNamesLTab.clear();

    FieldNamesTab = {"Ditches", "Thalwegs", "WaterCs"};
    FieldNamesLTab = {"DitchesL", "ThalwegsL", "WaterCsL"};

    for (int i = 0; i < FileNamesTab.size(); i++)
    {
      if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(FileNamesTab[i]).c_str()))
      {
        DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
        RS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );
        RSLayer = RS->GetLayer(0);
        for (int j = 0; j < RSLayer->GetFeatureCount(); j++)
        {
          RSFeature = RSLayer->GetFeature(j);
          RSGeometry = RSFeature->GetGeometryRef();
          m_SQLRequest = "SELECT ST_Length(ST_Intersection(geometry,(SELECT ST_Buffer(geometry,0.05) FROM " + m_OutputRSVectorFile.substr(0, m_OutputRSVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(RSFeature->GetFID()) + "))) AS Length "
              "FROM " + FileNamesTab[i].substr(0, FileNamesTab[i].find(".shp")) + " "
              "WHERE ST_Length(ST_Intersection(geometry,(SELECT ST_Buffer(geometry,0.05) FROM " + m_OutputRSVectorFile.substr(0, m_OutputRSVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(RSFeature->GetFID()) + "))) > 0.1";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() > 0)
          {
            Length = 0;
            for (int k = 0; k < SQLLayer->GetFeatureCount(); k++)
            {
              Length = Length + SQLLayer->GetFeature(k)->GetFieldAsDouble("Length");
            }
            RSFeature->SetField(FieldNamesLTab[i].c_str(), Length);
            RSLayer->SetFeature(RSFeature);
          }
          DataSource->ReleaseResultSet(SQLLayer);
        }
        RSLayer->SyncToDisk();
        OGRDataSource::DestroyDataSource(RS);
        OGRDataSource::DestroyDataSource(DataSource);
      }
    }

    RS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );
    RSLayer = RS->GetLayer(0);
    RSLayer->ResetReading();

    while ((RSFeature = RSLayer->GetNextFeature()) != nullptr)
    {
      for (int i = 0; i < FieldNamesLTab.size(); i++)
      {
        Length = RSFeature->GetFieldAsDouble("Length");
        LSLength = RSFeature->GetFieldAsDouble(FieldNamesLTab[i].c_str());
        if (LSLength > Length)
        {
          RSFeature->SetField(FieldNamesLTab[i].c_str(), Length);
          RSLayer->SetFeature(RSFeature);
        }
      }
    }

    RSLayer->SyncToDisk();
    OGRDataSource::DestroyDataSource(RS);

    RS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );
    RSLayer = RS->GetLayer(0);
    RSLayer->ResetReading();

    while ((RSFeature = RSLayer->GetNextFeature()) != nullptr)
    {
      for (int i = 0; i < FieldNamesLTab.size(); i++)
      {
        Length = RSFeature->GetFieldAsDouble("Length");
        LSLength = RSFeature->GetFieldAsDouble(FieldNamesLTab[i].c_str());
        if (LSLength != Length and LSLength < 0.05*6)
        {
          RSFeature->SetField(FieldNamesLTab[i].c_str(), 0);
          RSLayer->SetFeature(RSFeature);
        }
      }
    }

    RSLayer->SyncToDisk();
    OGRDataSource::DestroyDataSource(RS);

    // =====================================================================
    // Removing RS that do not have any linear structures attributed to them
    // =====================================================================

    RS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );
    RSLayer = RS->GetLayer(0);
    RSLayer->ResetReading();

    while ((RSFeature = RSLayer->GetNextFeature()) != nullptr)
    {
      FID = RSFeature->GetFID();
      if (RSFeature->GetFieldAsDouble("DitchesL") == 0 and
          RSFeature->GetFieldAsDouble("ThalwegsL") == 0 and RSFeature->GetFieldAsDouble("WaterCsL") == 0)
      {
        RSLayer->DeleteFeature(FID);
      }
    }

    m_SQLRequest = "REPACK " + m_OutputRSVectorFile.substr(0, m_OutputRSVectorFile.find("."));
    RS->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);

    OGRDataSource::DestroyDataSource(RS);

    // ====================================================================================
    // Setting attributes to contain the linear structure length to actual RS# length ratio
    // ====================================================================================

    RS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );
    RSLayer = RS->GetLayer(0);
    RSLayer->ResetReading();

    while ((RSFeature = RSLayer->GetNextFeature()) != nullptr)
    {
      for (int i = 0; i < FieldNamesLTab.size(); i++)
      {
        Length = RSFeature->GetFieldAsDouble("Length");
        LSLength = RSFeature->GetFieldAsDouble(FieldNamesLTab[i].c_str());
        RSFeature->SetField(FieldNamesTab[i].c_str(), LSLength/Length);
        RSLayer->SetFeature(RSFeature);
      }
      RSLayer->SyncToDisk();
    }

    FileNamesTab.clear();
    FieldNamesTab.clear();
    FieldNamesLTab.clear();

    RSLayer->SyncToDisk();
    OGRDataSource::DestroyDataSource(RS);
  }

  BLOCK_EXIT(__PRETTY_FUNCTION__)
}


// ====================================================================
// ====================================================================


void LandProcessor::createLI()
{
  BLOCK_ENTER(__PRETTY_FUNCTION__)

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int FID;

  double Value, LSLength, Length;

  std::string ID, IDTo;

  std::vector <std::string> FieldNamesTab, FieldNamesLTab, FileNamesTab = {m_InputHedgesFile, m_InputGrassBandFile, m_InputBenchesFile};

  OGRDataSource *DataSource, *LNR, *LI;
  OGRLayer *SQLLayer, *LNRLayer, *LILayer;
  OGRFeature *LNRFeature, *LIFeature;
  OGRGeometry *LIGeometry, *IntersectionGeometry;
  OGRFeatureDefn *FeatureDefn;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  // ======================================================================================
  // Checking if necessary vector files are present in the appropriate folder
  // ======================================================================================

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::createLI(): " + m_OutputLinearEntitiesVectorFile + ": no such file in the output vector directory");
  }

  if (!openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputHedgesFile).c_str()) and
      !openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputGrassBandFile).c_str()) and
      !openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputBenchesFile).c_str()))
  {
    std::cout << "LandProcessor::createLI(): There are no linear structure files in the input vector directory that could be used to create LI vector" << std::endl;
  }
  else
  {

    // ====================================================================
    // Creating the LI file
    // ====================================================================

    LNR = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str(), TRUE );

    if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLIVectorFile).c_str()))
    {
      mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputLIVectorFile).c_str());
    }

    LI = mp_VectorDriver->CopyDataSource(LNR, getOutputVectorPath(m_OutputLIVectorFile).c_str(), nullptr);

    OGRDataSource::DestroyDataSource(LI);
    OGRDataSource::DestroyDataSource(LNR);

    // ====================================================================================
    // Creating attributes and populating them
    // ====================================================================================

    LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );

    FieldNamesTab.clear();

    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    FieldNamesTab = {"FlowDist", "Hedges", "GrassBs", "Benches", "HedgesL", "GrassBsL", "BenchesL", "SurfToLen"};

    for (int i = 0; i < FieldNamesTab.size(); i++)
    {
      createField(LILayer, FieldNamesTab[i].c_str(), OFTReal);
    }

    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    while ((LIFeature = LILayer->GetNextFeature()) != nullptr)
    {
      for (int i = 0; i < FieldNamesTab.size(); i++)
      {
        LIFeature->SetField(FieldNamesTab[i].c_str(), 0);
        LILayer->SetFeature(LIFeature);
      }
    }

    FieldNamesTab.clear();

    OGRDataSource::DestroyDataSource(LI);

    // =================================================================================================
    // Deleting all LI that have 'None' as value of the attribute IDTo (not crossed by the surface flow)
    // =================================================================================================

    LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    while ((LIFeature = LILayer->GetNextFeature()) != nullptr)
    {
      FID = LIFeature->GetFID();
      IDTo = "None";
      if (LIFeature->GetFieldAsString("IDTo") == IDTo)
      {
        LILayer->DeleteFeature(FID);
      }
    }

    m_SQLRequest = "REPACK " + m_OutputLIVectorFile.substr(0, m_OutputLIVectorFile.find("."));
    LI->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
    OGRDataSource::DestroyDataSource(LI);

    // ====================================================================
    // Setting IDTo attribute in case when there is an RS with the same ID
    // ====================================================================

    if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputRSVectorFile).c_str()))
    {
      DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );

      LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
      LILayer = LI->GetLayer(0);

      while ((LIFeature = LILayer->GetNextFeature()) != nullptr)
      {
        ID = LIFeature->GetFieldAsString("ID");
        m_SQLRequest = "SELECT ID FROM " + m_OutputRSVectorFile.substr(0, m_OutputRSVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark;
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          LIFeature->SetField("IDTo", ID.c_str());
          LILayer->SetFeature(LIFeature);
        }
        DataSource->ReleaseResultSet(SQLLayer);
      }

      LILayer->SyncToDisk();
      OGRDataSource::DestroyDataSource(LI);
    }

    // =========================================================================================
    // Setting attributes that contain information about linear structures attributed to them
    // =========================================================================================

    FieldNamesTab.clear();
    FieldNamesLTab.clear();

    FieldNamesTab = {"Hedges", "GrassBs", "Benches"};
    FieldNamesLTab = {"HedgesL", "GrassBsL", "BenchesL"};

    for (int i = 0; i < FileNamesTab.size(); i++)
    {
      if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(FileNamesTab[i]).c_str()))
      {
        DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
        LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
        LILayer = LI->GetLayer(0);
        for (int j = 0; j < LILayer->GetFeatureCount(); j++)
        {
          LIFeature = LILayer->GetFeature(j);
          LIGeometry = LIFeature->GetGeometryRef();
          m_SQLRequest = "SELECT ST_Length(ST_Intersection(geometry,(SELECT ST_Buffer(geometry,0.05) FROM " + m_OutputLIVectorFile.substr(0, m_OutputLIVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(LIFeature->GetFID()) + "))) AS Length "
              "FROM " + FileNamesTab[i].substr(0, FileNamesTab[i].find(".shp")) + " "
              "WHERE ST_Length(ST_Intersection(geometry,(SELECT ST_Buffer(geometry,0.05) FROM " + m_OutputLIVectorFile.substr(0, m_OutputLIVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(LIFeature->GetFID()) + "))) > 0.1";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() > 0)
          {
            Length = 0;
            for (int k = 0; k < SQLLayer->GetFeatureCount(); k++)
            {
              Length = Length + SQLLayer->GetFeature(k)->GetFieldAsDouble("Length");
            }
            LIFeature->SetField(FieldNamesLTab[i].c_str(), Length);
            LILayer->SetFeature(LIFeature);
          }
          DataSource->ReleaseResultSet(SQLLayer);
        }
        LILayer->SyncToDisk();
        OGRDataSource::DestroyDataSource(LI);
        OGRDataSource::DestroyDataSource(DataSource);
      }
    }

    LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    while ((LIFeature = LILayer->GetNextFeature()) != nullptr)
    {
      for (int i = 0; i < FieldNamesLTab.size(); i++)
      {
        Length = LIFeature->GetFieldAsDouble("Length");
        LSLength = LIFeature->GetFieldAsDouble(FieldNamesLTab[i].c_str());
        if (LSLength > Length)
        {
          LIFeature->SetField(FieldNamesLTab[i].c_str(), Length);
          LILayer->SetFeature(LIFeature);
        }
      }
    }

    LILayer->SyncToDisk();
    OGRDataSource::DestroyDataSource(LI);

    LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    while ((LIFeature = LILayer->GetNextFeature()) != nullptr)
    {
      for (int i = 0; i < FieldNamesLTab.size(); i++)
      {
        Length = LIFeature->GetFieldAsDouble("Length");
        LSLength = LIFeature->GetFieldAsDouble(FieldNamesLTab[i].c_str());
        if (LSLength != Length and LSLength < 0.05*6)
        {
          LIFeature->SetField(FieldNamesLTab[i].c_str(), 0);
          LILayer->SetFeature(LIFeature);
        }
      }
    }

    LILayer->SyncToDisk();
    OGRDataSource::DestroyDataSource(LI);

    // =================================================================================
    // Removing those entities that do not have any linear structures attributed to them
    // =================================================================================

    LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    while ((LIFeature = LILayer->GetNextFeature()) != nullptr)
    {
      FID = LIFeature->GetFID();
      if (LIFeature->GetFieldAsDouble("HedgesL") == 0 and
          LIFeature->GetFieldAsDouble("GrassBsL") == 0 and LIFeature->GetFieldAsDouble("BenchesL") == 0)
      {
        LILayer->DeleteFeature(FID);
      }
    }

    m_SQLRequest = "REPACK " + m_OutputLIVectorFile.substr(0, m_OutputLIVectorFile.find("."));
    LI->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);

    OGRDataSource::DestroyDataSource(LI);

    // ===========================================================================================
    // Setting attributes to contain the ration of linear structure length to actual entity length
    // ===========================================================================================

    LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    while ((LIFeature = LILayer->GetNextFeature()) != nullptr)
    {
      for (int i = 0; i < FieldNamesLTab.size(); i++)
      {
        Length = LIFeature->GetFieldAsDouble("Length");
        LSLength = LIFeature->GetFieldAsDouble(FieldNamesLTab[i].c_str());
        LIFeature->SetField(FieldNamesTab[i].c_str(), LSLength/Length);
        LILayer->SetFeature(LIFeature);
      }
      LILayer->SyncToDisk();
    }

    FileNamesTab.clear();
    FieldNamesTab.clear();
    FieldNamesLTab.clear();

    LILayer->SyncToDisk();
    OGRDataSource::DestroyDataSource(LI);
  }

  BLOCK_EXIT(__PRETTY_FUNCTION__)
}


// ====================================================================
// ====================================================================


void LandProcessor::setSUParameters()
{
  BLOCK_ENTER(__PRETTY_FUNCTION__)

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int FID, N1, N2;

  std::string ID, IDTo, IDNew, IDToNew;

  std::vector <std::string> FieldNamesTab;

  OGRDataSource *DataSource, *SU;
  OGRLayer *SQLLayer, *SULayer;
  OGRFeature *SUFeature;
  OGRGeometry *SUGeometry;
  OGRFeatureDefn *FeatureDefn;

  // ======================================================================================
  // Checking if necessary vector files are present in the appropriate folder
  // ======================================================================================

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputSUVectorFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::setSUParameters(): " + m_OutputSUVectorFile + ": no such file in the output vector directory");
  }

  // ==================================================================================
  // Set IDTo parameter based on LI# and RS#
  // ==================================================================================

  SU = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSUVectorFile).c_str(), TRUE );
  SULayer = SU->GetLayer(0);
  SULayer->ResetReading();

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );

  while ((SUFeature = SULayer->GetNextFeature()) != nullptr)
  {
    ID = SUFeature->GetFieldAsString("ID");
    IDNew = "SU#" + ID;
    SUFeature->SetField("ID", IDNew.c_str());
    SULayer->SetFeature(SUFeature);
  }

  SULayer->SyncToDisk();
  SULayer = SU->GetLayer(0);
  SULayer->ResetReading();

  while ((SUFeature = SULayer->GetNextFeature()) != nullptr)
  {
    ID = SUFeature->GetFieldAsString("ID");
    IDTo = SUFeature->GetFieldAsString("IDTo");
    if (IDTo.find("-") != std::string::npos)
    {
      if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLIVectorFile).c_str()))
      {
        m_SQLRequest = "SELECT * FROM " + m_OutputLIVectorFile.substr(0, m_OutputLIVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark;
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() == 0)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputRSVectorFile).c_str()))
          {
            m_SQLRequest = "SELECT * FROM " + m_OutputRSVectorFile.substr(0, m_OutputRSVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            if (SQLLayer->GetFeatureCount() == 0)
            {
              DataSource->ReleaseResultSet(SQLLayer);
              IDToNew = "SU#" + IDTo.substr(IDTo.find("-")+1, IDTo.length()-1);
              SUFeature->SetField("IDTo", IDToNew.c_str());
              SULayer->SetFeature(SUFeature);
            }
            else
            {
              DataSource->ReleaseResultSet(SQLLayer);
              IDToNew = "RS#" + IDTo;
              SUFeature->SetField("IDTo", IDToNew.c_str());
              SULayer->SetFeature(SUFeature);
            }
          }
        }
        else
        {
          DataSource->ReleaseResultSet(SQLLayer);
          IDToNew = "LI#" + IDTo;
          SUFeature->SetField("IDTo", IDToNew.c_str());
          SULayer->SetFeature(SUFeature);
        }
      }
    }
    else
    {
      IDToNew = "SU#" + IDTo;
      SUFeature->SetField("IDTo", IDToNew.c_str());
      SULayer->SetFeature(SUFeature);
    }
  }

  OGRDataSource::DestroyDataSource(SU);
  OGRDataSource::DestroyDataSource(DataSource);

  // ==================================================================================
  // Set IDTo parameter in case there are no LI# or RS# that correspond to it
  // ==================================================================================

  SU = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSUVectorFile).c_str(), TRUE );
  SULayer = SU->GetLayer(0);
  SULayer->ResetReading();

  while ((SUFeature = SULayer->GetNextFeature()) != nullptr)
  {
    ID = SUFeature->GetFieldAsString("ID");
    IDTo = SUFeature->GetFieldAsString("IDTo");
    if (IDTo.find("-") != std::string::npos and
        IDTo.find("RS") == std::string::npos and
        IDTo.find("LI") == std::string::npos)
    {
      IDToNew = "SU#" + IDTo.substr(IDTo.find("-")+1, IDTo.length()-1);
      SUFeature->SetField("IDTo", IDToNew.c_str());
      SULayer->SetFeature(SUFeature);
    }
  }

  OGRDataSource::DestroyDataSource(SU);

  // ==================================================================================
  // Setting flow distance parameter for SU#
  // ==================================================================================

  SU = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSUVectorFile).c_str(), TRUE );
  SULayer = SU->GetLayer(0);
  SULayer->ResetReading();

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );

  while ((SUFeature = SULayer->GetNextFeature()) != nullptr)
  {
    ID = SUFeature->GetFieldAsString("ID");
    IDTo = SUFeature->GetFieldAsString("IDTo");
    if (IDTo.find("SU") != std::string::npos)
    {
      if (IDTo != "SU#0N1")
      {
        m_SQLRequest = "SELECT ST_Intersection(geometry, "
            "(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + ")) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark;
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeature(0)->GetGeometryRef()->getGeometryType() == 1)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + " AND "
              "ST_Within((SELECT ST_Centroid(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + "),"
              "(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + "))";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          N1 = SQLLayer->GetFeatureCount();
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + " AND "
              "ST_Within((SELECT ST_Centroid(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
              "(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          N2 = SQLLayer->GetFeatureCount();
          DataSource->ReleaseResultSet(SQLLayer);
          if (N1 == 1 and N2 == 1)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_Centroid(geometry),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "), geometry)))) + "
                "ST_Distance((SELECT ST_Centroid(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "geometry)))) AS Distance FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark; //TODO: simplify complicated SQL request
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
          else if (N1 == 0 and N2 == 0)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_PointOnSurface(geometry),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "), geometry)))) + "
                "ST_Distance((SELECT ST_PointOnSurface(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "geometry)))) AS Distance FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
          else if (N1 == 1 and N2 == 0)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_Centroid(geometry),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "), geometry)))) + "
                "ST_Distance((SELECT ST_PointOnSurface(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "geometry)))) AS Distance FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
          else if (N1 == 0 and N2 == 1)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_PointOnSurface(geometry),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "), geometry)))) + "
                "ST_Distance((SELECT ST_Centroid(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "geometry)))) AS Distance FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
        }
        else
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + " AND "
              "ST_Within((SELECT ST_Centroid(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + "),"
              "(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + "))";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          N1 = SQLLayer->GetFeatureCount();
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + " AND "
              "ST_Within((SELECT ST_Centroid(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
              "(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          N2 = SQLLayer->GetFeatureCount();
          DataSource->ReleaseResultSet(SQLLayer);
          if (N1 == 1 and N2 == 1)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_Centroid(geometry), ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),"
                "(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) + "
                "ST_Distance((SELECT ST_Centroid(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID ==" + m_QMark + IDTo + m_QMark + "))), 0.5)) "
                "AS Distance FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
          else if (N1 == 0 and N2 == 0)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_PointOnSurface(geometry), ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),"
                "(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) + "
                "ST_Distance((SELECT ST_PointOnSurface(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) "
                "AS Distance FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
          else if (N1 == 1 and N2 == 0)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_Centroid(geometry), ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),"
                "(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) + "
                "ST_Distance((SELECT ST_PointOnSurface(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) "
                "AS Distance FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
          else if (N1 == 0 and N2 == 1)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_PointOnSurface(geometry), ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),"
                "(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) + "
                "ST_Distance((SELECT ST_Centroid(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) "
                "AS Distance FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
        }
      }
    }
    else if (IDTo.find("LI") != std::string::npos)
    {
      IDTo = IDTo.substr(IDTo.find("#")+1, IDTo.length()-1);
      m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + " AND "
          "ST_Within((SELECT ST_Centroid(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + "),"
          "(SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + "))";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() == 1)
      {
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT geometry, ST_Distance(ST_Centroid(geometry), (SELECT geometry FROM " + m_OutputLIVectorFile.substr(0, m_OutputLIVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark +"))" //Used to be ST_Line_Interpolate_Point(ST_LineMerge(geometry), 0.5) in the attribution
            " AS Distance FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark;
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
        SULayer->SetFeature(SUFeature);
        DataSource->ReleaseResultSet(SQLLayer);
      }
      else if (SQLLayer->GetFeatureCount() == 0)
      {
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT geometry, ST_Distance(ST_PointOnSurface(geometry), (SELECT geometry FROM " + m_OutputLIVectorFile.substr(0, m_OutputLIVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark +"))"
            " AS Distance FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark; //16/12/2016/17:00 Simplified the selection here
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
        SULayer->SetFeature(SUFeature);
        DataSource->ReleaseResultSet(SQLLayer);
      }
      else
      {
        DataSource->ReleaseResultSet(SQLLayer);
      }
    }
    else if (IDTo.find("RS") != std::string::npos)
    {
      IDTo = IDTo.substr(IDTo.find("#")+1, IDTo.length()-1);
      m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + " AND "
          "ST_Within(ST_Centroid(geometry), (SELECT geometry FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark + "))"; //16/12/2016/17:00 Simplified the selection here
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() == 1)
      {
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT geometry, ST_Distance(ST_Centroid(geometry), (SELECT geometry FROM " + m_OutputRSVectorFile.substr(0, m_OutputRSVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark +"))" //TODO: Use ST_Line_Interpolate_Point(ST_LineMerge(geometry), 0.5) in the attribution
            " AS Distance FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark; //16/12/2016/17:00 Simplified the selection here
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
        SULayer->SetFeature(SUFeature);
        DataSource->ReleaseResultSet(SQLLayer);
      }
      else if (SQLLayer->GetFeatureCount() == 0)
      {
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT geometry, ST_Distance(ST_PointOnSurface(geometry), (SELECT geometry FROM " + m_OutputRSVectorFile.substr(0, m_OutputRSVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDTo + m_QMark +"))"
            " AS Distance FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark; //16/12/2016/17:00 Simplified the selection here
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
        SULayer->SetFeature(SUFeature);
        DataSource->ReleaseResultSet(SQLLayer);
      }
      else
      {
        DataSource->ReleaseResultSet(SQLLayer);
      }
    }
  }

  OGRDataSource::DestroyDataSource(DataSource);
  OGRDataSource::DestroyDataSource(SU);

  BLOCK_EXIT(__PRETTY_FUNCTION__)
}


// ====================================================================
// ====================================================================


void LandProcessor::setRSParameters()
{

  BLOCK_ENTER(__PRETTY_FUNCTION__)

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  double Value, Length, Surface;

  std::string ID, IDTo, IDFrom, IDNew;

  OGRDataSource *DataSource, *RS;
  OGRLayer *SQLLayer, *RSLayer;
  OGRFeature *RSFeature;

  // ======================================================================================
  // Checking if necessary vector files are present in the appropriate folder
  // ======================================================================================

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputRSVectorFile).c_str()))
  {
    std::cout << "LandProcessor::setRSParameters(): " + m_OutputRSVectorFile + ": no such file in the output vector directory" << std::endl;
  }
  else
  {

    // =====================================================================================
    // Calculate ration S(SU#)/L(RS#) (if possible). For that, we extract the identifier
    // of the SU# that drain into that RS#, look it up in the SU# file, extract its surface and
    // divide it by the length of the RS
    // =====================================================================================

    RS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );
    RSLayer = RS->GetLayer(0);
    RSLayer->ResetReading();

    DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );

    while((RSFeature = RSLayer->GetNextFeature()) != nullptr)
    {
      ID = RSFeature->GetFieldAsString("ID");
      IDFrom = "SU#" + ID.substr(0, ID.find("-"));
      m_SQLRequest = "SELECT Surface FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDFrom + m_QMark + " "
          "AND IDTo == " + m_QMark + "RS#" + ID + m_QMark;
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        Length = RSFeature->GetFieldAsDouble("Length");
        Surface = SQLLayer->GetFeature(0)->GetFieldAsDouble("Surface");
        Value = Surface/Length;
        RSFeature->SetField("SurfToLen", Value);
        RSLayer->SetFeature(RSFeature);
      }
      DataSource->ReleaseResultSet(SQLLayer);
    }

    OGRDataSource::DestroyDataSource(DataSource);
    OGRDataSource::DestroyDataSource(RS);

    // ================================================================================
    //  Set ID as RS# as following: RS#145N1-134N7
    // ================================================================================

    RS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );
    RSLayer = RS->GetLayer(0);
    RSLayer->ResetReading();

    while ((RSFeature = RSLayer->GetNextFeature()) != nullptr)
    {
      ID = RSFeature->GetFieldAsString("ID");
      IDNew = "RS#" + ID;
      RSFeature->SetField("ID", IDNew.c_str());
      RSLayer->SetFeature(RSFeature);
    }

    OGRDataSource::DestroyDataSource(RS);
  }

  BLOCK_EXIT(__PRETTY_FUNCTION__)
}


// ====================================================================
// ====================================================================


void LandProcessor::setLIParameters()
{
  BLOCK_ENTER(__PRETTY_FUNCTION__)

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  double Value, Length, Surface;

  std::string ID, IDTo, IDFrom, IDNew;

  OGRDataSource *DataSource, *LI;
  OGRLayer *SQLLayer, *LILayer;
  OGRFeature *LIFeature;
  OGRGeometry *LIGeometry;

  // ======================================================================================
  // Checking if necessary vector files are present in the appropriate folder
  // ======================================================================================

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLIVectorFile).c_str()))
  {
    std::cout << "LandProcessor::setLIParameters(): " + m_OutputLIVectorFile + ": no such file in the output vector directory" << std::endl;
  }
  else
  {
    // ==========================================================================================
    // Calculate ration S(SU#)/L(LI#) (if possible)
    // ==========================================================================================

    LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );

    while((LIFeature = LILayer->GetNextFeature()) != nullptr)
    {
      ID = LIFeature->GetFieldAsString("ID");
      IDFrom = "SU#" + ID.substr(0, ID.find("-"));
      m_SQLRequest = "SELECT Surface FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + IDFrom + m_QMark + " AND "
          "IDTo == " + m_QMark + "LI#" + ID + m_QMark;
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        Length = LIFeature->GetFieldAsDouble("Length");
        Surface = SQLLayer->GetFeature(0)->GetFieldAsDouble("Surface");
        Value = Surface/Length;
        LIFeature->SetField("SurfToLen", Value);
        LILayer->SetFeature(LIFeature);
      }
      DataSource->ReleaseResultSet(SQLLayer);
    }

    OGRDataSource::DestroyDataSource(DataSource);
    OGRDataSource::DestroyDataSource(LI);

    // ==================================================================================
    // Setting flow distance parameter for LI#
    // ==================================================================================

    LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );

    while ((LIFeature = LILayer->GetNextFeature()) != nullptr)
    {
      ID = LIFeature->GetFieldAsString("ID");
      IDTo = LIFeature->GetFieldAsString("IDTo");
      if (IDTo.find("-") == std::string::npos)
      {
        m_SQLRequest = "SELECT ROWID as RID FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + "SU#" + IDTo + m_QMark + " AND "
            "ST_Within(ST_Centroid(geometry), geometry)";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() == 1)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT geometry, ST_Distance((SELECT ST_Centroid(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + "SU#" + IDTo + m_QMark + "),"
              "geometry) AS Distance FROM " + m_OutputLIVectorFile.substr(0, m_OutputLIVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark; //16/12/2016/17:00 Simplified the selection here
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          LIFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
          LILayer->SetFeature(LIFeature);
          DataSource->ReleaseResultSet(SQLLayer);
        }
        else if (SQLLayer->GetFeatureCount() == 0)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT geometry, ST_Distance((SELECT ST_PointOnSurface(geometry) FROM " + m_OutputSUVectorFile.substr(0, m_OutputSUVectorFile.find(".")) + " WHERE ID == " + m_QMark + "SU#" + IDTo + m_QMark + "),"
              "geometry) AS Distance FROM " + m_OutputLIVectorFile.substr(0, m_OutputLIVectorFile.find(".")) + " WHERE ID == " + m_QMark + ID + m_QMark; //16/12/2016/17:00 Simplified the selection here
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          LIFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
          LILayer->SetFeature(LIFeature);
          DataSource->ReleaseResultSet(SQLLayer);
        }
        else
        {
          DataSource->ReleaseResultSet(SQLLayer);
        }
      }
    }

    OGRDataSource::DestroyDataSource(DataSource);
    OGRDataSource::DestroyDataSource(LI);

    // ================================================================================
    // Set ID as following: LI#145N1-134N7
    // ================================================================================

    LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    while ((LIFeature = LILayer->GetNextFeature()) != nullptr)
    {
      ID = LIFeature->GetFieldAsString("ID");
      IDNew = "LI#" + ID;
      LIFeature->SetField("ID", IDNew.c_str());
      LILayer->SetFeature(LIFeature);
    }

    OGRDataSource::DestroyDataSource(LI);

    // ================================================================================
    //  Set IDTo as following: in case of SU# for the receiving feature:
    //	SU#15N1, in case of RS# for the receiving feature: RS#15N1-164N3
    //  (and 15N1-164N3 must be the same as ID)
    // ================================================================================

    LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    while ((LIFeature = LILayer->GetNextFeature()) != nullptr)
    {
      IDTo = LIFeature->GetFieldAsString("IDTo");
      if (IDTo.find("-") == std::string::npos)
      {
        IDTo = "SU#" + IDTo;
        LIFeature->SetField("IDTo", IDTo.c_str());
        LILayer->SetFeature(LIFeature);
      }
      else
      {
        IDTo = "RS#" + IDTo;
        LIFeature->SetField("IDTo", IDTo.c_str());
        LILayer->SetFeature(LIFeature);
      }
    }

    OGRDataSource::DestroyDataSource(LI);
  }

  BLOCK_EXIT(__PRETTY_FUNCTION__)
}


// ====================================================================
// ====================================================================


void LandProcessor::extractPlotsLimits()
{
  BLOCK_ENTER(__PRETTY_FUNCTION__)


  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int FID, IDPlotA, IDPlotB;

  std::string IDPlotAB, IDPlotBA;

  std::vector <std::string> FieldNamesTab;
  std::vector <std::string> IDPlotTab;

  OGRDataSource *DataSource, *Plots, *PlotsLimits;
  OGRLayer *SQLLayer, *PlotsLayer, *PlotsLimitsLayer;
  OGRFeature *PlotsFeature, *PlotsLimitsFeature;
  OGRGeometry *NewGeometry;
  OGRFeatureDefn *FeatureDefn;

  // ======================================================================================
  // Checking if necessary vector files are present in the appropriate folder
  // ======================================================================================

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsVectorFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::extractPlotsLimits(): " + getOutputVectorPath(m_OutputPlotsVectorFile) + ": no such file in the output vector directory");
  }

  // ====================================================================
  // Extracting parcel limits in order to obtain the well-attributes LSs
  // ====================================================================

  Plots = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsVectorFile).c_str(), TRUE );
  PlotsLayer = Plots->GetLayer(0);
  PlotsLayer->ResetReading();
  mp_SRS = Plots->GetLayer(0)->GetSpatialRef();

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsLimitsVectorFile).c_str()))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputPlotsLimitsVectorFile).c_str());
  }

  PlotsLimits = mp_VectorDriver->CreateDataSource( getOutputVectorPath(m_OutputPlotsLimitsVectorFile).c_str(), nullptr);
  PlotsLimitsLayer = PlotsLimits->CreateLayer( "plotslimits", mp_SRS, wkbLineString, nullptr );

  FieldNamesTab.clear();
  FieldNamesTab = {"ID", "IDPlotA","IDPlotB"};

  for (int i = 0; i < FieldNamesTab.size(); i++)
  {
    createField(PlotsLimitsLayer, FieldNamesTab[i].c_str(), OFTInteger);
  }

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsVectorFile).c_str(), TRUE );
  PlotsLimitsFeature = OGRFeature::CreateFeature(PlotsLimitsLayer->GetLayerDefn());
  IDPlotTab.clear();

  while ((PlotsFeature = PlotsLayer->GetNextFeature()) != nullptr)
  {
    FID = PlotsFeature->GetFID();
    IDPlotA = PlotsFeature->GetFieldAsInteger(m_IDFieldName.c_str());
    m_SQLRequest = "SELECT ST_CollectionExtract(ST_Intersection(geometry, "
        "(SELECT geometry FROM " + m_OutputPlotsVectorFile.substr(0, m_OutputPlotsVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FID) + ")), 2), " + m_IDFieldName + " FROM " + m_OutputPlotsVectorFile.substr(0, m_OutputPlotsVectorFile.find(".")) + " WHERE "
        "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + m_OutputPlotsVectorFile.substr(0, m_OutputPlotsVectorFile.find(".")) + " WHERE ROWID == " + std::to_string(FID) + "))) > 0";
    SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
    if (SQLLayer->GetFeatureCount() != 0)
    {
      for (int i = 0; i < SQLLayer->GetFeatureCount(); i++)
      {
        IDPlotB = SQLLayer->GetFeature(i)->GetFieldAsInteger(m_IDFieldName.c_str());
        IDPlotBA = std::to_string(IDPlotB) + "-" + std::to_string(IDPlotA);
        if (std::find(IDPlotTab.begin(), IDPlotTab.end(), IDPlotBA) == IDPlotTab.end())
        {
          PlotsLimitsFeature->SetField(FieldNamesTab[0].c_str(), PlotsLimits->GetLayer(0)->GetFeatureCount());
          PlotsLimitsFeature->SetField(FieldNamesTab[1].c_str(), IDPlotA);
          PlotsLimitsFeature->SetField(FieldNamesTab[2].c_str(), IDPlotB);
          NewGeometry = PlotsFeature->GetGeometryRef()->Intersection(SQLLayer->GetFeature(i)->GetGeometryRef());
          PlotsLimitsFeature->SetGeometry(NewGeometry);
          PlotsLimitsLayer->CreateFeature(PlotsLimitsFeature);
          PlotsLimits->SyncToDisk();
          IDPlotAB = std::to_string(IDPlotA) + "-" + std::to_string(IDPlotB);;
          IDPlotTab.push_back(IDPlotAB);
        }
      }
    }
    DataSource->ReleaseResultSet(SQLLayer);
  }

  FieldNamesTab.clear();
  IDPlotTab.clear();
  OGRDataSource::DestroyDataSource(Plots);
  OGRDataSource::DestroyDataSource(DataSource);
  OGRDataSource::DestroyDataSource(PlotsLimits);

  BLOCK_EXIT(__PRETTY_FUNCTION__)
}


// ====================================================================
// ====================================================================


void LandProcessor::attributeLinearStructures()
{
  BLOCK_ENTER(__PRETTY_FUNCTION__)

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int IDPlot, IDLimit, FID, N;
  double Distance;
  std::string Options, FileName;

  std::vector <std::string> FieldNamesTab, FileNamesTab = {m_InputDitchesFile, m_InputHedgesFile, m_InputGrassBandFile, m_InputRivesrFile, m_InputBenchesFile};

  OGRDataSource *DataSource, *LS, *Plots, *PlotsLimits, *Intersection;
  OGRLayer *SQLLayer, *LSLayer, *PlotsLayer, *PlotsLimitsLayer, *IntersectionLayer;
  OGRFeature *LSFeature, *PlotsFeature, *NewFeature, *IntersectionFeature;
  OGRGeometry *NewGeometry, *StartPointGeometry, *EndPointGeometry;
  OGRFeatureDefn *FeatureDefn;

  OGRLineString *IntersectionLineString, *NewLine;
  OGRMultiLineString *IntersectionMultiLineString;
  OGRGeometryCollection *IntersectionGeometryCollection;
  OGRPoint IntersectionPoint;

  OGRSpatialReference *VectorSRS;
  const char *VectorDataEPSGCode;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  // ======================================================================================
  // Checking if necessary vector files are present in the appropriate folder
  // ======================================================================================

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsLimitsVectorFile).c_str()))
  {
    throw std::runtime_error("LandProcessor::attributeLinearStructures(): " + m_OutputPlotsLimitsVectorFile + ": no such file in the output vector directory");
  }

  for (int i = 0; i < FileNamesTab.size(); i++)
  {
    if (FileNamesTab[i].empty() or !openfluid::tools::Filesystem::isFile(getInputVectorPath(FileNamesTab[i]).c_str()))
    {
      std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": no such file in the input vector directory" << std::endl;
    }
    else
    {
      LS = OGRSFDriverRegistrar::Open( getInputVectorPath(FileNamesTab[i]).c_str(), TRUE );
      if (LS->GetDriver()->GetName() != m_VectorDriverName)
      {
        std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": vector file is not in ESRI Shapefile format" << std::endl;
        OGRDataSource::DestroyDataSource(LS);
      }
      else
      {
        LSLayer = LS->GetLayer(0);
        VectorSRS = LSLayer->GetSpatialRef();
        if (!VectorSRS)
        {
          std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": no spatial reference information is provided for the vector file" << std::endl;
          OGRDataSource::DestroyDataSource(LS);
        }
        else
        {
          VectorDataEPSGCode = VectorSRS->GetAttrValue("AUTHORITY",1);
          if (!VectorDataEPSGCode)
          {
            std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": no EPSG code provided for the vector data" << std::endl;
            OGRDataSource::DestroyDataSource(LS);
          }
          else
          {
            if (VectorDataEPSGCode != m_EPSGCode.substr(m_EPSGCode.find(":")+1, m_EPSGCode.length()-1))
            {
              std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": vector data EPSG code does not correspond to the default EPSG code" << std::endl;
              OGRDataSource::DestroyDataSource(LS);
            }
            else
            {

              m_SQLRequest = "REPACK " + FileNamesTab[i].substr(0, FileNamesTab[i].find("."));
              LS->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
              LSLayer = LS->GetLayer(0);
              LSLayer->ResetReading();

              N = 0;

              while ((LSFeature = LSLayer->GetNextFeature()) != 0)
              {
                if (LSFeature->GetGeometryRef()->IsValid() == 0 or LSFeature->GetGeometryRef()->IsSimple() == 0)
                {
                  N += 1;
                  break;
                }
              }

              if (N != 0)
              {
                std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": one or several features have invalid geometry" << std::endl;
                OGRDataSource::DestroyDataSource(LS);
              }
              else
              {
                LSLayer = LS->GetLayer(0);
                LSLayer->ResetReading();

                N = 0;

                while ((LSFeature = LSLayer->GetNextFeature()) != 0)
                {
                  if (LSFeature->GetGeometryRef()->getGeometryType() != 2)
                  {
                    N += 1;
                    break;
                  }
                }

                if (N != 0)
                {
                  std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": one or several features have geometry that are not of type 'LINESTRING'" << std::endl;
                  OGRDataSource::DestroyDataSource(LS);
                }
                else
                {

                  // ==================================================================
                  // Attribution of the original linear structures to the parcel limits
                  // ==================================================================
                  LSLayer = LS->GetLayer(0);
                  FeatureDefn = LSLayer->GetLayerDefn();
                  FieldNamesTab.clear();

                  for (int i = 0; i < FeatureDefn->GetFieldCount(); i++)
                  {
                    FieldNamesTab.push_back(FeatureDefn->GetFieldDefn(i)->GetNameRef());
                  }

                  if(std::find(FieldNamesTab.begin(), FieldNamesTab.end(), m_IDFieldName) != FieldNamesTab.end())
                  {
                    for (int i = 0; i < LSLayer->GetFeatureCount(); i++)
                    {
                      LSFeature = LSLayer->GetFeature(i);
                      LSFeature->SetField(m_IDFieldName.c_str(), i+1);
                      LSLayer->SetFeature(LSFeature);
                    }
                  }
                  else
                  {
                    createField(LSLayer, m_IDFieldName.c_str(), OFTInteger);
                    for (int i = 0; i < LSLayer->GetFeatureCount(); i++)
                    {
                      LSFeature = LSLayer->GetFeature(i);
                      LSFeature->SetField(m_IDFieldName.c_str(), i+1);
                      LSLayer->SetFeature(LSFeature);
                    }
                  }

                  FieldNamesTab.clear();
                  LSLayer->SyncToDisk();

                  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath("linearstructure.shp").c_str()))
                  {
                    mp_VectorDriver->DeleteDataSource(getOutputVectorPath("linearstructure.shp").c_str());
                  }

                  DataSource = mp_VectorDriver->CopyDataSource(LS, getOutputVectorPath("linearstructure.shp").c_str(), nullptr);
                  OGRDataSource::DestroyDataSource(LS);
                  OGRDataSource::DestroyDataSource(DataSource);
                  LS = OGRSFDriverRegistrar::Open( getOutputVectorPath("linearstructure.shp").c_str(), TRUE );
                  mp_SRS = LS->GetLayer(0)->GetSpatialRef();
                  LSLayer = LS->GetLayer(0);
                  FeatureDefn = LSLayer->GetLayerDefn();

                  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath("intersection.shp").c_str()))
                  {
                    mp_VectorDriver->DeleteDataSource(getOutputVectorPath("intersection.shp").c_str());
                  }

                  Intersection = mp_VectorDriver->CreateDataSource(getOutputVectorPath("intersection.shp").c_str());
                  IntersectionLayer = Intersection->CreateLayer("intersection", mp_SRS, wkbLineString);
                  FieldNamesTab.clear();
                  FieldNamesTab = {m_IDFieldName, "IDLS", "IDPlot", "IDLimit", "PointDist"};

                  for (int i = 0; i < FieldNamesTab.size(); i++)
                  {
                    if (i < 4)
                    {
                      createField(IntersectionLayer, FieldNamesTab[i].c_str(), OFTInteger);
                    }
                    else
                    {
                      createField(IntersectionLayer, FieldNamesTab[i].c_str(), OFTReal);
                    }
                  }

                  FieldNamesTab.clear();
                  IntersectionLayer->SyncToDisk();
                  LSLayer = LS->GetLayer(0);
                  Plots = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsVectorFile).c_str(), TRUE );
                  PlotsLayer = Plots->GetLayer(0);

                  for (int i = 0; i < LSLayer->GetFeatureCount(); i++)
                  {
                    for (int j = 0; j < PlotsLayer->GetFeatureCount(); j++)
                    {
                      if (LSLayer->GetFeature(i)->GetGeometryRef()->Intersection(PlotsLayer->GetFeature(j)->GetGeometryRef())->getGeometryType() == 5) //if intersection is a multi-line
                      {
                        IntersectionMultiLineString = (OGRMultiLineString*) LSLayer->GetFeature(i)->GetGeometryRef()->Intersection(PlotsLayer->GetFeature(j)->GetGeometryRef());
                        for (int k = 0; k < IntersectionMultiLineString->getNumGeometries(); k++)
                        {
                          IntersectionLineString = (OGRLineString*) IntersectionMultiLineString->getGeometryRef(k);
                          NewFeature = OGRFeature::CreateFeature(IntersectionLayer->GetLayerDefn());
                          NewFeature->SetField(m_IDFieldName.c_str(), IntersectionLayer->GetFeatureCount());
                          NewFeature->SetField("IDLS", LSLayer->GetFeature(i)->GetFieldAsInteger(m_IDFieldName.c_str()));
                          NewFeature->SetField("IDPlot", PlotsLayer->GetFeature(j)->GetFieldAsInteger(m_IDFieldName.c_str()));
                          NewFeature->SetGeometry(IntersectionLineString);
                          IntersectionLayer->CreateFeature(NewFeature);
                          IntersectionLayer->SyncToDisk();
                        }
                      }
                      else if (LSLayer->GetFeature(i)->GetGeometryRef()->Intersection(PlotsLayer->GetFeature(j)->GetGeometryRef())->getGeometryType() == 2)
                      {
                        IntersectionLineString = (OGRLineString*) LSLayer->GetFeature(i)->GetGeometryRef()->Intersection(PlotsLayer->GetFeature(j)->GetGeometryRef());
                        NewFeature = OGRFeature::CreateFeature(IntersectionLayer->GetLayerDefn());
                        NewFeature->SetField(m_IDFieldName.c_str(), IntersectionLayer->GetFeatureCount());
                        NewFeature->SetField("IDLS", LSLayer->GetFeature(i)->GetFieldAsInteger(m_IDFieldName.c_str()));
                        NewFeature->SetField("IDPlot", PlotsLayer->GetFeature(j)->GetFieldAsInteger(m_IDFieldName.c_str()));
                        NewFeature->SetGeometry(IntersectionLineString);
                        IntersectionLayer->CreateFeature(NewFeature);
                        IntersectionLayer->SyncToDisk();
                      }
                      else if (LSLayer->GetFeature(i)->GetGeometryRef()->Intersection(PlotsLayer->GetFeature(j)->GetGeometryRef())->getGeometryType() == 7)
                      {
                        IntersectionGeometryCollection = (OGRGeometryCollection*) LSLayer->GetFeature(i)->GetGeometryRef()->Intersection(PlotsLayer->GetFeature(j)->GetGeometryRef());
                        if (IntersectionGeometryCollection->getNumGeometries() != 0)
                        {
                          for (int k = 0; k < IntersectionGeometryCollection->getNumGeometries(); k++)
                          {
                            if (IntersectionGeometryCollection->getGeometryRef(k)->getGeometryType() == 2)
                            {
                              IntersectionLineString = (OGRLineString*) IntersectionGeometryCollection->getGeometryRef(k);
                              NewFeature = OGRFeature::CreateFeature(IntersectionLayer->GetLayerDefn());
                              NewFeature->SetField(m_IDFieldName.c_str(), IntersectionLayer->GetFeatureCount());
                              NewFeature->SetField("IDLS", LSLayer->GetFeature(i)->GetFieldAsInteger(m_IDFieldName.c_str()));
                              NewFeature->SetField("IDPlot", PlotsLayer->GetFeature(j)->GetFieldAsInteger(m_IDFieldName.c_str()));
                              NewFeature->SetGeometry(IntersectionLineString);
                              IntersectionLayer->CreateFeature(NewFeature);
                              IntersectionLayer->SyncToDisk();
                            }
                            else if (IntersectionGeometryCollection->getGeometryRef(k)->getGeometryType() == 5)
                            {
                              IntersectionMultiLineString = (OGRMultiLineString*) IntersectionGeometryCollection->getGeometryRef(k);
                              for (int l = 0; l < IntersectionMultiLineString->getNumGeometries(); l++)
                              {
                                IntersectionLineString = (OGRLineString*) IntersectionMultiLineString->getGeometryRef(l);
                                NewFeature = OGRFeature::CreateFeature(IntersectionLayer->GetLayerDefn());
                                NewFeature->SetField(m_IDFieldName.c_str(), IntersectionLayer->GetFeatureCount());
                                NewFeature->SetField("IDLS", LSLayer->GetFeature(i)->GetFieldAsInteger(m_IDFieldName.c_str()));
                                NewFeature->SetField("IDPlot", PlotsLayer->GetFeature(j)->GetFieldAsInteger(m_IDFieldName.c_str()));
                                NewFeature->SetGeometry(IntersectionLineString);
                                IntersectionLayer->CreateFeature(NewFeature);
                                IntersectionLayer->SyncToDisk();
                              }
                            }
                          }
                        }
                      }
                    }
                  }

                  OGRDataSource::DestroyDataSource(Intersection);
                  Intersection = OGRSFDriverRegistrar::Open( getOutputVectorPath("intersection.shp").c_str(), TRUE );
                  IntersectionLayer = Intersection->GetLayer(0);
                  IntersectionLayer->ResetReading();
                  PlotsLimits = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsLimitsVectorFile).c_str(), TRUE );
                  PlotsLimitsLayer = PlotsLimits->GetLayer(0);
                  DataSource =  OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );

                  while ((IntersectionFeature = IntersectionLayer->GetNextFeature()) != nullptr)
                  {
                    FID = IntersectionFeature->GetFID();
                    IDPlot = IntersectionFeature->GetFieldAsInteger("IDPlot");
                    m_SQLRequest = "SELECT *, geometry, ROWID AS RID"
                        " FROM " + m_OutputPlotsLimitsVectorFile.substr(0, m_OutputPlotsLimitsVectorFile.find(".")) + " WHERE IDPlotA == " + std::to_string(IDPlot) + " OR IDPlotB == " + std::to_string(IDPlot) + " ORDER BY "
                        "ST_Distance(geometry,(SELECT ST_Line_Interpolate_Point(ST_LineMerge(geometry), 0.5) FROM intersection WHERE ROWID == " + std::to_string(FID) + " )) ASC, "
                        "ST_Distance(geometry,(SELECT ST_EndPoint(geometry) FROM intersection WHERE ROWID == " + std::to_string(FID) + " )) ASC, "
                        "ST_Distance(geometry,(SELECT ST_StartPoint(geometry) FROM intersection WHERE ROWID == " + std::to_string(FID) + ")) ASC";
                    SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
                    IDLimit = SQLLayer->GetFeature(0)->GetFieldAsInteger("RID");
                    IntersectionFeature->SetField("IDLimit", SQLLayer->GetFeature(0)->GetFieldAsInteger("ID"));
                    IntersectionLayer->SetFeature(IntersectionFeature);
                    DataSource->ReleaseResultSet(SQLLayer);
                    IntersectionLineString = (OGRLineString*) IntersectionFeature->GetGeometryRef();
                    IntersectionLineString->segmentize(0.5);
                    Distance = 0;
                    for (int i = 0; i < IntersectionLineString->getNumPoints(); i++)
                    {
                      IntersectionLineString->getPoint(i,&IntersectionPoint);
                      if (PlotsLimitsLayer->GetFeature(IDLimit)->GetGeometryRef()->Distance(&IntersectionPoint) > Distance)
                      {
                        Distance = PlotsLimitsLayer->GetFeature(IDLimit)->GetGeometryRef()->Distance(&IntersectionPoint);
                      }
                    }
                    IntersectionFeature->SetField("PointDist", Distance);
                    IntersectionLayer->SetFeature(IntersectionFeature);
                  }

                  OGRDataSource::DestroyDataSource(DataSource);
                  OGRDataSource::DestroyDataSource(Intersection);
                  OGRDataSource::DestroyDataSource(PlotsLimits);

                  Intersection = OGRSFDriverRegistrar::Open( getOutputVectorPath("intersection.shp").c_str(), TRUE );
                  DataSource =  OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
                  IntersectionLayer = Intersection->GetLayer(0);

                  for (int i = 0; i < IntersectionLayer->GetFeatureCount(); i++)
                  {
                    IntersectionFeature = IntersectionLayer->GetFeature(i);
                    FID = IntersectionFeature->GetFID();
                    IDLimit = IntersectionFeature->GetFieldAsInteger("IDLimit");
                    Distance = IntersectionFeature->GetFieldAsDouble("PointDist");
                    if (Distance > 5.00000000e-07)
                    {
                      m_SQLRequest = "SELECT ST_CollectionExtract(ST_Intersection(geometry, "
                          "(SELECT ST_Buffer(geometry, PointDist*1.5) FROM intersection WHERE ROWID == " + std::to_string(FID) + ")), 2) "
                          "FROM " + m_OutputPlotsLimitsVectorFile.substr(0, m_OutputPlotsLimitsVectorFile.find(".")) + " WHERE ID == " + std::to_string(IDLimit);
                      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
                      IntersectionFeature->SetGeometryDirectly(SQLLayer->GetFeature(0)->GetGeometryRef());
                      IntersectionLayer->SetFeature(IntersectionFeature);
                      IntersectionLayer->SyncToDisk();
                      DataSource->ReleaseResultSet(SQLLayer);
                    }
                  }

                  IntersectionLayer->SyncToDisk();
                  IntersectionLayer->ResetReading();

                  while ((IntersectionFeature = IntersectionLayer->GetNextFeature()) != nullptr)
                  {
                    FID = IntersectionFeature->GetFID();
                    Distance = IntersectionFeature->GetFieldAsDouble("PointDist");
                    IntersectionLineString = (OGRLineString*) IntersectionFeature->GetGeometryRef();
                    if (IntersectionLineString->get_Length() <= Distance*1.5*3)
                    {
                      IntersectionLayer->DeleteFeature(FID);
                    }
                  }

                  m_SQLRequest = "REPACK intersection";
                  Intersection->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
                  OGRDataSource::DestroyDataSource(Intersection);

                  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(FileNamesTab[i]).c_str()))
                  {
                    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(FileNamesTab[i]).c_str());
                  }

                  Intersection = OGRSFDriverRegistrar::Open( getOutputVectorPath("intersection.shp").c_str(), TRUE );

                  if (Intersection->GetLayer(0)->GetFeatureCount() != 0)
                  {
                    m_SQLRequest = "SELECT ST_Union(geometry), " + m_IDFieldName + ", IDLS FROM intersection GROUP BY " + m_IDFieldName + ", IDLS, IDLimit"; // 15/12/16/14:34 Dobavila m_IDFieldName
                    SQLLayer = Intersection->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
                    FileName = FileNamesTab[i].substr(0, FileNamesTab[i].find(".shp"));
                    Intersection->CopyLayer(SQLLayer, FileName.c_str(), nullptr);
                    Intersection->ReleaseResultSet(SQLLayer);
                  }

                  OGRDataSource::DestroyDataSource(Intersection);
                  OGRDataSource::DestroyDataSource(LS);
                  OGRDataSource::DestroyDataSource(Plots);
                  OGRDataSource::DestroyDataSource(DataSource);

                  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath("linearstructure.shp").c_str()))
                  {
                    mp_VectorDriver->DeleteDataSource(getOutputVectorPath("linearstructure.shp").c_str());
                  }

                  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath("intersection.shp").c_str()))
                  {
                    mp_VectorDriver->DeleteDataSource(getOutputVectorPath("intersection.shp").c_str());
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  BLOCK_EXIT(__PRETTY_FUNCTION__)
}


// =====================================================================
// =====================================================================


int LandProcessor::extractFromRasterToPoint(GDALDataset *Dataset, unsigned int RasterBandIndex)
{

  int Value;
  int *ScanValue = (int*) CPLMalloc(sizeof(int));

  GDALRasterIO(getRasterBand(Dataset, RasterBandIndex), GF_Read,
               calculateOffset(getCoordinatesOfPoint().first,getCoordinatesOfPoint().second).first,
               calculateOffset(getCoordinatesOfPoint().first,getCoordinatesOfPoint().second).second,
               1, 1,
               ScanValue,
               1, 1, GDT_Int32,
               0, 0);

  Value = ScanValue[0];
  CPLFree(ScanValue);
  return Value;

}


// =====================================================================
// =====================================================================


std::pair <double, double> LandProcessor::getCoordinatesOfPoint()
{

  std::pair <double, double> Coord;
  Coord = std::make_pair(mp_Point->getX(), mp_Point->getY());
  return Coord;

}


// =====================================================================
// =====================================================================


std::pair <double, double> LandProcessor::calculateOffset(double XCoord, double YCoord)
{

  std::pair <double, double> Offset;

  Offset = std::make_pair(
      (XCoord - mp_GeoTransformVal[0])/mp_GeoTransformVal[1],
      (YCoord - mp_GeoTransformVal[3])/mp_GeoTransformVal[5]);

  return Offset;

}


// =====================================================================
// =====================================================================


void LandProcessor::getCentroidPoint(OGRGeometry *Geometry)
{

  mp_Point = new OGRPoint;
  Geometry->Centroid(mp_Point);

}


// =====================================================================
// =====================================================================


void LandProcessor::createField(OGRLayer *LayerName, std::string FieldName, OGRFieldType FieldType)
{

  OGRFieldDefn FieldDefn(FieldName.c_str(), FieldType);
  LayerName->CreateField(&FieldDefn);

}


// =====================================================================
// =====================================================================


GDALRasterBand* LandProcessor::getRasterBand(GDALDataset *Dataset, unsigned int RastreBandIndex)
{

  return Dataset->GetRasterBand(RastreBandIndex);

}

// =====================================================================
// =====================================================================

void LandProcessor::getGeoTransform(GDALDataset *Dataset)
{

  mp_GeoTransformVal = new double[6];
  GDALGetGeoTransform(Dataset, mp_GeoTransformVal);

}


// =====================================================================
// =====================================================================


void LandProcessor::checkLinearStructuresVectorData()
{

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  std::vector <std::string> FieldNamesTab, FileNamesTab = {m_InputDitchesFile, m_InputHedgesFile, m_InputGrassBandFile, m_InputRivesrFile, m_InputBenchesFile, m_InputThalwegsFile};

  OGRDataSource *LS;
  OGRLayer *LSLayer;
  OGRFeature *LSFeature;

  OGRSpatialReference *VectorSRS;
  const char *VectorDataEPSGCode;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  for (int i = 0; i < FileNamesTab.size(); i++)
  {
    if (FileNamesTab[i].empty() or !openfluid::tools::Filesystem::isFile(getInputVectorPath(FileNamesTab[i]).c_str()))
    {
      std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": no such file in the input vector directory" << std::endl;
    }
    else
    {
      LS = OGRSFDriverRegistrar::Open( getInputVectorPath(FileNamesTab[i]).c_str(), TRUE );
      if (LS->GetDriver()->GetName() != m_VectorDriverName)
      {
        std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": vector file is not in ESRI Shapefile format" << std::endl;
        OGRDataSource::DestroyDataSource(LS);
      }
      else
      {
        LSLayer = LS->GetLayer(0);
        VectorSRS = LSLayer->GetSpatialRef();
        if (!VectorSRS)
        {
          std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": no spatial reference information is provided for the vector file" << std::endl;
          OGRDataSource::DestroyDataSource(LS);
        }
        else
        {
          VectorDataEPSGCode = VectorSRS->GetAttrValue("AUTHORITY",1);
          if (!VectorDataEPSGCode)
          {
            std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": no EPSG code provided for the vector data" << std::endl;
            OGRDataSource::DestroyDataSource(LS);
          }
          else
          {
            if (VectorDataEPSGCode != m_EPSGCode.substr(m_EPSGCode.find(":")+1, m_EPSGCode.length()-1))
            {
              std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": vector data EPSG code does not correspond to the default EPSG code" << std::endl;
              OGRDataSource::DestroyDataSource(LS);
            }
            else
            {

              m_SQLRequest = "REPACK " + FileNamesTab[i].substr(0, FileNamesTab[i].find("."));
              LS->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
              LSLayer = LS->GetLayer(0);
              LSLayer->ResetReading();

              while ((LSFeature = LSLayer->GetNextFeature()) != 0)
              {
                if (LSFeature->GetGeometryRef()->IsValid() == 0)
                {
                  std::cout << "The feature with FID " << LSFeature->GetFID() << " of the vector file " << FileNamesTab[i] << " has invalid geometry" << std::endl;
                }
                if (LSFeature->GetGeometryRef()->IsSimple() == 0)
                {
                  std::cout << "The feature with FID " << LSFeature->GetFID() << " of the vector file " << FileNamesTab[i] << " has geometry that is not simple" << std::endl;
                }
                if (LSFeature->GetGeometryRef()->IsEmpty() == 1)
                {
                  std::cout << "The feature with FID " << LSFeature->GetFID() << " of the vector file " << FileNamesTab[i] << " has empty geometry" << std::endl;
                }
              }

              LSLayer = LS->GetLayer(0);
              LSLayer->ResetReading();

              while ((LSFeature = LSLayer->GetNextFeature()) != 0)
              {
                if (LSFeature->GetGeometryRef()->getGeometryType() != 2)
                {
                  std::cout << "The feature with FID " << LSFeature->GetFID() << " of the vector file " << FileNamesTab[i] << " has geometry that is not of LINESTRING type" << std::endl;
                }
              }
            }
          }
        }
      }
    }
  }
}


// =====================================================================
// =====================================================================


void LandProcessor::checkPolygonVectorData()
{

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  std::vector <std::string> FieldNamesTab, FileNamesTab = {m_InputPlotsFile};

  OGRDataSource *Polygon;
  OGRLayer *PolygonLayer;
  OGRFeature *PolygonFeature;

  OGRSpatialReference *VectorSRS;
  const char *VectorDataEPSGCode;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  for (int i = 0; i < FileNamesTab.size(); i++)
  {
    if (FileNamesTab[i].empty() or !openfluid::tools::Filesystem::isFile(getInputVectorPath(FileNamesTab[i]).c_str()))
    {
      std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": no such file in the input vector directory" << std::endl;
    }
    else
    {
      Polygon = OGRSFDriverRegistrar::Open( getInputVectorPath(FileNamesTab[i]).c_str(), TRUE );
      if (Polygon->GetDriver()->GetName() != m_VectorDriverName)
      {
        std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": vector file is not in ESRI Shapefile format" << std::endl;
        OGRDataSource::DestroyDataSource(Polygon);
      }
      else
      {
        PolygonLayer = Polygon->GetLayer(0);
        VectorSRS = PolygonLayer->GetSpatialRef();
        if (!VectorSRS)
        {
          std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": no spatial reference information is provided for the vector file" << std::endl;
          OGRDataSource::DestroyDataSource(Polygon);
        }
        else
        {
          VectorDataEPSGCode = VectorSRS->GetAttrValue("AUTHORITY",1);
          if (!VectorDataEPSGCode)
          {
            std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": no EPSG code provided for the vector data" << std::endl;
            OGRDataSource::DestroyDataSource(Polygon);
          }
          else
          {
            if (VectorDataEPSGCode != m_EPSGCode.substr(m_EPSGCode.find(":")+1, m_EPSGCode.length()-1))
            {
              std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTab[i] << ": vector data EPSG code does not correspond to the default EPSG code" << std::endl;
              OGRDataSource::DestroyDataSource(Polygon);
            }
            else
            {

              m_SQLRequest = "REPACK " + FileNamesTab[i].substr(0, FileNamesTab[i].find("."));
              Polygon->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
              PolygonLayer = Polygon->GetLayer(0);
              PolygonLayer->ResetReading();

              while ((PolygonFeature = PolygonLayer->GetNextFeature()) != 0)
              {
                if (PolygonFeature->GetGeometryRef()->IsValid() == 0)
                {
                  std::cout << "The feature with FID " << PolygonFeature->GetFID() << " of the vector file " << FileNamesTab[i] << " has invalid geometry" << std::endl;
                }
                if (PolygonFeature->GetGeometryRef()->Boundary()->IsSimple() == 0)
                {
                  std::cout << "The feature with FID " << PolygonFeature->GetFID() << " of the vector file " << FileNamesTab[i] << " has geometry that is not simple" << std::endl;
                }
                if (PolygonFeature->GetGeometryRef()->IsEmpty() == 1)
                {
                  std::cout << "The feature with FID " << PolygonFeature->GetFID() << " of the vector file " << FileNamesTab[i] << " has empty geometry" << std::endl;
                }
              }

              PolygonLayer = Polygon->GetLayer(0);
              PolygonLayer->ResetReading();

              while ((PolygonFeature = PolygonLayer->GetNextFeature()) != 0)
              {
                if (PolygonFeature->GetGeometryRef()->getGeometryType() != 3)
                {
                  std::cout << "The feature with FID " << PolygonFeature->GetFID() << " of the vector file " << FileNamesTab[i] << " has  geometry that is not of POLYGON type" << std::endl;
                }
              }
            }
          }
        }
      }
    }
  }
}


