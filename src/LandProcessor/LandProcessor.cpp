/*
    LandProcessor
    Copyright (C) 2016- INRA

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


/**
  @file LandProcessor.cpp
  @author Ekaterina Zadonina <ekaterina.zadonina@inra.fr>
  @author Jean-Christophe FABRE <jean-christophe.fabre@inra.fr>
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

#include <LandProcessor/LandProcessor.hpp>
#include <LandProcessor/Helpers.hpp>


LandProcessor::LandProcessor(const std::string& InputPath,
                             const std::string& OutputPath,
                             const std::string& ReleasePath):
  m_InputPath(InputPath), m_OutputPath(OutputPath), m_ReleasePath(ReleasePath)
{
  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);


  openfluid::tools::Filesystem::makeDirectory(OutputPath);
  openfluid::tools::Filesystem::makeDirectory(getOutputVectorPath());
  openfluid::tools::Filesystem::makeDirectory(getOutputRasterPath());

  openfluid::tools::Filesystem::makeDirectory(ReleasePath);
  openfluid::tools::Filesystem::makeDirectory(getReleaseVectorPath());
  openfluid::tools::Filesystem::makeDirectory(getReleaseRasterPath());

  m_IsReady = m_IsReady && openfluid::tools::Filesystem::isDirectory(getInputVectorPath());
  m_IsReady = m_IsReady && openfluid::tools::Filesystem::isDirectory(getInputRasterPath());

  openfluid::tools::Filesystem::makeDirectory(openfluid::base::Environment::getTempDir());

  openfluid::tools::Filesystem::makeDirectory("/tmp/landprocessor-grass");
  openfluid::utils::GrassGISProxy GRASS("/tmp/landprocessor-grass","temp");

  GRASS.createLocation(m_EPSGCode.c_str());

  OGRRegisterAll();
  GDALAllRegister();
}


// =====================================================================
// =====================================================================


LandProcessor::~LandProcessor()
{
  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);
}


// =====================================================================
// =====================================================================


std::string LandProcessor::getLayerNameFromFilename(const std::string& Filename)
{
  return Filename.substr(0,Filename.find("."));
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


std::string LandProcessor::getInputRasterPath(const std::string& Filename) const
{
  std::string Path = m_InputPath+"/"+m_RasterDir;

  if (!Filename.empty())
    Path += "/"+Filename;

  return Path;

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


std::string LandProcessor::getOutputRasterPath(const std::string& Filename) const
{
  std::string Path = m_OutputPath+"/"+m_RasterDir;

  if (!Filename.empty())
    Path += "/"+Filename;

  return Path;

}


// =====================================================================
// =====================================================================


std::string LandProcessor::getReleaseVectorPath(const std::string& Filename) const
{
  std::string Path = m_ReleasePath+"/"+m_VectorDir;

  if (!Filename.empty())
    Path += "/"+Filename;

  return Path;
}


// =====================================================================
// =====================================================================


std::string LandProcessor::getReleaseRasterPath(const std::string& Filename) const
{
  std::string Path = m_ReleasePath+"/"+m_RasterDir;

  if (!Filename.empty())
    Path += "/"+Filename;

  return Path;
}


// =====================================================================
// =====================================================================


/**
  @internal
  <hr>
  Files
    - m_InputPlotsVectorFile
    - m_InputDEMFile
    - m_OutputPlotsVectorFile
*/
void LandProcessor::preprocessVectorData()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  std::string FIDAB, FIDBA;

  int FIDA, FIDB;

  std::vector<std::string> FieldNamesTable, ClipTable, ClipedTable;
  OGRDataSource *DataSource, *Plots;
  OGRLayer *SQLLayer, *PlotsLayer;
  OGRFeature *PlotsFeature;
  OGRGeometry *PlotsGeometry, *NewGeometry;
  OGRFeatureDefn *FeatureDefn;

  GDALDataset *DEM;

  OGRSpatialReference *VectorSRS, RasterSRS;
  const char *VectorDataEPSGCode, *RasterDataEPSGCode;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  if (!openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputPlotsVectorFile)))
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsVectorFile + ": there is no such file input vector directory");
  }

  Plots = OGRSFDriverRegistrar::Open( getInputVectorPath(m_InputPlotsVectorFile).c_str(), TRUE );

  if (Plots->GetDriver()->GetName() != m_VectorDriverName)
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsVectorFile + ": vector file is not in ESRI Shapefile format");
  }

  PlotsLayer = Plots->GetLayer(0);
  VectorSRS = PlotsLayer->GetSpatialRef();

  if (!VectorSRS)
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsVectorFile + ": no spatial reference information is provided for the vector file");
  }

  VectorDataEPSGCode = VectorSRS->GetAttrValue("AUTHORITY",1);

  if (!VectorDataEPSGCode)
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsVectorFile + ": no EPSG code provided for the vector data");
  }

  if (VectorDataEPSGCode != m_EPSGCode.substr(m_EPSGCode.find(":")+1, m_EPSGCode.length()-1))
  {
    throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsVectorFile + ":  vector data EPSG code does not correspond to the default EPSG code");
  }

  if (!openfluid::tools::Filesystem::isFile(getInputRasterPath(m_InputDEMFile)))
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

  const std::string InputPlotsLayerName = getLayerNameFromFilename(m_InputPlotsVectorFile);

  m_SQLRequest = "REPACK " + InputPlotsLayerName;
  Plots->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
  PlotsLayer = Plots->GetLayer(0);
  PlotsLayer->ResetReading();

  while ((PlotsFeature = PlotsLayer->GetNextFeature()) != 0)
  {
    if (PlotsFeature->GetGeometryRef()->IsValid() == 0)
    {
      throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsVectorFile + ": one or several features have invalid geometry");
    }
  }

  PlotsLayer = Plots->GetLayer(0);
  PlotsLayer->ResetReading();

  while ((PlotsFeature = PlotsLayer->GetNextFeature()) != 0)
  {
    if (PlotsFeature->GetGeometryRef()->getGeometryType() != 3)
    {
      throw std::runtime_error("LandProcessor::preprocessVectorData(): " + m_InputPlotsVectorFile + ": one or several features have geometry that are not of type 'POLYGON'");
    }
  }

  // =====================================================================
  // Creation or calculation of the unique ID field
  // =====================================================================

  PlotsLayer = Plots->GetLayer(0);
  PlotsLayer->ResetReading();
  FeatureDefn = PlotsLayer->GetLayerDefn();

  FieldNamesTable.clear();

  for (int i = 0; i < FeatureDefn->GetFieldCount(); i++)
  {
    FieldNamesTable.push_back(PlotsLayer->GetLayerDefn()->GetFieldDefn(i)->GetNameRef());
  }

  if(std::find(FieldNamesTable.begin(), FieldNamesTable.end(), m_IDFieldName) != FieldNamesTable.end())
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

  FieldNamesTable.clear();
  OGRDataSource::DestroyDataSource(Plots);

  // =====================================================================
  // Correcting overlaps and gaps between polygons, if at all possible
  // =====================================================================

  Plots = OGRSFDriverRegistrar::Open( getInputVectorPath(m_InputPlotsVectorFile).c_str(), TRUE );

  PlotsLayer = Plots->GetLayer(0);
  PlotsLayer->ResetReading();

  DataSource = OGRSFDriverRegistrar::Open( getInputVectorPath().c_str(), TRUE );

  ClipTable.push_back("0");

  while(ClipTable.size() != ClipedTable.size())
  {
    ClipTable.clear();
    PlotsLayer = Plots->GetLayer(0);
    PlotsLayer->ResetReading();
    while ((PlotsFeature = PlotsLayer->GetNextFeature()) != nullptr)
    {
      PlotsGeometry = PlotsFeature->GetGeometryRef();
      DataSource = OGRSFDriverRegistrar::Open( getInputVectorPath().c_str(), TRUE );
      FIDA = PlotsFeature->GetFID();
      m_SQLRequest = "SELECT *, ROWID AS RID FROM " + InputPlotsLayerName + " WHERE "
          "ST_Overlaps(geometry, (SELECT geometry FROM " + InputPlotsLayerName + " WHERE ROWID == " + std::to_string(PlotsFeature->GetFID()) + ")) AND ROWID != " + std::to_string(PlotsFeature->GetFID()) + " AND "
          " ST_GeometryType(ST_Intersection(geometry, (SELECT geometry FROM " + InputPlotsLayerName + " WHERE ROWID == " + std::to_string(PlotsFeature->GetFID()) + "))) != " + m_QMark + "LINESTRING" + m_QMark + " AND "
          "ST_GeometryType(ST_Intersection(geometry, (SELECT geometry FROM " + InputPlotsLayerName + " WHERE ROWID == " + std::to_string(PlotsFeature->GetFID()) + "))) != " + m_QMark + "MULTILINESTRING" + m_QMark;
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        NewGeometry = nullptr;
        for (int i = 0; i < SQLLayer->GetFeatureCount(); i++)
        {
          FIDB = SQLLayer->GetFeature(i)->GetFieldAsInteger("RID");
          FIDBA = std::to_string(FIDB) + "/" + std::to_string(FIDA);
          if (std::find(ClipTable.begin(), ClipTable.end(), FIDBA) == ClipTable.end())
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
            ClipTable.push_back(FIDAB);
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
    ClipedTable.clear();
    while ((PlotsFeature = PlotsLayer->GetNextFeature()) != nullptr)
    {
      PlotsGeometry = PlotsFeature->GetGeometryRef();
      DataSource = OGRSFDriverRegistrar::Open( getInputVectorPath().c_str(), TRUE );
      FIDA = PlotsFeature->GetFID();
      m_SQLRequest = "SELECT *, ROWID AS RID FROM " + InputPlotsLayerName + " WHERE "
          "ST_Overlaps(geometry, (SELECT geometry FROM " + InputPlotsLayerName + " WHERE ROWID == " + std::to_string(PlotsFeature->GetFID()) + ")) AND ROWID != " + std::to_string(PlotsFeature->GetFID()) + " AND "
          " ST_GeometryType(ST_Intersection(geometry, (SELECT geometry FROM " + InputPlotsLayerName + " WHERE ROWID == " + std::to_string(PlotsFeature->GetFID()) + "))) != " + m_QMark + "LINESTRING" + m_QMark + " AND "
          "ST_GeometryType(ST_Intersection(geometry, (SELECT geometry FROM " + InputPlotsLayerName + " WHERE ROWID == " + std::to_string(PlotsFeature->GetFID()) + "))) != " + m_QMark + "MULTILINESTRING" + m_QMark;
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        for (int i = 0; i < SQLLayer->GetFeatureCount(); i++)
        {
          FIDB = SQLLayer->GetFeature(i)->GetFieldAsInteger("RID");
          FIDBA = std::to_string(FIDB) + "/" + std::to_string(FIDA);
          if (std::find(ClipedTable.begin(), ClipedTable.end(), FIDBA) == ClipedTable.end())
          {
            FIDAB = std::to_string(FIDA) + "/" + std::to_string(FIDB);
            ClipedTable.push_back(FIDAB);
          }
        }
      }
      DataSource->ReleaseResultSet(SQLLayer);
      OGRDataSource::DestroyDataSource(DataSource);
    }
  }

  ClipTable.clear();
  ClipedTable.clear();
  OGRDataSource::DestroyDataSource(Plots);

  openfluid::utils::GrassGISProxy GRASS("/tmp/bvservice-grass","temp");
  GRASS.setOutputFile("/tmp/bvservice-grass/procesvectordata.out");
  GRASS.setErrorFile("/tmp/bvservice-grass/processvectordata.err");

  GRASS.appendTask("v.in.ogr", {{"input", QString::fromStdString(getInputVectorPath(m_InputPlotsVectorFile))},
                                {"output", "plotssnapped"}, {"snap", QString::fromStdString(m_SnapDistance)}},
                   {"--o"});

  GRASS.appendTask("r.in.gdal",{{"input",QString::fromStdString(getInputRasterPath(m_InputDEMFile))},
                                {"output","dem"}},
                   {"--o"});

  GRASS.appendTask("g.region", {{"raster", "dem"}});

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsVectorFile)))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputPlotsVectorFile).c_str());
  }

  GRASS.appendTask("v.out.ogr", {{"input", "plotssnapped"}, {"type", "area"},
                                 {"output", QString::fromStdString(getOutputVectorPath(m_OutputPlotsVectorFile))}},
                   {"-s"});

  if (GRASS.runJob() != 0)
  {
    exit(-1);
  }

  // =====================================================================
  // Removing attributes that will not be used in the processing
  // =====================================================================

  Plots = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsVectorFile).c_str(), TRUE );
  PlotsLayer = Plots->GetLayer(0);
  FeatureDefn = PlotsLayer->GetLayerDefn();

  FieldNamesTable.clear();

  for (int i = 0; i < FeatureDefn->GetFieldCount(); i++)
  {
    if (FeatureDefn->GetFieldDefn(i)->GetNameRef() != m_IDFieldName)
    {
      FieldNamesTable.push_back(FeatureDefn->GetFieldDefn(i)->GetNameRef());
    }
  }

  if (FieldNamesTable.size() > 0)
  {
    for (int i = 0; i < FieldNamesTable.size(); i++)
    {
      int FieldIndex = FeatureDefn->GetFieldIndex(FieldNamesTable.at(i).c_str());
      PlotsLayer->DeleteField(FieldIndex);
    }
  }

  FieldNamesTable.clear();

  OGRDataSource::DestroyDataSource(Plots);
  GDALClose( DEM );

}


// =====================================================================
// =====================================================================


/**
  @internal
  <hr>
  Files
    - m_InputDEMFile
    - m_OutputPlotsVectorFile
    - m_OutputPlotsRasterVectorFile
    - m_OutputDEMFile
    - m_OutputCatchmentsRasterFile
    - m_OutputDownslopeRasterFile
    - m_OutputPlotsRasterFile
    - m_OutputSlopeRasterFile
    - m_OutputDrainageRasterFile
    - m_OutputIDRasterFile
    - m_OutputOutletsRasterFile
    - m_OutputReceiversRasterFile
*/
void LandProcessor::preprocessRasterData()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

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

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsVectorFile)))
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

  UnionGeometry->Buffer(m_GeoTransformVal[1])->getEnvelope(Envelope);

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
                                {"nsres", QString::fromStdString(std::to_string(m_GeoTransformVal[1]))},
                                {"ewres", QString::fromStdString(std::to_string(m_GeoTransformVal[1]))}});

  delete Envelope;

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

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsRasterVectorFile)))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputPlotsRasterVectorFile).c_str());
  }

  GRASS.appendTask("v.out.ogr", {{"input", QString::fromStdString("plotsrastervector")},
                                 {"type", QString::fromStdString("area")},
                                 {"output", QString::fromStdString(getOutputVectorPath(m_OutputPlotsRasterVectorFile))},
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



  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputDEMFile)))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputDEMFile).c_str());
  }

  GRASS.appendTask("r.out.gdal", {{"input",QString::fromStdString("demclean")},
                                  {"output", QString::fromStdString(getOutputRasterPath(m_OutputDEMFile))},
                                  {"format", QString::fromStdString("GTiff")}},
                   {"--o"});

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputSlopeRasterFile)))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputSlopeRasterFile).c_str());
  }

  GRASS.appendTask("r.out.gdal", {{"input",QString::fromStdString("slope")},
                                  {"output", QString::fromStdString(getOutputRasterPath(m_OutputSlopeRasterFile))},
                                  {"format", QString::fromStdString("GTiff")}},
                   {"--o"});

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputDrainageRasterFile)))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputDrainageRasterFile).c_str());
  }

  GRASS.appendTask("r.out.gdal", {{"input",QString::fromStdString("drainageraster")},
                                  {"output", QString::fromStdString(getOutputRasterPath(m_OutputDrainageRasterFile))},
                                  {"format", QString::fromStdString("GTiff")}},
                   {"--o"});

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputPlotsRasterFile)))
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
  // Creation of cell ID raster that will be used to set outlets, receivers and catchments IDs
  // =========================================================================================

  DEM = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputDEMFile).c_str(), GA_ReadOnly );
  DEMBand = DEM->GetRasterBand(1);
  getGeoTransform(DEM);

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputIDRasterFile)))
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
  RasterID->SetGeoTransform(m_GeoTransformVal);
  RasterID->SetProjection(DEM->GetProjectionRef());
  GDALClose( RasterID );

  // =====================================================================
  // Outlets raster creation
  // =====================================================================

  PlotsR = (GDALDataset *) GDALOpen(getOutputRasterPath(m_OutputPlotsRasterFile).c_str(), GA_ReadOnly );
  PlotsRBand = PlotsR->GetRasterBand(1);

  Drainage = (GDALDataset *) GDALOpen(getOutputRasterPath(m_OutputDrainageRasterFile).c_str(), GA_ReadOnly );
  DrainageBand = Drainage->GetRasterBand(1);

  RasterID = (GDALDataset *) GDALOpen(getOutputRasterPath(m_OutputIDRasterFile).c_str(), GA_ReadOnly );
  RasterIDBand = RasterID->GetRasterBand(1);

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputOutletsRasterFile)))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputOutletsRasterFile).c_str());
  }

  Outlets = mp_RasterDriver->Create(getOutputRasterPath(m_OutputOutletsRasterFile).c_str(), DEM->GetRasterXSize(), DEM->GetRasterYSize(), 1, GDT_Int32, Options );
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
  Outlets->SetGeoTransform(m_GeoTransformVal);
  Outlets->SetProjection(DEM->GetProjectionRef());
  GDALClose( Outlets );

  // =====================================================================
  // Receivers raster creation
  // =====================================================================

  Outlets = (GDALDataset *) GDALOpen( getOutputRasterPath(m_OutputOutletsRasterFile).c_str(), GA_ReadOnly );
  OutletsBand = Outlets->GetRasterBand(1);

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputReceiversRasterFile)))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputReceiversRasterFile).c_str());
  }

  Receivers = mp_RasterDriver->Create( getOutputRasterPath(m_OutputReceiversRasterFile).c_str(), DEM->GetRasterXSize(), DEM->GetRasterYSize(), 1, GDT_Int32, Options );
  ReceiversBand = Receivers->GetRasterBand(1);

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputDownslopeRasterFile)))
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
  Receivers->SetGeoTransform(m_GeoTransformVal);
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
  Downslope->SetGeoTransform(m_GeoTransformVal);
  Downslope->SetProjection(DEM->GetProjectionRef());
  GDALClose( Downslope );

  // ============================================================================================
  // Catchments raster creation
  // First, we copy the outlet points to the new raster
  // Then we fill the cells that belong to the same parcel working outwards using drainage raster
  // ============================================================================================

  if (openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputCatchmentsRasterFile)))
  {
    mp_RasterDriver->Delete(getOutputRasterPath(m_OutputCatchmentsRasterFile).c_str());
  }

  Catchments = mp_RasterDriver->Create(getOutputRasterPath(m_OutputCatchmentsRasterFile).c_str(),
                                       DEM->GetRasterXSize(), DEM->GetRasterYSize(), 1, GDT_Int32, Options );
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
        OutletsBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, OutletsMiddleRow, DEM->GetRasterXSize(),
                               1, GDT_Int32, 0, 0 );
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
    CatchmentsBand->RasterIO( GF_Write, 0, i, DEM->GetRasterXSize(), 1, Row, DEM->GetRasterXSize(),
                              1, GDT_Int32, 0, 0 );
  }

  CPLFree(OutletsMiddleRow);
  CPLFree(Row);

  CatchmentsBand->FlushCache();
  CatchmentsBand->SetNoDataValue(-9999);
  Catchments->SetGeoTransform(m_GeoTransformVal);
  Catchments->SetProjection(DEM->GetProjectionRef());
  GDALClose( Catchments );

  Catchments = (GDALDataset *) GDALOpen(getOutputRasterPath(m_OutputCatchmentsRasterFile).c_str(), GA_Update );
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
    PlotsRBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, PlotsMiddleRow, DEM->GetRasterXSize(),
                          1, GDT_Int32, 0, 0 );
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
    CatchmentsBand->RasterIO( GF_Read, 0, i, DEM->GetRasterXSize(), 1, CatchmentsMiddleRow, DEM->GetRasterXSize(),
                              1, GDT_Int32, 0, 0 );
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
  Catchments->SetGeoTransform(m_GeoTransformVal);
  Catchments->SetProjection(DEM->GetProjectionRef());
  GDALClose( Catchments );
  GDALClose( Outlets );
  GDALClose( PlotsR );
  GDALClose( Drainage );
  GDALClose( RasterID );
  GDALClose( DEM );

}


// =====================================================================
// =====================================================================


/**
  @internal
  <hr>
  Files
    - m_InputPlotsVectorFile
    - m_OutputCatchmentsRasterFile
    - m_OutputDownslopeRasterFile
    - m_OutputPlotsRasterFile
    - m_OutputOutletsRasterFile
    - m_OutputReceiversRasterFile
    - m_OutputCatchmentsGroupedVectorFile
    - m_OutputCatchmentsVectorFile
    - m_OutputEntitiesGroupedVectorFile
    - m_OutputEntitiesVectorFile
    - m_OutputLinearEntitiesVectorFile
    - m_OutputSurfaceEntitiesVectorFile
    - m_OutputPlotsVectorFile
    - m_OutputOutletsVectorFile
    - m_OutputReceiversVectorFile
    - m_OutputPlotsAndEntitiesUnionVectorFile
*/
void LandProcessor::createSRFandLNR()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);


  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int i, FieldIndex, Value, IDCat, IDCat2, IDCat3, IDN, ID2N, ID3N, IDCatR, IDPRO, IDPR, IDPR2,
  IDPlotO, IDPlot, IDPlot2, IDPlot3, IDPlotNew, IDPlot2New, IDPlot3New,
  FID, IDField, IDMax, N, Beginning, End, pos, IDPlotAbove;

  std::string FieldName, FileName, ID, ID2, ID3, IDOld, IDNew, ID2New,
  ID3New, IDAbove, IDTo, FaceA, FaceB, FaceAB, FaceBA;

  std::vector<std::string> FieldNamesTable, IDNewTable, ID2NewTable, ID3NewTable, IDOldTable, FacesTable;

  std::vector<int> FIDTable, FIDToRemoveTable, IDPlotTable, IDPlotNewTable, IDPlot2NewTable, IDPlot3NewTable,
  IDPlotLargeTable, MissingIDTable, FIDSmallTable, FeaturesToDeleteTable,
  IDPlotAllSmallTable, FIDSliverTable, FIDConnectionProblemTable;

  std::vector < std::vector <int>> OutletsTable, ReceiversTable;

  std::vector< std::vector<int>> Table(5,std::vector<int >(0)), BigTable(5,std::vector<int >(0));
  std::vector< std::vector<int>> CatchmentsTable(10,std::vector<int >(0));
  std::vector< std::vector<int>> IDPlotIDMaxTable(2,std::vector<int >(0));

  GDALDataset *Dataset;

  OGRDataSource *Plots, *Catchments, *Outlets, *Receivers, *DataSource, *Entities, *Union, *Final, *LNR, *SRF;
  OGRLayer *PlotsLayer, *CatchmentsLayer, *SQLLayer, *OutletsLayer, *ReceiversLayer, *EntitiesLayer,
  *UnionLayer, *FinalLayer, *LNRLayer, *SRFLayer;
  OGRFeature *CatchmentsFeature, *OutletsFeature, *ReceiversFeature,
  *SQLFeature, *EntitiesFeature, *EntitiesFeatureToUpdate, *NewFeature, *AggregatingFeature,
  *UnionFeature, *UnionFeatureToUpdate, *FinalFeature,
  *LNRFeature, *SRFFeature;
  OGRGeometry *CatchmentsGeometry, *OutletsGeometry, *ReceiversGeometry, *SQLGeometry,
  *EntitiesGeometry, *NewGeometry, *AggregatingGeometry, *UnionGeometry;
  OGRFeatureDefn *FeatureDefn;

  CPLErr ERR;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  if (!openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputCatchmentsRasterFile)))
  {
    throw std::runtime_error("LandProcessor::createSRFandLNR(): " + m_OutputCatchmentsRasterFile + ": no such file in the output raster directory");
  }

  if (!openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputOutletsRasterFile)))
  {
    throw std::runtime_error("LandProcessor::createSRFandLNR(): " + m_OutputOutletsRasterFile + ": no such file in the output raster directory");
  }

  if (!openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputReceiversRasterFile)))
  {
    throw std::runtime_error("LandProcessor::createSRFandLNR(): " + m_OutputReceiversRasterFile + ": no such file in the output raster directory");
  }

  if (!openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputDownslopeRasterFile)))
  {
    throw std::runtime_error("LandProcessor::createSRFandLNR(): " + m_OutputDownslopeRasterFile + ": no such file in the output raster directory");
  }

  // ================================================================================================
  // Creating a catchment vector file based on the catchment raster (IDCat is the catchment ID)
  // ================================================================================================

  Plots = OGRSFDriverRegistrar::Open( getInputVectorPath(m_InputPlotsVectorFile).c_str(), TRUE );

  mp_SRS = Plots->GetLayer(0)->GetSpatialRef();

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputCatchmentsVectorFile)))
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

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputOutletsVectorFile)))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputOutletsVectorFile).c_str());
  }

  Outlets = mp_VectorDriver->CreateDataSource(getOutputVectorPath(m_OutputOutletsVectorFile).c_str());
  OutletsLayer = Outlets->CreateLayer("outlets", mp_SRS, wkbPolygon);

  FieldNamesTable.clear();
  FieldNamesTable = {"IDCatO", "IDPlotO", "IDPRO"};

  for (int i = 0; i < FieldNamesTable.size(); i++)
  {
    createField(OutletsLayer, FieldNamesTable[i].c_str(), OFTInteger);
  }

  GDALClose( Dataset );

  Dataset = (GDALDataset *) GDALOpen(getOutputRasterPath(m_OutputOutletsRasterFile).c_str(), GA_ReadOnly );
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

  Dataset = (GDALDataset *) GDALOpen(getOutputRasterPath(m_OutputPlotsRasterFile).c_str(), GA_ReadOnly );
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

  Dataset = (GDALDataset *) GDALOpen(getOutputRasterPath(m_OutputDownslopeRasterFile).c_str(), GA_ReadOnly );
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

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputReceiversVectorFile)))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputReceiversVectorFile).c_str());
  }

  Receivers = mp_VectorDriver->CreateDataSource(getOutputVectorPath(m_OutputReceiversVectorFile).c_str());
  ReceiversLayer = Receivers->CreateLayer("receivers", mp_SRS, wkbPolygon);

  FieldNamesTable.clear();
  FieldNamesTable = {"IDPRR", "IDPlotR", "IDCatR"};

  for (int i = 0; i < FieldNamesTable.size(); i++)
  {
    createField(ReceiversLayer, FieldNamesTable.at(i).c_str(), OFTInteger);
  }

  FieldNamesTable.clear();

  GDALClose( Dataset );

  Dataset = (GDALDataset *) GDALOpen(getOutputRasterPath(m_OutputReceiversRasterFile).c_str(), GA_ReadOnly );
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

  Dataset = (GDALDataset *) GDALOpen(getOutputRasterPath(m_OutputPlotsRasterFile).c_str(), GA_ReadOnly );
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

  Dataset = (GDALDataset *) GDALOpen(getOutputRasterPath(m_OutputCatchmentsRasterFile).c_str(), GA_ReadOnly );
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

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputCatchmentsGroupedVectorFile)))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputCatchmentsGroupedVectorFile).c_str());
  }

  const std::string CatchmentsLayerName = getLayerNameFromFilename(m_OutputCatchmentsVectorFile);
  const std::string CatchmentsGroupedLayerName = getLayerNameFromFilename(m_OutputCatchmentsGroupedVectorFile);

  Catchments = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputCatchmentsVectorFile).c_str(), TRUE );
  m_SQLRequest = "SELECT ST_Union(geometry), IDCat FROM " +CatchmentsLayerName + " GROUP BY IDCat";
  FileName = CatchmentsGroupedLayerName;
  SQLLayer = Catchments->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
  Catchments->CopyLayer(SQLLayer, FileName.c_str(), nullptr);
  Catchments->ReleaseResultSet(SQLLayer);
  OGRDataSource::DestroyDataSource(Catchments);

  // ===============================================================================
  // Populating catchments vector attribute table
  // using outlets and receivers vectors
  // The attributes obtain during this operation will be used for further regrouping
  // ===============================================================================

  Outlets = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputOutletsVectorFile).c_str(), TRUE );
  OutletsLayer = Outlets->GetLayer(0);

  Receivers = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputReceiversVectorFile).c_str(), TRUE );
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

  Catchments = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputCatchmentsGroupedVectorFile).c_str(), TRUE );
  CatchmentsLayer = Catchments->GetLayer(0);

  FieldNamesTable.clear();
  FieldNamesTable = {"IDPlot", "IDPR", "IDCat2", "IDPlot2", "IDPR2", "IDCat3", "IDPlot3", "IDPR3"};

  for (int i = 0; i < FieldNamesTable.size(); i++)
  {
    createField(CatchmentsLayer, FieldNamesTable[i].c_str(), OFTInteger);
  }

  FieldNamesTable.clear();
  FieldNamesTable = {"ID","ID2","ID3"};

  for (int i = 0; i < FieldNamesTable.size(); i++)
  {
    createField(CatchmentsLayer, FieldNamesTable[i].c_str(),  OFTString);
  }

  FieldNamesTable.clear();

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

  // ======================================================================
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
  // ======================================================================

  Catchments = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputCatchmentsGroupedVectorFile).c_str(), TRUE );
  CatchmentsLayer = Catchments->GetLayer(0);
  CatchmentsLayer->ResetReading();

  while ((CatchmentsFeature = CatchmentsLayer->GetNextFeature()) != nullptr)
  {
    CatchmentsTable[0].push_back(CatchmentsFeature->GetFID());
    CatchmentsTable[1].push_back(CatchmentsFeature->GetFieldAsInteger("IDCat"));
    CatchmentsTable[2].push_back(CatchmentsFeature->GetFieldAsInteger("IDPlot"));
    CatchmentsTable[3].push_back(CatchmentsFeature->GetFieldAsInteger("IDCat2"));
    CatchmentsTable[4].push_back(CatchmentsFeature->GetFieldAsInteger("IDPlot2"));
    CatchmentsTable[5].push_back(CatchmentsFeature->GetFieldAsInteger("IDCat3"));
    CatchmentsTable[6].push_back(CatchmentsFeature->GetFieldAsInteger("IDPlot3"));
    CatchmentsTable[7].push_back(0);
    CatchmentsTable[8].push_back(0);
    CatchmentsTable[9].push_back(0);
  }

  OGRDataSource::DestroyDataSource(Catchments);

  Table.clear();

  for (int i = 0; i < CatchmentsTable[0].size(); i++)
  {
    FID = CatchmentsTable[0][i];
    IDPlot = CatchmentsTable[2][i];
    IDPlot2 = CatchmentsTable[4][i];
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
          CatchmentsTable[7][i] = 1;
          CatchmentsTable[8][i] = 1;
          CatchmentsTable[9][i] = 1;
        }
        else
        {
          Table[0].push_back(0);
          Table[1].push_back(1);
          Table[2].push_back(0);
          Table[3].push_back(1);
          Table[4].push_back(FID);
          CatchmentsTable[7][i] = 1;
          CatchmentsTable[8][i] = 1;
          CatchmentsTable[9][i] = 1;;
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

  DataSource = OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );

  while (Table[0].size() != 0)
  {
    i = 0;
    while (i < Table[0].size())
    {
      if (Table[0][i] != 0)
      {
        FID = CatchmentsTable[0][Table[4][i]];
        IDPlot = CatchmentsTable[2][Table[4][i]];
        IDN = CatchmentsTable[7][Table[4][i]];
        ID2N = CatchmentsTable[8][Table[4][i]];
        ID3N = CatchmentsTable[9][Table[4][i]];
        m_SQLRequest = "SELECT ROWID as RID FROM " + CatchmentsGroupedLayerName + " WHERE IDPlot2 == 0 "
            "AND IDPlot == " + std::to_string(IDPlot) + " AND "
            "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + CatchmentsGroupedLayerName + " WHERE ROWID == " + std::to_string(FID) +"))) > 0";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          FIDTable.clear();
          for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
          {
            FIDTable.push_back(SQLLayer->GetFeature(j)->GetFieldAsInteger("RID"));
          }
          DataSource->ReleaseResultSet(SQLLayer);
          for (int k = 0; k < FIDTable.size(); k++)
          {
            if (CatchmentsTable[7][FIDTable[k]] == 0)
            {
              CatchmentsTable[7][FIDTable[k]] = IDN;
              CatchmentsTable[8][FIDTable[k]] = ID2N;
              CatchmentsTable[9][FIDTable[k]] = ID3N;
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

          }
          FIDTable.clear();
        }
        else
        {
          DataSource->ReleaseResultSet(SQLLayer);
        }
      }
      i += 1;
    }
    IDMax = Table[1].back();
    Table[0].clear();
    Table[1].clear();
    Table[2].clear();
    Table[3].clear();
    Table[4].clear();
    for (int j = 0; j < CatchmentsTable[0].size(); j++)
    {
      FID = CatchmentsTable[0][j];
      IDN = CatchmentsTable[7][j];
      IDPlot = CatchmentsTable[2][j];
      IDPlot2 = CatchmentsTable[4][j];
      if (IDPlot2 == 0 and IDN == 0)
      {
        if (std::find(Table[0].begin(), Table[0].end(), IDPlot) == Table[0].end())
        {
          IDN = IDMax+1;
          CatchmentsTable[7][j] = IDN;
          CatchmentsTable[8][j] = 1;
          CatchmentsTable[9][j] = 1;
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

  for (int i = 0; i < CatchmentsTable[0].size(); i++)
  {
    IDCat2 = CatchmentsTable[3][i];
    ID2N = CatchmentsTable[8][i];
    if (ID2N == 0)
    {
      for (int j = 0; j < CatchmentsTable[0].size(); j++)
      {
        if (CatchmentsTable[1][j] == IDCat2 and CatchmentsTable[7][j] != 0)
        {
          CatchmentsTable[8][i] = CatchmentsTable[7][j];
          CatchmentsTable[9][i] = CatchmentsTable[8][j];
        }
      }
    }
  }

  // =======================================================================================
  // Setting the rest of the catchments with ID, ID2 and ID3 using the same principal as above
  // =======================================================================================

  N = 0;

  for (int i = 0; i < CatchmentsTable[0].size(); i++)
  {
    if (CatchmentsTable[7][i] == 0)
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
    for (int i = 0; i < CatchmentsTable[0].size(); i++)
    {
      FID = CatchmentsTable[0][i];
      IDPlot = CatchmentsTable[2][i];
      IDN = CatchmentsTable[7][i];
      IDPlot2 = CatchmentsTable[4][i];
      ID2N = CatchmentsTable[8][i];
      if (IDN == 0 and ID2N != 0)
      {
        if (std::find(Table[0].begin(), Table[0].end(), IDPlot) == Table[0].end() and
            std::find(BigTable[0].begin(), BigTable[0].end(), IDPlot) == BigTable[0].end())
        {
          CatchmentsTable[7][i] = 1;
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
          CatchmentsTable[7][i] = IDMax+1;
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
          FID = CatchmentsTable[0][Table[4][i]];
          IDPlot = CatchmentsTable[2][Table[4][i]];
          IDN = CatchmentsTable[7][Table[4][i]];
          IDPlot2 = CatchmentsTable[4][Table[4][i]];
          ID2N = CatchmentsTable[8][Table[4][i]];
          ID3N = CatchmentsTable[9][Table[4][i]];
          m_SQLRequest = "SELECT ROWID as RID FROM " + CatchmentsGroupedLayerName + " "
              "WHERE IDPlot == " + std::to_string(IDPlot) + " AND IDPlot2 == " + std::to_string(IDPlot2) + " AND "
              "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + CatchmentsGroupedLayerName + " WHERE ROWID == " + std::to_string(FID) + "))) > 0";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() != 0)
          {
            for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
            {
              FIDTable.push_back(SQLLayer->GetFeature(j)->GetFieldAsInteger("RID"));
            }
            DataSource->ReleaseResultSet(SQLLayer);
            for (int k = 0; k < FIDTable.size(); k++)
            {
              if (CatchmentsTable[7][FIDTable[k]] == 0 and
                  CatchmentsTable[8][FIDTable[k]] == ID2N)
              {
                CatchmentsTable[7][FIDTable[k]] = IDN;
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
            }
            FIDTable.clear();
          }
          else
          {
            DataSource->ReleaseResultSet(SQLLayer);
          }
        }
        i += 1;
      }
      Table[0].clear();
      Table[1].clear();
      Table[2].clear();
      Table[3].clear();
      Table[4].clear();
      for (int j = 0; j < CatchmentsTable[0].size(); j++)
      {
        FID = CatchmentsTable[0][j];
        IDPlot = CatchmentsTable[2][j];
        IDN = CatchmentsTable[7][j];
        IDPlot2 = CatchmentsTable[4][j];
        ID2N = CatchmentsTable[8][j];
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
            CatchmentsTable[7][j] = IDMax+1;
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
      for (int i = 0; i < CatchmentsTable[0].size(); i++)
      {
        IDCat2 = CatchmentsTable[3][i];
        ID2N = CatchmentsTable[8][i];
        if (ID2N == 0)
        {
          for (int j = 0; j < CatchmentsTable[0].size(); j++)
          {
            if (CatchmentsTable[1][j] == IDCat2 and CatchmentsTable[7][j] != 0)
            {
              CatchmentsTable[8][i] = CatchmentsTable[7][j];
              CatchmentsTable[9][i] = CatchmentsTable[8][j];
            }
          }
        }
      }
    }
    N = 0;
    for (int i = 0; i < CatchmentsTable[0].size(); i++)
    {
      if (CatchmentsTable[7][i] == 0)
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

  for (int i = 0; i < CatchmentsTable[0].size(); i++)
  {
    IDPlot = CatchmentsTable[2][i];
    IDPlot2 = CatchmentsTable[4][i];
    IDPlot3 = CatchmentsTable[6][i];
    IDN = CatchmentsTable[7][i];
    ID2N = CatchmentsTable[8][i];
    ID3N = CatchmentsTable[9][i];
    CatchmentsFeature = CatchmentsLayer->GetFeature(CatchmentsTable[0][i]);
    IDNew = std::to_string(IDPlot)+"N"+std::to_string(IDN);
    ID2New = std::to_string(IDPlot2)+"N"+std::to_string(ID2N);
    ID3New = std::to_string(IDPlot3)+"N"+std::to_string(ID3N);
    CatchmentsFeature->SetField("ID", IDNew.c_str());
    CatchmentsFeature->SetField("ID2", ID2New.c_str());
    CatchmentsFeature->SetField("ID3", ID3New.c_str());
    CatchmentsLayer->SetFeature(CatchmentsFeature);
  }

  OGRDataSource::DestroyDataSource(Catchments);

  CatchmentsTable[0].clear();
  CatchmentsTable[1].clear();
  CatchmentsTable[2].clear();
  CatchmentsTable[3].clear();
  CatchmentsTable[4].clear();
  CatchmentsTable[5].clear();
  CatchmentsTable[6].clear();
  CatchmentsTable[7].clear();
  CatchmentsTable[8].clear();
  CatchmentsTable[9].clear();

  // ======================================================================
  //  Regrouping catchments that have the same attributes into new entities
  // ======================================================================

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputEntitiesVectorFile)))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str());
  }

  const std::string EntitiesLayerName = getLayerNameFromFilename(m_OutputEntitiesVectorFile);

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
  m_SQLRequest = "SELECT ST_Union(geometry), IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3 FROM " + CatchmentsGroupedLayerName + " GROUP BY IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3";
  SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
  FileName = EntitiesLayerName;
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
  m_SQLRequest = "SELECT ROWID as RID, IDPlot FROM " + EntitiesLayerName + " WHERE ST_Area(geometry) > " + std::to_string(m_MinEntSize);
  SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());

  for(int i = 0; i < SQLLayer->GetFeatureCount(); i++)
  {
    if (std::find(IDPlotLargeTable.begin(), IDPlotLargeTable.end(), SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot")) == IDPlotLargeTable.end())
    {
      IDPlotLargeTable.push_back(SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot"));
    }
  }

  DataSource->ReleaseResultSet(SQLLayer);
  m_SQLRequest = "SELECT ROWID as RID, IDPlot FROM " + EntitiesLayerName + " WHERE ST_Area(geometry) <= " + std::to_string(m_MinEntSize);
  SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());

  for(int i = 0; i < SQLLayer->GetFeatureCount(); i++)
  {
    if (std::find(IDPlotLargeTable.begin(), IDPlotLargeTable.end(), SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot")) == IDPlotLargeTable.end() and std::find(IDPlotAllSmallTable.begin(), IDPlotAllSmallTable.end(), SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot")) == IDPlotAllSmallTable.end())
    {
      IDPlotAllSmallTable.push_back(SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot"));
    }
  }

  DataSource->ReleaseResultSet(SQLLayer);

  IDNewTable.clear();
  ID2NewTable.clear();
  ID3NewTable.clear();
  IDPlot2NewTable.clear();
  IDPlot3NewTable.clear();

  if (IDPlotAllSmallTable.size() != 0)
  {
    for (int i = 0; i < IDPlotAllSmallTable.size(); i++)
    {
      m_SQLRequest = "SELECT IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3 FROM " + EntitiesLayerName + " WHERE IDPlot == " + std::to_string(IDPlotAllSmallTable[i]) + " "
          "AND IDPlot != IDPlot3 ORDER BY ST_Area(geometry) DESC";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        IDNewTable.push_back(SQLLayer->GetFeature(0)->GetFieldAsString("ID"));
        IDPlot2NewTable.push_back(SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot2"));
        ID2NewTable.push_back(SQLLayer->GetFeature(0)->GetFieldAsString("ID2"));
        IDPlot3NewTable.push_back(SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot3"));
        ID3NewTable.push_back(SQLLayer->GetFeature(0)->GetFieldAsString("ID3"));
      }
      DataSource->ReleaseResultSet(SQLLayer);
    }
    Entities = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
    EntitiesLayer = Entities->GetLayer(0);
    FeatureDefn = EntitiesLayer->GetLayerDefn();
    FIDToRemoveTable.clear();
    for (int i = 0; i < IDPlotAllSmallTable.size(); i++)
    {
      NewGeometry = nullptr;
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + EntitiesLayerName + " WHERE IDPlot == " + std::to_string(IDPlotAllSmallTable[i]) + " ORDER BY ST_Area(geometry) DESC";
      SQLLayer = Entities->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      FIDToRemoveTable.clear();
      for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
      {
        SQLFeature = SQLLayer->GetFeature(j);
        SQLGeometry = SQLFeature->GetGeometryRef();
        FIDToRemoveTable.push_back(SQLFeature->GetFieldAsInteger("RID"));
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
      for (int j = 0; j < FIDToRemoveTable.size(); j++)
      {
        EntitiesLayer->DeleteFeature(FIDToRemoveTable[j]);
      }
      m_SQLRequest = "REPACK " + EntitiesLayerName;
      Entities->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
      EntitiesLayer->SyncToDisk();
      FIDToRemoveTable.clear();
      NewFeature = OGRFeature::CreateFeature(FeatureDefn);
      NewFeature->SetField("IDPlot", IDPlotAllSmallTable[i]);
      NewFeature->SetField("ID", IDNewTable[i].c_str());
      NewFeature->SetField("IDPlot2", IDPlot2NewTable[i]);
      NewFeature->SetField("ID2", ID2NewTable[i].c_str());
      NewFeature->SetField("IDPlot3", IDPlot3NewTable[i]);
      NewFeature->SetField("ID3", ID3NewTable[i].c_str());
      NewFeature->SetGeometry(NewGeometry);
      EntitiesLayer->CreateFeature(NewFeature);
      EntitiesLayer->SyncToDisk();
    }
    OGRDataSource::DestroyDataSource(Entities);
    Entities = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
    EntitiesLayer = Entities->GetLayer(0);
    EntitiesLayer->ResetReading();
    while ((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
    {
      if (std::find(IDPlotAllSmallTable.begin(), IDPlotAllSmallTable.end(), EntitiesFeature->GetFieldAsInteger("IDPlot2")) != IDPlotAllSmallTable.end())
      {
        pos = std::find(IDPlotAllSmallTable.begin(), IDPlotAllSmallTable.end(), EntitiesFeature->GetFieldAsInteger("IDPlot2")) - IDPlotAllSmallTable.begin();
        EntitiesFeature->SetField("ID2", IDNewTable[pos].c_str());
        EntitiesFeature->SetField("IDPlot3", IDPlot2NewTable[pos]);
        EntitiesFeature->SetField("ID3", ID2NewTable[pos].c_str());
        EntitiesLayer->SetFeature(EntitiesFeature);
      }
    }
    EntitiesLayer->SyncToDisk();
    EntitiesLayer = Entities->GetLayer(0);
    EntitiesLayer->ResetReading();
    while ((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
    {
      if (std::find(IDPlotAllSmallTable.begin(), IDPlotAllSmallTable.end(), EntitiesFeature->GetFieldAsInteger("IDPlot3")) != IDPlotAllSmallTable.end())
      {
        pos = std::find(IDPlotAllSmallTable.begin(), IDPlotAllSmallTable.end(), EntitiesFeature->GetFieldAsInteger("IDPlot3")) - IDPlotAllSmallTable.begin();
        EntitiesFeature->SetField("ID3", IDNewTable[pos].c_str());
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
  m_SQLRequest = "SELECT ROWID as RID, IDPlot FROM " + EntitiesLayerName + " WHERE ST_Area(geometry) <= " + std::to_string(m_MinEntSize);
  SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());

  for(int i = 0; i < SQLLayer->GetFeatureCount(); i++)
  {
    if (std::find(IDPlotLargeTable.begin(), IDPlotLargeTable.end(), SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot")) != IDPlotLargeTable.end() and
        std::find(IDPlotAllSmallTable.begin(), IDPlotAllSmallTable.end(), SQLLayer->GetFeature(i)->GetFieldAsInteger("IDPlot")) == IDPlotAllSmallTable.end())
    {
      FIDSmallTable.push_back(SQLLayer->GetFeature(i)->GetFieldAsInteger("RID"));
    }
  }

  DataSource->ReleaseResultSet(SQLLayer);
  OGRDataSource::DestroyDataSource(DataSource);

  IDOldTable.clear();
  IDNewTable.clear();
  ID2NewTable.clear();
  IDPlotNewTable.clear();
  IDPlot2NewTable.clear();
  IDPlot3NewTable.clear();

  FeaturesToDeleteTable = FIDSmallTable;
  Entities = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
  EntitiesLayer = Entities->GetLayer(0);

  while(FIDSmallTable.size() != 0)
  {
    Beginning = FIDSmallTable.size();
    for (int i = 0; i < FIDSmallTable.size(); i++)
    {
      DataSource = OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );
      EntitiesFeature = EntitiesLayer->GetFeature(FIDSmallTable[i]);
      ID = EntitiesFeature->GetFieldAsString("ID");
      EntitiesGeometry = EntitiesFeature->GetGeometryRef();
      IDPlot = EntitiesFeature->GetFieldAsInteger("IDPlot");
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + EntitiesLayerName + " WHERE "
          "IDPlot == " + std::to_string(IDPlot) + " AND ID3 != " + m_QMark + ID + m_QMark + " AND "
          "ST_Area(geometry) > " + std::to_string(m_MinEntSize) + " AND "
          "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + EntitiesLayerName + " WHERE ROWID == " + std::to_string(FIDSmallTable[i])+ "))) > 0 "
          "ORDER BY ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + EntitiesLayerName + " WHERE ROWID == " + std::to_string(FIDSmallTable[i])+ "))) DESC, ST_Area(geometry) DESC";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        SQLFeature = SQLLayer->GetFeature(0);
        FID = SQLFeature->GetFieldAsInteger("RID");
        IDOldTable.push_back(EntitiesFeature->GetFieldAsString("ID"));
        IDNewTable.push_back(SQLFeature->GetFieldAsString("ID"));
        ID2NewTable.push_back(SQLFeature->GetFieldAsString("ID2"));
        IDPlotNewTable.push_back(SQLFeature->GetFieldAsInteger("IDPlot"));
        IDPlot2NewTable.push_back(SQLFeature->GetFieldAsInteger("IDPlot2"));
        AggregatingFeature = EntitiesLayer->GetFeature(FID);
        AggregatingGeometry = AggregatingFeature->GetGeometryRef();
        AggregatingGeometry = AggregatingGeometry->Union(EntitiesGeometry);
        AggregatingFeature->SetGeometry(AggregatingGeometry);
        EntitiesLayer->SetFeature(AggregatingFeature);
        EntitiesLayer->SyncToDisk();
        FIDToRemoveTable.push_back(FIDSmallTable[i]);
        EntitiesFeature = nullptr;
      }
      DataSource->ReleaseResultSet(SQLLayer);
      OGRDataSource::DestroyDataSource(DataSource);
    }
    for (int k = 0; k < FIDToRemoveTable.size(); k++)
    {
      if (std::find(FIDSmallTable.begin(), FIDSmallTable.end(), FIDToRemoveTable[k]) != FIDSmallTable.end())
      {
        pos = std::find(FIDSmallTable.begin(), FIDSmallTable.end(), FIDToRemoveTable[k]) - FIDSmallTable.begin();
        FIDSmallTable.erase(FIDSmallTable.begin()+pos);
      }
    }
    FIDToRemoveTable.clear();
    End = FIDSmallTable.size();
    if (Beginning == End)
    {
      break;
    }
  }

  for (int i = 0; i < FeaturesToDeleteTable.size(); i++)
  {
    if (std::find(FIDSmallTable.begin(), FIDSmallTable.end(), FeaturesToDeleteTable[i]) == FIDSmallTable.end())
    {
      EntitiesLayer->DeleteFeature(FeaturesToDeleteTable[i]);
    }
  }

  m_SQLRequest = "REPACK " + EntitiesLayerName;
  Entities->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
  EntitiesLayer->SyncToDisk();
  EntitiesLayer = Entities->GetLayer(0);
  EntitiesLayer->ResetReading();

  while ((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
  {
    ID2 = EntitiesFeature->GetFieldAsString("ID2");
    if (std::find(IDOldTable.begin(), IDOldTable.end(), ID2.c_str()) != IDOldTable.end())
    {
      pos = std::find(IDOldTable.begin(), IDOldTable.end(), ID2.c_str()) - IDOldTable.begin();
      EntitiesFeature->SetField("ID2", IDNewTable[pos].c_str());
      EntitiesFeature->SetField("IDPlot3", IDPlot2NewTable[pos]);
      EntitiesFeature->SetField("ID3", ID2NewTable[pos].c_str());
      EntitiesLayer->SetFeature(EntitiesFeature);
    }
  }

  EntitiesLayer->SyncToDisk();
  EntitiesLayer = Entities->GetLayer(0);
  EntitiesLayer->ResetReading();

  while ((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
  {
    ID3 = EntitiesFeature->GetFieldAsString("ID3");
    if (std::find(IDOldTable.begin(), IDOldTable.end(), ID3.c_str()) != IDOldTable.end())
    {
      pos = std::find(IDOldTable.begin(), IDOldTable.end(), ID3.c_str()) - IDOldTable.begin();
      EntitiesFeature->SetField("ID3", IDNewTable[pos].c_str());
      EntitiesLayer->SetFeature(EntitiesFeature);
    }
  }

  OGRDataSource::DestroyDataSource(Entities);

  IDNewTable.clear();
  ID2NewTable.clear();
  IDPlotNewTable.clear();
  IDPlot2NewTable.clear();
  IDOldTable.clear();

  //=================================================================
  // Regrouping entities that belong to the same plot, drain into the
  // same entity and share the common border
  //=================================================================

  N = 1;

  while (N != 0)
  {
    FIDToRemoveTable.clear();
    IDOldTable.clear();
    IDNewTable.clear();
    IDPlotNewTable.clear();
    Entities = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
    EntitiesLayer = Entities->GetLayer(0);
    EntitiesLayer->ResetReading();
    DataSource  = OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );
    NewGeometry = nullptr;
    while((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
    {
      FID = EntitiesFeature->GetFID();
      IDPlot = EntitiesFeature->GetFieldAsInteger("IDPlot");
      ID = EntitiesFeature->GetFieldAsString("ID");
      ID2 = EntitiesFeature->GetFieldAsString("ID2");
      if (ID != "0N1" and std::find(FIDToRemoveTable.begin(), FIDToRemoveTable.end(), FID) == FIDToRemoveTable.end()
          and std::find(IDPlotNewTable.begin(), IDPlotNewTable.end(), IDPlot) == IDPlotNewTable.end())
      {
        m_SQLRequest = "SELECT geometry, ROWID as RID, * FROM " + EntitiesLayerName + " WHERE "
            "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + EntitiesLayerName + " WHERE ROWID == " + std::to_string(FID) + "))) > 0 AND "
            "IDPlot == " + std::to_string(IDPlot) + " AND ID2 == " + m_QMark + ID2 + m_QMark + " AND "
            "ROWID != " + std::to_string(FID) + " AND "
            "ST_Touches("
            "ST_Intersection(geometry,(SELECT geometry FROM " + EntitiesLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")), "
            "ST_Intersection((SELECT geometry FROM " + EntitiesLayerName + " WHERE ROWID == " + std::to_string(FID) + "),"
            "(SELECT geometry FROM " + EntitiesLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + "))) "
            "ORDER BY "
            "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM "
            "" + EntitiesLayerName + " WHERE "
            "ROWID == " + std::to_string(FID) + ")))";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          NewGeometry = nullptr;
          for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
          {
            if (std::find(FIDToRemoveTable.begin(), FIDToRemoveTable.end(), SQLLayer->GetFeature(j)->GetFieldAsInteger("RID")) == FIDToRemoveTable.end())
            {
              FIDToRemoveTable.push_back(SQLLayer->GetFeature(j)->GetFieldAsInteger("RID"));
              IDOldTable.push_back(SQLLayer->GetFeature(j)->GetFieldAsString("ID"));
              IDNewTable.push_back(ID.c_str());
              IDPlotNewTable.push_back(IDPlot);
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
            EntitiesGeometry = EntitiesFeature->GetGeometryRef();
            EntitiesGeometry = EntitiesGeometry->Union(NewGeometry);
            EntitiesFeature->SetGeometry(EntitiesGeometry);
            EntitiesLayer->SetFeature(EntitiesFeature);
          }
        }
        DataSource->ReleaseResultSet(SQLLayer);
      }
    }
    OGRDataSource::DestroyDataSource(DataSource);
    EntitiesLayer = Entities->GetLayer(0);
    EntitiesLayer->ResetReading();
    if (FIDToRemoveTable.size() != 0)
    {
      for (int i = 0; i < FIDToRemoveTable.size(); i++)
      {
        EntitiesLayer->DeleteFeature(FIDToRemoveTable[i]);
      }
    }
    FIDToRemoveTable.clear();
    m_SQLRequest = "REPACK " + EntitiesLayerName;
    Entities->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
    Entities->SyncToDisk();
    OGRDataSource::DestroyDataSource(Entities);
    Entities  = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
    EntitiesLayer = Entities->GetLayer(0);
    EntitiesLayer->ResetReading();
    while((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
    {
      if (std::find(IDOldTable.begin(), IDOldTable.end(), EntitiesFeature->GetFieldAsString("ID2")) != IDOldTable.end())
      {
        pos = std::find(IDOldTable.begin(), IDOldTable.end(), EntitiesFeature->GetFieldAsString("ID2")) - IDOldTable.begin();
        EntitiesFeature->SetField("ID2", IDNewTable[pos].c_str());
        EntitiesFeature->SetField("IDPlot2", IDPlotNewTable[pos]);
        EntitiesLayer->SetFeature(EntitiesFeature);
      }
    }
    EntitiesLayer->SyncToDisk();
    Entities->SyncToDisk();
    EntitiesLayer = Entities->GetLayer(0);
    EntitiesLayer->ResetReading();
    while((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
    {
      if (std::find(IDOldTable.begin(), IDOldTable.end(), EntitiesFeature->GetFieldAsString("ID3")) != IDOldTable.end())
      {
        pos = std::find(IDOldTable.begin(), IDOldTable.end(), EntitiesFeature->GetFieldAsString("ID3")) - IDOldTable.begin();
        EntitiesFeature->SetField("ID3", IDNewTable[pos].c_str());
        EntitiesFeature->SetField("IDPlot3", IDPlotNewTable[pos]);
        EntitiesLayer->SetFeature(EntitiesFeature);
      }
    }
    EntitiesLayer->SyncToDisk();
    Entities->SyncToDisk();
    OGRDataSource::DestroyDataSource(Entities);
    Entities = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
    EntitiesLayer = Entities->GetLayer(0);
    EntitiesLayer->ResetReading();
    DataSource  = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
    N = 0;
    while((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
    {
      FID = EntitiesFeature->GetFID();
      IDPlot = EntitiesFeature->GetFieldAsInteger("IDPlot");
      ID = EntitiesFeature->GetFieldAsString("ID");
      ID2 = EntitiesFeature->GetFieldAsString("ID2");
      if (ID != "0N1")
      {
        m_SQLRequest = "SELECT geometry, ROWID as RID, * FROM " + EntitiesLayerName + " WHERE "
            "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + EntitiesLayerName + " WHERE ROWID == " + std::to_string(FID) + "))) > 0 AND "
            "IDPlot == " + std::to_string(IDPlot) + " AND ID2 == " + m_QMark + ID2 + m_QMark + " AND "
            "ROWID != " + std::to_string(FID) + " AND "
            "ST_Touches("
            "ST_Intersection(geometry,(SELECT geometry FROM " + EntitiesLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")), "
            "ST_Intersection((SELECT geometry FROM " + EntitiesLayerName + " WHERE ROWID == " + std::to_string(FID) + "),"
            "(SELECT geometry FROM " + EntitiesLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + "))) "
            "ORDER BY "
            "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM "
            "" + EntitiesLayerName + " WHERE "
            "ROWID == " + std::to_string(FID) + ")))";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          N += 1;
        }
        DataSource->ReleaseResultSet(SQLLayer);
      }
    }
    OGRDataSource::DestroyDataSource(DataSource);
    OGRDataSource::DestroyDataSource(Entities);
  }

  FIDToRemoveTable.clear();
  IDOldTable.clear();
  IDNewTable.clear();
  IDPlotNewTable.clear();

  // =========================================================================================
  // Finding and aggregating entities that are "locked" inside a parcel
  // =========================================================================================

  FIDToRemoveTable.clear();

  N = 1;

  while (N != 0)
  {
    Entities = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
    EntitiesLayer = Entities->GetLayer(0);
    EntitiesLayer->ResetReading();
    DataSource  = OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );
    while((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
    {
      FID = EntitiesFeature->GetFID();
      IDPlot = EntitiesFeature->GetFieldAsInteger("IDPlot");
      ID = EntitiesFeature->GetFieldAsString("ID");
      ID2 = EntitiesFeature->GetFieldAsString("ID2");
      if (ID != "0N1" and ID2 != "0N1")
      {
        m_SQLRequest = "SELECT ROWID as RID, * FROM " + EntitiesLayerName + " WHERE "
            "ST_Touches(geometry,(SELECT geometry FROM " + EntitiesLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
            "ROWID != " + std::to_string(FID);
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() <= 2)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID, * FROM " + EntitiesLayerName + " WHERE "
              "ST_Touches(geometry,(SELECT geometry FROM " + EntitiesLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
              "IDPlot == " + std::to_string(IDPlot) + " AND ROWID != " + std::to_string(FID);
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() != 0)
          {
            EntitiesGeometry = EntitiesFeature->GetGeometryRef();
            AggregatingFeature = EntitiesLayer->GetFeature(SQLLayer->GetFeature(0)->GetFieldAsInteger("RID"));
            AggregatingGeometry = AggregatingFeature->GetGeometryRef();
            AggregatingGeometry = AggregatingGeometry->Union(EntitiesGeometry);
            AggregatingFeature->SetGeometry(AggregatingGeometry);
            EntitiesLayer->SetFeature(AggregatingFeature);
            FIDToRemoveTable.push_back(FID);
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
    OGRDataSource::DestroyDataSource(Entities);
    Entities = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
    EntitiesLayer = Entities->GetLayer(0);
    EntitiesLayer->ResetReading();
    if (FIDToRemoveTable.size() != 0)
    {
      for (int i = 0; i < FIDToRemoveTable.size(); i++)
      {
        EntitiesLayer->DeleteFeature(FIDToRemoveTable[i]);
      }
    }
    FIDToRemoveTable.clear();
    m_SQLRequest = "REPACK " + EntitiesLayerName;
    Entities->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
    OGRDataSource::DestroyDataSource(Entities);
    Entities = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
    EntitiesLayer = Entities->GetLayer(0);
    EntitiesLayer->ResetReading();
    DataSource = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesVectorFile).c_str(), TRUE );
    N = 0;
    while((EntitiesFeature = EntitiesLayer->GetNextFeature()) != nullptr)
    {
      FID = EntitiesFeature->GetFID();
      IDPlot = EntitiesFeature->GetFieldAsInteger("IDPlot");
      ID = EntitiesFeature->GetFieldAsString("ID");
      ID2 = EntitiesFeature->GetFieldAsString("ID2");
      if (ID != "0N1" and ID2 != "0N1")
      {
        m_SQLRequest = "SELECT ROWID as RID, * FROM " + EntitiesLayerName + " WHERE "
            "ST_Touches(geometry,(SELECT geometry FROM " + EntitiesLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
            "ROWID != " + std::to_string(FID);
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() <= 2)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID, * FROM " + EntitiesLayerName + " WHERE "
              "ST_Touches(geometry,(SELECT geometry FROM " + EntitiesLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
              "IDPlot == " + std::to_string(IDPlot) + " AND ROWID != " + std::to_string(FID);
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() != 0)
          {
            N += 1;
          }
          DataSource->ReleaseResultSet(SQLLayer);
        }
        else
        {
          DataSource->ReleaseResultSet(SQLLayer);
        }
      }
    }
    OGRDataSource::DestroyDataSource(Entities);
    OGRDataSource::DestroyDataSource(DataSource);
  }

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

  FieldNamesTable.clear();
  FIDTable.clear();

  FieldNamesTable = {m_IDFieldName,"IDPlot","ID","IDPlot2","ID2","IDPlot3","ID3"};

  GRASS.appendTask("v.db.renamecolumn", {{"map", "union"},
                                         {"column", QString::fromStdString("a_" + FieldNamesTable[0] + "," + FieldNamesTable[0])}});

  for (int i = 1; i < FieldNamesTable.size(); i++)
  {
    GRASS.appendTask("v.db.renamecolumn", {{"map", "union"},
                                           {"column", QString::fromStdString("b_" + FieldNamesTable[i] + "," + FieldNamesTable[i])}});
  }

  GRASS.appendTask("v.db.dropcolumn", {{"map", "union"},
                                       {"column", QString::fromStdString("a_cat,b_cat")}});

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile)))
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

  Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
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

  const std::string PlotsAndEntitiesUnionLayerName = getLayerNameFromFilename(m_OutputPlotsAndEntitiesUnionVectorFile);

  m_SQLRequest = "REPACK " + PlotsAndEntitiesUnionLayerName;
  Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
  OGRDataSource::DestroyDataSource(Union);

  Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
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

  Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  UnionFeature = nullptr;
  UnionGeometry = nullptr;

  IDPlotTable.clear();
  MissingIDTable.clear();
  FIDSliverTable.clear();

  while ((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    if (std::find(IDPlotTable.begin(), IDPlotTable.end(), UnionFeature->GetFieldAsInteger("IDPlot")) == IDPlotTable.end())
    {
      IDPlotTable.push_back(UnionFeature->GetFieldAsInteger("IDPlot"));
    }
  }

  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  while ((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    if (std::find(IDPlotTable.begin(), IDPlotTable.end(), UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) == IDPlotTable.end() and
        std::find(MissingIDTable.begin(), MissingIDTable.end(), UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) == MissingIDTable.end())
    {
      MissingIDTable.push_back(UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str()));
    }
  }

  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  UnionFeature = nullptr;

  while ((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    if (UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str()) != UnionFeature->GetFieldAsInteger("IDPlot") and
        std::find(MissingIDTable.begin(), MissingIDTable.end(), UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) == MissingIDTable.end())
    {
      FIDSliverTable.push_back(UnionFeature->GetFID());
    }
  }

  OGRDataSource::DestroyDataSource(Union);

  // ======================================================================================
  // Regrouping entities that have at least one other feature that belong to the same plot,
  // have the same original parcel ID as the entities in question and whose
  // original parcel ID attribute is the same as IDPlot
  // ======================================================================================

  FIDToRemoveTable.clear();

  while (FIDSliverTable.size() != 0)
  {
    Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    for (int i = 0; i < FIDSliverTable.size(); i++)
    {
      UnionFeature = UnionLayer->GetFeature(FIDSliverTable[i]);
      UnionGeometry = UnionFeature->GetGeometryRef();
      DataSource = OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );
      AggregatingFeature = nullptr;
      AggregatingGeometry = nullptr;
      m_SQLRequest = "SELECT ROWID as RID FROM " + PlotsAndEntitiesUnionLayerName + " WHERE " + m_IDFieldName + " == IDPlot "
          "AND " + m_IDFieldName + " == " + std::to_string(UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) + "  AND "
          "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FIDSliverTable[i]) + "))) > 0 "
          "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FIDSliverTable[i]) + "))) DESC";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() > 1)
      {
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT ROWID as RID FROM " + PlotsAndEntitiesUnionLayerName + " WHERE " + m_IDFieldName + " == "
            "IDPlot AND " + m_IDFieldName + " == " + std::to_string(UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) + " AND "
            "ID2 == " + m_QMark + UnionFeature->GetFieldAsString("ID") + m_QMark + " AND "
            "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FIDSliverTable[i]) + "))) > 0 "
            "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FIDSliverTable[i]) + "))) DESC";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() == 0)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + PlotsAndEntitiesUnionLayerName + " WHERE " + m_IDFieldName + " == IDPlot "
              "AND " + m_IDFieldName + " == " + std::to_string(UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) + " AND "
              "ID == " + m_QMark + UnionFeature->GetFieldAsString("ID2") + m_QMark + " AND "
              "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FIDSliverTable[i]) + "))) > 0 "
              "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FIDSliverTable[i]) + "))) DESC";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() == 0)
          {
            DataSource->ReleaseResultSet(SQLLayer);
            m_SQLRequest = "SELECT ROWID as RID FROM " + PlotsAndEntitiesUnionLayerName + " WHERE " + m_IDFieldName + " == IDPlot "
                "AND " + m_IDFieldName + " == " + std::to_string(UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) + " AND "
                "ID2 == " + m_QMark + UnionFeature->GetFieldAsString("ID2") + m_QMark + " AND "
                "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FIDSliverTable[i]) + "))) > 0 "
                "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FIDSliverTable[i]) + "))) DESC";
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            if (SQLLayer->GetFeatureCount() == 0)
            {
              DataSource->ReleaseResultSet(SQLLayer);
              m_SQLRequest = "SELECT ROWID as RID FROM " + PlotsAndEntitiesUnionLayerName + " WHERE " + m_IDFieldName + " == IDPlot "
                  "AND " + m_IDFieldName + " == " + std::to_string(UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) + " AND "
                  "ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FIDSliverTable[i]) + "))) > 0 "
                  "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FIDSliverTable[i]) + "))) DESC";
              SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
              FID = SQLLayer->GetFeature(0)->GetFieldAsInteger("RID");
              DataSource->ReleaseResultSet(SQLLayer);
              AggregatingFeature = UnionLayer->GetFeature(FID);
              AggregatingGeometry = AggregatingFeature->GetGeometryRef();
              AggregatingFeature->SetGeometry(AggregatingGeometry->Union(UnionGeometry));
              UnionLayer->SetFeature(AggregatingFeature);
              UnionLayer->SyncToDisk();
              FIDToRemoveTable.push_back(FIDSliverTable[i]);
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
              FIDToRemoveTable.push_back(FIDSliverTable[i]);
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
            FIDToRemoveTable.push_back(FIDSliverTable[i]);
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
          FIDToRemoveTable.push_back(FIDSliverTable[i]);
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
        FIDToRemoveTable.push_back(FIDSliverTable[i]);
      }
      else
      {
        DataSource->ReleaseResultSet(SQLLayer);
      }
      OGRDataSource::DestroyDataSource(DataSource);
    }
    if (FIDToRemoveTable.size() != 0)
    {
      for (int i = 0; i < FIDToRemoveTable.size(); i++)
      {
        UnionLayer->DeleteFeature(FIDToRemoveTable[i]);
      }
    }
    FIDToRemoveTable.clear();
    m_SQLRequest = "REPACK " + PlotsAndEntitiesUnionLayerName;
    Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
    OGRDataSource::DestroyDataSource(Union);
    Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    UnionLayer->ResetReading();
    UnionFeature = nullptr;
    FIDSliverTable.clear();
    while ((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
    {
      if (UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str()) != UnionFeature->GetFieldAsInteger("IDPlot") and std::find(MissingIDTable.begin(), MissingIDTable.end(), UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str())) == MissingIDTable.end())
      {
        FIDSliverTable.push_back(UnionFeature->GetFID());
      }
    }
    OGRDataSource::DestroyDataSource(Union);
  }

  // ==================================================================================
  // Regrouping entities that do not have any entities that they can be aggregated with
  // ==================================================================================

  Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  FeatureDefn = UnionLayer->GetLayerDefn();
  NewFeature = nullptr;

  if (MissingIDTable.size() != 0)
  {
    for (int i = 0; i < MissingIDTable.size(); i++)
    {
      SQLFeature = nullptr;
      SQLGeometry = nullptr;
      NewGeometry = nullptr;
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + PlotsAndEntitiesUnionLayerName + " WHERE " + m_IDFieldName + " == " + std::to_string(MissingIDTable[i]) + " ORDER BY ST_Area(geometry) DESC";
      SQLLayer = Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      IDField = SQLLayer->GetFeature(0)->GetFieldAsInteger(m_IDFieldName.c_str());
      IDPlot2 = SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot2");
      ID2 = SQLLayer->GetFeature(0)->GetFieldAsString("ID2");
      IDPlot3 = SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot3");
      ID3 = SQLLayer->GetFeature(0)->GetFieldAsString("ID3");
      FIDToRemoveTable.clear();
      for(int j = 0; j < SQLLayer->GetFeatureCount(); j++)
      {
        SQLFeature = SQLLayer->GetFeature(j);
        SQLGeometry = SQLFeature->GetGeometryRef();
        FIDToRemoveTable.push_back(SQLFeature->GetFieldAsInteger("RID"));
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
      if (NewGeometry->getGeometryType() == 6)
      {
        for (int k = 0; k < ((OGRMultiPolygon*) NewGeometry)->getNumGeometries(); k++)
        {
          IDNew = std::to_string(IDField) + "N" + std::to_string(k+1);
          NewFeature = OGRFeature::CreateFeature(FeatureDefn);
          NewFeature->SetField("OFLD_ID", IDField);
          NewFeature->SetField("IDPlot", IDField);
          NewFeature->SetField("ID", IDNew.c_str());
          NewFeature->SetField("IDPlot2", IDPlot2);
          NewFeature->SetField("ID2", ID2.c_str());
          NewFeature->SetField("IDPlot3", IDPlot3);
          NewFeature->SetField("ID3", ID3.c_str());
          NewFeature->SetGeometry(((OGRMultiPolygon*) NewGeometry)->getGeometryRef(k));
          UnionLayer->CreateFeature(NewFeature);
          UnionLayer->SyncToDisk();
        }
      }
      else if (NewGeometry->getGeometryType() == 3)
      {
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
      }
      for (int l = 0; l < FIDToRemoveTable.size(); l++)
      {
        UnionLayer->DeleteFeature(FIDToRemoveTable[l]);
      }
      m_SQLRequest = "REPACK " + PlotsAndEntitiesUnionLayerName;
      Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
      Union->SyncToDisk();
    }
  }

  OGRDataSource::DestroyDataSource(Union);

  IDPlotTable.clear();
  MissingIDTable.clear();
  FIDSliverTable.clear();
  FIDToRemoveTable.clear();

  // ====================================================================================
  // Regrouping of those features that have the same attributes and share a common border
  // ====================================================================================

  Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  UnionFeature = nullptr;
  UnionGeometry = nullptr;

  DataSource = OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );
  SQLLayer = nullptr;
  SQLFeature = nullptr;
  SQLGeometry = nullptr;

  FIDToRemoveTable.clear();

  while ((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    UnionGeometry = UnionFeature->GetGeometryRef();
    NewGeometry = nullptr;
    FID = UnionFeature->GetFID();
    IDField = UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str());
    IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
    ID = UnionFeature->GetFieldAsString("ID");
    if (std::find(FIDToRemoveTable.begin(), FIDToRemoveTable.end(), FID) == FIDToRemoveTable.end())
    {
      m_SQLRequest = "SELECT ROWID as RID, geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE " + m_IDFieldName + " == " + std::to_string(IDField) + " AND "
          "IDPlot == " + std::to_string(IDPlot) + " AND ID == " + m_QMark + ID + m_QMark + " AND "
          "RID != " + std::to_string(FID) + " ORDER BY ST_Area(geometry) DESC";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        for (int i = 0; i < SQLLayer->GetFeatureCount(); i++)
        {
          SQLFeature = SQLLayer->GetFeature(i);
          SQLGeometry = SQLFeature->GetGeometryRef();
          FIDToRemoveTable.push_back(SQLFeature->GetFieldAsInteger("RID"));
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

  if (FIDToRemoveTable.size() != 0)
  {
    for (int i = 0; i < FIDToRemoveTable.size(); i++)
    {
      UnionLayer->DeleteFeature(FIDToRemoveTable[i]);
    }
  }

  FIDToRemoveTable.clear();
  m_SQLRequest = "REPACK " + PlotsAndEntitiesUnionLayerName;
  Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
  OGRDataSource::DestroyDataSource(Union);

  // =======================================================================================
  // In case 0N1 is of MULTIPOLYGON type, we separate it into POLYGONES, create new features
  // for them and (except for the largest one) set attributes as those of the feature that it
  // shares the longest boundary with
  // =======================================================================================

  Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  UnionFeature = nullptr;
  UnionGeometry = nullptr;

  NewGeometry = nullptr;

  FIDToRemoveTable.clear();

  while ((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    FID = UnionFeature->GetFID();
    ID = UnionFeature->GetFieldAsString("ID");
    UnionGeometry = UnionFeature->GetGeometryRef();
    if (ID == "0N1" and UnionGeometry->getGeometryType() == 6)
    {
      for (int i = 0; i < ((OGRMultiPolygon *) UnionGeometry)->getNumGeometries(); i++)
      {
        NewGeometry = ((OGRMultiPolygon*) UnionGeometry)->getGeometryRef(i);
        NewFeature = OGRFeature::CreateFeature(UnionLayer->GetLayerDefn());
        NewFeature->SetField(m_IDFieldName.c_str(), 0);
        NewFeature->SetField("IDPlot", 0);
        NewFeature->SetField("ID", ID.c_str());
        NewFeature->SetField("IDPlot2", 0);
        NewFeature->SetField("ID2", ID.c_str());
        NewFeature->SetField("IDPlot3", 0);
        NewFeature->SetField("ID3", ID.c_str());
        NewFeature->SetGeometry(NewGeometry);
        UnionLayer->CreateFeature(NewFeature);
      }
      FIDToRemoveTable.push_back(FID);
    }
  }

  UnionLayer->SyncToDisk();
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  if (FIDToRemoveTable.size() != 0)
  {
    for (int i = 0; i < FIDToRemoveTable.size(); i++)
    {
      UnionLayer->DeleteFeature(FIDToRemoveTable[i]);
    }
  }

  FIDToRemoveTable.clear();

  m_SQLRequest = "REPACK " + PlotsAndEntitiesUnionLayerName;
  Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
  OGRDataSource::DestroyDataSource(Union);

  Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  DataSource = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );

  m_SQLRequest = "SELECT ROWID as RID FROM "
      "" + PlotsAndEntitiesUnionLayerName + " "
      "WHERE ID == " + m_QMark + "0N1" + m_QMark + " ORDER BY ST_Area(geometry) DESC";
  SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());

  if (SQLLayer->GetFeatureCount() > 1)
  {
    for (int i = 1; i < SQLLayer->GetFeatureCount(); i++)
    {
      FIDTable.push_back(SQLLayer->GetFeature(i)->GetFieldAsInteger("RID"));
    }
    DataSource->ReleaseResultSet(SQLLayer);
    for (int j = 0; j < FIDTable.size(); j++)
    {
      m_SQLRequest = "SELECT * FROM " + PlotsAndEntitiesUnionLayerName + " "
          "WHERE ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FIDTable[j]) + "))) > 0 "
          "AND ID != " + m_QMark + "0N1" + m_QMark + " "
          "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FIDTable[j]) + "))) DESC";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        UnionFeature = UnionLayer->GetFeature(FIDTable[j]);
        UnionFeature->SetField(m_IDFieldName.c_str(), SQLLayer->GetFeature(0)->GetFieldAsInteger(m_IDFieldName.c_str()));
        UnionFeature->SetField("IDPlot", SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot"));
        UnionFeature->SetField("ID", SQLLayer->GetFeature(0)->GetFieldAsString("ID"));
        UnionFeature->SetField("IDPlot2", SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot2"));
        UnionFeature->SetField("ID2", SQLLayer->GetFeature(0)->GetFieldAsString("ID2"));
        UnionFeature->SetField("IDPlot3", SQLLayer->GetFeature(0)->GetFieldAsInteger("IDPlot3"));
        UnionFeature->SetField("ID3", SQLLayer->GetFeature(0)->GetFieldAsString("ID3"));
        UnionLayer->SetFeature(UnionFeature);
      }
      DataSource->ReleaseResultSet(SQLLayer);
    }

  }
  else
  {
    DataSource->ReleaseResultSet(SQLLayer);
  }

  OGRDataSource::DestroyDataSource(Union);
  OGRDataSource::DestroyDataSource(DataSource);

  // =======================================================================================
  // If entities of the MULTIPOLYGON type that do not touch each other were formed, we split
  // them into separate entities of the POLYGON type
  // =======================================================================================

  Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  IDPlotIDMaxTable[0].clear();
  IDPlotIDMaxTable[1].clear();

  while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
    ID = UnionFeature->GetFieldAsString("ID");
    IDN = std::stoi(ID.substr(ID.find("N")+1, ID.length()-1));
    if (std::find(IDPlotIDMaxTable[0].begin(), IDPlotIDMaxTable[0].end(), IDPlot) == IDPlotIDMaxTable[0].end())
    {
      IDPlotIDMaxTable[0].push_back(IDPlot);
      IDPlotIDMaxTable[1].push_back(IDN);
    }
    else
    {
      pos = std::find(IDPlotIDMaxTable[0].begin(), IDPlotIDMaxTable[0].end(), IDPlot) - IDPlotIDMaxTable[0].begin();
      if (IDPlotIDMaxTable[1][pos] < IDN)
      {
        IDPlotIDMaxTable[1][pos] = IDN;
      }
    }
  }

  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  FIDToRemoveTable.clear();

  while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    if (UnionFeature->GetGeometryRef()->getGeometryType() == 6)
    {
      for (int i = 0; i < ((OGRMultiPolygon*)UnionFeature->GetGeometryRef())->getNumGeometries(); i++)
      {
        NewGeometry = ((OGRMultiPolygon*)UnionFeature->GetGeometryRef())->getGeometryRef(i);
        UnionGeometry = nullptr;
        for (int j = 0; j < ((OGRMultiPolygon*)UnionFeature->GetGeometryRef())->getNumGeometries(); j++)
        {
          if (j != i)
          {
            if (UnionGeometry == nullptr)
            {
              UnionGeometry = ((OGRMultiPolygon*)UnionFeature->GetGeometryRef())->getGeometryRef(j);
            }
            else
            {
              UnionGeometry = UnionGeometry->Union(((OGRMultiPolygon*)UnionFeature->GetGeometryRef())->getGeometryRef(j));
            }
          }
        }
        if (!NewGeometry->Touches(UnionGeometry) and ((OGRPolygon*)NewGeometry)->get_Area() < ((OGRPolygon*)UnionGeometry)->get_Area())
        {
          FIDToRemoveTable.push_back(UnionFeature->GetFID());
          OGRFeature *NewFeature = OGRFeature::CreateFeature(UnionLayer->GetLayerDefn());
          NewFeature->SetField(m_IDFieldName.c_str(), UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str()));
          NewFeature->SetField("IDPlot", UnionFeature->GetFieldAsInteger("IDPlot"));
          NewFeature->SetField("ID", UnionFeature->GetFieldAsString("ID"));
          NewFeature->SetField("IDPlot2", UnionFeature->GetFieldAsInteger("IDPlot2"));
          NewFeature->SetField("ID2", UnionFeature->GetFieldAsString("ID2"));
          NewFeature->SetField("IDPlot3", UnionFeature->GetFieldAsInteger("IDPlot3"));
          NewFeature->SetField("ID3", UnionFeature->GetFieldAsString("ID3"));
          NewFeature->SetGeometry(UnionGeometry);
          UnionLayer->CreateFeature(NewFeature);
          UnionLayer->SyncToDisk();
          Union->SyncToDisk();
          pos = std::find(IDPlotIDMaxTable[0].begin(), IDPlotIDMaxTable[0].end(), UnionFeature->GetFieldAsInteger("IDPlot")) - IDPlotIDMaxTable[0].begin();
          IDMax = IDPlotIDMaxTable[1][pos];
          IDPlotIDMaxTable[1][pos] = IDMax+1;
          IDNew = std::to_string(UnionFeature->GetFieldAsInteger("IDPlot"))+"N"+std::to_string(IDMax+1);
          NewFeature = OGRFeature::CreateFeature(UnionLayer->GetLayerDefn());
          NewFeature->SetField(m_IDFieldName.c_str(), UnionFeature->GetFieldAsInteger(m_IDFieldName.c_str()));
          NewFeature->SetField("IDPlot", UnionFeature->GetFieldAsInteger("IDPlot"));
          NewFeature->SetField("ID", IDNew.c_str());
          NewFeature->SetField("IDPlot2", UnionFeature->GetFieldAsInteger("IDPlot2"));
          NewFeature->SetField("ID2", UnionFeature->GetFieldAsString("ID2"));
          NewFeature->SetField("IDPlot3", UnionFeature->GetFieldAsInteger("IDPlot3"));
          NewFeature->SetField("ID3", UnionFeature->GetFieldAsString("ID3"));
          NewFeature->SetGeometry(NewGeometry);
          UnionLayer->CreateFeature(NewFeature);
          UnionLayer->SyncToDisk();
          Union->SyncToDisk();
        }
      }
    }
  }

  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  for (int i = 0; i < FIDToRemoveTable.size(); i++)
  {
    UnionLayer->DeleteFeature(FIDToRemoveTable[i]);
  }

  m_SQLRequest = "REPACK " + PlotsAndEntitiesUnionLayerName;
  Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);

  OGRDataSource::DestroyDataSource(Union);

  FIDToRemoveTable.clear();
  IDPlotIDMaxTable[0].clear();
  IDPlotIDMaxTable[1].clear();

  // ===============================================================================
  // Checking for missing connections between entities
  // ===============================================================================

  Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  DataSource  = OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );

  FIDConnectionProblemTable.clear();

  while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    ID = UnionFeature->GetFieldAsString("ID");
    ID2 = UnionFeature->GetFieldAsString("ID2");
    if (ID != "0N1")
    {
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + " AND "
          "ST_Touches(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + "))";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() == 0)
      {
        FIDConnectionProblemTable.push_back(UnionFeature->GetFID());
      }
      DataSource->ReleaseResultSet(SQLLayer);
    }
  }

  OGRDataSource::DestroyDataSource(DataSource);
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  if (FIDConnectionProblemTable.size() != 0)
  {
    while (FIDConnectionProblemTable.size() != 0)
    {
      Beginning = FIDConnectionProblemTable.size();
      FIDToRemoveTable.clear();
      DataSource  = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
      for (int i = 0; i < FIDConnectionProblemTable.size(); i++)
      {
        UnionFeature = UnionLayer->GetFeature(FIDConnectionProblemTable[i]);
        IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
        ID = UnionFeature->GetFieldAsString("ID");
        ID2 = UnionFeature->GetFieldAsString("ID2");
        m_SQLRequest = "SELECT IDPlot, ID, IDPlot2, ID2 FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID2 == " + m_QMark + ID2 + m_QMark + " AND "
            "IDPlot != " + std::to_string(IDPlot) + " AND "
            "ST_Touches(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
            "ST_Touches(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")) "
            "ORDER BY ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + "))) DESC,"
            " ST_Length(ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + "))) DESC";
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
          FIDToRemoveTable.push_back(UnionFeature->GetFID());
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID2 == " + m_QMark + ID + m_QMark;
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
      for (int i = 0; i < FIDToRemoveTable.size(); i++)
      {
        if (std::find(FIDConnectionProblemTable.begin(), FIDConnectionProblemTable.end(), FIDToRemoveTable[i]) != FIDConnectionProblemTable.end())
        {
          pos = std::find(FIDConnectionProblemTable.begin(), FIDConnectionProblemTable.end(), FIDToRemoveTable[i]) - FIDConnectionProblemTable.begin();
          FIDConnectionProblemTable.erase(FIDConnectionProblemTable.begin() + pos);
        }
      }
      End = FIDConnectionProblemTable.size();
      if (Beginning == End)
      {
        break;
      }
    }
  }

  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  if (FIDConnectionProblemTable.size() != 0)
  {
    while (FIDConnectionProblemTable.size() != 0)
    {
      Beginning = FIDConnectionProblemTable.size();
      FIDToRemoveTable.clear();
      DataSource  = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
      for (int i = 0; i < FIDConnectionProblemTable.size(); i++)
      {
        UnionFeature = UnionLayer->GetFeature(FIDConnectionProblemTable[i]);
        IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
        ID = UnionFeature->GetFieldAsString("ID");
        ID2 = UnionFeature->GetFieldAsString("ID2");
        m_SQLRequest = "SELECT IDPlot, ID, IDPlot2, ID2 FROM " + PlotsAndEntitiesUnionLayerName + " WHERE IDPlot != " + std::to_string(IDPlot) + " AND "
            "ST_Touches(geometry, (SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
            "ID2 != " + m_QMark + ID + m_QMark + " AND ID3 != " + m_QMark + ID + m_QMark + " ORDER BY "
            "ST_Distance(geometry, (SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")) ASC, "
            "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + "))) DESC";
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
          FIDToRemoveTable.push_back(UnionFeature->GetFID());
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID2 == " + m_QMark + ID + m_QMark;
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
      for (int i = 0; i < FIDToRemoveTable.size(); i++)
      {
        if (std::find(FIDConnectionProblemTable.begin(), FIDConnectionProblemTable.end(), FIDToRemoveTable[i]) != FIDConnectionProblemTable.end())
        {
          pos = std::find(FIDConnectionProblemTable.begin(), FIDConnectionProblemTable.end(), FIDToRemoveTable[i]) - FIDConnectionProblemTable.begin();
          FIDConnectionProblemTable.erase(FIDConnectionProblemTable.begin() + pos);
        }
      }
      OGRDataSource::DestroyDataSource(DataSource);
      End = FIDConnectionProblemTable.size();
      if (Beginning == End)
      {
        break;
      }
    }
  }

  FIDConnectionProblemTable.clear();
  OGRDataSource::DestroyDataSource(Union);

  // =========================================================================================
  // Finding and aggregating entities that are "locked" inside a parcel
  // (could happen as a result of sliver re-attribution)
  // =========================================================================================

  Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();
  DataSource  = OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );

  FIDConnectionProblemTable.clear();

  while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
  {
    IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
    ID = UnionFeature->GetFieldAsString("ID");
    ID2 = UnionFeature->GetFieldAsString("ID2");
    if (ID != "0N1")
    {
      m_SQLRequest = "SELECT ROWID as RID, * FROM " + PlotsAndEntitiesUnionLayerName + " WHERE "
          "ST_Touches(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND IDPlot != " + std::to_string(IDPlot);
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() == 0)
      {
        FIDConnectionProblemTable.push_back(UnionFeature->GetFID());
      }
      DataSource->ReleaseResultSet(SQLLayer);
    }
  }

  UnionLayer = Union->GetLayer(0);
  UnionLayer->ResetReading();

  FIDToRemoveTable.clear();

  if (FIDConnectionProblemTable.size() != 0)
  {
    DataSource  = OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );
    for (int i = 0; i < FIDConnectionProblemTable.size(); i++)
    {
      UnionFeature = UnionLayer->GetFeature(FIDConnectionProblemTable[i]);
      IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
      ID = UnionFeature->GetFieldAsString("ID");
      ID2 = UnionFeature->GetFieldAsString("ID2");
      m_SQLRequest = "SELECT ROWID AS RID, ST_Union(geometry, (SELECT geometry "
          "FROM " + PlotsAndEntitiesUnionLayerName + " "
          "WHERE ROWID == " + std::to_string(FIDConnectionProblemTable[i]) + ")), ROWID AS RID, IDPlot, "
          "ID, IDPlot2, ID2, IDPlot3, ID3 FROM "
          "" + PlotsAndEntitiesUnionLayerName + " WHERE"
          " IDPlot == " + std::to_string(IDPlot) + " AND "
          "ST_Touches(geometry, (SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE "
          "ROWID == " + std::to_string(FIDConnectionProblemTable[i]) + ")) AND "
          "ID2 != " + m_QMark + ID + m_QMark + " AND ID3 != " + m_QMark + ID + m_QMark + " ORDER BY "
          "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " "
          "WHERE ROWID == " + std::to_string(FIDConnectionProblemTable[i]) + "))) DESC";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
        {
          if (std::find(FIDConnectionProblemTable.begin(), FIDConnectionProblemTable.end(), SQLLayer->GetFeature(j)->GetFieldAsInteger("RID")) == FIDConnectionProblemTable.end())
          {
            IDNew = SQLLayer->GetFeature(j)->GetFieldAsString("ID");
            IDPlot2New = SQLLayer->GetFeature(j)->GetFieldAsInteger("IDPlot2");
            ID2New = SQLLayer->GetFeature(j)->GetFieldAsString("ID2");
            IDPlot3New = SQLLayer->GetFeature(j)->GetFieldAsInteger("IDPlot3");
            ID3New = SQLLayer->GetFeature(j)->GetFieldAsString("ID3");
            UnionGeometry = UnionFeature->GetGeometryRef();
            NewFeature = OGRFeature::CreateFeature(UnionLayer->GetLayerDefn());
            NewFeature->SetField(m_IDFieldName.c_str(), IDPlot);
            NewFeature->SetField("IDPlot", IDPlot);
            NewFeature->SetField("ID", IDNew.c_str());
            NewFeature->SetField("IDPlot2", IDPlot2New);
            NewFeature->SetField("ID2", ID2New.c_str());
            NewFeature->SetField("IDPlot3", IDPlot3New);
            NewFeature->SetField("ID3", ID3New.c_str());
            NewFeature->SetGeometry(SQLLayer->GetFeature(j)->GetGeometryRef());
            UnionLayer->CreateFeature(NewFeature);
            FIDToRemoveTable.push_back(SQLLayer->GetFeature(j)->GetFieldAsInteger("RID"));
            FIDToRemoveTable.push_back(FIDConnectionProblemTable[i]);
            DataSource->ReleaseResultSet(SQLLayer);
            m_SQLRequest = "SELECT ROWID as RID FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID2 == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            if (SQLLayer->GetFeatureCount() != 0)
            {
              for (int k = 0; k < SQLLayer->GetFeatureCount(); k++)
              {
                FID = SQLLayer->GetFeature(k)->GetFieldAsInteger("RID");
                UnionFeatureToUpdate = UnionLayer->GetFeature(FID);
                UnionFeatureToUpdate->SetField("ID2", IDNew.c_str());
                UnionFeatureToUpdate->SetField("IDPlot3", IDPlot2New);
                UnionFeatureToUpdate->SetField("ID3", ID2New.c_str());
                UnionLayer->SetFeature(UnionFeatureToUpdate);
              }
            }
            DataSource->ReleaseResultSet(SQLLayer);
            m_SQLRequest = "SELECT ROWID as RID FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID3 == " + m_QMark + ID + m_QMark;
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
            break;
          }
        }
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

  if (FIDToRemoveTable.size() != 0)
  {
    for (int i = 0; i < FIDToRemoveTable.size(); i++)
    {
      UnionLayer->DeleteFeature(FIDToRemoveTable[i]);
    }
  }

  FIDToRemoveTable.clear();
  m_SQLRequest = "REPACK " + PlotsAndEntitiesUnionLayerName;
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
    FIDToRemoveTable.clear();
    IDOldTable.clear();
    IDNewTable.clear();
    IDPlotNewTable.clear();
    Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    UnionLayer->ResetReading();
    DataSource  = OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );
    NewGeometry = nullptr;
    while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
    {
      FID = UnionFeature->GetFID();
      IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
      ID = UnionFeature->GetFieldAsString("ID");
      ID2 = UnionFeature->GetFieldAsString("ID2");
      if (ID != "0N1" and std::find(FIDToRemoveTable.begin(), FIDToRemoveTable.end(), FID) == FIDToRemoveTable.end()
          and std::find(IDPlotNewTable.begin(), IDPlotNewTable.end(), IDPlot) == IDPlotNewTable.end())
      {
        m_SQLRequest = "SELECT geometry, ROWID as RID, * FROM " + PlotsAndEntitiesUnionLayerName + " WHERE "
            "ST_Touches(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FID) + ")) AND "
            "IDPlot == " + std::to_string(IDPlot) + " AND ID2 == " + m_QMark + ID2 + m_QMark + " AND "
            "ROWID != " + std::to_string(FID) + " AND "
            "ST_Touches("
            "ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")), "
            "ST_Intersection((SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FID) + "),"
            "(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")))";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() != 0)
        {
          NewGeometry = nullptr;
          for (int j = 0; j < SQLLayer->GetFeatureCount(); j++)
          {
            if (std::find(FIDToRemoveTable.begin(), FIDToRemoveTable.end(), SQLLayer->GetFeature(j)->GetFieldAsInteger("RID")) == FIDToRemoveTable.end())
            {
              FIDToRemoveTable.push_back(SQLLayer->GetFeature(j)->GetFieldAsInteger("RID"));
              IDOldTable.push_back(SQLLayer->GetFeature(j)->GetFieldAsString("ID"));
              IDNewTable.push_back(ID.c_str());
              IDPlotNewTable.push_back(IDPlot);
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
    if (FIDToRemoveTable.size() != 0)
    {
      for (int i = 0; i < FIDToRemoveTable.size(); i++)
      {
        UnionLayer->DeleteFeature(FIDToRemoveTable[i]);
      }
    }
    FIDToRemoveTable.clear();
    m_SQLRequest = "REPACK " + PlotsAndEntitiesUnionLayerName;
    Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
    Union->SyncToDisk();
    OGRDataSource::DestroyDataSource(Union);
    Union  = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    UnionLayer->ResetReading();
    while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
    {
      if (std::find(IDOldTable.begin(), IDOldTable.end(), UnionFeature->GetFieldAsString("ID2")) != IDOldTable.end())
      {
        pos = std::find(IDOldTable.begin(), IDOldTable.end(), UnionFeature->GetFieldAsString("ID2")) - IDOldTable.begin();
        UnionFeature->SetField("ID2", IDNewTable[pos].c_str());
        UnionFeature->SetField("IDPlot2", IDPlotNewTable[pos]);
        UnionLayer->SetFeature(UnionFeature);
      }
    }
    UnionLayer->SyncToDisk();
    Union->SyncToDisk();
    UnionLayer = Union->GetLayer(0);
    UnionLayer->ResetReading();
    while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
    {
      if (std::find(IDOldTable.begin(), IDOldTable.end(), UnionFeature->GetFieldAsString("ID3")) != IDOldTable.end())
      {
        pos = std::find(IDOldTable.begin(), IDOldTable.end(), UnionFeature->GetFieldAsString("ID3")) - IDOldTable.begin();
        UnionFeature->SetField("ID3", IDNewTable[pos].c_str());
        UnionFeature->SetField("IDPlot3", IDPlotNewTable[pos]);
        UnionLayer->SetFeature(UnionFeature);
      }
    }
    UnionLayer->SyncToDisk();
    Union->SyncToDisk();
    OGRDataSource::DestroyDataSource(Union);
    Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    UnionLayer->ResetReading();
    DataSource  = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    N = 0;
    while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
    {
      FID = UnionFeature->GetFID();
      IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
      ID = UnionFeature->GetFieldAsString("ID");
      ID2 = UnionFeature->GetFieldAsString("ID2");
      if (ID != "0N1")
      {
        m_SQLRequest = "SELECT ROWID as RID FROM " + PlotsAndEntitiesUnionLayerName + " WHERE "
            "ST_Touches(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FID) + ")) AND "
            "IDPlot == " + std::to_string(IDPlot) + " AND ID2 == " + m_QMark + ID2 + m_QMark + " AND "
            "ROWID != " + std::to_string(FID) + " AND "
            "ST_Touches("
            "ST_Intersection(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")), "
            "ST_Intersection((SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ROWID == " + std::to_string(FID) + "),"
            "(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + ")))";
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

  FIDToRemoveTable.clear();
  IDOldTable.clear();
  IDNewTable.clear();
  IDPlotNewTable.clear();

  // ==========================================================================
  // Checking for entities that have only two other entities as neighbors and
  // one of them belong to the same plot. In this case, the entity in question
  // will be, if possible, aggregated with the entity of the same plot
  // that it shears the border with in order to avoid fracturing the future LNR
  // ==========================================================================

  FIDToRemoveTable.clear();

  N = 1;

  while (N != 0)
  {
    Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    UnionLayer->ResetReading();
    DataSource  = OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );
    while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
    {
      FID = UnionFeature->GetFID();
      IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
      ID = UnionFeature->GetFieldAsString("ID");
      ID2 = UnionFeature->GetFieldAsString("ID2");
      if (ID != "0N1" and ID2 != "0N1")
      {
        m_SQLRequest = "SELECT ROWID as RID, * FROM " + PlotsAndEntitiesUnionLayerName + " WHERE "
            "ST_Touches(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
            "ROWID != " + std::to_string(FID);
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() <= 2)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID, * FROM " + PlotsAndEntitiesUnionLayerName + " WHERE "
              "ST_Touches(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
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
            FIDToRemoveTable.push_back(FID);
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
    OGRDataSource::DestroyDataSource(Union);
    Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    UnionLayer->ResetReading();
    if (FIDToRemoveTable.size() != 0)
    {
      for (int i = 0; i < FIDToRemoveTable.size(); i++)
      {
        UnionLayer->DeleteFeature(FIDToRemoveTable[i]);
      }
    }
    FIDToRemoveTable.clear();
    m_SQLRequest = "REPACK " + PlotsAndEntitiesUnionLayerName;
    Union->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
    OGRDataSource::DestroyDataSource(Union);
    Union = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    UnionLayer = Union->GetLayer(0);
    UnionLayer->ResetReading();
    DataSource = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str(), TRUE );
    N = 0;
    while((UnionFeature = UnionLayer->GetNextFeature()) != nullptr)
    {
      FID = UnionFeature->GetFID();
      IDPlot = UnionFeature->GetFieldAsInteger("IDPlot");
      ID = UnionFeature->GetFieldAsString("ID");
      ID2 = UnionFeature->GetFieldAsString("ID2");
      if (ID != "0N1" and ID2 != "0N1")
      {
        m_SQLRequest = "SELECT ROWID as RID, * FROM " + PlotsAndEntitiesUnionLayerName + " WHERE "
            "ST_Touches(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
            "ROWID != " + std::to_string(FID);
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() <= 2)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID, * FROM " + PlotsAndEntitiesUnionLayerName + " WHERE "
              "ST_Touches(geometry,(SELECT geometry FROM " + PlotsAndEntitiesUnionLayerName + " WHERE ID == " + m_QMark + ID + m_QMark + ")) AND "
              "IDPlot == " + std::to_string(IDPlot) + " AND ROWID != " + std::to_string(FID);
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() != 0)
          {
            N += 1;
          }
          DataSource->ReleaseResultSet(SQLLayer);
        }
        else
        {
          DataSource->ReleaseResultSet(SQLLayer);
        }
      }
    }
    OGRDataSource::DestroyDataSource(Union);
    OGRDataSource::DestroyDataSource(DataSource);
  }

  // ============================================================================
  // Final regrouping
  // ============================================================================

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputEntitiesGroupedVectorFile)))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputEntitiesGroupedVectorFile).c_str());
  }

  const std::string EntitiesGroupedLayerName = getLayerNameFromFilename(m_OutputEntitiesGroupedVectorFile);

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
  m_SQLRequest = "SELECT ST_Union(geometry), ID, ID2, ID3 FROM " + PlotsAndEntitiesUnionLayerName + " GROUP BY ID, ID2, ID3 ORDER BY IDPlot ASC, ID ASC";
  SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
  FileName = EntitiesGroupedLayerName;
  DataSource->CopyLayer(SQLLayer, FileName.c_str(), nullptr);
  DataSource->ReleaseResultSet(SQLLayer);
  OGRDataSource::DestroyDataSource(DataSource);

  // =========================================================================================
  // Making sure that ID are consecutive (i.e., that if there are only three entities that
  // belong to the same plot number 5, their IDs will be 5N1, 5N2 and 5N3, not 5N2, 5N9, 5N15)
  // =========================================================================================

  Final = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesGroupedVectorFile).c_str(), TRUE );
  FinalLayer = Final->GetLayer(0);
  FinalLayer->ResetReading();

  IDOldTable.clear();
  IDNewTable.clear();

  while ((FinalFeature = FinalLayer->GetNextFeature()) != nullptr)
  {
    FID = FinalFeature->GetFID();
    ID = FinalFeature->GetFieldAsString("ID");
    IDPlot = std::stoi(ID.substr(0, ID.find("N")));
    if (FID == 0)
    {
      if (ID.substr(ID.find("N")+1, ID.length()-1) != "1")
      {
        IDOldTable.push_back(ID.c_str());
        IDNew = std::to_string(IDPlot) + "N1";
        IDNewTable.push_back(IDNew);
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
          IDOldTable.push_back(ID.c_str());
          IDNew = std::to_string(IDPlot) + "N" + std::to_string(std::stoi(IDAbove.substr(IDAbove.find("N")+1, IDAbove.length()-1))+1);
          IDNewTable.push_back(IDNew.c_str());
          FinalFeature->SetField("ID", IDNew.c_str());
          FinalLayer->SetFeature(FinalFeature);
        }
      }
      else
      {
        if (ID.substr(ID.find("N")+1, ID.length()-1) != "1")
        {
          IDOldTable.push_back(ID.c_str());
          IDNew = std::to_string(IDPlot) + "N1";
          IDNewTable.push_back(IDNew);
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
    if (std::find(IDOldTable.begin(), IDOldTable.end(), ID2.c_str()) != IDOldTable.end())
    {
      pos = std::find(IDOldTable.begin(), IDOldTable.end(), ID2.c_str()) - IDOldTable.begin();
      FinalFeature->SetField("ID2", IDNewTable[pos].c_str());
      FinalLayer->SetFeature(FinalFeature);
    }
    if (std::find(IDOldTable.begin(), IDOldTable.end(), ID3.c_str()) != IDOldTable.end())
    {
      pos = std::find(IDOldTable.begin(), IDOldTable.end(), ID3.c_str()) - IDOldTable.begin();
      FinalFeature->SetField("ID3", IDNewTable[pos].c_str());
      FinalLayer->SetFeature(FinalFeature);
    }
  }

  IDOldTable.clear();
  IDNewTable.clear();
  OGRDataSource::DestroyDataSource(Final);

  // ================================================================================
  // Creating SRF.shp and LNR.shp files
  // ================================================================================

  // ================================================================================
  // Creating linear entities vector file - LNR.shp
  // ================================================================================

  Final = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesGroupedVectorFile).c_str(), TRUE );
  FinalLayer = Final->GetLayer(0);
  FinalLayer->ResetReading();
  mp_SRS = Final->GetLayer(0)->GetSpatialRef();

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLinearEntitiesVectorFile)))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str());
  }

  LNR = mp_VectorDriver->CreateDataSource( getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str(), nullptr); // New boundaries file
  LNRLayer = LNR->CreateLayer( "LNR", mp_SRS, wkbLineString, nullptr );

  FieldNamesTable.clear();
  FieldNamesTable = {"FaceA","FaceB","ID","IDTo"};

  for (int i = 0; i < FieldNamesTable.size(); i++)
  {
    createField(LNRLayer, FieldNamesTable[i].c_str(), OFTString);
  }

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputEntitiesGroupedVectorFile).c_str(), TRUE );
  LNRFeature = OGRFeature::CreateFeature(LNRLayer->GetLayerDefn());
  FacesTable.clear();

  while ((FinalFeature = FinalLayer->GetNextFeature()) != nullptr)
  {
    FID = FinalFeature->GetFID();
    FaceA = FinalFeature->GetFieldAsString("ID");
    IDPlot = std::stoi(FaceA.substr(0, FaceA.find("N")));
    if (IDPlot != 0)
    {
      m_SQLRequest = "SELECT ST_CollectionExtract(ST_Intersection(geometry, (SELECT geometry FROM " + EntitiesGroupedLayerName + " WHERE ROWID == " + std::to_string(FID) + ")),2), ID "
          "FROM " + EntitiesGroupedLayerName + " WHERE "
          "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + EntitiesGroupedLayerName + " WHERE ROWID == " + std::to_string(FID) + "))) > 0";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() != 0)
      {
        for (int i = 0; i < SQLLayer->GetFeatureCount(); i++)
        {
          FaceB = SQLLayer->GetFeature(i)->GetFieldAsString("ID");
          FaceBA = FaceB + "-" + FaceA;
          if (std::find(FacesTable.begin(), FacesTable.end(), FaceBA.c_str()) == FacesTable.end())
          {
            LNRFeature->SetField(FieldNamesTable[0].c_str(), FaceA.c_str());
            LNRFeature->SetField(FieldNamesTable[1].c_str(), FaceB.c_str());
            NewGeometry = FinalFeature->GetGeometryRef()->Intersection(SQLLayer->GetFeature(i)->GetGeometryRef());
            LNRFeature->SetGeometry(NewGeometry);
            LNRLayer->CreateFeature(LNRFeature);
            LNR->SyncToDisk();
            FaceAB = FaceA + "-" + FaceB;
            FacesTable.push_back(FaceAB.c_str());
          }
        }
      }
      DataSource->ReleaseResultSet(SQLLayer);
    }
  }

  FieldNamesTable.clear();
  FacesTable.clear();
  OGRFeature::DestroyFeature(LNRFeature);
  OGRDataSource::DestroyDataSource(LNR);
  OGRDataSource::DestroyDataSource(DataSource);
  OGRDataSource::DestroyDataSource(Final);

  // ================================================================================
  // Creating surface entities vector file - SRF.shp
  // ================================================================================

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile)))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str());
  }

  const std::string SurfaceEntitiesLayerName = getLayerNameFromFilename(m_OutputSurfaceEntitiesVectorFile);

  Final = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputEntitiesGroupedVectorFile).c_str(), TRUE );
  m_SQLRequest = "SELECT ST_Union(geometry), ID, ID2 FROM " + EntitiesGroupedLayerName + " GROUP BY ID, ID2";
  FileName = SurfaceEntitiesLayerName;
  SQLLayer = Final->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
  Final->CopyLayer(SQLLayer, FileName.c_str(), nullptr);
  Final->ReleaseResultSet(SQLLayer);
  OGRDataSource::DestroyDataSource(Final);

  FieldNamesTable.clear();
  FieldNamesTable = {"IDTo"};

  SRF = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );
  SRFLayer = SRF->GetLayer(0);
  SRFLayer->ResetReading();

  for (int i = 0; i < FieldNamesTable.size(); i++)
  {
    createField(SRFLayer, FieldNamesTable[i].c_str(), OFTString);
  }

  FieldNamesTable.clear();
  OGRDataSource::DestroyDataSource(SRF);

  // ================================================================================
  // Filling IDTo attribute for the linear entities
  // ================================================================================

  LNR = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str(), TRUE );
  LNRLayer = LNR->GetLayer(0);
  LNRLayer->ResetReading();

  DataSource = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );

  while ((LNRFeature = LNRLayer->GetNextFeature()) != nullptr)
  {
    FaceA = LNRFeature->GetFieldAsString("FaceA");
    FaceB = LNRFeature->GetFieldAsString("FaceB");
    m_SQLRequest = "SELECT ID, ID2 FROM " + SurfaceEntitiesLayerName + " WHERE "
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

  FieldNamesTable.clear();
  FieldNamesTable = {"FaceA", "FaceB"};

  FeatureDefn = LNRLayer->GetLayerDefn();

  if (FieldNamesTable.size() > 0)
  {
    for (int i = 0; i < FieldNamesTable.size(); i++)
    {
      FieldIndex = FeatureDefn->GetFieldIndex(FieldNamesTable[i].c_str());
      LNRLayer->DeleteField(FieldIndex);
    }
  }

  LNR->SyncToDisk();
  FieldNamesTable.clear();

  LNRLayer = LNR->GetLayer(0);

  FieldNamesTable = {"Length"};

  for (int i = 0; i < FieldNamesTable.size(); i++)
  {
    createField(LNRLayer, FieldNamesTable[i].c_str(), OFTReal);
  }

  FieldNamesTable.clear();

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
    m_SQLRequest = "SELECT * FROM " + SurfaceEntitiesLayerName + " WHERE ID == " + m_QMark + ID2 + m_QMark + " AND "
        "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + SurfaceEntitiesLayerName + " WHERE ROWID == " + std::to_string(FID) + "))) > 0";
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

  FieldNamesTable.clear();
  FieldNamesTable = {"ID2"};

  FeatureDefn = SRFLayer->GetLayerDefn();

  if (FieldNamesTable.size() > 0)
  {
    for (int i = 0; i < FieldNamesTable.size(); i++)
    {
      FieldIndex = FeatureDefn->GetFieldIndex(FieldNamesTable[i].c_str());
      SRFLayer->DeleteField(FieldIndex);
    }
  }

  SRF->SyncToDisk();
  FieldNamesTable.clear();

  SRFLayer = SRF->GetLayer(0);
  SRFLayer->ResetReading();

  FieldNamesTable.clear();
  FieldNamesTable = {"LandUse", "Surface",  "SlopeMin", "SlopeMax", "SlopeMean"};

  for (int i = 0; i < FieldNamesTable.size(); i++)
  {
    if (FieldNamesTable[i] == "LandUse")
    {
      createField(SRFLayer, FieldNamesTable[i].c_str(), OFTString);
    }
    else
    {
      createField(SRFLayer, FieldNamesTable[i].c_str(), OFTReal);
    }
  }

  FieldNamesTable.clear();
  OGRDataSource::DestroyDataSource(SRF);
  OGRDataSource::DestroyDataSource(DataSource);

  m_VectorFilesToRelease.push_back(m_OutputSurfaceEntitiesVectorFile);
  m_VectorFilesToRelease.push_back(m_OutputLinearEntitiesVectorFile);
}


// ====================================================================
// ====================================================================


/**
  @internal
  <hr>
  Files
    - m_InputPlotsVectorFile
    - m_OutputSlopeRasterFile
    - m_OutputSurfaceEntitiesVectorFile
*/
void LandProcessor::setSRFParameters()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int IDPlot, XOff, YOff, XCount, YCount, NPoints;
  double Surface, XOrigin, YOrigin, XMin, XMax, YMin, YMax, PixelWidth, PixelHeight, SlopeMin, SlopeMax, SlopeMean;

  std::string ID, IDTo, LandUseValue;

  std::vector <std::string> FieldNamesTable;
  std::vector <double> XPointsTable, YPointsTable, SlopeValuesTable;

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

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile)))
  {
    throw std::runtime_error("LandProcessor::setSRFParameters(): " + m_OutputSurfaceEntitiesVectorFile + ": no such file in the output vector directory");
  }

  if (!openfluid::tools::Filesystem::isFile(getOutputRasterPath(m_OutputSlopeRasterFile)))
  {
    throw std::runtime_error("LandProcessor::setSRFParameters(): " + m_OutputSlopeRasterFile + ": no such file in the output raster directory");
  }

  // ==================================================================================
  // Setting surface parameter for SRF (surface entities)
  // ==================================================================================

  SRF = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );

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

  Plots = OGRSFDriverRegistrar::Open( getInputVectorPath(m_InputPlotsVectorFile).c_str(), TRUE );
  PlotsLayer = Plots->GetLayer(0);
  FeatureDefn = PlotsLayer->GetLayerDefn();

  const std::string InputPlotsLayerName = getLayerNameFromFilename(m_InputPlotsVectorFile);

  FieldNamesTable.clear();

  for (int i = 0; i < FeatureDefn->GetFieldCount(); i++)
  {
    FieldNamesTable.push_back(PlotsLayer->GetLayerDefn()->GetFieldDefn(i)->GetNameRef());
  }

  if(std::find(FieldNamesTable.begin(), FieldNamesTable.end(), m_LandUseFieldName) != FieldNamesTable.end())
  {
    while ((SRFFeature = SRFLayer->GetNextFeature()) != nullptr)
    {
      ID = SRFFeature->GetFieldAsString("ID");
      IDPlot = std::stoi(ID.substr(0, ID.find("N")));
      if (IDPlot != 0)
      {
        m_SQLRequest = "SELECT LandUse FROM " + InputPlotsLayerName + " WHERE " + m_IDFieldName + " == "  + std::to_string(IDPlot);
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

  FieldNamesTable.clear();
  OGRDataSource::DestroyDataSource(Plots);
  OGRDataSource::DestroyDataSource(SRF);

  // ================================================================================
  // Minimum, maximum and mean slope calculations for SRF
  // ================================================================================

  SRF = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );
  SRFLayer = SRF->GetLayer(0);

  Slope = (GDALDataset *) GDALOpen(getOutputRasterPath(m_OutputSlopeRasterFile).c_str(), GA_ReadOnly );
  SlopeBand = Slope->GetRasterBand(1);

  mp_RasterDriver = GetGDALDriverManager()->GetDriverByName("MEM");

  getGeoTransform(Slope);

  XOrigin = m_GeoTransformVal[0];
  YOrigin = m_GeoTransformVal[3];
  PixelWidth = m_GeoTransformVal[1];
  PixelHeight = m_GeoTransformVal[5];

  SRFLayer = SRF->GetLayer(0);
  SRFLayer->ResetReading();

  while((SRFFeature = SRFLayer->GetNextFeature()) != nullptr)
  {
    ID = SRFFeature->GetFieldAsString("ID");
    IDPlot = std::stoi(ID.substr(0, ID.find("N")));
    if (IDPlot != 0)
    {
      SRFGeometry = SRFFeature->GetGeometryRef();
      XPointsTable.clear();
      YPointsTable.clear();
      if (SRFGeometry->getGeometryType() == 3)
      {
        SRFPolygon = (OGRPolygon*) SRFGeometry;
        NPoints = SRFPolygon->getExteriorRing()->getNumPoints();
        for (int i = 0; i < NPoints; i++)
        {
          SRFPolygon->getExteriorRing()->getPoint(i,&SRFPoint);
          XPointsTable.push_back(SRFPoint.getX());
          YPointsTable.push_back(SRFPoint.getY());
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
            XPointsTable.push_back(SRFPoint.getX());
            YPointsTable.push_back(SRFPoint.getY());
          }
        }
      }
      XMin = XPointsTable.at(std::distance(std::begin(XPointsTable), std::min_element(std::begin(XPointsTable), std::end(XPointsTable))));
      XMax = XPointsTable.at(std::distance(std::begin(XPointsTable), std::max_element(std::begin(XPointsTable), std::end(XPointsTable))));
      YMin = YPointsTable.at(std::distance(std::begin(YPointsTable), std::min_element(std::begin(YPointsTable), std::end(YPointsTable))));
      YMax = YPointsTable.at(std::distance(std::begin(YPointsTable), std::max_element(std::begin(YPointsTable), std::end(YPointsTable))));
      XPointsTable.clear();
      YPointsTable.clear();
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
      SlopeValuesTable.clear();
      for (int i = 0; i < Target->GetRasterYSize(); i++)
      {
        for (int j = 0; j < Target->GetRasterXSize(); j++)
        {
          SlopeBand->RasterIO( GF_Read, XOff, YOff+i, Target->GetRasterXSize(), 1, SlopeRow, Target->GetRasterXSize(), 1, GDT_Float32, 0, 0 );
          TargetBand->RasterIO( GF_Read, 0, i, Target->GetRasterXSize(), 1, TargetRow, Target->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
          if (TargetRow[j] != 0)
          {
            SlopeValuesTable.push_back(SlopeRow[j]);
          }
        }
      }
      SlopeMin = SlopeValuesTable.at(std::distance(std::begin(SlopeValuesTable), std::min_element(std::begin(SlopeValuesTable), std::end(SlopeValuesTable))));
      SlopeMax = SlopeValuesTable.at(std::distance(std::begin(SlopeValuesTable), std::max_element(std::begin(SlopeValuesTable), std::end(SlopeValuesTable))));
      SlopeMean = 0;
      for (int k = 0; k < SlopeValuesTable.size(); k++)
      {
        SlopeMean += SlopeValuesTable[k];
      }
      SlopeMean = SlopeMean/SlopeValuesTable.size();
      SlopeValuesTable.clear();
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

}


// =====================================================================
// =====================================================================


/**
  @internal
  <hr>
  Files
    - m_OutputLinearEntitiesVectorFile
*/
void LandProcessor::setLNRParameters()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  double Length;

  OGRDataSource *LNR;
  OGRLayer *LNRLayer;
  OGRFeature *LNRFeature;

  OGRMultiLineString* LNRMultiLineString;

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLinearEntitiesVectorFile)))
  {
    throw std::runtime_error("LandProcessor::setLNRParameters(): " + m_OutputLinearEntitiesVectorFile + ": no such file in the output vector directory");
  }

  // ================================================================================
  //  Calculate length for LNR (linear entities)
  // ================================================================================

  LNR = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str(), TRUE );
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

}


// =====================================================================
// =====================================================================

/**
  @internal
  <hr>
  Files
    - m_OutputSurfaceEntitiesVectorFile
    - m_OutputSUVectorFile
*/
void LandProcessor::createSU()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  std::string FieldName, FileName;

  OGRDataSource *SRF, *SU;
  OGRLayer *SQLLayer, *SRFLayer, *SULayer;
  OGRFeature *SRFFeature, *SUFeature;
  OGRGeometry *SUGeometry;
  OGRFeatureDefn *FeatureDefn;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile)))
  {
    throw std::runtime_error("LandProcessor::createSU(): " + m_OutputSurfaceEntitiesVectorFile + ": no such file in the output vector directory");
  }

  // ====================================================================
  // Creating the SU file
  // ====================================================================

  SRF = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputSurfaceEntitiesVectorFile).c_str(), TRUE );

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputSUVectorFile)))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputSUVectorFile).c_str());
  }

  const std::string SurfaceEntitiesLayerName = getLayerNameFromFilename(m_OutputSurfaceEntitiesVectorFile);
  const std::string SULayerName = getLayerNameFromFilename(m_OutputSUVectorFile);

  m_SQLRequest = "SELECT * FROM " + SurfaceEntitiesLayerName + " WHERE ID != " + m_QMark + "0N1" + m_QMark;
  FileName = SULayerName;
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

  m_VectorFilesToRelease.push_back(m_OutputSUVectorFile);

}


// =====================================================================
// =====================================================================

/**
  @internal
  <hr>
  Files
    - m_InputDitchesVectorFile
    - m_InputThalwegsVectorFile
    - m_InputRiversVectorFile
    - m_OutputLinearEntitiesVectorFile
    - m_OutputRSVectorFile
*/
void LandProcessor::createRS()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int FID;

  double Value, LSLength, Length;

  std::string ID, IDTo;

  std::vector <std::string> FieldNamesTable, FieldNamesLTable, FileNamesTab =
    {m_InputDitchesVectorFile, m_InputThalwegsVectorFile, m_InputRiversVectorFile};

  OGRDataSource *DataSource, *LNR, *RS;
  OGRLayer *SQLLayer, *LNRLayer, *RSLayer;
  OGRFeature *LNRFeature, *RSFeature;
  OGRGeometry *RSGeometry, *IntersectionGeometry;
  OGRFeatureDefn *FeatureDefn;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLinearEntitiesVectorFile)))
  {
    throw std::runtime_error("LandProcessor::createRS(): " + m_OutputLinearEntitiesVectorFile + ": no such file in the output vector directory");
  }


  if (!openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputDitchesVectorFile)) and
      !openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputThalwegsVectorFile)) and
      !openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputRiversVectorFile)))
  {
    std::cout << "LandProcessor::createRS(): There are no linear structure files in the input vector directory that could be used to create RS vector" << std::endl;
  }
  else
  {

    // ====================================================================
    // Creating the RS file
    // ====================================================================

    LNR = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str(), TRUE );

    if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputRSVectorFile)))
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

    FieldNamesTable.clear();

    RSLayer = RS->GetLayer(0);
    RSLayer->ResetReading();

    FieldNamesTable = {"FlowDist", "Ditches", "Thalwegs", "WaterCs", "DitchesL", "ThalwegsL", "WaterCsL", "SurfToLen"};

    for (int i = 0; i < FieldNamesTable.size(); i++)
    {
      createField(RSLayer, FieldNamesTable[i].c_str(), OFTReal);
    }

    RSLayer = RS->GetLayer(0);
    RSLayer->ResetReading();

    while ((RSFeature = RSLayer->GetNextFeature()) != nullptr)
    {
      for (int i = 0; i < FieldNamesTable.size(); i++)
      {
        RSFeature->SetField(FieldNamesTable[i].c_str(), 0);
        RSLayer->SetFeature(RSFeature);
      }
    }

    FieldNamesTable.clear();

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

    const std::string RSLayerName = getLayerNameFromFilename(m_OutputRSVectorFile);

    m_SQLRequest = "REPACK " + RSLayerName;
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

    FieldNamesTable.clear();
    FieldNamesLTable.clear();

    FieldNamesTable = {"Ditches", "Thalwegs", "WaterCs"};
    FieldNamesLTable = {"DitchesL", "ThalwegsL", "WaterCsL"};

    for (int i = 0; i < FileNamesTab.size(); i++)
    {
      if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(FileNamesTab[i])))
      {
        DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
        RS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );
        RSLayer = RS->GetLayer(0);
        const std::string FileLayerName = FileNamesTab[i].substr(0, FileNamesTab[i].find(".shp"));

        for (int j = 0; j < RSLayer->GetFeatureCount(); j++)
        {
          RSFeature = RSLayer->GetFeature(j);
          RSGeometry = RSFeature->GetGeometryRef();
          m_SQLRequest = "SELECT ST_Length(ST_Intersection(geometry,(SELECT ST_Buffer(geometry,0.05) FROM " + RSLayerName + " WHERE ROWID == " + std::to_string(RSFeature->GetFID()) + "))) AS Length "
              "FROM " + FileLayerName + " "
              "WHERE ST_Length(ST_Intersection(geometry,(SELECT ST_Buffer(geometry,0.05) FROM " + RSLayerName + " WHERE ROWID == " + std::to_string(RSFeature->GetFID()) + "))) > 0.1";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() > 0)
          {
            Length = 0;
            for (int k = 0; k < SQLLayer->GetFeatureCount(); k++)
            {
              Length = Length + SQLLayer->GetFeature(k)->GetFieldAsDouble("Length");
            }
            RSFeature->SetField(FieldNamesLTable[i].c_str(), Length);
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
      for (int i = 0; i < FieldNamesLTable.size(); i++)
      {
        Length = RSFeature->GetFieldAsDouble("Length");
        LSLength = RSFeature->GetFieldAsDouble(FieldNamesLTable[i].c_str());
        if (LSLength > Length)
        {
          RSFeature->SetField(FieldNamesLTable[i].c_str(), Length);
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
      for (int i = 0; i < FieldNamesLTable.size(); i++)
      {
        Length = RSFeature->GetFieldAsDouble("Length");
        LSLength = RSFeature->GetFieldAsDouble(FieldNamesLTable[i].c_str());
        if (LSLength != Length and LSLength < 0.05*6)
        {
          RSFeature->SetField(FieldNamesLTable[i].c_str(), 0);
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

    m_SQLRequest = "REPACK " + RSLayerName;
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
      for (int i = 0; i < FieldNamesLTable.size(); i++)
      {
        Length = RSFeature->GetFieldAsDouble("Length");
        LSLength = RSFeature->GetFieldAsDouble(FieldNamesLTable[i].c_str());
        RSFeature->SetField(FieldNamesTable[i].c_str(), LSLength/Length);
        RSLayer->SetFeature(RSFeature);
      }
      RSLayer->SyncToDisk();
    }

    FileNamesTab.clear();
    FieldNamesTable.clear();
    FieldNamesLTable.clear();

    RSLayer->SyncToDisk();
    OGRDataSource::DestroyDataSource(RS);
  }

  m_VectorFilesToRelease.push_back(m_OutputRSVectorFile);

}


// ====================================================================
// ====================================================================


/**
  @internal
  <hr>
  Files
    - m_InputHedgesVectorFile
    - m_InputGrassBandVectorFile
    - m_InputBenchesVectorFile
    - m_OutputLinearEntitiesVectorFile
    - m_OutputRSVectorFile
    - m_OutputLIVectorFile
*/
void LandProcessor::createLI()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int FID;

  double Value, LSLength, Length;

  std::string ID, IDTo;

  std::vector <std::string> FieldNamesTable, FieldNamesLTable, FileNamesTab =
    {m_InputHedgesVectorFile, m_InputGrassBandVectorFile, m_InputBenchesVectorFile};

  OGRDataSource *DataSource, *LNR, *LI;
  OGRLayer *SQLLayer, *LNRLayer, *LILayer;
  OGRFeature *LNRFeature, *LIFeature;
  OGRGeometry *LIGeometry, *IntersectionGeometry;
  OGRFeatureDefn *FeatureDefn;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLinearEntitiesVectorFile)))
  {
    throw std::runtime_error("LandProcessor::createLI(): " + m_OutputLinearEntitiesVectorFile + ": no such file in the output vector directory");
  }

  if (!openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputHedgesVectorFile)) and
      !openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputGrassBandVectorFile)) and
      !openfluid::tools::Filesystem::isFile(getInputVectorPath(m_InputBenchesVectorFile)))
  {
    std::cout << "LandProcessor::createLI(): There are no linear structure files in the input vector directory that could be used to create LI vector" << std::endl;
  }
  else
  {

    // ====================================================================
    // Creating the LI file
    // ====================================================================

    LNR = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLinearEntitiesVectorFile).c_str(), TRUE );

    if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLIVectorFile)))
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

    FieldNamesTable.clear();

    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    FieldNamesTable = {"FlowDist", "Hedges", "GrassBs", "Benches", "HedgesL", "GrassBsL", "BenchesL", "SurfToLen"};

    for (int i = 0; i < FieldNamesTable.size(); i++)
    {
      createField(LILayer, FieldNamesTable[i].c_str(), OFTReal);
    }

    LILayer = LI->GetLayer(0);
    LILayer->ResetReading();

    while ((LIFeature = LILayer->GetNextFeature()) != nullptr)
    {
      for (int i = 0; i < FieldNamesTable.size(); i++)
      {
        LIFeature->SetField(FieldNamesTable[i].c_str(), 0);
        LILayer->SetFeature(LIFeature);
      }
    }

    FieldNamesTable.clear();

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

    const std::string LILayerName = getLayerNameFromFilename(m_OutputLIVectorFile);
    const std::string RSLayerName = getLayerNameFromFilename(m_OutputRSVectorFile);

    m_SQLRequest = "REPACK " + LILayerName;
    LI->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
    OGRDataSource::DestroyDataSource(LI);

    // ====================================================================
    // Setting IDTo attribute in case when there is an RS with the same ID
    // ====================================================================

    if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputRSVectorFile)))
    {
      DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputRSVectorFile).c_str(), TRUE );

      LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
      LILayer = LI->GetLayer(0);

      while ((LIFeature = LILayer->GetNextFeature()) != nullptr)
      {
        ID = LIFeature->GetFieldAsString("ID");
        m_SQLRequest = "SELECT ID FROM " + RSLayerName + " WHERE ID == " + m_QMark + ID + m_QMark;
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

    FieldNamesTable.clear();
    FieldNamesLTable.clear();

    FieldNamesTable = {"Hedges", "GrassBs", "Benches"};
    FieldNamesLTable = {"HedgesL", "GrassBsL", "BenchesL"};

    for (int i = 0; i < FileNamesTab.size(); i++)
    {
      if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(FileNamesTab[i])))
      {
        DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );
        LI = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLIVectorFile).c_str(), TRUE );
        LILayer = LI->GetLayer(0);
        for (int j = 0; j < LILayer->GetFeatureCount(); j++)
        {
          LIFeature = LILayer->GetFeature(j);
          LIGeometry = LIFeature->GetGeometryRef();
          m_SQLRequest = "SELECT ST_Length(ST_Intersection(geometry,(SELECT ST_Buffer(geometry,0.05) FROM " + LILayerName + " WHERE ROWID == " + std::to_string(LIFeature->GetFID()) + "))) AS Length "
              "FROM " + FileNamesTab[i].substr(0, FileNamesTab[i].find(".shp")) + " "
              "WHERE ST_Length(ST_Intersection(geometry,(SELECT ST_Buffer(geometry,0.05) FROM " + LILayerName + " WHERE ROWID == " + std::to_string(LIFeature->GetFID()) + "))) > 0.1";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          if (SQLLayer->GetFeatureCount() > 0)
          {
            Length = 0;
            for (int k = 0; k < SQLLayer->GetFeatureCount(); k++)
            {
              Length = Length + SQLLayer->GetFeature(k)->GetFieldAsDouble("Length");
            }
            LIFeature->SetField(FieldNamesLTable[i].c_str(), Length);
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
      for (int i = 0; i < FieldNamesLTable.size(); i++)
      {
        Length = LIFeature->GetFieldAsDouble("Length");
        LSLength = LIFeature->GetFieldAsDouble(FieldNamesLTable[i].c_str());
        if (LSLength > Length)
        {
          LIFeature->SetField(FieldNamesLTable[i].c_str(), Length);
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
      for (int i = 0; i < FieldNamesLTable.size(); i++)
      {
        Length = LIFeature->GetFieldAsDouble("Length");
        LSLength = LIFeature->GetFieldAsDouble(FieldNamesLTable[i].c_str());
        if (LSLength != Length and LSLength < 0.05*6)
        {
          LIFeature->SetField(FieldNamesLTable[i].c_str(), 0);
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

    m_SQLRequest = "REPACK " + LILayerName;
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
      for (int i = 0; i < FieldNamesLTable.size(); i++)
      {
        Length = LIFeature->GetFieldAsDouble("Length");
        LSLength = LIFeature->GetFieldAsDouble(FieldNamesLTable[i].c_str());
        LIFeature->SetField(FieldNamesTable[i].c_str(), LSLength/Length);
        LILayer->SetFeature(LIFeature);
      }
      LILayer->SyncToDisk();
    }

    FileNamesTab.clear();
    FieldNamesTable.clear();
    FieldNamesLTable.clear();

    LILayer->SyncToDisk();
    OGRDataSource::DestroyDataSource(LI);
  }

  m_VectorFilesToRelease.push_back(m_OutputLIVectorFile);

}


// ====================================================================
// ====================================================================


/**
  @internal
  <hr>
  Files
    - m_OutputRSVectorFile
    - m_OutputLIVectorFile
    - m_OutputSUVectorFile
*/
void LandProcessor::setSUParameters()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int FID, N1, N2;

  std::string ID, IDTo, IDNew, IDToNew;

  std::vector <std::string> FieldNamesTable;

  OGRDataSource *DataSource, *SU;
  OGRLayer *SQLLayer, *SULayer;
  OGRFeature *SUFeature;
  OGRGeometry *SUGeometry;
  OGRFeatureDefn *FeatureDefn;

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputSUVectorFile)))
  {
    throw std::runtime_error("LandProcessor::setSUParameters(): " + m_OutputSUVectorFile + ": no such file in the output vector directory");
  }

  // ==================================================================================
  // Set IDTo parameter based on LI# and RS#
  // ==================================================================================

  const std::string RSLayerName = getLayerNameFromFilename(m_OutputRSVectorFile);
  const std::string SULayerName = getLayerNameFromFilename(m_OutputSUVectorFile);
  const std::string LILayerName = getLayerNameFromFilename(m_OutputLIVectorFile);

  SU = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputSUVectorFile).c_str(), TRUE );
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
      if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLIVectorFile)))
      {
        m_SQLRequest = "SELECT * FROM " + LILayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark;
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() == 0)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputRSVectorFile)))
          {
            m_SQLRequest = "SELECT * FROM " + RSLayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark;
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

  SU = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputSUVectorFile).c_str(), TRUE );
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

  DataSource = OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );

  while ((SUFeature = SULayer->GetNextFeature()) != nullptr)
  {
    ID = SUFeature->GetFieldAsString("ID");
    IDTo = SUFeature->GetFieldAsString("IDTo");
    if (IDTo.find("SU") != std::string::npos)
    {
      if (IDTo != "SU#0N1")
      {
        m_SQLRequest = "SELECT ST_Intersection(geometry, "
            "(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + ")) FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark;
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeature(0)->GetGeometryRef()->getGeometryType() == 1)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark + " AND "
              "ST_Within((SELECT ST_Centroid(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark + "),"
              "(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark + "))";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          N1 = SQLLayer->GetFeatureCount();
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + " AND "
              "ST_Within((SELECT ST_Centroid(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
              "(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          N2 = SQLLayer->GetFeatureCount();
          DataSource->ReleaseResultSet(SQLLayer);
          if (N1 == 1 and N2 == 1)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_Centroid(geometry),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "), geometry)))) + "
                "ST_Distance((SELECT ST_Centroid(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "geometry)))) AS Distance FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark; //TODO: simplify complicated SQL request
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
          else if (N1 == 0 and N2 == 0)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_PointOnSurface(geometry),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "), geometry)))) + "
                "ST_Distance((SELECT ST_PointOnSurface(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "geometry)))) AS Distance FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
          else if (N1 == 1 and N2 == 0)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_Centroid(geometry),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "), geometry)))) + "
                "ST_Distance((SELECT ST_PointOnSurface(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "geometry)))) AS Distance FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
          else if (N1 == 0 and N2 == 1)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_PointOnSurface(geometry),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "), geometry)))) + "
                "ST_Distance((SELECT ST_Centroid(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "((SELECT ST_Intersection((SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "geometry)))) AS Distance FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
        }
        else
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark + " AND "
              "ST_Within((SELECT ST_Centroid(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark + "),"
              "(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark + "))";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          N1 = SQLLayer->GetFeatureCount();
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT ROWID as RID FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + " AND "
              "ST_Within((SELECT ST_Centroid(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
              "(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))";
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          N2 = SQLLayer->GetFeatureCount();
          DataSource->ReleaseResultSet(SQLLayer);
          if (N1 == 1 and N2 == 1)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_Centroid(geometry), ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),"
                "(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) + "
                "ST_Distance((SELECT ST_Centroid(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),(SELECT geometry FROM " + SULayerName + " WHERE ID ==" + m_QMark + IDTo + m_QMark + "))), 0.5)) "
                "AS Distance FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
          else if (N1 == 0 and N2 == 0)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_PointOnSurface(geometry), ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),"
                "(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) + "
                "ST_Distance((SELECT ST_PointOnSurface(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) "
                "AS Distance FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
          else if (N1 == 1 and N2 == 0)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_Centroid(geometry), ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),"
                "(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) + "
                "ST_Distance((SELECT ST_PointOnSurface(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) "
                "AS Distance FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark;
            SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
            SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
            SULayer->SetFeature(SUFeature);
            DataSource->ReleaseResultSet(SQLLayer);
          }
          else if (N1 == 0 and N2 == 1)
          {
            m_SQLRequest = "SELECT geometry, ST_Distance(ST_PointOnSurface(geometry), ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),"
                "(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) + "
                "ST_Distance((SELECT ST_Centroid(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "),"
                "ST_Line_Interpolate_Point(ST_LineMerge(ST_Intersection((geometry),(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark + "))), 0.5)) "
                "AS Distance FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark;
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
      m_SQLRequest = "SELECT ROWID as RID FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark + " AND "
          "ST_Within((SELECT ST_Centroid(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark + "),"
          "(SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark + "))";
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() == 1)
      {
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT geometry, ST_Distance(ST_Centroid(geometry), (SELECT geometry FROM " + LILayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark +"))" //Used to be ST_Line_Interpolate_Point(ST_LineMerge(geometry), 0.5) in the attribution
            " AS Distance FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark;
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
        SULayer->SetFeature(SUFeature);
        DataSource->ReleaseResultSet(SQLLayer);
      }
      else if (SQLLayer->GetFeatureCount() == 0)
      {
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT geometry, ST_Distance(ST_PointOnSurface(geometry), (SELECT geometry FROM " + LILayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark +"))"
            " AS Distance FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark; //16/12/2016/17:00 Simplified the selection here
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
      m_SQLRequest = "SELECT ROWID as RID FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark + " AND "
          "ST_Within(ST_Centroid(geometry), (SELECT geometry FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark + "))"; //16/12/2016/17:00 Simplified the selection here
      SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
      if (SQLLayer->GetFeatureCount() == 1)
      {
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT geometry, ST_Distance(ST_Centroid(geometry), (SELECT geometry FROM " + RSLayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark +"))" //TODO: Use ST_Line_Interpolate_Point(ST_LineMerge(geometry), 0.5) in the attribution
            " AS Distance FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark; //16/12/2016/17:00 Simplified the selection here
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        SUFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
        SULayer->SetFeature(SUFeature);
        DataSource->ReleaseResultSet(SQLLayer);
      }
      else if (SQLLayer->GetFeatureCount() == 0)
      {
        DataSource->ReleaseResultSet(SQLLayer);
        m_SQLRequest = "SELECT geometry, ST_Distance(ST_PointOnSurface(geometry), (SELECT geometry FROM " + RSLayerName + " WHERE ID == " + m_QMark + IDTo + m_QMark +"))"
            " AS Distance FROM " + SULayerName + " WHERE ID == " + m_QMark + ID + m_QMark; //16/12/2016/17:00 Simplified the selection here
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

}


// ====================================================================
// ====================================================================


/**
  @internal
  <hr>
  Files
    - m_OutputRSVectorFile
*/
void LandProcessor::setRSParameters()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  double Value, Length, Surface;

  std::string ID, IDTo, IDFrom, IDNew;

  OGRDataSource *DataSource, *RS;
  OGRLayer *SQLLayer, *RSLayer;
  OGRFeature *RSFeature;

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputRSVectorFile)))
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

    const std::string SULayerName = getLayerNameFromFilename(m_OutputSUVectorFile);

    DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );

    while((RSFeature = RSLayer->GetNextFeature()) != nullptr)
    {
      ID = RSFeature->GetFieldAsString("ID");
      IDFrom = "SU#" + ID.substr(0, ID.find("-"));
      m_SQLRequest = "SELECT Surface FROM " + SULayerName + " WHERE ID == " + m_QMark + IDFrom + m_QMark + " "
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
}


// ====================================================================
// ====================================================================


/**
  @internal
  <hr>
  Files
    - m_OutputLIVectorFile
*/
void LandProcessor::setLIParameters()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  double Value, Length, Surface;

  std::string ID, IDTo, IDFrom, IDNew;

  OGRDataSource *DataSource, *LI;
  OGRLayer *SQLLayer, *LILayer;
  OGRFeature *LIFeature;
  OGRGeometry *LIGeometry;

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLIVectorFile)))
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

    const std::string LILayerName = getLayerNameFromFilename(m_OutputLIVectorFile);
    const std::string SULayerName = getLayerNameFromFilename(m_OutputSUVectorFile);

    DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath().c_str(), TRUE );

    while((LIFeature = LILayer->GetNextFeature()) != nullptr)
    {
      ID = LIFeature->GetFieldAsString("ID");
      IDFrom = "SU#" + ID.substr(0, ID.find("-"));
      m_SQLRequest = "SELECT Surface FROM " + SULayerName + " WHERE ID == " + m_QMark + IDFrom + m_QMark + " AND "
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
        m_SQLRequest = "SELECT ROWID as RID FROM " + SULayerName + " WHERE ID == " + m_QMark + "SU#" + IDTo + m_QMark + " AND "
            "ST_Within(ST_Centroid(geometry), geometry)";
        SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
        if (SQLLayer->GetFeatureCount() == 1)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT geometry, ST_Distance((SELECT ST_Centroid(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + "SU#" + IDTo + m_QMark + "),"
              "geometry) AS Distance FROM " + LILayerName + " WHERE ID == " + m_QMark + ID + m_QMark; //16/12/2016/17:00 Simplified the selection here
          SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
          LIFeature->SetField("FlowDist", SQLLayer->GetFeature(0)->GetFieldAsDouble("Distance"));
          LILayer->SetFeature(LIFeature);
          DataSource->ReleaseResultSet(SQLLayer);
        }
        else if (SQLLayer->GetFeatureCount() == 0)
        {
          DataSource->ReleaseResultSet(SQLLayer);
          m_SQLRequest = "SELECT geometry, ST_Distance((SELECT ST_PointOnSurface(geometry) FROM " + SULayerName + " WHERE ID == " + m_QMark + "SU#" + IDTo + m_QMark + "),"
              "geometry) AS Distance FROM " + LILayerName + " WHERE ID == " + m_QMark + ID + m_QMark; //16/12/2016/17:00 Simplified the selection here
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
}


// ====================================================================
// ====================================================================


void LandProcessor::extractPlotsLimits()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int FID, IDPlotA, IDPlotB;

  std::string IDPlotAB, IDPlotBA;

  std::vector <std::string> FieldNamesTable;
  std::vector <std::string> IDPlotTable;

  OGRDataSource *DataSource, *Plots, *PlotsLimits;
  OGRLayer *SQLLayer, *PlotsLayer, *PlotsLimitsLayer;
  OGRFeature *PlotsFeature, *PlotsLimitsFeature;
  OGRGeometry *NewGeometry;
  OGRFeatureDefn *FeatureDefn;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsVectorFile)))
  {
    throw std::runtime_error("LandProcessor::extractPlotsLimits(): " + m_OutputPlotsVectorFile + ": no such file in the output vector directory");
  }

  // ====================================================================
  // Extracting parcel limits in order to obtain the well-attributes LSs
  // ====================================================================

  Plots = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsVectorFile).c_str(), TRUE );
  PlotsLayer = Plots->GetLayer(0);
  PlotsLayer->ResetReading();
  mp_SRS = Plots->GetLayer(0)->GetSpatialRef();

  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsLimitsVectorFile)))
  {
    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputPlotsLimitsVectorFile).c_str());
  }

  PlotsLimits = mp_VectorDriver->CreateDataSource(getOutputVectorPath(m_OutputPlotsLimitsVectorFile).c_str(), nullptr);
  PlotsLimitsLayer = PlotsLimits->CreateLayer("plotslimits", mp_SRS, wkbLineString, nullptr );

  FieldNamesTable.clear();
  FieldNamesTable = {"ID", "IDPlotA","IDPlotB"};

  for (int i = 0; i < FieldNamesTable.size(); i++)
  {
    createField(PlotsLimitsLayer, FieldNamesTable[i].c_str(), OFTInteger);
  }

  DataSource = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputPlotsVectorFile).c_str(), TRUE );
  PlotsLimitsFeature = OGRFeature::CreateFeature(PlotsLimitsLayer->GetLayerDefn());
  IDPlotTable.clear();

  std::string WorkPlotsLayerName = m_OutputPlotsVectorFile.substr(0, m_OutputPlotsVectorFile.find("."));

  while ((PlotsFeature = PlotsLayer->GetNextFeature()) != nullptr)
  {
    FID = PlotsFeature->GetFID();
    IDPlotA = PlotsFeature->GetFieldAsInteger(m_IDFieldName.c_str());
    m_SQLRequest = "SELECT ST_CollectionExtract(ST_Intersection(geometry, "
        "(SELECT geometry FROM " + WorkPlotsLayerName + " WHERE ROWID == " + std::to_string(FID) + ")), 2), " + m_IDFieldName + " FROM " + WorkPlotsLayerName + " WHERE "
        "ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + WorkPlotsLayerName + " WHERE ROWID == " + std::to_string(FID) + "))) > 0";
    SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
    if (SQLLayer->GetFeatureCount() != 0)
    {
      for (int i = 0; i < SQLLayer->GetFeatureCount(); i++)
      {
        IDPlotB = SQLLayer->GetFeature(i)->GetFieldAsInteger(m_IDFieldName.c_str());
        IDPlotBA = std::to_string(IDPlotB) + "-" + std::to_string(IDPlotA);
        if (std::find(IDPlotTable.begin(), IDPlotTable.end(), IDPlotBA) == IDPlotTable.end())
        {
          PlotsLimitsFeature->SetField(FieldNamesTable[0].c_str(), PlotsLimits->GetLayer(0)->GetFeatureCount());
          PlotsLimitsFeature->SetField(FieldNamesTable[1].c_str(), IDPlotA);
          PlotsLimitsFeature->SetField(FieldNamesTable[2].c_str(), IDPlotB);
          NewGeometry = PlotsFeature->GetGeometryRef()->Intersection(SQLLayer->GetFeature(i)->GetGeometryRef());
          PlotsLimitsFeature->SetGeometry(NewGeometry);
          PlotsLimitsLayer->CreateFeature(PlotsLimitsFeature);
          PlotsLimits->SyncToDisk();
          IDPlotAB = std::to_string(IDPlotA) + "-" + std::to_string(IDPlotB);;
          IDPlotTable.push_back(IDPlotAB);
        }
      }
    }
    DataSource->ReleaseResultSet(SQLLayer);
  }

  FieldNamesTable.clear();
  IDPlotTable.clear();
  OGRDataSource::DestroyDataSource(Plots);
  OGRDataSource::DestroyDataSource(DataSource);
  OGRDataSource::DestroyDataSource(PlotsLimits);

}


// ====================================================================
// ====================================================================


/**
  @internal
  <hr>
  Files
    - m_InputDitchesVectorFile
    - m_InputHedgesVectorFile
    - m_InputGrassBandVectorFile
    - m_InputRiversVectorFile
    - m_InputThalwegsVectorFile
    - m_InputBenchesVectorFile
    - m_OutputPlotsVectorFile
    - m_OutputLinearStructureVectorFile
    - m_OutputIntersectionVectorFile
    - m_OutputPlotsLimitsVectorFile
*/
void LandProcessor::attributeLinearStructures()
{

  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  int IDPlot, IDLimit, FID, N;
  double Distance;
  std::string Options, FileName;

  std::vector <std::string> FieldNamesTable, FileNamesTable = {m_InputDitchesVectorFile, m_InputHedgesVectorFile, m_InputGrassBandVectorFile, m_InputRiversVectorFile, m_InputBenchesVectorFile};

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

  // ==================================================================
  // Attribution of the original linear structures to the parcel limits
  // ==================================================================

  const std::string WorkPlotsLimitsLayerName = getLayerNameFromFilename(m_OutputPlotsLimitsVectorFile);

  if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsLimitsVectorFile)))
  {
    throw std::runtime_error("LandProcessor::attributeLinearStructures(): " + m_OutputPlotsLimitsVectorFile + ": no such file in the output vector directory");
  }

  for (int i = 0; i < FileNamesTable.size(); i++)
  {
    if (FileNamesTable[i].empty() or !openfluid::tools::Filesystem::isFile(getInputVectorPath(FileNamesTable[i])))
    {
      std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": no such file in the input vector directory" << std::endl;
    }
    else
    {
      LS = OGRSFDriverRegistrar::Open( getInputVectorPath(FileNamesTable[i]).c_str(), TRUE );
      if (LS->GetDriver()->GetName() != m_VectorDriverName)
      {
        std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": vector file is not in ESRI Shapefile format" << std::endl;
        OGRDataSource::DestroyDataSource(LS);
      }
      else
      {
        LSLayer = LS->GetLayer(0);
        VectorSRS = LSLayer->GetSpatialRef();
        if (!VectorSRS)
        {
          std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": no spatial reference information is provided for the vector file" << std::endl;
          OGRDataSource::DestroyDataSource(LS);
        }
        else
        {
          VectorDataEPSGCode = VectorSRS->GetAttrValue("AUTHORITY",1);
          if (!VectorDataEPSGCode)
          {
            std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": no EPSG code provided for the vector data" << std::endl;
            OGRDataSource::DestroyDataSource(LS);
          }
          else
          {
            if (VectorDataEPSGCode != m_EPSGCode.substr(m_EPSGCode.find(":")+1, m_EPSGCode.length()-1))
            {
              std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": vector data EPSG code does not correspond to the default EPSG code" << std::endl;
              OGRDataSource::DestroyDataSource(LS);
            }
            else
            {

              m_SQLRequest = "REPACK " + FileNamesTable[i].substr(0, FileNamesTable[i].find("."));
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
                std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": one or several features have invalid geometry" << std::endl;
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
                  std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": one or several features have geometry that are not of type 'LINESTRING'" << std::endl;
                  OGRDataSource::DestroyDataSource(LS);
                }
                else
                {
                  LSLayer = LS->GetLayer(0);
                  FeatureDefn = LSLayer->GetLayerDefn();
                  FieldNamesTable.clear();

                  for (int i = 0; i < FeatureDefn->GetFieldCount(); i++)
                  {
                    FieldNamesTable.push_back(FeatureDefn->GetFieldDefn(i)->GetNameRef());
                  }

                  if(std::find(FieldNamesTable.begin(), FieldNamesTable.end(), m_IDFieldName) != FieldNamesTable.end())
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

                  FieldNamesTable.clear();
                  LSLayer->SyncToDisk();

                  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLinearStructureVectorFile)))
                  {
                    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputLinearStructureVectorFile).c_str());
                  }

                  DataSource = mp_VectorDriver->CopyDataSource(LS, getOutputVectorPath(m_OutputLinearStructureVectorFile).c_str(), nullptr);
                  OGRDataSource::DestroyDataSource(LS);
                  OGRDataSource::DestroyDataSource(DataSource);
                  LS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputLinearStructureVectorFile).c_str(), TRUE );
                  mp_SRS = LS->GetLayer(0)->GetSpatialRef();
                  LSLayer = LS->GetLayer(0);
                  FeatureDefn = LSLayer->GetLayerDefn();

                  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputIntersectionVectorFile)))
                  {
                    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputIntersectionVectorFile).c_str());
                  }

                  Intersection = mp_VectorDriver->CreateDataSource(getOutputVectorPath(m_OutputIntersectionVectorFile).c_str());
                  IntersectionLayer = Intersection->CreateLayer("intersection", mp_SRS, wkbLineString);
                  FieldNamesTable.clear();
                  FieldNamesTable = {m_IDFieldName, "IDLS", "IDPlot", "IDLimit", "PointDist"};

                  for (int i = 0; i < FieldNamesTable.size(); i++)
                  {
                    if (i < 4)
                    {
                      createField(IntersectionLayer, FieldNamesTable[i].c_str(), OFTInteger);
                    }
                    else
                    {
                      createField(IntersectionLayer, FieldNamesTable[i].c_str(), OFTReal);
                    }
                  }

                  FieldNamesTable.clear();
                  IntersectionLayer->SyncToDisk();
                  LSLayer = LS->GetLayer(0);
                  Plots = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsVectorFile).c_str(), TRUE );
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
                  Intersection = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputIntersectionVectorFile).c_str(), TRUE );
                  IntersectionLayer = Intersection->GetLayer(0);
                  IntersectionLayer->ResetReading();
                  PlotsLimits = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputPlotsLimitsVectorFile).c_str(), TRUE );
                  PlotsLimitsLayer = PlotsLimits->GetLayer(0);
                  DataSource =  OGRSFDriverRegistrar::Open(getOutputVectorPath().c_str(), TRUE );

                  while ((IntersectionFeature = IntersectionLayer->GetNextFeature()) != nullptr)
                  {
                    FID = IntersectionFeature->GetFID();
                    IDPlot = IntersectionFeature->GetFieldAsInteger("IDPlot");
                    m_SQLRequest = "SELECT *, geometry, ROWID AS RID" // changes here LineInter and EndPointDist, StartPointDist removed
                        " FROM " + WorkPlotsLimitsLayerName + " WHERE IDPlotA == " + std::to_string(IDPlot) + " OR IDPlotB == " + std::to_string(IDPlot) + " ORDER BY "
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

                  Intersection = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputIntersectionVectorFile).c_str(), TRUE );
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
                          "FROM " + WorkPlotsLimitsLayerName + " WHERE ID == " + std::to_string(IDLimit);
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

                  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(FileNamesTable[i])))
                  {
                    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(FileNamesTable[i]).c_str());
                  }

                  Intersection = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputIntersectionVectorFile).c_str(), TRUE );

                  if (Intersection->GetLayer(0)->GetFeatureCount() != 0)
                  {
                    m_SQLRequest = "SELECT ST_Union(geometry), " + m_IDFieldName + ", IDLS FROM intersection GROUP BY " + m_IDFieldName + ", IDLS, IDLimit"; // 15/12/16/14:34 Dobavila m_IDFieldName
                    SQLLayer = Intersection->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
                    FileName = FileNamesTable[i].substr(0, FileNamesTable[i].find(".shp"));
                    Intersection->CopyLayer(SQLLayer, FileName.c_str(), nullptr);
                    Intersection->ReleaseResultSet(SQLLayer);
                  }

                  OGRDataSource::DestroyDataSource(Intersection);
                  OGRDataSource::DestroyDataSource(LS);
                  OGRDataSource::DestroyDataSource(Plots);
                  OGRDataSource::DestroyDataSource(DataSource);

                  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLinearStructureVectorFile)))
                  {
                    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputLinearStructureVectorFile).c_str());
                  }

                  if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputIntersectionVectorFile)))
                  {
                    mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputIntersectionVectorFile).c_str());
                  }
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


void LandProcessor::releaseVectorFile(const std::string& Filename)
{
  VERBOSE_MESSAGE(2,"Releasing vector file : " << Filename);

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());


  // Removing destination file if exists

  OGRDataSource* TestDest = OGRSFDriverRegistrar::Open(getReleaseVectorPath(Filename).c_str(),false);

  if (TestDest)
  {
    OGRDataSource::DestroyDataSource(TestDest);
    mp_VectorDriver->DeleteDataSource(getReleaseVectorPath(Filename).c_str());
  }


  // Perform copy

  OGRDataSource* Source = OGRSFDriverRegistrar::Open(getOutputVectorPath(Filename).c_str(),false);

  if (!Source)
  {
    throw std::runtime_error("LandProcessor::releasVectorFile(): " + Filename + ": no such file in the output vector directory");
  }

  OGRDataSource* Copy = mp_VectorDriver->CopyDataSource(Source,getReleaseVectorPath(Filename).c_str(),nullptr);

  OGRDataSource::DestroyDataSource(Copy);

  OGRDataSource::DestroyDataSource(Source);
}


// =====================================================================
// =====================================================================


void LandProcessor::releaseRasterFile(const std::string& Filename) const
{
  VERBOSE_MESSAGE(2,"Releasing raster file : " << Filename);

  throw std::runtime_error("not implemented");
}


// =====================================================================
// =====================================================================


void LandProcessor::releaseFiles()
{
  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

  for (auto& ReleaseFile : m_VectorFilesToRelease)
  {
    releaseVectorFile(ReleaseFile);
  }

  for (auto& ReleaseFile : m_RasterFilesToRelease)
  {
    releaseRasterFile(ReleaseFile);
  }
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
      (XCoord - m_GeoTransformVal[0])/m_GeoTransformVal[1],
      (YCoord - m_GeoTransformVal[3])/m_GeoTransformVal[5]);

  return Offset;

}


// =====================================================================
// =====================================================================


void LandProcessor::getCentroidPoint(OGRGeometry *Geometry)
{
  delete mp_Point;
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
  GDALGetGeoTransform(Dataset, m_GeoTransformVal);
}


// =====================================================================
// =====================================================================


void LandProcessor::checkLinearVectorData()
{

  // =====================================================================
  // Variables that are used in this block
  // =====================================================================

  std::vector <std::string> FieldNamesTable, FileNamesTable = {m_InputDitchesVectorFile, m_InputHedgesVectorFile, m_InputGrassBandVectorFile, m_InputRiversVectorFile, m_InputBenchesVectorFile, m_InputThalwegsVectorFile};

  OGRDataSource *LS;
  OGRLayer *LSLayer;
  OGRFeature *LSFeature;

  OGRSpatialReference *VectorSRS;
  const char *VectorDataEPSGCode;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  for (int i = 0; i < FileNamesTable.size(); i++)
  {
    if (FileNamesTable[i].empty() or !openfluid::tools::Filesystem::isFile(getInputVectorPath(FileNamesTable[i])))
    {
      std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": no such file in the input vector directory" << std::endl;
    }
    else
    {
      LS = OGRSFDriverRegistrar::Open( getInputVectorPath(FileNamesTable[i]).c_str(), TRUE );
      if (LS->GetDriver()->GetName() != m_VectorDriverName)
      {
        std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": vector file is not in ESRI Shapefile format" << std::endl;
        OGRDataSource::DestroyDataSource(LS);
      }
      else
      {
        LSLayer = LS->GetLayer(0);
        VectorSRS = LSLayer->GetSpatialRef();
        if (!VectorSRS)
        {
          std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": no spatial reference information is provided for the vector file" << std::endl;
          OGRDataSource::DestroyDataSource(LS);
        }
        else
        {
          VectorDataEPSGCode = VectorSRS->GetAttrValue("AUTHORITY",1);
          if (!VectorDataEPSGCode)
          {
            std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": no EPSG code provided for the vector data" << std::endl;
            OGRDataSource::DestroyDataSource(LS);
          }
          else
          {
            if (VectorDataEPSGCode != m_EPSGCode.substr(m_EPSGCode.find(":")+1, m_EPSGCode.length()-1))
            {
              std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": vector data EPSG code does not correspond to the default EPSG code" << std::endl;
              OGRDataSource::DestroyDataSource(LS);
            }
            else
            {

              m_SQLRequest = "REPACK " + FileNamesTable[i].substr(0, FileNamesTable[i].find("."));
              LS->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
              LSLayer = LS->GetLayer(0);
              LSLayer->ResetReading();

              while ((LSFeature = LSLayer->GetNextFeature()) != 0)
              {
                if (LSFeature->GetGeometryRef()->IsValid() == 0)
                {
                  std::cout << "The feature with FID " << LSFeature->GetFID() << " of the vector file " << FileNamesTable[i] << " has invalid geometry" << std::endl;
                }
                if (LSFeature->GetGeometryRef()->IsSimple() == 0)
                {
                  std::cout << "The feature with FID " << LSFeature->GetFID() << " of the vector file " << FileNamesTable[i] << " has geometry that is not simple" << std::endl;
                }
                if (LSFeature->GetGeometryRef()->IsEmpty() == 1)
                {
                  std::cout << "The feature with FID " << LSFeature->GetFID() << " of the vector file " << FileNamesTable[i] << " has empty geometry" << std::endl;
                }
              }

              LSLayer = LS->GetLayer(0);
              LSLayer->ResetReading();

              while ((LSFeature = LSLayer->GetNextFeature()) != 0)
              {
                if (LSFeature->GetGeometryRef()->getGeometryType() != 2)
                {
                  std::cout << "The feature with FID " << LSFeature->GetFID() << " of the vector file " << FileNamesTable[i] << " has geometry that is not of LINESTRING type" << std::endl;
                }
              }
              OGRDataSource::DestroyDataSource(LS);
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

  std::vector <std::string> FieldNamesTable, FileNamesTable = {m_InputPlotsVectorFile};

  OGRDataSource *Polygon;
  OGRLayer *PolygonLayer;
  OGRFeature *PolygonFeature;

  OGRSpatialReference *VectorSRS;
  const char *VectorDataEPSGCode;

  mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

  for (int i = 0; i < FileNamesTable.size(); i++)
  {
    if (FileNamesTable[i].empty() or !openfluid::tools::Filesystem::isFile(getInputVectorPath(FileNamesTable[i])))
    {
      std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": no such file in the input vector directory" << std::endl;
    }
    else
    {
      Polygon = OGRSFDriverRegistrar::Open( getInputVectorPath(FileNamesTable[i]).c_str(), TRUE );
      if (Polygon->GetDriver()->GetName() != m_VectorDriverName)
      {
        std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": vector file is not in ESRI Shapefile format" << std::endl;
        OGRDataSource::DestroyDataSource(Polygon);
      }
      else
      {
        PolygonLayer = Polygon->GetLayer(0);
        VectorSRS = PolygonLayer->GetSpatialRef();
        if (!VectorSRS)
        {
          std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": no spatial reference information is provided for the vector file" << std::endl;
          OGRDataSource::DestroyDataSource(Polygon);
        }
        else
        {
          VectorDataEPSGCode = VectorSRS->GetAttrValue("AUTHORITY",1);
          if (!VectorDataEPSGCode)
          {
            std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": no EPSG code provided for the vector data" << std::endl;
            OGRDataSource::DestroyDataSource(Polygon);
          }
          else
          {
            if (VectorDataEPSGCode != m_EPSGCode.substr(m_EPSGCode.find(":")+1, m_EPSGCode.length()-1))
            {
              std::cout << "LandProcessor::attributeLinearStructures(): " << FileNamesTable[i] << ": vector data EPSG code does not correspond to the default EPSG code" << std::endl;
              OGRDataSource::DestroyDataSource(Polygon);
            }
            else
            {

              m_SQLRequest = "REPACK " + FileNamesTable[i].substr(0, FileNamesTable[i].find("."));
              Polygon->ExecuteSQL(m_SQLRequest.c_str(), nullptr, nullptr);
              PolygonLayer = Polygon->GetLayer(0);
              PolygonLayer->ResetReading();

              while ((PolygonFeature = PolygonLayer->GetNextFeature()) != 0)
              {
                if (PolygonFeature->GetGeometryRef()->IsValid() == 0)
                {
                  std::cout << "The feature with FID " << PolygonFeature->GetFID() << " of the vector file " << FileNamesTable[i] << " has invalid geometry" << std::endl;
                }
                if (PolygonFeature->GetGeometryRef()->Boundary()->IsSimple() == 0)
                {
                  std::cout << "The feature with FID " << PolygonFeature->GetFID() << " of the vector file " << FileNamesTable[i] << " has geometry that is not simple" << std::endl;
                }
                if (PolygonFeature->GetGeometryRef()->IsEmpty() == 1)
                {
                  std::cout << "The feature with FID " << PolygonFeature->GetFID() << " of the vector file " << FileNamesTable[i] << " has empty geometry" << std::endl;
                }
              }

              PolygonLayer = Polygon->GetLayer(0);
              PolygonLayer->ResetReading();

              while ((PolygonFeature = PolygonLayer->GetNextFeature()) != 0)
              {
                if (PolygonFeature->GetGeometryRef()->getGeometryType() != 3)
                {
                  std::cout << "The feature with FID " << PolygonFeature->GetFID() << " of the vector file " << FileNamesTable[i] << " has  geometry that is not of POLYGON type" << std::endl;
                }
              }
              OGRDataSource::DestroyDataSource(Polygon);
            }
          }
        }
      }
    }
  }
}