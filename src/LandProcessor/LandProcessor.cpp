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

//
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
#include <memory>
#include <cstring>
#include <stdexcept>
#include <clocale>
#include <stdio.h>
#include <gdal/ogr_core.h>
#include <gdal/gdal_alg.h>
#include <gdal/gdal_priv.h>
#include <gdal/ogr_spatialref.h>
#include <gdal/ogr_geometry.h>
#include <gdal/cpl_conv.h>

#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/LineString.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/operation/valid/IsValidOp.h>
#include <geos/geom/Point.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/Coordinate.h>
#include <geos/operation/overlay/snap/GeometrySnapper.h>

#include <openfluid/tools/Filesystem.hpp>
#include <openfluid/utils/GrassGISProxy.hpp>
#include <openfluid/base/Environment.hpp>
#include <openfluid/tools/DataHelpers.hpp>
#include <openfluid/landr/GEOSHelpers.hpp>

#include <LandProcessor/LandProcessor.hpp>
#include <LandProcessor/VectorProcessing.hpp>
#include <LandProcessor/RasterProcessing.hpp>
#include <LandProcessor/VectorRegrouping.hpp>
#include <LandProcessor/LinearEntitiesProcessing.hpp>
#include <LandProcessor/ArealEntitiesProcessing.hpp>
#include <LandProcessor/Helpers.hpp>


LandProcessor::LandProcessor(const std::string& InputPath,
							//const std::string& AtworkPath,
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

  if (!openfluid::tools::Filesystem::isDirectory(getInputVectorPath()))
    throw std::runtime_error("LandProcessor::LandProcessor(): input vector directory does not exists (" +
                             getInputVectorPath()+")");

  if (!openfluid::tools::Filesystem::isDirectory(getInputRasterPath()))
    throw std::runtime_error("LandProcessor::LandProcessor(): input raster directory does not exists (" +
                             getInputRasterPath()+")");

  openfluid::tools::Filesystem::makeDirectory(openfluid::base::Environment::getTempDir());

  m_TmpPath = openfluid::tools::Filesystem::makeUniqueSubdirectory(openfluid::base::Environment::getTempDir(),
                                                                   "landprocessor");
  m_GrassTmpPath = m_TmpPath + "/grass";

  openfluid::tools::Filesystem::makeDirectory(m_GrassTmpPath);
  openfluid::utils::GrassGISProxy GRASS(QString::fromStdString(m_GrassTmpPath),
                                        QString::fromStdString(m_GrassLocation));

  GRASS.createLocation(m_EPSGCode.c_str());

  OGRRegisterAll();
  GDALAllRegister();

  std::setlocale(LC_NUMERIC, "C");

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
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	return Filename.substr(0,Filename.find("."));
}


// =====================================================================
// =====================================================================


std::string LandProcessor::getInputVectorPath(const std::string& Filename) const
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::string Path = m_InputPath+"/"+m_VectorDir;

	if (!Filename.empty())
		Path += "/"+Filename;

	return Path;

}


// =====================================================================
// =====================================================================


std::string LandProcessor::getInputRasterPath(const std::string& Filename) const
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::string Path = m_InputPath+"/"+m_RasterDir;

	if (!Filename.empty())
		Path += "/"+Filename;

	return Path;

}


// =====================================================================
// =====================================================================


std::string LandProcessor::getOutputVectorPath(const std::string& Filename) const
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::string Path = m_OutputPath+"/"+m_VectorDir;

	if (!Filename.empty())
		Path += "/"+Filename;

	return Path;

}


// =====================================================================
// =====================================================================


std::string LandProcessor::getOutputRasterPath(const std::string& Filename) const
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::string Path = m_OutputPath+"/"+m_RasterDir;

	if (!Filename.empty())
		Path += "/"+Filename;

	return Path;

}


// =====================================================================
// =====================================================================


std::string LandProcessor::getReleaseVectorPath(const std::string& Filename) const
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::string Path = m_ReleasePath+"/"+m_VectorDir;

	if (!Filename.empty())
		Path += "/"+Filename;

	return Path;
}


// =====================================================================
// =====================================================================


std::string LandProcessor::getReleaseRasterPath(const std::string& Filename) const
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

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

	std::vector <int> Types = {wkbPolygon,wkbMultiPolygon};

	{
		VectorProcessing VP(getInputVectorPath(), m_InputPlotsVectorFile);

		VP.checkEPSG();

		VP.checkGeometryType(Types);

		VP.checkGeometries();

		VP.copyFile(getOutputVectorPath());
	}

	std::vector<std::string> FieldsToKeep;
	FieldsToKeep.push_back(m_LandUseFieldName);
	FieldsToKeep.push_back(m_IDFieldName);

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputPlotsVectorFile);

		VP.correctGeometry();

		VP.createIDField();

		VP.deleteUnnecessaryFields(FieldsToKeep);
	}

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

	openfluid::utils::GrassGISProxy GRASS(QString::fromStdString(m_GrassTmpPath),
                                        QString::fromStdString(m_GrassLocation));

	GRASS.setOutputFile(QString::fromStdString(m_TmpPath)+"/procesrasterdata.out");

	GRASS.setErrorFile(QString::fromStdString(m_TmpPath)+"/processrasterdata.err");

	{
		RasterProcessing RP(getInputRasterPath(), m_InputDEMFile, GA_ReadOnly);

		RP.checkEPSG();

		std::vector<std::pair<double,double>> Extent;

		std::vector<double> GeoTransformValues;

		VectorProcessing VP(getOutputVectorPath(), m_OutputPlotsVectorFile);

		Extent = VP.getLayerExtent();

		GeoTransformValues = RP.getGeoTransformValues();

		GRASS.appendTask("v.in.ogr", {{"input", QString::fromStdString(getOutputVectorPath(m_OutputPlotsVectorFile))},
			{"output", "plotsvector"}, {"snap", QString::fromStdString(m_SnapDistance)}}, {"--o"});

		GRASS.appendTask("r.in.gdal",{{"input",QString::fromStdString(getInputRasterPath(m_InputDEMFile))},
			{"output","dem"}}, {"--o"});

		GRASS.appendTask("g.region", {{"align", QString::fromStdString("dem")},
			{"n", QString::fromStdString(std::to_string(Extent[0].first+GeoTransformValues[1]))},
			{"s", QString::fromStdString(std::to_string(Extent[0].second-GeoTransformValues[1]))},
			{"e", QString::fromStdString(std::to_string(Extent[1].first+GeoTransformValues[1]))},
			{"w", QString::fromStdString(std::to_string(Extent[1].second-GeoTransformValues[1]))},
			{"nsres", QString::fromStdString(std::to_string(GeoTransformValues[1]))},
			{"ewres", QString::fromStdString(std::to_string(GeoTransformValues[1]))}});

		GRASS.appendTask("v.to.rast", {{"input", QString::fromStdString("plotsvector")}, {"output", QString::fromStdString("plotsraster")},
			{"type", QString::fromStdString("area")}, {"use", QString::fromStdString("attr")},
			{"attribute_column", QString::fromStdString(m_IDFieldName)}}, {"--o"});

		GRASS.appendTask("r.to.vect", {{"input", QString::fromStdString("plotsraster")}, {"output", QString::fromStdString("plotsrastervector")},
			{"type", QString::fromStdString("area")}, {"column", QString::fromStdString("IDPlot")}}, {"--o"});

		GRASS.appendTask("v.db.dropcolumn", {{"map", QString::fromStdString("plotsrastervector")}, {"column", QString::fromStdString("label")}});

		VP.deleteFile(getOutputVectorPath(m_OutputPlotsRasterVectorFile));

		GRASS.appendTask("v.out.ogr", {{"input", QString::fromStdString("plotsrastervector")}, {"type", QString::fromStdString("area")},
			{"output", QString::fromStdString(getOutputVectorPath(m_OutputPlotsRasterVectorFile))}, {"format", QString::fromStdString("ESRI_Shapefile")}}, {"-s"});

		GRASS.appendTask("g.remove", {{"type", QString::fromStdString("raster")}, {"name", QString::fromStdString("drainageraster")}}, {"-f"});

		GRASS.appendTask("r.fill.dir", {{"input",QString::fromStdString("dem")}, {"output",QString::fromStdString("demclean")},
			{"direction",QString::fromStdString("direction")}, {"areas",QString::fromStdString("areas")}}, {"--o"});

		RP.deleteFile(getOutputRasterPath(m_OutputProblemAreasFile));

		GRASS.appendTask("r.out.gdal", {{"input",QString::fromStdString("areas")}, {"output", QString::fromStdString(getOutputRasterPath(m_OutputProblemAreasFile))},
			{"format", QString::fromStdString("GTiff")}}, {"--o"});

		if (GRASS.runJob() != 0)
		{
			throw std::runtime_error("LandProcessor::preprocessRasterData() : unable to run GRASS job (see file " + m_TmpPath + "/processrasterdata.err)");
		}
	}

	{
		std::pair <double, double> MinMaxRasterValues;

		RasterProcessing RP(getOutputRasterPath(), m_OutputProblemAreasFile, GA_ReadOnly);

		int Count = 0;

		while (RP.getMinMaxValues().first != 0 and RP.getMinMaxValues().second != 0)
		{

			RP.deleteFile();

			GRASS.appendTask("g.remove", {{"type",QString::fromStdString("raster")}, {"name",QString::fromStdString("dem")}}, {"-f"});

			GRASS.appendTask("g.rename", {{"raster",QString::fromStdString("demclean,dem")}});

			GRASS.appendTask("r.fill.dir", {{"input",QString::fromStdString("dem")}, {"output",QString::fromStdString("demclean")},
				{"direction",QString::fromStdString("direction")}, {"areas",QString::fromStdString("areas")}}, {"--o"});

			GRASS.appendTask("r.out.gdal", {{"input",QString::fromStdString("areas")}, {"output", QString::fromStdString(getOutputRasterPath(m_OutputProblemAreasFile))},
				{"format", QString::fromStdString("GTiff")}}, {"--o"});

			if (GRASS.runJob() != 0)
			{
				throw std::runtime_error("LandProcessor::preprocessRasterData() : unable to run GRASS job (see file " + m_TmpPath + "/processrasterdata.err)");
			}

			RP.getMinMaxValues();

			Count++;

			if (Count >= 5)
			{
				break;
			}
		}

		GRASS.appendTask("r.slope.aspect", {{"elevation", QString::fromStdString("demclean")}, {"slope", QString::fromStdString("slope")}}, {"--o"});

		GRASS.appendTask("r.watershed", {{"elevation", QString::fromStdString("demclean")}, {"drainage", QString::fromStdString("drainageraster")}}, {"-sb"});

		RP.deleteFile(getOutputRasterPath(m_OutputDEMFile));

		GRASS.appendTask("r.out.gdal", {{"input",QString::fromStdString("demclean")}, {"output", QString::fromStdString(getOutputRasterPath(m_OutputDEMFile))},
			{"format", QString::fromStdString("GTiff")}}, {"--o"});

		RP.deleteFile(getOutputRasterPath(m_OutputSlopeRasterFile));

		GRASS.appendTask("r.out.gdal", {{"input",QString::fromStdString("slope")}, {"output", QString::fromStdString(getOutputRasterPath(m_OutputSlopeRasterFile))},
			{"format", QString::fromStdString("GTiff")}}, {"--o"});

		RP.deleteFile(getOutputRasterPath(m_OutputDrainageRasterFile));

		GRASS.appendTask("r.out.gdal", {{"input",QString::fromStdString("drainageraster")}, {"output", QString::fromStdString(getOutputRasterPath(m_OutputDrainageRasterFile))},
			{"format", QString::fromStdString("GTiff")}}, {"--o"});

		RP.deleteFile(getOutputRasterPath(m_OutputPlotsRasterFile));

		GRASS.appendTask("r.out.gdal", {{"input", QString::fromStdString("plotsraster")}, {"output", QString::fromStdString(getOutputRasterPath(m_OutputPlotsRasterFile))},
			{"format", QString::fromStdString("GTiff")}, {"type", QString::fromStdString("Int16")},
			{"nodata", QString::fromStdString("-9999")}}, {"-f"});

		if (GRASS.runJob() != 0)
		{
			throw std::runtime_error("LandProcessor::preprocessRasterData() : unable to run GRASS job (see file " + m_TmpPath + "/processrasterdata.err)");
		}
	}

	{
		RasterProcessing RP(getOutputRasterPath(), m_OutputDEMFile, GA_ReadOnly);

		RP.createIDRaster(m_OutputIDRasterFile);

		RP.createOutletRaster(m_OutputPlotsRasterFile, m_OutputDrainageRasterFile, m_OutputIDRasterFile, m_OutputOutletsRasterFile);

		RP.createReceiversRaster(m_OutputOutletsRasterFile, m_OutputDrainageRasterFile, m_OutputIDRasterFile, m_OutputReceiversRasterFile);

		RP.createConnectingRaster(m_OutputOutletsRasterFile, m_OutputDrainageRasterFile, m_OutputIDRasterFile, m_OutputDownslopeRasterFile);

		RP.createCatchmentsRaster(m_OutputOutletsRasterFile, m_OutputDrainageRasterFile, m_OutputPlotsRasterFile, m_OutputCatchmentsRasterFile);

	}
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
void LandProcessor::createCatchmentsVector()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	// ================================================================================================
	// Creating a catchment vector file based on the catchment raster (IDCat is the catchment ID)
	// ================================================================================================

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputPlotsVectorFile);

		VP.createDataSource(wkbPolygon, getOutputVectorPath(), m_OutputCatchmentsVectorFile);
		VP.createDataSource(wkbPolygon, getOutputVectorPath(), m_OutputOutletsVectorFile);
		VP.createDataSource(wkbPolygon, getOutputVectorPath(), m_OutputReceiversVectorFile);

	}

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputCatchmentsVectorFile);

		std::pair<OGRFieldType,std::string> FieldTypesAndNames;
		FieldTypesAndNames = std::make_pair(OFTInteger,"IDCat");

		VP.polygonizeRaster(FieldTypesAndNames, getOutputRasterPath(), m_OutputCatchmentsRasterFile);

	}

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputOutletsVectorFile);

		std::pair<OGRFieldType,std::string> FieldTypesAndNames;

		FieldTypesAndNames = std::make_pair(OFTInteger,"IDCatO");

		VP.polygonizeRaster(FieldTypesAndNames, getOutputRasterPath(), m_OutputOutletsRasterFile);

		VP.createField(OFTInteger,"IDPlotO");
		VP.createField(OFTInteger,"IDPRO");

		VP.fillFieldFromRaster("IDPlotO", getOutputRasterPath(), m_OutputPlotsRasterFile);
		VP.fillFieldFromRaster("IDPRO", getOutputRasterPath(), m_OutputDownslopeRasterFile);

	}

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputReceiversVectorFile);

		std::pair<OGRFieldType,std::string> FieldTypesAndNames;

		FieldTypesAndNames = std::make_pair(OFTInteger,"IDPRR");

		VP.polygonizeRaster(FieldTypesAndNames, getOutputRasterPath(), m_OutputReceiversRasterFile);

		VP.createField(OFTInteger,"IDPlotR");
		VP.createField(OFTInteger,"IDCatR");

		VP.fillFieldFromRaster("IDPlotR", getOutputRasterPath(), m_OutputPlotsRasterFile);
		VP.fillFieldFromRaster("IDCatR", getOutputRasterPath(), m_OutputCatchmentsRasterFile);

	}

	std::vector<std::string> FilesToProcess = {m_OutputOutletsVectorFile, m_OutputReceiversVectorFile};

	for(unsigned int i = 0; i < FilesToProcess.size(); i++)
	{

		std::vector <unsigned int> FIDList;

		{
			VectorProcessing VP(getOutputVectorPath(), FilesToProcess[i]);

			OGRLayer *Layer = VP.getLayer();
			OGRFeature *Feature = nullptr;

			Layer->ResetReading();

			while ((Feature = Layer->GetNextFeature()) != nullptr)
			{
				if (Feature->GetFieldAsInteger(0) == 0)
				{
					FIDList.push_back(Feature->GetFID());
			}
				OGRFeature::DestroyFeature(Feature);
			}
		}

		if(!FIDList.empty())
		{
			{
				VectorProcessing VP(getOutputVectorPath(), FilesToProcess[i]);

				VP.deleteFeatures(FIDList);

				VP.repack();
			}
		}
	}

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputCatchmentsVectorFile);

		m_SQLRequest = "SELECT ST_Union(geometry), IDCat FROM " + VP.getLayerNameFromFileName(m_OutputCatchmentsVectorFile) + " GROUP BY IDCat";

		VP.executeSQLRequest(m_SQLRequest, m_OutputCatchmentsGroupedVectorFile);

	}

	  // ===============================================================================
	  // Populating catchments vector attribute table
	  // using outlets and receivers vectors
	  // The attributes obtain during this operation will be used for further regrouping
	  // ===============================================================================

	  std::vector < std::vector <int>> OutletsTable, ReceiversTable;

	  {
		  VectorProcessing VP(getOutputVectorPath(), m_OutputOutletsVectorFile);

		  OutletsTable = VP.transformToVector();
	  }

	  {
		  VectorProcessing VP(getOutputVectorPath(), m_OutputReceiversVectorFile);

		  ReceiversTable = VP.transformToVector();
	  }

	  // ==============================================================================
	  // Creating necessary fields in the attribute table and populating them
	  // ("IDPlot", "IDPR", "IDCat2", "IDPlot2", "IDPR2", "IDCat3", "IDPlot3", "IDPR3")
	  // based on previously created OutletsTable and ReceiversTable.
	  // ==============================================================================

	  {
		  VectorProcessing VP(getOutputVectorPath(), m_OutputCatchmentsGroupedVectorFile);

		  std::vector<std::pair<OGRFieldType, std::string>> FieldTypesAndNames;

		  std::vector<std::string> FieldNames = {"IDPlot", "IDPR", "IDCat2", "IDPlot2", "IDPR2", "IDCat3", "IDPlot3", "IDPR3", "ID", "ID2", "ID3"};

		  for (unsigned int i = 0; i < FieldNames.size(); i++)
		  {
			  if(i <= 7)
			  {
				  FieldTypesAndNames.push_back(std::make_pair(OFTInteger, FieldNames[i]));
			  }
			  else
			  {
				  FieldTypesAndNames.push_back(std::make_pair(OFTString, FieldNames[i]));
			  }
		  }

		  for (unsigned int i = 0; i < FieldTypesAndNames.size(); i++)
		  {
			  VP.createField(FieldTypesAndNames[i].first, FieldTypesAndNames[i].second);
		  }
	  }

	  {
		  VectorProcessing VP(getOutputVectorPath(), m_OutputCatchmentsGroupedVectorFile);

		  OGRLayer *Layer = VP.getLayer();
		  OGRFeature *Feature = nullptr;

		  Layer->ResetReading();

		  while ((Feature = Layer->GetNextFeature()) != nullptr)
		  {
			  int IDCat = 0, IDPlot = 0, IDPR = 0, IDCat2 = 0, IDPlot2 = 0, IDPR2 = 0, IDCat3 = 0, IDPlot3 = 0, IDPR3 = 0;
			  IDCat = Feature->GetFieldAsInteger("IDCat");
			  if (IDCat != 0)
			  {
				  if (std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat) != OutletsTable[0].end())
				  {
					  IDPlot =  OutletsTable[1].at(std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat) - OutletsTable[0].begin());
					  IDPR = OutletsTable[2].at(std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat) - OutletsTable[0].begin());
				  }
				  if(IDPR != 0)
				  {
					  if (std::find(ReceiversTable[0].begin(), ReceiversTable[0].end(), IDPR) != ReceiversTable[0].end())
					  {
						  IDCat2 = ReceiversTable[2].at(std::find(ReceiversTable[0].begin(), ReceiversTable[0].end(), IDPR) - ReceiversTable[0].begin());
					  }
					  if(IDCat2 != 0)
					  {
						  if (std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat2) != OutletsTable[0].end())
						  {
							  IDPlot2 = OutletsTable[1].at(std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat2) - OutletsTable[0].begin());
							  IDPR2 = OutletsTable[2].at(std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat2) - OutletsTable[0].begin());
						  }
						  if(IDPR2 != 0)
						  {
							  if (std::find(ReceiversTable[0].begin(), ReceiversTable[0].end(), IDPR2) != ReceiversTable[0].end())
							  {
								  IDCat3 = ReceiversTable[2].at(std::find(ReceiversTable[0].begin(), ReceiversTable[0].end(), IDPR2) - ReceiversTable[0].begin());
							  }
							  if(IDCat3 != 0)
							  {
								  if (std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat3) != OutletsTable[0].end())
								  {
									  IDPlot3 = OutletsTable[1].at(std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat3) - OutletsTable[0].begin());
									  IDPR3 = OutletsTable[2].at(std::find(OutletsTable[0].begin(), OutletsTable[0].end(), IDCat3) - OutletsTable[0].begin());
								  }
							  }
						  }
					  }
				  }
			  }
			  Feature->SetField("IDPlot", IDPlot);
			  Feature->SetField("IDPR", IDPR);
			  Feature->SetField("IDCat2", IDCat2);
			  Feature->SetField("IDPlot2", IDPlot2);
			  Feature->SetField("IDPR2", IDPR2);
			  Feature->SetField("IDCat3", IDCat3);
			  Feature->SetField("IDPlot3", IDPlot3);
			  Feature->SetField("IDPR3", IDPR3);
			  Layer->SetFeature(Feature);
			  OGRFeature::DestroyFeature(Feature);
		  }
	  }
}


// =====================================================================
// =====================================================================


void LandProcessor::labelCatchments()
{

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

	std::vector<std::pair<std::vector< std::vector<int>>, std::vector<int>>> Catchments;

	std::vector< std::vector<int>> Table(5,std::vector<int >(0)), MainTable(5,std::vector<int >(0));

	{
		std::vector< std::vector<int>> CatchmentsFields(10,std::vector<int >(1));

		std::vector<int> CatchmentsGeometry;

		OGRDataSource *DS = OGRSFDriverRegistrar::Open(getOutputVectorPath(m_OutputCatchmentsGroupedVectorFile).c_str(), TRUE );
		OGRLayer *Layer = DS->GetLayer(0);
		OGRFeature *Feature = nullptr;

		Layer->ResetReading();

		while ((Feature = Layer->GetNextFeature()) != nullptr)
		{
			for(unsigned int i = 0; i < CatchmentsFields.size(); i++)
			{
				CatchmentsFields[i].clear();
			}
			CatchmentsGeometry.clear();

			OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();

			CatchmentsFields[0].push_back(Feature->GetFID());
			CatchmentsFields[1].push_back(Feature->GetFieldAsInteger("IDCat"));
			CatchmentsFields[2].push_back(Feature->GetFieldAsInteger("IDPlot"));
			CatchmentsFields[3].push_back(Feature->GetFieldAsInteger("IDCat2"));
			CatchmentsFields[4].push_back(Feature->GetFieldAsInteger("IDPlot2"));
			CatchmentsFields[5].push_back(Feature->GetFieldAsInteger("IDCat3"));
			CatchmentsFields[6].push_back(Feature->GetFieldAsInteger("IDPlot3"));
			CatchmentsFields[7].push_back(0);
			CatchmentsFields[8].push_back(0);
			CatchmentsFields[9].push_back(0);
			for(unsigned int i = 0; i < Layer->GetFeatureCount(); i++)
			{
				OGRFeature *FeatureToTest = Layer->GetFeature(i);
				OGRGeometry *GeometryToTest = FeatureToTest->GetGeometryRef()->clone();
				if(i != Feature->GetFID() and Feature->GetFieldAsInteger("IDPlot") == FeatureToTest->GetFieldAsInteger("IDPlot") and
						Feature->GetFieldAsInteger("IDPlot2") == FeatureToTest->GetFieldAsInteger("IDPlot2"))
				{
					OGRGeometry *Intersection = GeometryToTest->Intersection(Geometry);
					if(Intersection->getGeometryType() == wkbLineString or Intersection->getGeometryType() == wkbMultiLineString)
					{
						CatchmentsGeometry.push_back(i);
					}
					delete(Intersection);
				}

				OGRFeature::DestroyFeature(FeatureToTest);
				delete(GeometryToTest);
			}
			OGRFeature::DestroyFeature(Feature);
			delete(Geometry);
			Catchments.push_back(std::make_pair(CatchmentsFields,CatchmentsGeometry));
		}

		OGRDataSource::DestroyDataSource(DS);

	}

	{
		// =======================================================================================
		// Setting catchments with ID, ID2 and ID3
		// =======================================================================================

		for (unsigned int i = 0; i < Catchments.size(); i++)
		{
			int IDPlot = Catchments[i].first.at(2)[0];
			int IDPlot2 = Catchments[i].first.at(4)[0];
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
						Table[4].push_back(i);
						Catchments[i].first.at(7)[0] = 1;
						Catchments[i].first.at(8)[0] = 1;
						Catchments[i].first.at(9)[0] = 1;
					}
					else
					{
						Table[0].push_back(0);
						Table[1].push_back(1);
						Table[2].push_back(0);
						Table[3].push_back(1);
						Table[4].push_back(i);
						Catchments[i].first.at(7)[0] = 1;
						Catchments[i].first.at(8)[0] = 1;
						Catchments[i].first.at(9)[0] = 1;
					}
				}
			}
		}

		for(unsigned int i = 0; i < Table.size(); i++)
		{
			MainTable[i] = Table[i];
		}

		if (Table[0].size() != 0)
		{
			while (Table[0].size() != 0)
			{
				int i = 0;

				while (i < Table[0].size())
				{
					if (Table[0][i] != 0)
					{
						int FID = Table[4][i];
						int IDPlot = Catchments[FID].first.at(2)[0];
						int IDN = Catchments[FID].first.at(7)[0];
						int ID2N = Catchments[FID].first.at(8)[0];
						int ID3N = Catchments[FID].first.at(9)[0];

						if(Catchments[FID].second.size() != 0)
						{
							for (unsigned int j = 0; j < Catchments[FID].second.size(); j++)
							{
								int NewFID = Catchments[FID].second.at(j);

								if (Catchments[NewFID].first.at(7)[0] == 0)
								{
									Catchments[NewFID].first.at(7)[0] = IDN;
									Catchments[NewFID].first.at(8)[0] = ID2N;
									Catchments[NewFID].first.at(9)[0] = ID3N;
									Table[0].push_back(IDPlot);
									Table[1].push_back(IDN);
									Table[2].push_back(0);
									Table[3].push_back(1);
									Table[4].push_back(NewFID);
									MainTable[0].push_back(IDPlot);
									MainTable[1].push_back(IDN);
									MainTable[2].push_back(0);
									MainTable[3].push_back(1);
									MainTable[4].push_back(NewFID);
								}
							}
						}
					  }
					  i += 1;
				  }

				  int IDMax = Table[1].back();

				  for(unsigned int i = 0; i < Table.size(); i++)
				  {
					  Table[i].clear();
				  }

				  for (unsigned int i = 0; i < Catchments.size(); i++)
				  {
					  int IDN = Catchments[i].first.at(7)[0];
					  int IDPlot = Catchments[i].first.at(2)[0];
					  int IDPlot2 = Catchments[i].first.at(4)[0];
					  if (IDPlot2 == 0 and IDN == 0)
					  {
						  if (std::find(Table[0].begin(), Table[0].end(), IDPlot) == Table[0].end())
						  {
							  IDN = IDMax+1;
							  Catchments[i].first.at(7)[0] = IDN;
							  Catchments[i].first.at(8)[0] = 1;
							  Catchments[i].first.at(9)[0] = 1;
							  Table[0].push_back(IDPlot);
							  Table[1].push_back(IDMax+1);
							  Table[2].push_back(0);
							  Table[3].push_back(1);
							  Table[4].push_back(i);
							  MainTable[0].push_back(IDPlot);
							  MainTable[1].push_back(IDMax+1);
							  MainTable[2].push_back(0);
							  MainTable[3].push_back(1);
							  MainTable[4].push_back(i);
						  }
					  }
				  }
			  }
		}

		for (unsigned int i = 0; i < Catchments.size(); i++)
		{
			int IDCat2 = Catchments[i].first.at(3)[0];
			int ID2N = Catchments[i].first.at(8)[0];
			if (ID2N == 0)
			{
				for (unsigned int j = 0; j < Catchments.size(); j++)
				{
					if (Catchments[j].first.at(1)[0] == IDCat2 and Catchments[j].first.at(7)[0] != 0)
					{
						Catchments[i].first.at(8)[0] = Catchments[j].first.at(7)[0];
						Catchments[i].first.at(9)[0] = Catchments[j].first.at(8)[0];
					}
				}
			}
		}

		int N = 0;

		for (unsigned int i = 0; i < Catchments.size(); i++)
		{
			if (Catchments[i].first.at(7)[0] == 0)
			{
				N += 1;
			}
		}

		if(N != 0)
		{
			while (N > 0)
			{
				for(unsigned int i = 0; i < Table.size(); i++)
				{
					Table[i].clear();
				}

				for (unsigned int i = 0; i < Catchments.size(); i++)
				{
					int IDPlot = Catchments[i].first.at(2)[0];
					int IDN = Catchments[i].first.at(7)[0];
					int IDPlot2 = Catchments[i].first.at(4)[0];
					int ID2N = Catchments[i].first.at(8)[0];

					if (IDN == 0 and ID2N != 0)
					{
						if (std::find(Table[0].begin(), Table[0].end(), IDPlot) == Table[0].end() and
								std::find(MainTable[0].begin(), MainTable[0].end(), IDPlot) == MainTable[0].end())
						{
							Catchments[i].first.at(7)[0] = 1;
							Table[0].push_back(IDPlot);
							Table[1].push_back(1);
							Table[2].push_back(IDPlot2);
							Table[3].push_back(ID2N);
							Table[4].push_back(i);
							MainTable[0].push_back(IDPlot);
							MainTable[1].push_back(1);
							MainTable[2].push_back(IDPlot2);
							MainTable[3].push_back(ID2N);
							MainTable[4].push_back(i);
						}
						else if (std::find(Table[0].begin(), Table[0].end(), IDPlot) == Table[0].end() and
								std::find(MainTable[0].begin(), MainTable[0].end(), IDPlot) != MainTable[0].end())
						{
							int IDMax = 0;

							for (unsigned int j = 0; j < MainTable[0].size(); j++)
							{
								if (MainTable[0][j] == IDPlot and MainTable[1][j] > IDMax)
								{
									IDMax = MainTable[1][j];
								}
							}

							Catchments[i].first.at(7)[0] = IDMax+1;
							Table[0].push_back(IDPlot);
							Table[1].push_back(IDMax+1);
							Table[2].push_back(IDPlot2);
							Table[3].push_back(ID2N);
							Table[4].push_back(i);
							MainTable[0].push_back(IDPlot);
							MainTable[1].push_back(IDMax+1);
							MainTable[2].push_back(IDPlot2);
							MainTable[3].push_back(ID2N);
							MainTable[4].push_back(i);
						}
					}
				}
				while (Table[0].size() != 0)
				{
					int i = 0;

					while (i < Table[0].size())
					{
						if (Table[0][i] != 0)
						{
							int FID = Table[4][i];
							int IDPlot = Catchments[FID].first.at(2)[0];
							int IDN = Catchments[FID].first.at(7)[0];
							int IDPlot2 = Catchments[FID].first.at(4)[0];
							int ID2N = Catchments[FID].first.at(8)[0];
							int ID3N = Catchments[FID].first.at(9)[0];

							if(Catchments[FID].second.size() != 0)
							{
								for (unsigned int j = 0; j < Catchments[FID].second.size(); j++)
								{
									int NewFID = Catchments[FID].second.at(j);

									if (Catchments[NewFID].first.at(7)[0] == 0 and Catchments[NewFID].first.at(8)[0] == ID2N)
									{
										Catchments[NewFID].first.at(7)[0] = IDN;
										Table[0].push_back(IDPlot);
										Table[1].push_back(IDN);
										Table[2].push_back(IDPlot2);
										Table[3].push_back(ID2N);
										Table[4].push_back(NewFID);
										MainTable[0].push_back(IDPlot);
										MainTable[1].push_back(IDN);
										MainTable[2].push_back(IDPlot2);
										MainTable[3].push_back(ID2N);
										MainTable[4].push_back(NewFID);
									}
								}
							}
						}
						i += 1;
					}

					for(unsigned int i = 0; i < Table.size(); i++)
					{
						Table[i].clear();
					}

					for (unsigned int i = 0; i < Catchments.size(); i++)
					{
						int IDPlot = Catchments[i].first.at(2)[0];
						int IDN = Catchments[i].first.at(7)[0];
						int IDPlot2 = Catchments[i].first.at(4)[0];
						int ID2N = Catchments[i].first.at(8)[0];

						if (IDN == 0 and ID2N != 0)
						{
							if (std::find(Table[0].begin(), Table[0].end(), IDPlot) == Table[0].end())
							{
								int IDMax = 0;

								for (unsigned int j = 0; j < MainTable[0].size(); j++)
								{
									if (MainTable[0][j] == IDPlot and MainTable[1][j] > IDMax)
									{
										IDMax = MainTable[1][j];
									}
								}

								Catchments[i].first.at(7)[0] = IDMax+1;
								Table[0].push_back(IDPlot);
								Table[1].push_back(IDMax+1);
								Table[2].push_back(IDPlot2);
								Table[3].push_back(ID2N);
								Table[4].push_back(i);
								MainTable[0].push_back(IDPlot);
								MainTable[1].push_back(IDMax+1);
								MainTable[2].push_back(IDPlot2);
								MainTable[3].push_back(ID2N);
								MainTable[4].push_back(i);
							}
						}
					}

					for (unsigned int i = 0; i < Catchments.size(); i++)
					{
						int IDCat2 = Catchments[i].first.at(3)[0];
						int ID2N = Catchments[i].first.at(8)[0];

						if (ID2N == 0)
						{
							for (unsigned int j = 0; j < Catchments.size(); j++)
							{
								if (Catchments[j].first.at(1)[0] == IDCat2 and Catchments[j].first.at(7)[0] != 0)
								{
									Catchments[i].first.at(8)[0] = Catchments[j].first.at(7)[0];
									Catchments[i].first.at(9)[0] = Catchments[j].first.at(8)[0];
								}
							}
						}
					}
				}

				N = 0;

				for (unsigned int i = 0; i < Catchments.size(); i++)
				{
					if (Catchments[i].first.at(7)[0] == 0)
					{
						N += 1;
					}
				}
			}
		}

		OGRDataSource *DS = OGRSFDriverRegistrar::Open( getOutputVectorPath(m_OutputCatchmentsGroupedVectorFile).c_str(), TRUE );
		OGRLayer *Layer = DS->GetLayer(0);

		for (unsigned int i = 0; i < Catchments.size(); i++)
		{
			int IDPlot = Catchments[i].first.at(2)[0];
			int IDPlot2 = Catchments[i].first.at(4)[0];
			int IDPlot3 = Catchments[i].first.at(6)[0];
			int IDN = Catchments[i].first.at(7)[0];
			int ID2N = Catchments[i].first.at(8)[0];
			int ID3N = Catchments[i].first.at(9)[0];
			OGRFeature *Feature = Layer->GetFeature(i);
			std::string IDNew = std::to_string(IDPlot)+"N"+std::to_string(IDN);
			std::string ID2New = std::to_string(IDPlot2)+"N"+std::to_string(ID2N);
			std::string ID3New = std::to_string(IDPlot3)+"N"+std::to_string(ID3N);
			Feature->SetField("ID", IDNew.c_str());
			Feature->SetField("ID2", ID2New.c_str());
			Feature->SetField("ID3", ID3New.c_str());
			Layer->SetFeature(Feature);
			OGRFeature::DestroyFeature(Feature);
		}

		OGRDataSource::DestroyDataSource(DS);
	}

}


// =====================================================================
// =====================================================================


void LandProcessor::createEntitiesVector()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputCatchmentsGroupedVectorFile);

		m_SQLRequest = "SELECT ST_Union(geometry), IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3 FROM " + VP.getLayerNameFromFileName(m_OutputCatchmentsGroupedVectorFile) + " GROUP BY IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3";

		VP.executeSQLRequest(m_SQLRequest, m_OutputEntitiesVectorFile);
	}
}


// =====================================================================
// =====================================================================


void LandProcessor::regroupEntitiesVector()
{
	  VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	  {
		  VectorRegrouping VR(getOutputVectorPath(), m_OutputEntitiesVectorFile);

		  VR.regroupSmallFeaturesThatConstituteSinglePlot(m_MinEntSize);

		  VR.regroupSmallFeaturesWithLargerFeatures(m_MinEntSize);
	  }

	  {
		  VectorRegrouping VR(getOutputVectorPath(), m_OutputEntitiesVectorFile);

		  VR.regroupFeaturesThatShareAttributesAndLimits();
	  }

	  {
		  VectorRegrouping VR(getOutputVectorPath(), m_OutputEntitiesVectorFile);

		  VR.regroupLockedFeatures();
	  }
}


// =====================================================================
// =====================================================================


void LandProcessor::createUnionVector()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

	openfluid::utils::GrassGISProxy GRASS(QString::fromStdString(m_GrassTmpPath),
	                                        QString::fromStdString(m_GrassLocation));

	GRASS.setOutputFile(QString::fromStdString(m_TmpPath)+"/createsrfandlnr5.out");
	GRASS.setErrorFile(QString::fromStdString(m_TmpPath)+"/createsrfandlnr5.err");

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
	                                 {"output", "unionplotsentities"},
	                                 {"snap", QString::fromStdString(m_SnapDistance)}},
	                   {"--o"});

	std::vector <std::string> FieldNamesList = {m_IDFieldName,"IDPlot","ID","IDPlot2","ID2","IDPlot3","ID3"};

	GRASS.appendTask("v.db.renamecolumn", {{"map", "unionplotsentities"},
	                                         {"column", QString::fromStdString("a_" + FieldNamesList[0] + "," + FieldNamesList[0])}});

	if(!m_LandUseFieldName.empty())
	{
		GRASS.appendTask("v.db.dropcolumn", {{"map", "unionplotsentities"},
				{"column", QString::fromStdString("a_" + m_LandUseFieldName)}});
	}

	for (unsigned int i = 1; i < FieldNamesList.size(); i++)
	{
		GRASS.appendTask("v.db.renamecolumn", {{"map", "unionplotsentities"},
	                                           {"column", QString::fromStdString("b_" + FieldNamesList[i] + "," + FieldNamesList[i])}});
	}

	GRASS.appendTask("v.db.dropcolumn", {{"map", "unionplotsentities"},
	                                       {"column", QString::fromStdString("a_cat,b_cat")}});

	if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile)))
	{
		mp_VectorDriver->DeleteDataSource(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str());
	}

	GRASS.appendTask("v.out.ogr", {{"input", "unionplotsentities"},
	                                 {"type", "area"},
	                                 {"output", QString::fromStdString(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile))}},
	                   {"-s"});

	if (GRASS.runJob() != 0)
	{
		throw std::runtime_error("LandProcessor::createSRFandLNR5() : unable to run GRASS job (see file " + m_TmpPath + "/createsrfandlnr5.err)");
	}

	VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

	VR.createOuterPolygon();

}


// =====================================================================
// =====================================================================


void LandProcessor::regroupUnionVector()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> FIDSliverList;

	// ======================================================================================
	// Regrouping sliver feature
	// ======================================================================================

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

		FIDSliverList = VR.getSliverFIDList();

		while (FIDSliverList.size() != 0)
		{

			FIDSliverList = VR.getSliverFIDList();

			unsigned int Start = FIDSliverList.size();

			// ======================================================================================
			// Regrouping features that have only one non-sliver polygon as neighbor
			// ======================================================================================

			VR.regroupSliverFeaturesWithSingleNonSliverNeighbor();

			// ======================================================================================
			// Regrouping features based on attributes and geometry
			// ======================================================================================

			VR.regroupSliverFeaturesBasedOnGeometryAndAttributes();

			FIDSliverList = VR.getSliverFIDList();

			unsigned int End = FIDSliverList.size();

			if(Start == End)
			{
				break;
			}
		}

		// ======================================================================================
		// Regrouping sliver features based on geometry (longest shared edge)
		// ======================================================================================

		FIDSliverList = VR.getSliverFIDList();

		while (FIDSliverList.size() != 0)
		{
			FIDSliverList = VR.getSliverFIDList();

			unsigned int Start = FIDSliverList.size();

			VR.regroupSliverFeaturesBasedOnGeometry();

			FIDSliverList = VR.getSliverFIDList();

			unsigned int End = FIDSliverList.size();

			if(Start == End)
			{
				break;
			}
		}

		// ======================================================================================
		// Regrouping sliver features based on attributes (neighbor feature that belongs to the same parcel)
		// ======================================================================================

		FIDSliverList = VR.getSliverFIDList();

		while (FIDSliverList.size() != 0)
		{
			FIDSliverList = VR.getSliverFIDList();

			unsigned int Start = FIDSliverList.size();

			VR.regroupSliverFeaturesBasedOnAttributes();

			FIDSliverList = VR.getSliverFIDList();

			unsigned int End = FIDSliverList.size();

			if(Start == End)
			{
				break;
			}
		}

		// ======================================================================================
		// Regrouping features with parcel ID that is missing from catchments vector
		// ======================================================================================

		VR.regroupFeaturesWithMissingPlotID();
	}

	// ======================================================================================
	// Regrouping of union features
	// ======================================================================================

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

		m_SQLRequest = "SELECT ST_Union(geometry), " + m_IDFieldName + ", IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3 FROM " + VR.getLayerNameFromFileName(m_OutputPlotsAndEntitiesUnionVectorFile) + " GROUP BY " + m_IDFieldName + ", IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3";

		VR.executeSQLRequest(m_SQLRequest, "tmp.shp");
	}

	{
		VectorRegrouping VR(getOutputVectorPath(), "tmp.shp");

		VR.deleteFile(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile).c_str());

		VR.copyDataSourceToDataSource(m_OutputPlotsAndEntitiesUnionVectorFile, getOutputVectorPath());
	}

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

		VR.deleteFile(getOutputVectorPath("tmp.shp"));

		VR.correctOuterPolygon();

		VR.correctMultiparts();
	}

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

		m_SQLRequest = "SELECT ST_Union(geometry), IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3 FROM " + VR.getLayerNameFromFileName(m_OutputPlotsAndEntitiesUnionVectorFile) + " GROUP BY IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3";

		VR.executeSQLRequest(m_SQLRequest, "tmp.shp");
	}

	{
		VectorRegrouping VR(getOutputVectorPath(), "tmp.shp");

		VR.deleteFile(getOutputVectorPath(m_OutputPlotsAndEntitiesUnionVectorFile));

		VR.copyDataSourceToDataSource(m_OutputPlotsAndEntitiesUnionVectorFile, getOutputVectorPath());
	}

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

		VR.deleteFile(getOutputVectorPath("tmp.shp"));
	}

	// ======================================================================================
	// Check for missing connections and set new connections if necessary
	// ======================================================================================

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

		std::vector<unsigned int> FIDWithMissingConnectionList = VR.getFIDWithMissingConnectionList();

		while (FIDWithMissingConnectionList.size() != 0)
		{
			unsigned int Start = VR.getFIDWithMissingConnectionList().size();

			VR.setNewConnectionBasedOnGeometry();

			unsigned int End = VR.getFIDWithMissingConnectionList().size();

			if (Start == End)
			{
				break;
			}
		}
	}

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

		std::vector<unsigned int> FIDWithMissingConnectionList = VR.getFIDWithMissingConnectionList();

		while (FIDWithMissingConnectionList.size() != 0)
		{
			unsigned int Start = VR.getFIDWithMissingConnectionList().size();

			VR.setNewConnectionBasedOnDistance();

			unsigned int End = VR.getFIDWithMissingConnectionList().size();

			if (Start == End)
			{
				break;
			}
		}
	}

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

		std::vector<unsigned int> FIDWithMissingConnectionList = VR.getFIDWithMissingConnectionList();

		while (FIDWithMissingConnectionList.size() != 0)
		{
			unsigned int Start = VR.getFIDWithMissingConnectionList().size();

			VR.dissolveInnerFeatures();

			unsigned int End = VR.getFIDWithMissingConnectionList().size();

			if (Start == End)
			{
				break;
			}
		}
	}

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

	    while(!VR.areThereFeaturesWithCommonReceivingFeature())
	    {
	    	VR.regroupFeaturesBasedOnCommonReceivingFeature();

	    	if(!VR.areThereFeaturesWithCommonReceivingFeature())
	    	{
	    		break;
	    	}
	    }
	}

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

		VR.regroupFeaturesThatShareAttributesAndLimits();
	}

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

		VR.regroupLockedFeatures();
	}
}


// =====================================================================
// =====================================================================


void LandProcessor::createGroupedEntitiesVector()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	// ============================================================================
	// Final regrouping
	// ============================================================================

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputPlotsAndEntitiesUnionVectorFile);

		m_SQLRequest = "SELECT ST_Union(geometry), ID, ID2, ID3 FROM " + VR.getLayerNameFromFileName(m_OutputPlotsAndEntitiesUnionVectorFile) + " GROUP BY ID, ID2, ID3 ORDER BY IDPlot ASC, ID ASC";

		VR.getLayerNameFromFileName(m_OutputEntitiesGroupedVectorFile);

		VR.executeSQLRequest(m_SQLRequest, m_OutputEntitiesGroupedVectorFile);

	}

	{
		VectorRegrouping VR(getOutputVectorPath(), m_OutputEntitiesGroupedVectorFile);

		VR.setConsecutiveID();
	}

}


// =====================================================================
// =====================================================================


void LandProcessor::createLNRVector()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	// ================================================================================
	// Creating linear entities vector LNR.shp files
	// ================================================================================

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputEntitiesGroupedVectorFile);

		VP.createDataSource(wkbLineString, getOutputVectorPath(), m_OutputLinearEntitiesVectorFile);
	}

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputLinearEntitiesVectorFile);

		std::vector<std::pair<OGRFieldType, std::string>> FieldTypesAndNames;

		std::vector<std::string> FieldNames = {"FaceA","FaceB","ID","IDTo"};

		for (unsigned int i = 0; i < FieldNames.size(); i++)
		{
			FieldTypesAndNames.push_back(std::make_pair(OFTString, FieldNames[i]));
		}

		for (unsigned int i = 0; i < FieldTypesAndNames.size(); i++)
		{
			VP.createField(FieldTypesAndNames[i].first, FieldTypesAndNames[i].second);
		}
	}

	{
		LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputLinearEntitiesVectorFile);

		LEP.setLinearEntitiesFaces(m_OutputEntitiesGroupedVectorFile);

		m_SQLRequest = "SELECT ST_Union(geometry), FaceA, FaceB, ID, IDTo FROM " + LEP.getLayerNameFromFileName(m_OutputLinearEntitiesVectorFile) + " GROUP BY FaceA, FaceB, ID, IDTo";

		LEP.executeSQLRequest(m_SQLRequest, "tmp.shp");
	}

	{
		LinearEntitiesProcessing LEP(getOutputVectorPath(), "tmp.shp");

		LEP.deleteFile(getOutputVectorPath(m_OutputLinearEntitiesVectorFile));

		LEP.copyDataSourceToDataSource(m_OutputLinearEntitiesVectorFile, getOutputVectorPath());
	}

	{
		LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputLinearEntitiesVectorFile);

		LEP.deleteFile(getOutputVectorPath("tmp.shp"));
	}

}


// =====================================================================
// =====================================================================


void LandProcessor::createARLVector()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	// ================================================================================
	// Creating surface entities vector file - ARL.shp
	// ================================================================================

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputEntitiesGroupedVectorFile);

		m_SQLRequest = "SELECT ST_Union(geometry), ID, ID2 FROM " + VP.getLayerNameFromFileName(m_OutputEntitiesGroupedVectorFile) + " GROUP BY ID, ID2";

		VP.executeSQLRequest(m_SQLRequest, m_OutputArealEntitiesVectorFile);
	}

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputArealEntitiesVectorFile);

		std::vector<std::pair<OGRFieldType, std::string>> FieldTypesAndNames;

		std::vector<std::string> FieldNames = {"IDTo"};

		for (unsigned int i = 0; i < FieldNames.size(); i++)
		{
			FieldTypesAndNames.push_back(std::make_pair(OFTString, FieldNames[i]));
		}

		for (unsigned int i = 0; i < FieldTypesAndNames.size(); i++)
		{
			VP.createField(FieldTypesAndNames[i].first, FieldTypesAndNames[i].second);
		}
	}
}


// =====================================================================
// =====================================================================


void LandProcessor::setLNRIDs()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputLinearEntitiesVectorFile);

	LEP.setLinearEntitiesIDs(m_OutputArealEntitiesVectorFile);

	std::vector<std::string> AttributesToKeep = {"ID", "IDTo"};

	LEP.deleteUnnecessaryFields(AttributesToKeep);
}


// =====================================================================
// =====================================================================


void LandProcessor::setARLIDs()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	ArealEntitiesProcessing AEP(getOutputVectorPath(), m_OutputArealEntitiesVectorFile);

	AEP.setArealEntitiesIDs();

	std::vector<std::string> AttributesToKeep = {"ID", "IDTo"};

	AEP.deleteUnnecessaryFields(AttributesToKeep);
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
void LandProcessor::setARLAttributes()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	ArealEntitiesProcessing AEP(getOutputVectorPath(), m_OutputArealEntitiesVectorFile);

	AEP.createArealEntitiesAttributes();

	AEP.setArealEntitiesSurface();

	AEP.setArealEntitiesLanduse(m_LandUseFieldName, m_OutputPlotsVectorFile);

	AEP.setArealEntitiesSlope(m_OutputSlopeRasterFile, getOutputRasterPath());
}


// =====================================================================
// =====================================================================


/**
  @internal
  <hr>
  Files
    - m_OutputLinearEntitiesVectorFile
*/
void LandProcessor::setLNRAttributes()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputLinearEntitiesVectorFile);

	LEP.createLinearEntitiesAttributes();

	LEP.setLinearEntitiesLength();
}


// =====================================================================
// =====================================================================


void LandProcessor::releaseARLAndLNRVectors()
{
	 VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	 m_VectorFilesToRelease.push_back(m_OutputArealEntitiesVectorFile);

	 m_VectorFilesToRelease.push_back(m_OutputLinearEntitiesVectorFile);
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
void LandProcessor::createSUVector()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputArealEntitiesVectorFile);

		m_SQLRequest = "SELECT * FROM " + VP.getLayerNameFromFileName(m_OutputArealEntitiesVectorFile) + " WHERE ID != " + "'" + "0N1" + "'";

		VP.executeSQLRequest(m_SQLRequest, m_OutputSUVectorFile);
	}

	// ====================================================================
	// Adding FlowDist attribute that will contain flow distance
	// ====================================================================

	{
		VectorProcessing VP(getOutputVectorPath(), m_OutputSUVectorFile);

		std::vector<std::pair<OGRFieldType, std::string>> FieldTypesAndNames;

		std::vector<std::string> FieldNames = {"FlowDist"};

		for (unsigned int i = 0; i < FieldNames.size(); i++)
		{
			FieldTypesAndNames.push_back(std::make_pair(OFTReal, FieldNames[i]));
		}

		for (unsigned int i = 0; i < FieldTypesAndNames.size(); i++)
		{
			VP.createField(FieldTypesAndNames[i].first, FieldTypesAndNames[i].second);
		}

		for (unsigned int i = 0; i < FieldNames.size(); i++)
		{
			VP.setDefaultFieldValue(FieldNames[i]);
		}
	}
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
void LandProcessor::createRSVector()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	{
		std::vector<std::string> FileNamesList = {m_InputDitchesVectorFile, m_InputThalwegsVectorFile, m_InputRiversVectorFile};

		VectorProcessing VP(getOutputVectorPath(), m_OutputLinearEntitiesVectorFile);

		if (!openfluid::tools::Filesystem::isFile(getInputVectorPath(FileNamesList[0])) and
				!openfluid::tools::Filesystem::isFile(getInputVectorPath(FileNamesList[1])) and
				!openfluid::tools::Filesystem::isFile(getInputVectorPath(FileNamesList[2])))
		{
			std::cout << "LandProcessor::createRSVector(): There are no linear structure files in the input vector directory that could be used to create RS vector" << std::endl;
		}
		else
		{
			VP.copyFile(getOutputVectorPath(), m_OutputRSVectorFile);

			{
				{
					{
						LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputRSVectorFile);

						std::vector<std::pair<OGRFieldType, std::string>> FieldTypesAndNames;

						std::vector<std::string> FieldNames = {"FlowDist", "Ditches", "Thalwegs", "WaterCs", "DitchesL", "ThalwegsL", "WaterCsL", "SurfToLen"};

						for (unsigned int i = 0; i < FieldNames.size(); i++)
						{
							FieldTypesAndNames.push_back(std::make_pair(OFTReal, FieldNames[i]));
						}

						for (unsigned int i = 0; i < FieldTypesAndNames.size(); i++)
						{
							LEP.createField(FieldTypesAndNames[i].first, FieldTypesAndNames[i].second);
						}

						for (unsigned int i = 0; i < FieldNames.size(); i++)
						{
							LEP.setDefaultFieldValue(FieldNames[i]);
						}
					}

					std::vector<unsigned int> FIDList;

					{
						LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputRSVectorFile);

						OGRLayer *Layer = LEP.getLayer();
						OGRFeature *Feature = nullptr;

						Layer->ResetReading();

						while ((Feature = Layer->GetNextFeature()) != nullptr)
						{
							std::string IDTo = "None";
							if (Feature->GetFieldAsString("IDTo") == IDTo)
							{
								FIDList.push_back(Feature->GetFID());
							}
							OGRFeature::DestroyFeature(Feature);
						}
					}

					if(!FIDList.empty())
					{
						LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputRSVectorFile);

						LEP.deleteFeatures(FIDList);

						LEP.repack();
					}

					{
						LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputRSVectorFile);

						OGRLayer *Layer = LEP.getLayer();
						OGRFeature *Feature = nullptr;
						Layer->ResetReading();

						while ((Feature = Layer->GetNextFeature()) != nullptr)
						{
							std::string IDTo = "None";
							Feature->SetField("IDTo", IDTo.c_str());
							Layer->SetFeature(Feature);
							OGRFeature::DestroyFeature(Feature);
						}
						Layer->SyncToDisk();
					}
				}

				std::vector<std::string> FieldNamesList = {"Ditches", "Thalwegs", "WaterCs"};
				std::vector<std::string> FieldNamesLList = {"DitchesL", "ThalwegsL", "WaterCsL"};

				std::vector <int> GeometryTypes = {2,5};

				LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputRSVectorFile);

				for (unsigned int i = 0; i < FileNamesList.size(); i++)
				{
					if(LEP.doesDataExist(FileNamesList[i], getInputVectorPath()) and LEP.canDataBeOpened(FileNamesList[i], getInputVectorPath()) and
							LEP.isDriverSupported(FileNamesList[i], getInputVectorPath()) and
							LEP.isSRSInfoProvided(FileNamesList[i], getInputVectorPath()) and
							LEP.isOfDeclaredGeometryType(GeometryTypes, FileNamesList[i], getInputVectorPath()) and
							LEP.areGeometriesConform(FileNamesList[i], getInputVectorPath()))
					{
						{
							{
								VectorProcessing VP(getInputVectorPath(), FileNamesList[i]);

								VP.copyFile(getOutputVectorPath(), FileNamesList[i]);
							}

							{
								VectorProcessing VP(getOutputVectorPath(), FileNamesList[i]);

								VP.createIDField();
							}
						}

						{
							{
								LinearEntitiesProcessing LinearStructure(getOutputVectorPath(), FileNamesList[i]);

								LinearStructure.creatLinearStructuresIntersectionVector(m_OutputIntersectionVectorFile);
							}

							{
								LinearEntitiesProcessing Intersection(getOutputVectorPath(), m_OutputIntersectionVectorFile);

								Intersection.setLinearStructuresIntersectionVector(FileNamesList[i], m_OutputPlotsVectorFile);

								Intersection.attributeLinearStructures(m_OutputPlotsVectorFile);
							}
						}

						if(!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputIntersectionVectorFile)))
						{
							std::cout << "LandProcessor::createRSVector(): " << FileNamesList[i] << ": reattribution of linear structures produced no results" << std::endl;
						}
						else
						{
							{
								LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputRSVectorFile);

								LEP.setAttributedLinearStructureLength(FieldNamesLList[i], FieldNamesList[i], m_OutputIntersectionVectorFile);

								LEP.deleteFile(getOutputVectorPath(m_OutputIntersectionVectorFile));

								LEP.deleteFile(getOutputVectorPath(FileNamesList[i]));
							}
						}
					}
				}

				{
					std::vector<unsigned int> FIDList;

					{
						LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputRSVectorFile);

						OGRLayer *Layer = LEP.getLayer();
						OGRFeature *Feature = nullptr;
						Layer->ResetReading();

						while ((Feature = Layer->GetNextFeature()) != nullptr)
						{
							unsigned int N = 0;
							for(unsigned int i = 0; i < FieldNamesLList.size(); i++)
							{
								if(Feature->GetFieldAsDouble(FieldNamesLList[i].c_str()) > 0)
								{
									N++;
								}
							}
							if(N == 0)
							{
								FIDList.push_back(Feature->GetFID());
							}
							OGRFeature::DestroyFeature(Feature);
						}
					}

					if(!FIDList.empty())
					{
						{
							LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputRSVectorFile);

							LEP.deleteFeatures(FIDList);

							LEP.repack();
						}
					}
				}
				m_VectorFilesToRelease.push_back(m_OutputRSVectorFile);
			}
		}
	}
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
void LandProcessor::createLIVector()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	{
		std::vector<std::string> FileNamesList = {m_InputHedgesVectorFile, m_InputGrassBandVectorFile, m_InputBenchesVectorFile};

		VectorProcessing VP(getOutputVectorPath(), m_OutputLinearEntitiesVectorFile);

		if (!openfluid::tools::Filesystem::isFile(getInputVectorPath(FileNamesList[0])) and
				!openfluid::tools::Filesystem::isFile(getInputVectorPath(FileNamesList[1])) and
				!openfluid::tools::Filesystem::isFile(getInputVectorPath(FileNamesList[2])))
		{
			std::cout << "LandProcessor::createLIVector(): There are no linear structure files in the input vector directory that could be used to create RS vector" << std::endl;
		}
		else
		{
			VP.copyFile(getOutputVectorPath(), m_OutputLIVectorFile);

			{
				{
					LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputLIVectorFile);

					std::vector<std::pair<OGRFieldType, std::string>> FieldTypesAndNames;

					std::vector<std::string> FieldNames = {"FlowDist", "Hedges", "GrassBs", "Benches", "HedgesL", "GrassBsL", "BenchesL", "SurfToLen"};

					for (unsigned int i = 0; i < FieldNames.size(); i++)
					{
						FieldTypesAndNames.push_back(std::make_pair(OFTReal, FieldNames[i]));
					}

					for (unsigned int i = 0; i < FieldTypesAndNames.size(); i++)
					{
						LEP.createField(FieldTypesAndNames[i].first, FieldTypesAndNames[i].second);
					}

					for (unsigned int i = 0; i < FieldNames.size(); i++)
					{
						LEP.setDefaultFieldValue(FieldNames[i]);
					}

					std::vector<unsigned int> FIDList;

					OGRLayer *Layer = LEP.getLayer(), *RSLayer = nullptr;
					OGRFeature *Feature = nullptr, *RSFeature = nullptr;

					Layer->ResetReading();

					if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputRSVectorFile)))
					{
						{
							VectorProcessing VP(getOutputVectorPath(), m_OutputRSVectorFile);

							RSLayer = VP.getLayer();

							Layer->ResetReading();

							while ((Feature = Layer->GetNextFeature()) != nullptr)
							{
								std::string ID = Feature->GetFieldAsString("ID");
								RSLayer->ResetReading();
								while ((RSFeature = RSLayer->GetNextFeature()) != nullptr)
								{
									std::string IDRS = RSFeature->GetFieldAsString("ID");
									if(ID == IDRS)
									{
										Feature->SetField("IDTo",  IDRS.c_str());
										Layer->SetFeature(Feature);
									}
									OGRFeature::DestroyFeature(RSFeature);
								}
								OGRFeature::DestroyFeature(Feature);
							}
						}
					}

					Layer->SyncToDisk();

				}

				std::vector<std::string> FieldNamesList = {"Hedges", "GrassBs", "Benches"};
				std::vector<std::string> FieldNamesLList = {"HedgesL", "GrassBsL", "BenchesL"};

				std::vector <int> GeometryTypes = {2,5};

				LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputLIVectorFile);

				for (unsigned int i = 0; i < FileNamesList.size(); i++)
				{
					if(LEP.doesDataExist(FileNamesList[i], getInputVectorPath()) and LEP.canDataBeOpened(FileNamesList[i], getInputVectorPath()) and
							LEP.isDriverSupported(FileNamesList[i], getInputVectorPath()) and
							LEP.isSRSInfoProvided(FileNamesList[i], getInputVectorPath()) and
							LEP.isOfDeclaredGeometryType(GeometryTypes, FileNamesList[i], getInputVectorPath()) and
							LEP.areGeometriesConform(FileNamesList[i], getInputVectorPath()))
					{
						{
							{
								VectorProcessing VP(getInputVectorPath(), FileNamesList[i]);

								VP.copyFile(getOutputVectorPath(), FileNamesList[i]);
							}

							{
								VectorProcessing VP(getOutputVectorPath(), FileNamesList[i]);

								VP.createIDField();
							}
						}

						{
							{
								LinearEntitiesProcessing LinearStructure(getOutputVectorPath(), FileNamesList[i]);

								LinearStructure.creatLinearStructuresIntersectionVector(m_OutputIntersectionVectorFile);
							}

							{
								LinearEntitiesProcessing Intersection(getOutputVectorPath(), m_OutputIntersectionVectorFile);

								Intersection.setLinearStructuresIntersectionVector(FileNamesList[i], m_OutputPlotsVectorFile);

								Intersection.attributeLinearStructures(m_OutputPlotsVectorFile);
							}
						}

						if(!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputIntersectionVectorFile)))
						{
							std::cout << "LandProcessor::createLIVector(): " << FileNamesList[i] << ": reattribution of linear structures produced no results" << std::endl;
						}
						else
						{
							{
								LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputLIVectorFile);

								LEP.setAttributedLinearStructureLength(FieldNamesLList[i], FieldNamesList[i], m_OutputIntersectionVectorFile);

								LEP.deleteFile(getOutputVectorPath(m_OutputIntersectionVectorFile));

								LEP.deleteFile(getOutputVectorPath(FileNamesList[i]));
							}
						}
					}
				}
			}
		}
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
void LandProcessor::setSUAttributes()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	// ==================================================================================
	// Set IDTo parameter based on LI# and RS#
	// ==================================================================================

	{
		ArealEntitiesProcessing AEP(getOutputVectorPath(), m_OutputSUVectorFile);

		OGRLayer *Layer = AEP.getLayer();
		OGRFeature *Feature = nullptr;
		Layer->ResetReading();

		while ((Feature = Layer->GetNextFeature()) != nullptr)
		{
			std::string ID = Feature->GetFieldAsString("ID");
			std::string  IDNew = "SU#" + ID;
			Feature->SetField("ID", IDNew.c_str());
			Layer->SetFeature(Feature);
			OGRFeature::DestroyFeature(Feature);
		}

		LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputLIVectorFile);

		OGRLayer *LinearLayer = LEP.getLayer();
		OGRFeature *LinearFeature = nullptr;

		Layer->ResetReading();

		while ((Feature = Layer->GetNextFeature()) != nullptr)
		{
			std::string ID = Feature->GetFieldAsString("ID");
			std::string IDTo = Feature->GetFieldAsString("IDTo");
			if (IDTo.find("-") != std::string::npos)
			{
				if (openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLIVectorFile)))
				{
					LinearLayer->ResetReading();

					while ((LinearFeature = LinearLayer->GetNextFeature()) != nullptr)
					{
						std::string LinearID = LinearFeature->GetFieldAsString("ID");
						{
							if(LinearID == IDTo)
							{
								std::string IDToNew = "LI#" + IDTo;
								Feature->SetField("IDTo", IDToNew.c_str());
								Layer->SetFeature(Feature);
							}
						}
						OGRFeature::DestroyFeature(LinearFeature);
					}
				}
			}
			else
			{
				std::string IDToNew = "SU#" + IDTo;
				Feature->SetField("IDTo", IDToNew.c_str());
				Layer->SetFeature(Feature);
			}
			OGRFeature::DestroyFeature(Feature);
		}

		Layer->ResetReading();

		while ((Feature = Layer->GetNextFeature()) != nullptr)
		{
			std::string ID = Feature->GetFieldAsString("ID");
			std::string IDTo = Feature->GetFieldAsString("IDTo");
			if (IDTo.find("-") != std::string::npos)
			{
				if(IDTo.find("LI#") == std::string::npos)
				{
					std::string IDToNew = "SU#" + IDTo.substr(IDTo.find("-")+1, IDTo.length()-1);
					Feature->SetField("IDTo", IDToNew.c_str());
					Layer->SetFeature(Feature);
				}
			}
			OGRFeature::DestroyFeature(Feature);
		}
		Layer->SyncToDisk();
	}

	// ==================================================================================
	// Setting flow distance parameter for SU#
	// ==================================================================================

	{
		ArealEntitiesProcessing AEP(getOutputVectorPath(), m_OutputSUVectorFile);

		AEP.setArealEntitiesFlowDistance(m_OutputLIVectorFile);
	}
}


// ====================================================================
// ====================================================================


/**
  @internal
  <hr>
  Files
    - m_OutputRSVectorFile
*/
void LandProcessor::setRSAttributes()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputRSVectorFile)))
	{
		std::cout << "LandProcessor::setRSAttributes(): " + m_OutputRSVectorFile + ": no such file in the output vector directory" << std::endl;
	}
	else
	{
	  // ================================================================================
	  //  Set ID as RS# as following: RS#145N1-134N7
	  // ================================================================================

		LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputRSVectorFile);

		OGRLayer *Layer = LEP.getLayer();
		OGRFeature *Feature = nullptr;

		Layer->ResetReading();

		while ((Feature = Layer->GetNextFeature()) != nullptr)
		{
			std::string ID = Feature->GetFieldAsString("ID");
			std::string IDNew = "RS#" + ID;
			Feature->SetField("ID", IDNew.c_str());
			Layer->SetFeature(Feature);
			OGRFeature::DestroyFeature(Feature);
		}
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
void LandProcessor::setLIAttributes()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputSUVectorFile)))
	{
		throw std::runtime_error("LandProcessor::setLIAttributes(): " + m_OutputSUVectorFile + ": no such file in the output vector directory");
	}

	if (!openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputLIVectorFile)))
	{
		throw std::runtime_error("LandProcessor::setLIAttributes(): " + m_OutputLIVectorFile + ": no such file in the output vector directory");
	}

	{
		{
			// ==========================================================================================
			// Calculate ration S(SU#)/L(LI#) (if possible)
			// ==========================================================================================

			LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputLIVectorFile);

			ArealEntitiesProcessing AEP(getOutputVectorPath(), m_OutputSUVectorFile);

			OGRLayer *LinearLayer = LEP.getLayer();
			OGRLayer *ArealLayer = AEP.getLayer();

			OGRFeature *LinearFeature = nullptr, *ArealFeature = nullptr;

			LinearLayer->ResetReading();

			while((LinearFeature = LinearLayer->GetNextFeature()) != nullptr)
			{
				std::string ID = LinearFeature->GetFieldAsString("ID");
				std::string IDFrom = "SU#" + ID.substr(0, ID.find("-"));

				ArealLayer->ResetReading();

				while((ArealFeature = ArealLayer->GetNextFeature()) != nullptr)
				{
					if((std::string) ArealFeature->GetFieldAsString("ID") == IDFrom and (std::string) ArealFeature->GetFieldAsString("IDTo") == "LI#" + ID)
					{
						double Length = LinearFeature->GetFieldAsDouble("Length");
						double Surface = ArealFeature->GetFieldAsDouble("Surface");
						double Value = Surface/Length;
						LinearFeature->SetField("SurfToLen", Value);
						LinearLayer->SetFeature(LinearFeature);
					}
					OGRFeature::DestroyFeature(ArealFeature);
				}
				OGRFeature::DestroyFeature(LinearFeature);
			}
		}

		{
			// ==================================================================================
			// Setting flow distance parameter for LI#
			// ==================================================================================

			LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputLIVectorFile);

			LEP.setLinearEntitiesFlowDistance(m_OutputSUVectorFile);
		}

		{
			// ================================================================================
			// Set ID as following: LI#145N1-134N7
			// ================================================================================

			LinearEntitiesProcessing LEP(getOutputVectorPath(), m_OutputLIVectorFile);

			OGRLayer *LinearLayer = LEP.getLayer();
			OGRFeature *LinearFeature = nullptr;

			LinearLayer->ResetReading();

			while ((LinearFeature = LinearLayer->GetNextFeature()) != nullptr)
			{
				std::string ID = LinearFeature->GetFieldAsString("ID");
				std::string IDNew = "LI#" + ID;
				LinearFeature->SetField("ID", IDNew.c_str());
				LinearLayer->SetFeature(LinearFeature);
				OGRFeature::DestroyFeature(LinearFeature);
			}

			// ================================================================================
			//  Set IDTo as following: in case of SU# for the receiving feature:
			//	SU#15N1, in case of RS# for the receiving feature: RS#15N1-164N3
			//  (and 15N1-164N3 must be the same as ID)
			// ================================================================================

			LinearLayer->ResetReading();

			while ((LinearFeature = LinearLayer->GetNextFeature()) != nullptr)
			{
				std::string IDTo = LinearFeature->GetFieldAsString("IDTo");
				if (IDTo.find("None") == std::string::npos)
				{
					if (IDTo.find("-") == std::string::npos)
					{
						IDTo = "SU#" + IDTo;
						LinearFeature->SetField("IDTo", IDTo.c_str());
						LinearLayer->SetFeature(LinearFeature);
					}
					else
					{
						IDTo = "RS#" + IDTo;
						LinearFeature->SetField("IDTo", IDTo.c_str());
						LinearLayer->SetFeature(LinearFeature);
					}
				}
				OGRFeature::DestroyFeature(LinearFeature);
			}
		}
	}
}


// ====================================================================
// ====================================================================


void LandProcessor::releaseSURSLIVectors()
{
	 VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	 m_VectorFilesToRelease.push_back(m_OutputSUVectorFile);

	 if(openfluid::tools::Filesystem::isFile(getOutputVectorPath(m_OutputRSVectorFile)))
	 {
		 m_VectorFilesToRelease.push_back(m_OutputRSVectorFile);
	 }

	 m_VectorFilesToRelease.push_back(m_OutputLIVectorFile);
}


// ====================================================================
// ====================================================================


void LandProcessor::releaseVectorFile(const std::string& Filename)
{

	VERBOSE_MESSAGE(2,"Releasing vector file : " << Filename);

	mp_VectorDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_VectorDriverName.c_str());

	//=========================================================================
	// Removing destination file if exists
	//=========================================================================

	OGRDataSource* TestDest = OGRSFDriverRegistrar::Open(getReleaseVectorPath(Filename).c_str(),false);

	if (TestDest)
	{
		OGRDataSource::DestroyDataSource(TestDest);
		mp_VectorDriver->DeleteDataSource(getReleaseVectorPath(Filename).c_str());
	}

  //=========================================================================
  // Perform copy
  //=========================================================================

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


void LandProcessor::setLandUseFieldName(const std::string& LandUseFieldName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	VERBOSE_MESSAGE(1,"Land use filed name set to \"" << LandUseFieldName << "\"");

	m_LandUseFieldName = LandUseFieldName;
}

