/*
 * ArealEntitiesProcessing.cpp
 *
 *  Created on: 6 oct. 2017
 *      Author: zadonina
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

#include <LandProcessor/RasterProcessing.hpp>
#include <LandProcessor/VectorProcessing.hpp>
#include <LandProcessor/LinearEntitiesProcessing.hpp>
#include <LandProcessor/ArealEntitiesProcessing.hpp>
#include <LandProcessor/Helpers.hpp>


ArealEntitiesProcessing::ArealEntitiesProcessing(const std::string &FilePath, const std::string &FileName):
				VectorProcessing(FilePath, FileName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	OGRRegisterAll();
	GDALAllRegister();

	std::setlocale(LC_NUMERIC, "C");

}


// =====================================================================
// =====================================================================


ArealEntitiesProcessing::~ArealEntitiesProcessing()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);
}


// =====================================================================
// =====================================================================


void ArealEntitiesProcessing::setArealEntitiesIDs()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	OGRDataSource *DS = OGRSFDriverRegistrar::Open( (m_FullFilePath).c_str(), TRUE );
	OGRLayer *Layer = DS->GetLayer(0);
	OGRFeature *Feature = nullptr;

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		std::string ID = mp_Feature->GetFieldAsString("ID");
		std::string ID2 = mp_Feature->GetFieldAsString("ID2");

		if(ID == "0N1" and ID2 == "0N1")
		{
			mp_Feature->SetField("IDTo", ID2.c_str());
			mp_Layer->SetFeature(mp_Feature);
		}
		else
		{
			mp_Geometry = mp_Feature->GetGeometryRef()->clone();

			Layer->ResetReading();

			while ((Feature = Layer->GetNextFeature()) != nullptr)
			{
				if((std::string)Feature->GetFieldAsString("ID") == ID2)
				{
					OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();
					OGRGeometry *IntersectionGeometry = mp_Geometry->Intersection(Geometry);
					if(IntersectionGeometry->getGeometryType() == 1)
					{
						mp_Feature->SetField("IDTo", ID2.c_str());
						mp_Layer->SetFeature(mp_Feature);
					}
					else
					{
						std::string IDTo = ID + "-" + ID2;
						mp_Feature->SetField("IDTo", IDTo.c_str());
						mp_Layer->SetFeature(mp_Feature);
					}
					delete(IntersectionGeometry);
					delete(Geometry);
				}
				OGRFeature::DestroyFeature(Feature);
			}
			delete(mp_Geometry);
		}
		OGRFeature::DestroyFeature(mp_Feature);
	}

	OGRDataSource::DestroyDataSource(DS);
}


// =====================================================================
// =====================================================================


void ArealEntitiesProcessing::createArealEntitiesAttributes()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<std::pair<OGRFieldType, std::string>> AttributeTypesAndNames;

	std::vector<std::string> AttributeNames = {"LandUse", "Surface",  "SlopeMin", "SlopeMax", "SlopeMean"};

	AttributeTypesAndNames.push_back(std::make_pair(OFTString, AttributeNames[0]));

	for (unsigned int i = 1; i < AttributeNames.size(); i++)
	{
		AttributeTypesAndNames.push_back(std::make_pair(OFTReal, AttributeNames[i]));
	}

	for (unsigned int i = 0; i < AttributeTypesAndNames.size(); i++)
	{
		createField(AttributeTypesAndNames[i].first, AttributeTypesAndNames[i].second);
	}
}


// =====================================================================
// =====================================================================


void ArealEntitiesProcessing::setArealEntitiesSurface()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		mp_Geometry = mp_Feature->GetGeometryRef()->clone();
		if (mp_Geometry->getGeometryType() == 3)
		{
			mp_Feature->SetField("Surface", ((OGRPolygon *) mp_Feature->GetGeometryRef())->get_Area());
			mp_Layer->SetFeature(mp_Feature);
			OGRFeature::DestroyFeature(mp_Feature);
			delete(mp_Geometry);
		}
		else if (mp_Geometry->getGeometryType() == 6)
		{
			double Surface = 0;
			OGRMultiPolygon* SRFMultiPolygon = (OGRMultiPolygon*) mp_Geometry;
			for (unsigned int i = 0; i < SRFMultiPolygon->getNumGeometries(); i++)
			{
				Surface += ((OGRPolygon*) SRFMultiPolygon->getGeometryRef(i))->get_Area();
			}
			mp_Feature->SetField("Surface", Surface);
			mp_Layer->SetFeature(mp_Feature);
			OGRFeature::DestroyFeature(mp_Feature);
			delete(mp_Geometry);
		}
		else
		{
			OGRFeature::DestroyFeature(mp_Feature);
			delete(mp_Geometry);
		}
	}
}


// =====================================================================
// =====================================================================


void ArealEntitiesProcessing::setArealEntitiesLanduse(const std::string &LandUseFieldName, const std::string &ParcelVectorFileName, std::string OutputVectorFilePath)
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(OutputVectorFilePath.empty())
	{
		OutputVectorFilePath = m_FilePath;
	}

	std::vector<std::pair<unsigned int, std::string>> ParcelIDAndLandUse;

	{
		VectorProcessing ParcelVP(OutputVectorFilePath, ParcelVectorFileName);

		if(ParcelVP.doesFieldExist(LandUseFieldName))
		{
			OGRLayer *ParcelLayer = ParcelVP.getLayer();
			OGRFeature *ParcelFeature = nullptr;

			while ((ParcelFeature = ParcelLayer->GetNextFeature()) != nullptr)
			{
				ParcelIDAndLandUse.push_back(std::make_pair(ParcelFeature->GetFieldAsInteger(m_IDFieldName.c_str()), ParcelFeature->GetFieldAsString(LandUseFieldName.c_str())));
				OGRFeature::DestroyFeature(ParcelFeature);
			}
		}
	}

	if(!ParcelIDAndLandUse.empty())
	{
		mp_Layer->ResetReading();

		while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			std::string ID = mp_Feature->GetFieldAsString("ID");
			unsigned int IDPlot = std::stoi(ID.substr(0, ID.find("N")));

			if (IDPlot != 0)
			{
				for(unsigned int i = 0; i < ParcelIDAndLandUse.size(); i++)
				{
					if(ParcelIDAndLandUse[i].first == IDPlot)
					{
						mp_Feature->SetField("LandUse", (ParcelIDAndLandUse[i].second).c_str());
						mp_Layer->SetFeature(mp_Feature);
					}
				}
			}
			else
			{
				mp_Feature->SetField("LandUse", "None");
				mp_Layer->SetFeature(mp_Feature);
			}
			OGRFeature::DestroyFeature(mp_Feature);
		}
	}
	else
	{
		mp_Layer->ResetReading();
		while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			mp_Feature->SetField("LandUse", "None");
			mp_Layer->SetFeature(mp_Feature);
			OGRFeature::DestroyFeature(mp_Feature);
		}
	}
}


// =====================================================================
// =====================================================================


void ArealEntitiesProcessing::setArealEntitiesSlope(const std::string &SlopeRasterFileName, const std::string &SlopeRasterFilePath)
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_Layer->ResetReading();

	RasterProcessing SlopeRP(SlopeRasterFilePath, SlopeRasterFileName, GA_ReadOnly);

	double XOrigin = SlopeRP.m_GeoTransform[0], YOrigin = SlopeRP.m_GeoTransform[3], PixelWidth = SlopeRP.m_GeoTransform[1], PixelHeight = SlopeRP.m_GeoTransform[5];

	OGRPoint Point;

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		std::string ID = mp_Feature->GetFieldAsString("ID");
		unsigned int IDPlot = std::stoi(ID.substr(0, ID.find("N")));
		if (IDPlot != 0)
		{
			mp_Geometry = mp_Feature->GetGeometryRef()->clone();
			std::vector <double> XPointsList;
			std::vector <double> YPointsList;
			if (mp_Geometry->getGeometryType() == 3)
			{
				OGRPolygon *Polygon = (OGRPolygon *) mp_Geometry;
				unsigned int NumberOfPoints = Polygon->getExteriorRing()->getNumPoints();
				for (unsigned int i = 0; i < NumberOfPoints; i++)
				{
					Polygon->getExteriorRing()->getPoint(i, &Point);
					XPointsList.push_back(Point.getX());
					YPointsList.push_back(Point.getY());
				}
			}
			else if (mp_Geometry->getGeometryType() == 6)
			{
				OGRMultiPolygon *MultiPolygon = (OGRMultiPolygon *) mp_Geometry;
				for (unsigned int i = 0; i < MultiPolygon->getNumGeometries(); i++)
				{
					OGRPolygon *Polygon = (OGRPolygon*) MultiPolygon->getGeometryRef(i);
					unsigned int NumberOfPoints = Polygon->getExteriorRing()->getNumPoints();
					for (unsigned int j = 0; j < NumberOfPoints; j++)
					{
						Polygon->getExteriorRing()->getPoint(j, &Point);
						XPointsList.push_back(Point.getX());
						YPointsList.push_back(Point.getY());
					}
				}
			}
			double XMin = XPointsList.at(std::distance(std::begin(XPointsList), std::min_element(std::begin(XPointsList), std::end(XPointsList))));
			double XMax = XPointsList.at(std::distance(std::begin(XPointsList), std::max_element(std::begin(XPointsList), std::end(XPointsList))));
			double YMin = YPointsList.at(std::distance(std::begin(YPointsList), std::min_element(std::begin(YPointsList), std::end(YPointsList))));
			double YMax = YPointsList.at(std::distance(std::begin(YPointsList), std::max_element(std::begin(YPointsList), std::end(YPointsList))));
			XPointsList.clear();
			YPointsList.clear();
			int XOff = int((XMin - XOrigin)/PixelWidth);
			int YOff = int((YOrigin - YMax)/PixelWidth);
			int XCount = int((XMax - XMin)/PixelWidth)+1;
			int YCount = int((YMax - YMin)/PixelWidth)+1;
			GDALDriver *Driver = GetGDALDriverManager()->GetDriverByName("MEM");
			GDALDataset *Target = Driver->Create("", XCount+1, YCount+1, 1, GDT_Byte, nullptr);
			double NewGeoTransformVal[6] = { XOrigin+(XOff*PixelWidth), PixelWidth, 0, YOrigin-(YOff*PixelWidth), 0, PixelHeight };
			Target->SetGeoTransform(NewGeoTransformVal);
			Target->SetProjection(SlopeRP.m_ProjectionRef);
			int BandList[1] = {1};
			double GeometryBurnValue[1] = {1};
			char **NewOptions = nullptr;
			NewOptions = CSLSetNameValue(NewOptions, "ALL_TOUCHED", "TRUE");
			GDALProgressFunc ProgressFunc = nullptr;
			GDALTransformerFunc   TransformerFunc = nullptr;
			CPLErr ERR;
			OGRGeometryH mp_GeometryH = (OGRGeometryH) mp_Geometry;
			std::vector <OGRGeometryH> GeometryList;
			GeometryList.push_back(mp_GeometryH);
			ERR = GDALRasterizeGeometries(Target,
					1,
					BandList,
					1,
					(OGRGeometryH*)&GeometryList[0],
					TransformerFunc,
					nullptr,
					GeometryBurnValue,
					NewOptions,
					ProgressFunc,
					nullptr);
			GDALRasterBand *TargetBand = Target->GetRasterBand(1);
			float *SlopeRow = (float*) CPLMalloc(sizeof(float) *Target->GetRasterXSize());
			int *TargetRow = (int*) CPLMalloc(sizeof(int) *Target->GetRasterXSize());
			std::vector <double> SlopeValuesTable;
			GeometryList.clear();
			for (unsigned int i = 0; i < Target->GetRasterYSize(); i++)
			{
				for (unsigned int j = 0; j < Target->GetRasterXSize(); j++)
				{
					SlopeRP.mp_Band->RasterIO( GF_Read, XOff, YOff+i, Target->GetRasterXSize(), 1, SlopeRow, Target->GetRasterXSize(), 1, GDT_Float32, 0, 0 );
					TargetBand->RasterIO( GF_Read, 0, i, Target->GetRasterXSize(), 1, TargetRow, Target->GetRasterXSize(), 1, GDT_Int32, 0, 0 );
					if (TargetRow[j] != 0)
					{
						SlopeValuesTable.push_back(SlopeRow[j]);
					}
				}
			}
			double SlopeMin = SlopeValuesTable.at(std::distance(std::begin(SlopeValuesTable), std::min_element(std::begin(SlopeValuesTable), std::end(SlopeValuesTable))));
			double SlopeMax = SlopeValuesTable.at(std::distance(std::begin(SlopeValuesTable), std::max_element(std::begin(SlopeValuesTable), std::end(SlopeValuesTable))));
			double SlopeMean = 0;
			for (unsigned int k = 0; k < SlopeValuesTable.size(); k++)
			{
				SlopeMean += SlopeValuesTable[k];
			}
			SlopeMean = SlopeMean/SlopeValuesTable.size();
			SlopeValuesTable.clear();
			mp_Feature->SetField("SlopeMin", SlopeMin);
			mp_Feature->SetField("SlopeMax", SlopeMax);
			mp_Feature->SetField("SlopeMean", SlopeMean);
			mp_Layer->SetFeature(mp_Feature);
			CSLDestroy(NewOptions);
			CPLFree(SlopeRow);
			CPLFree(TargetRow);
			GDALClose( Target );
			delete(mp_Geometry);
		}
		else
		{
			mp_Feature->SetField("SlopeMin", 0);
			mp_Feature->SetField("SlopeMax", 0);
			mp_Feature->SetField("SlopeMean", 0);
			mp_Layer->SetFeature(mp_Feature);
		}
		OGRFeature::DestroyFeature(mp_Feature);
	}
}


// =====================================================================
// =====================================================================


void ArealEntitiesProcessing::setArealEntitiesFlowDistance(const std::string &LinearEntityVectorFileName, std::string OutputVectorFilePath)
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(OutputVectorFilePath.empty())
	{
		OutputVectorFilePath = m_FilePath;
	}

	OGRPoint *AECentroid, *AEPointOnSurf, *Centroid, *PointOnSurf;

	mp_Layer->ResetReading();

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		std::string AEIDTo = mp_Feature->GetFieldAsString("IDTo");

		mp_Geometry = mp_Feature->GetGeometryRef()->clone();

		AECentroid = new OGRPoint;

		mp_Geometry->Centroid(AECentroid);

		geos::geom::Geometry *AEGeometryGEOS = (geos::geom::Geometry*) openfluid::landr::convertOGRGeometryToGEOS(mp_Geometry);

		geos::geom::Point * AEGeometryPointOnSurfGEOS = AEGeometryGEOS->getInteriorPoint();

		OGRGeometry *AEGeometryPointOnSurf = openfluid::landr::convertGEOSGeometryToOGR((GEOSGeom) ((geos::geom::Geometry *) AEGeometryPointOnSurfGEOS));
		AEPointOnSurf = (OGRPoint*) AEGeometryPointOnSurf->clone();

		delete AEGeometryGEOS;
		delete AEGeometryPointOnSurfGEOS;
		delete(AEGeometryPointOnSurf);

		if (AEIDTo.find("SU") != std::string::npos)
		{
			if (AEIDTo != "SU#0N1")
			{
				OGRLayer *Layer = getLayer();
				OGRFeature *Feature = nullptr;

				Layer->ResetReading();

				while((Feature = Layer->GetNextFeature()) != nullptr)
				{
					std::string ID = Feature->GetFieldAsString("ID");

					if(ID == AEIDTo)
					{
						OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();

						Centroid = new OGRPoint;

						Geometry->Centroid(Centroid);

						geos::geom::Geometry *GeometryGEOS = (geos::geom::Geometry*) openfluid::landr::convertOGRGeometryToGEOS(Geometry);

						geos::geom::Point * GeometryPointOnSurfGEOS = GeometryGEOS->getInteriorPoint();

						OGRGeometry *GeometryPointOnSurf = openfluid::landr::convertGEOSGeometryToOGR((GEOSGeom) ((geos::geom::Geometry *) GeometryPointOnSurfGEOS));
						PointOnSurf = (OGRPoint*) GeometryPointOnSurf->clone();

						delete GeometryGEOS;
						delete GeometryPointOnSurfGEOS;
						delete(GeometryPointOnSurf);

						OGRGeometry *IntersectionGeometry = Geometry->Intersection(mp_Geometry);

						if(AECentroid->Within(mp_Geometry))
						{
							if(Centroid->Within(Geometry))
							{
								double Distance = AECentroid->Distance(IntersectionGeometry) + IntersectionGeometry->Distance(Centroid);
								mp_Feature->SetField("FlowDist", Distance);
								mp_Layer->SetFeature(mp_Feature);
							}
							else
							{
								double Distance = AECentroid->Distance(IntersectionGeometry) + IntersectionGeometry->Distance(PointOnSurf);
								mp_Feature->SetField("FlowDist", Distance);
								mp_Layer->SetFeature(mp_Feature);
							}
						}
						else
						{
							if(Centroid->Within(Geometry))
							{
								double Distance = AEPointOnSurf->Distance(IntersectionGeometry) + IntersectionGeometry->Distance(Centroid);
								mp_Feature->SetField("FlowDist", Distance);
								mp_Layer->SetFeature(mp_Feature);
							}
							else
							{
								double Distance = AEPointOnSurf->Distance(IntersectionGeometry) + IntersectionGeometry->Distance(PointOnSurf);
								mp_Feature->SetField("FlowDist", Distance);
								mp_Layer->SetFeature(mp_Feature);
							}
						}

						delete Centroid;
						delete PointOnSurf;
						delete(Geometry);
						delete(IntersectionGeometry);
					}
					OGRFeature::DestroyFeature(Feature);
				}
			}
			else
			{
				mp_Feature->SetField("FlowDist", 0);
				mp_Layer->SetFeature(mp_Feature);
			}
		}
		else if (AEIDTo.find("LI") != std::string::npos)
		{
			AEIDTo = AEIDTo.substr(AEIDTo.find("#")+1, AEIDTo.length()-1);

			LinearEntitiesProcessing LEP(OutputVectorFilePath, LinearEntityVectorFileName);

			LEP.mp_Layer->ResetReading();

			while((LEP.mp_Feature = LEP.mp_Layer->GetNextFeature()) != nullptr)
			{
				std::string ID = LEP.mp_Feature->GetFieldAsString("ID");

				if(ID == AEIDTo)
				{
					LEP.mp_Geometry = LEP.mp_Feature->GetGeometryRef()->clone();

					if(AECentroid->Within(mp_Geometry))
					{
						double Distance = AECentroid->Distance(LEP.mp_Geometry);
						mp_Feature->SetField("FlowDist", Distance);
						mp_Layer->SetFeature(mp_Feature);
					}
					else
					{
						double Distance = AEPointOnSurf->Distance(LEP.mp_Geometry);
						mp_Feature->SetField("FlowDist", Distance);
						mp_Layer->SetFeature(mp_Feature);
					}
					delete(LEP.mp_Geometry);
				}
				OGRFeature::DestroyFeature(LEP.mp_Feature);
			}
		}
		OGRFeature::DestroyFeature(mp_Feature);
		delete AECentroid;
		delete AEPointOnSurf;
		delete(mp_Geometry);
	}
}
