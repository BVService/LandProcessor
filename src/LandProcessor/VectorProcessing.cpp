/*
 * VectorProcessing.cpp
 *
 *  Created on: 26 juil. 2017
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
#include <gdal/ogr_srs_api.h>
#include <gdal/ogr_geometry.h>
#include <gdal/cpl_string.h>

#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/LineString.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/operation/valid/IsValidOp.h>
#include <geos/geom/Point.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/Coordinate.h>
#include <geos/operation/overlay/snap/GeometrySnapper.h>

#include <openfluid/tools/FileHelpers.hpp>
#include <openfluid/tools/Filesystem.hpp>
#include <openfluid/utils/GrassGISProxy.hpp>
#include <openfluid/base/Environment.hpp>
#include <openfluid/tools/DataHelpers.hpp>
#include <openfluid/landr/GEOSHelpers.hpp>

#include <LandProcessor/VectorProcessing.hpp>
#include <LandProcessor/RasterProcessing.hpp>
#include <LandProcessor/Helpers.hpp>

#include <iostream>

VectorProcessing::VectorProcessing(const std::string &FilePath, const std::string &FileName):
					m_FilePath(FilePath), m_FileName(FileName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	OGRRegisterAll();
	GDALAllRegister();

	std::setlocale(LC_NUMERIC, "C");

	if (!openfluid::tools::Filesystem::isDirectory(m_FilePath))
	{
		throw std::runtime_error("VectorProcessing::VectorProcessing(): input vector directory " + m_FilePath + " does not exists");
	}

	if (!openfluid::tools::Filesystem::isFile(m_FilePath + "/" + m_FileName))
	{
		throw std::runtime_error("VectorProcessing::VectorProcessing(): input vector file " + m_FileName + " does not exists in the "
				"input directory" + m_FilePath);
	}

	m_FullFilePath = m_FilePath + "/" + m_FileName;

	mp_Driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(m_DriverName.c_str());

	mp_DataSource = OGRSFDriverRegistrar::Open((m_FullFilePath).c_str(), TRUE);

	if(mp_DataSource == nullptr)
	{
		throw std::runtime_error("VectorProcessing::VectorProcessing(): input vector file " + m_FileName + " could not be opened");
	}

	if (mp_DataSource->GetDriver() != mp_Driver)
	{
		OGRDataSource::DestroyDataSource(mp_DataSource);
		throw std::runtime_error("VectorProcessing::VectorProcessing(): vector driver " + ((std::string) mp_DataSource->GetDriver()->GetName()) + " is not supported");
	}

	mp_Layer = mp_DataSource->GetLayer(0);

	if(mp_Layer == nullptr)
	{
		throw std::runtime_error("VectorProcessing::VectorProcessing(): input vector file " + m_FileName + " could not be opened");
	}

	m_LayerName = mp_Layer->GetName();

	mp_SRS = mp_Layer->GetSpatialRef();

	if (mp_SRS == nullptr)
	{
		OGRDataSource::DestroyDataSource(mp_DataSource);
		throw std::runtime_error("VectorProcessing::VectorProcessing(): " + m_FullFilePath + ": no spatial reference information is provided for the vector file");
	}
}


// =====================================================================
// =====================================================================


VectorProcessing::~VectorProcessing()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(mp_DataSource != nullptr)
	{
		OGRDataSource::DestroyDataSource(mp_DataSource);
	}
}


// =====================================================================
// =====================================================================


std::string VectorProcessing::getLayerName(OGRLayer *Layer)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::string LayerName;

	if(!Layer)
	{
		LayerName = mp_Layer->GetName();
	}
	else
	{
		LayerName = Layer->GetName();
	}

	return LayerName;
}


// =====================================================================
// =====================================================================


OGRLayer* VectorProcessing::getLayer()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	return mp_Layer;
}


// =====================================================================
// =====================================================================


bool VectorProcessing::doesDataExist(const std::string &FileName, std::string FilePath)
{
	if(FilePath.empty())
	{
		FilePath = m_FilePath;
	}

	if (!openfluid::tools::Filesystem::isFile(FilePath + "/" + FileName))
	{
		std::cout << "VectorProcessing::doesDataExist(): input vector directory " << FileName << " does not exists in the directory " << FilePath << std::endl;
		return false;
	}
	else
	{
		return true;
	}
}


// =====================================================================
// =====================================================================


bool VectorProcessing::canDataBeOpened(const std::string &FileName, std::string FilePath)
{
	if(FilePath.empty())
	{
		FilePath = m_FilePath;
	}

	OGRDataSource *DS = OGRSFDriverRegistrar::Open((FilePath + "/" + FileName).c_str(), TRUE);

	if(!DS)
	{
		std::cout << "VectorProcessing::canDataBeOpened(): input vector file " << FileName << " could not be opened" << std::endl;
		OGRDataSource::DestroyDataSource(DS);
		return false;
	}
	else
	{
		OGRDataSource::DestroyDataSource(DS);
		return true;
	}
}


// =====================================================================
// =====================================================================


bool VectorProcessing::isDriverSupported(const std::string &FileName, std::string FilePath)
{
	if(FilePath.empty())
	{
		FilePath = m_FilePath;
	}

	OGRDataSource *DS = OGRSFDriverRegistrar::Open((FilePath + "/" + FileName).c_str(), TRUE);

	if (DS->GetDriver() != mp_Driver)
	{
		OGRDataSource::DestroyDataSource(DS);
		std::cout << "VectorProcessing::isDriverSupported(): vector driver " << ((std::string) DS->GetDriver()->GetName()) << " is not supported" << std::endl;
		return false;
	}
	else
	{
		OGRDataSource::DestroyDataSource(DS);
		return true;
	}
}


// =====================================================================
// =====================================================================


bool VectorProcessing::isSRSInfoProvided(const std::string &FileName, std::string FilePath)
{
	if(FilePath.empty())
	{
		FilePath = m_FilePath;
	}

	OGRDataSource *DS = OGRSFDriverRegistrar::Open((FilePath + "/" + FileName).c_str(), TRUE);

	OGRLayer *Layer = DS->GetLayer(0);

	OGRSpatialReference *SRS = Layer->GetSpatialRef();

	if ((!SRS) or (!SRS->GetAttrValue("AUTHORITY",1)) or ((std::string) (SRS->GetAttrValue("AUTHORITY",1)) != m_EPSGCode.substr(m_EPSGCode.find(":")+1, m_EPSGCode.length()-1)))
	{
		OGRDataSource::DestroyDataSource(DS);
		std::cout << "VectorProcessing::isSRSInfoProvided(): " << FilePath << "/" << FileName << ": "
				"no spatial reference information is provided or it's parameters do not correspond to the default parameters" << std::endl;
		return false;
	}
	else
	{
		OGRDataSource::DestroyDataSource(DS);
		return true;
	}
}


// =====================================================================
// =====================================================================


bool VectorProcessing::isOfDeclaredGeometryType(std::vector <int> GeometryTypes, const std::string &FileName, std::string FilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FilePath.empty())
	{
		FilePath = m_FilePath;
	}

	OGRDataSource *DS = OGRSFDriverRegistrar::Open((FilePath + "/" + FileName).c_str(), TRUE);

	OGRLayer *Layer = DS->GetLayer(0);

	if(std::find(GeometryTypes.begin(), GeometryTypes.end(), (int) Layer->GetGeomType()) == GeometryTypes.end())
	{
		OGRDataSource::DestroyDataSource(mp_DataSource);
		std::cout << "VectorProcessing::isOfDeclaredGeometryType(): layer geometry type of the vector file " << FilePath  << "/" << FileName << " does not correspond to the"
				" declared geometry type(s)" << std::endl;
		return false;
	}
	else
	{
		OGRDataSource::DestroyDataSource(DS);
		return true;
	}
}


// =====================================================================
// =====================================================================


bool VectorProcessing::areGeometriesConform(const std::string &FileName, std::string FilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FilePath.empty())
	{
		FilePath = m_FilePath;
	}

	unsigned int N = 0;

	OGRDataSource *DS = OGRSFDriverRegistrar::Open((FilePath + "/" + FileName).c_str(), TRUE);

	OGRLayer *Layer = DS->GetLayer(0);

	OGRFeature *Feature = nullptr;

	Layer->ResetReading();

	while ((Feature = Layer->GetNextFeature()) != nullptr)
	{
		OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();

		if (!Geometry->IsValid())
		{
			std::cout << "VectorProcessing::areGeometriesConform(): " << FilePath  << "/" << FileName << ": feature " << Feature->GetFID() << " has invalid geometry" << std::endl;
			N++;
		}
		if (!Geometry->IsSimple())
		{
			std::cout << "VectorProcessing::areGeometriesConform(): " << FilePath  << "/" << FileName << ": feature " << Feature->GetFID() << " has geometry that is not simple" << std::endl;
			N++;
		}
		if (Geometry->IsEmpty())
		{
			std::cout << "VectorProcessing::areGeometriesConform(): " << FilePath  << "/" << FileName << ": feature " << Feature->GetFID() << " has empty geometry" << std::endl;
			N++;
		}
		OGRFeature::DestroyFeature(Feature);
		delete(Geometry);
	}

	if(N != 0)
	{
		OGRDataSource::DestroyDataSource(DS);
		std::cout << "VectorProcessing::areGeometriesConform(): one or several geometries of the vector file " << FilePath  << "/" << FileName << " "
				"could not be used for calculations" << std::endl;
		return false;
	}
	else
	{
		OGRDataSource::DestroyDataSource(DS);
		return true;
	}
}


// =====================================================================
// =====================================================================


void VectorProcessing::checkEPSG()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if (!mp_SRS->GetAttrValue("AUTHORITY",1))
	{
		throw std::runtime_error("VectorProcessing::checkEPSG(): " + m_FullFilePath + ": no EPSG code provided for the vector data");
	}

	if ((std::string) (mp_SRS->GetAttrValue("AUTHORITY",1)) != m_EPSGCode.substr(m_EPSGCode.find(":")+1, m_EPSGCode.length()-1))
	{
		throw std::runtime_error("VectorProcessing::checkEPSG(): " + m_FullFilePath + ": vector data EPSG code does not correspond to the default EPSG code");
	}
}


// =====================================================================
// =====================================================================


void VectorProcessing::checkGeometryType(std::vector <int> GeometryTypes)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(std::find(GeometryTypes.begin(), GeometryTypes.end(), (int) mp_Layer->GetGeomType()) == GeometryTypes.end())
	{
		throw std::runtime_error("VectorProcessing::checkGeometryType(): layer geometry type of the vector file " + m_FullFilePath + " does not correspond to the"
				" declared geometry type(s)");
	}
}


// =====================================================================
// =====================================================================


void VectorProcessing::checkGeometries()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	unsigned int N = 0;

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		mp_Geometry = mp_Feature->GetGeometryRef()->clone();

		if (!mp_Geometry->IsValid())
		{
			std::cout << "VectorProcessing::checkGeometries(): " << m_FullFilePath.c_str() << ": feature " << mp_Feature->GetFID() << " has invalid geometry" << std::endl;
			N++;
		}
		if (!mp_Geometry->IsSimple())
		{
			std::cout << "VectorProcessing::checkGeometries(): " << m_FullFilePath.c_str() << ": feature " << mp_Feature->GetFID() << " has geometry that is not simple" << std::endl;
			N++;
		}
		if (mp_Geometry->IsEmpty())
		{
			std::cout << "VectorProcessing::checkGeometries(): " << m_FullFilePath.c_str() << ": feature " << mp_Feature->GetFID() << " has empty geometry" << std::endl;
			N++;
		}
		OGRFeature::DestroyFeature(mp_Feature);
		delete(mp_Geometry);
	}

	if(N != 0)
	{
		throw std::runtime_error("VectorProcessing::checkGeometries(): one or several geometries of the vector file " + m_FullFilePath + " could not be used for calculations");
	}
}


// =====================================================================
// =====================================================================


std::string VectorProcessing::getLayerNameFromFileName(std::string FileName)
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::string LayerName;

	if(FileName.empty())
	{
		FileName = m_FileName;
	}

	LayerName = FileName.substr(0,FileName.find("."));

	return LayerName;

}


// =====================================================================
// =====================================================================


std::string VectorProcessing::getLayerNameFromFilePath(std::string FullFilePath)
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::string LayerName;

	if(FullFilePath.empty())
	{
		FullFilePath = m_FullFilePath;
	}

	LayerName = FullFilePath.substr(FullFilePath.rfind("/")+1, FullFilePath.rfind(".")-(FullFilePath.rfind("/")+1));

	return LayerName;
}


// =====================================================================
// =====================================================================


void VectorProcessing::copyFile(std::string NewFilePath, std::string NewFileName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if (NewFileName.empty())
	{
		NewFileName = m_FileName;
	}

	std::string FullNewPathName = NewFilePath + "/" + NewFileName;

	OGRDataSource *NewDS = nullptr;

	deleteFile(FullNewPathName);

	NewDS = mp_Driver->CreateDataSource((FullNewPathName).c_str());

	NewDS->CopyLayer(mp_Layer, (getLayerNameFromFileName(NewFileName)).c_str(), nullptr);

	OGRDataSource::DestroyDataSource(NewDS);
}


// =====================================================================
// =====================================================================


void VectorProcessing::deleteFile(std::string FullFilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FullFilePath.empty())
	{
		FullFilePath = m_FullFilePath;
	}

	if(openfluid::tools::Filesystem::isFile(FullFilePath))
	{
		OGRDataSource *DS = OGRSFDriverRegistrar::Open((FullFilePath).c_str(), TRUE);

		std::string DriverName = DS->GetDriver()->GetName();

		OGRDataSource::DestroyDataSource(DS);

		if(DriverName != m_DriverName)
		{
			std::cout << "VectorProcessing::deleteFile(): vector file " << FullFilePath << " could not be deleted because vector driver " << DriverName << " "
					"is not supported" << std::endl;
		}
		else
		{
			mp_Driver->DeleteDataSource(FullFilePath.c_str());
		}
	}
}


// =====================================================================
// =====================================================================


void VectorProcessing::correctGeometry()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	correctEmpty();
	correctMultyparts();
	correctDuplicates();
    correctOverlaps();
    correctHoles();
	correctWithinPolygons();
	correctMultyparts();
}


// =====================================================================
// =====================================================================


std::vector<unsigned int> VectorProcessing::findWithinPolygons()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	m_FIDlist.clear();

	OGRFeature *Feature = nullptr;

	for(unsigned int i = 0; i < mp_Layer->GetFeatureCount(); i++)
	{
		mp_Feature = mp_Layer->GetFeature(i);
		mp_Geometry = mp_Feature->GetGeometryRef()->clone();

		mp_Layer->ResetReading();

		while((Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			OGRGeometry* Geometry = Feature->GetGeometryRef()->clone();
			if((mp_Geometry->Within(Geometry)) and (i != (Feature->GetFID())) and (!Geometry->Equals(mp_Geometry)))
			{
				m_FIDlist.push_back(i);
			}
			OGRFeature::DestroyFeature(Feature);
			delete(Geometry);
		}
		OGRFeature::DestroyFeature(mp_Feature);
		delete(mp_Geometry);
	}

	return m_FIDlist;
}


// =====================================================================
// =====================================================================


void VectorProcessing::correctWithinPolygons()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	findWithinPolygons();

	if(m_FIDlist.size() != 0)
	{
		for(unsigned int i = 0; i < m_FIDlist.size(); i++)
		{
			mp_Layer->DeleteFeature(m_FIDlist[i]);
		}
	}

	repack();

}


// =====================================================================
// =====================================================================


std::vector<unsigned int> VectorProcessing::findOverlaps(unsigned int FID)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	m_FIDlist.clear();

	mp_Feature = mp_Layer->GetFeature(FID);
	mp_Geometry = mp_Feature->GetGeometryRef()->clone();

	OGRFeature *Feature = nullptr;
	OGRGeometry *Geometry = nullptr;

	mp_Layer->ResetReading();

	while((Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		Geometry = Feature->GetGeometryRef()->clone();
		int FeatureFID = Feature->GetFID();
		if (FID != FeatureFID)
		{
			if (mp_Geometry->Overlaps(Geometry))
			{
				m_FIDlist.push_back(FeatureFID);
			}
		}
		OGRFeature::DestroyFeature(Feature);
		delete(Geometry);
	}

	OGRFeature::DestroyFeature(mp_Feature);
	delete(mp_Geometry);

	return m_FIDlist;

}


// =====================================================================
// =====================================================================


void VectorProcessing::correctOverlaps()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	unsigned int N = mp_Layer->GetFeatureCount();

	for (unsigned int i = 0; i < N; i++)
	{
		findOverlaps(i);
		if(m_FIDlist.size() != 0)
		{
			for (unsigned int j = 0; j < m_FIDlist.size(); j++)
			{
				mp_Feature = mp_Layer->GetFeature(i);
				mp_Geometry = mp_Feature->GetGeometryRef()->clone();
				OGRFeature *Feature = mp_Layer->GetFeature(m_FIDlist[j]);
				OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();
				geos::geom::Geometry *GeometryToCorrectGEOS = (geos::geom::Geometry*) openfluid::landr::convertOGRGeometryToGEOS(mp_Geometry);
				geos::geom::Geometry *GeometryGEOS = (geos::geom::Geometry*) openfluid::landr::convertOGRGeometryToGEOS(Geometry);
				geos::geom::Geometry *DifferenceGeometryGEOS = GeometryToCorrectGEOS->difference(GeometryGEOS);
				OGRGeometry *NewGeometry = openfluid::landr::convertGEOSGeometryToOGR((GEOSGeom) DifferenceGeometryGEOS);
				if(NewGeometry->getGeometryType() == 3)
				{
					mp_Feature->SetGeometry(NewGeometry);
					mp_Layer->SetFeature(mp_Feature);
					mp_Layer->SyncToDisk();
				}
				else if(NewGeometry->getGeometryType() == 6)
				{
					mp_Feature->SetGeometry(((OGRMultiPolygon *) NewGeometry)->getGeometryRef(0));
					mp_Layer->SetFeature(mp_Feature);

					for(unsigned int i = 1; i < ((OGRMultiPolygon *) NewGeometry)->getNumGeometries(); i++)
					{
						OGRFeature *NewFeature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());
						NewFeature->SetFrom(mp_Feature);
						NewFeature->SetGeometry(((OGRMultiPolygon *) NewGeometry)->getGeometryRef(i));
						mp_Layer->CreateFeature(NewFeature);
						mp_Layer->SyncToDisk();
						OGRFeature::DestroyFeature(NewFeature);
					}

					mp_Layer->SyncToDisk();
				}
				else if(NewGeometry->getGeometryType() == 7)
				{
					std::vector<int> PolygonList;

					for(unsigned int i = 0; i < ((OGRGeometryCollection *) NewGeometry)->getNumGeometries(); i++)
					{
						if(((OGRGeometryCollection *) NewGeometry)->getGeometryRef(i)->getGeometryType() == 3)
						{
							PolygonList.push_back(i);
						}
					}

					if(PolygonList.size() != 0)
					{
						mp_Feature->SetGeometry(((OGRGeometryCollection *) NewGeometry)->getGeometryRef(PolygonList[0]));
						mp_Layer->SetFeature(mp_Feature);

						for(unsigned int j = 1; j < PolygonList.size(); j++)
						{
							OGRFeature *NewFeature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());
							NewFeature->SetFrom(mp_Feature);
							NewFeature->SetGeometry(((OGRGeometryCollection *) NewGeometry)->getGeometryRef(PolygonList[j]));
							mp_Layer->CreateFeature(NewFeature);
							mp_Layer->SyncToDisk();
							OGRFeature::DestroyFeature(NewFeature);
						}
					}

					PolygonList.clear();
					mp_Layer->SyncToDisk();
				}
				else
				{
					mp_Layer->SyncToDisk();
				}

				delete(NewGeometry);
				delete DifferenceGeometryGEOS;
				delete GeometryToCorrectGEOS;
				delete GeometryGEOS;
				delete(mp_Geometry);
				delete(Geometry);
				OGRFeature::DestroyFeature(Feature);
				OGRFeature::DestroyFeature(mp_Feature);
			}
		}
	}
}


// =====================================================================
// =====================================================================


void VectorProcessing::correctMultyparts()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		mp_Geometry = mp_Feature->GetGeometryRef()->clone();
		if (mp_Geometry->getGeometryType() == wkbMultiPolygon)
		{
			OGRMultiPolygon *MultiPolygon = (OGRMultiPolygon *) mp_Geometry;

			for (unsigned int i = 1; i < MultiPolygon->getNumGeometries(); i++)
			{
				OGRFeature *NewFeature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());
				NewFeature->SetFrom(mp_Feature);
				NewFeature->SetGeometry(MultiPolygon->getGeometryRef(i));
				mp_Layer->CreateFeature(NewFeature);
				mp_Layer->SyncToDisk();
				OGRFeature::DestroyFeature(NewFeature);
			}

			mp_Feature->SetGeometry(MultiPolygon->getGeometryRef(0));
			mp_Layer->SetFeature(mp_Feature);
			mp_Layer->SyncToDisk();
		}
		OGRFeature::DestroyFeature(mp_Feature);
		delete(mp_Geometry);
	}

}


// =====================================================================
// =====================================================================


void VectorProcessing::correctHoles()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_Layer->ResetReading();

	OGRLinearRing *LinearRing = nullptr;
	OGRPolygon *Polygon = nullptr, *NewPolygon = nullptr;
	OGRPoint *Point;

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		Polygon = (OGRPolygon *) mp_Feature->GetGeometryRef();
		if(Polygon->getNumInteriorRings() != 0)
		{
			LinearRing = new OGRLinearRing;
			NewPolygon = new OGRPolygon;
			for (int j = 0; j < Polygon->getExteriorRing()->getNumPoints(); j++)
			{
				Point = new OGRPoint;
				Polygon->getExteriorRing()->getPoint(j,Point);
				LinearRing->addPoint(Point->getX(),Point->getY());
				delete Point;
			}
			NewPolygon->addRing(LinearRing);
			OGRGeometry *NewGeometry = (OGRGeometry *) NewPolygon;
			mp_Feature->SetGeometry(NewGeometry);
			mp_Layer->SetFeature(mp_Feature);
			delete LinearRing;
			delete NewPolygon;
		}
		OGRFeature::DestroyFeature(mp_Feature);
	}

	repack();

}


// =====================================================================
// =====================================================================


std::vector <unsigned int> VectorProcessing::findDuplicates(unsigned int FID)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	OGRFeature *FeatureInQuestion = mp_Layer->GetFeature(FID);
	OGRGeometry *GeometryInQuestion = FeatureInQuestion->GetGeometryRef()->clone();

	m_FIDlist.clear();

	mp_Layer->ResetReading();

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		mp_Geometry = mp_Feature->GetGeometryRef()->clone();
		if (GeometryInQuestion->Equals(mp_Geometry) and FID != mp_Feature->GetFID())
		{
			m_FIDlist.push_back(mp_Feature->GetFID());
		}
		OGRFeature::DestroyFeature(mp_Feature);
		delete(mp_Geometry);
	}

	OGRFeature::DestroyFeature(FeatureInQuestion);
	delete(GeometryInQuestion);

	return m_FIDlist;
}


// =====================================================================
// =====================================================================


void VectorProcessing::correctDuplicates()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	unsigned int N = mp_Layer->GetFeatureCount();

	mp_Layer->ResetReading();

	for(unsigned int i = 0; i < N; i++)
	{
		findDuplicates(i);
		if(m_FIDlist.size() != 0)
		{
			mp_Layer->DeleteFeature(i);
		}
	}

	repack();
}


// =====================================================================
// =====================================================================


void VectorProcessing::correctEmpty()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_Layer->ResetReading();

	for(unsigned int i = 0; i < mp_Layer->GetFeatureCount(); i++)
	{
		mp_Feature = mp_Layer->GetFeature(i);
		mp_Geometry = mp_Feature->GetGeometryRef();
		if(mp_Geometry->getGeometryType() == wkbPolygon)
		{
			if (((OGRPolygon*) (mp_Geometry))->get_Area() == 0)
			{
				mp_Layer->DeleteFeature(i);
			}
		}
		else if(mp_Geometry->getGeometryType() == wkbMultiPolygon)
		{
			double Area = 0;
			for(unsigned int j = 0; j < ((OGRMultiPolygon*) (mp_Geometry))->getNumGeometries(); j++)
			{
				Area += ((OGRPolygon*)(((OGRMultiPolygon*) (mp_Geometry))->getGeometryRef(j)))->get_Area();
			}
			if (Area == 0)
			{
				mp_Layer->DeleteFeature(i);
			}
		}

		OGRFeature::DestroyFeature(mp_Feature);
	}

	repack();
}


// =====================================================================
// =====================================================================


void VectorProcessing::createIDField()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(doesFieldExist(m_IDFieldName))
	{
		deleteField(m_IDFieldName);
	}

	createField();

	for (int i = 0; i < mp_Layer->GetFeatureCount(); i++)
	{
		mp_Feature = mp_Layer->GetFeature(i);
		mp_Feature->SetField(m_IDFieldName.c_str(), i+1);
		mp_Layer->SetFeature(mp_Feature);
		OGRFeature::DestroyFeature(mp_Feature);
	}
}


// =====================================================================
// =====================================================================


bool VectorProcessing::doesFieldExist(std::string FieldName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FieldName.empty())
	{
		FieldName = m_IDFieldName;
	}

	unsigned int N = 0;

	for (unsigned int i = 0; i < mp_Layer->GetLayerDefn()->GetFieldCount(); i++)
	{
		if((std::string) mp_Layer->GetLayerDefn()->GetFieldDefn(i)->GetNameRef() != FieldName)
		{
			N++;
		}
	}

	if(N == mp_Layer->GetLayerDefn()->GetFieldCount())
	{
		return false;
	}
	else
	{
		return true;
	}
}


// =====================================================================
// =====================================================================


void VectorProcessing::createField(OGRFieldType FieldType, std::string FieldName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FieldName.empty())
	{
		FieldName = m_IDFieldName;
	}

	if(doesFieldExist(FieldName))
	{
		int FieldIndex = mp_Layer->GetLayerDefn()->GetFieldIndex(FieldName.c_str());
		mp_Layer->DeleteField(FieldIndex);
	}

	OGRFieldDefn FieldDefn(FieldName.c_str(), FieldType);

	mp_Layer->CreateField(&FieldDefn);

}


// =====================================================================
// =====================================================================


void VectorProcessing::deleteField(std::string FieldName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FieldName.empty())
	{
		FieldName = m_IDFieldName;
	}

	if(doesFieldExist(FieldName))
	{
		unsigned int FieldIndex = mp_Layer->GetLayerDefn()->GetFieldIndex(FieldName.c_str());
		mp_Layer->DeleteField(FieldIndex);
	}
}


// =====================================================================
// =====================================================================


void VectorProcessing::deleteUnnecessaryFields(std::vector<std::string> FieldNamesToKeep)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<std::string> FieldNameList;

	for (unsigned int i = 0; i < mp_Layer->GetLayerDefn()->GetFieldCount(); i++)
	{
		if(std::find(FieldNamesToKeep.begin(), FieldNamesToKeep.end(), (std::string) mp_Layer->GetLayerDefn()->GetFieldDefn(i)->GetNameRef()) == FieldNamesToKeep.end())
		{
			FieldNameList.push_back(mp_Layer->GetLayerDefn()->GetFieldDefn(i)->GetNameRef());
		}
	}

	if (FieldNameList.size() != 0)
	{
		for (unsigned int i = 0; i < FieldNameList.size(); i++)
		{
			deleteField(FieldNameList[i]);
		}
	}

}


// =====================================================================
// =====================================================================


void VectorProcessing::repack(const std::string &FullFilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FullFilePath.empty())
	{
		mp_DataSource->ExecuteSQL(("REPACK " + m_LayerName).c_str(), nullptr, nullptr);

		mp_DataSource->SyncToDisk();
	}
	else
	{
		OGRDataSource *NewDS = OGRSFDriverRegistrar::Open((FullFilePath).c_str(), TRUE);

		NewDS->ExecuteSQL(("REPACK " + getLayerNameFromFilePath(FullFilePath)).c_str(), nullptr, nullptr);

		NewDS->SyncToDisk();

		OGRDataSource::DestroyDataSource(NewDS);
	}
}


// =====================================================================
// =====================================================================


std::vector<std::pair<double,double>> VectorProcessing::getLayerExtent()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<std::pair<double,double>> Table;

	Table.clear();

	OGREnvelope* Envelope = new OGREnvelope();

	mp_Layer->GetExtent(Envelope);

	double MaxY, MinY, MaxX, MinX;

	MaxY = Envelope->MaxY;
	MinY = Envelope->MinY;
	MaxX = Envelope->MaxX;
	MinX = Envelope->MinX;

	Table.push_back(std::make_pair(MaxY,MinY));
	Table.push_back(std::make_pair(MaxX,MinX));

	delete Envelope;

	return Table;

}


// =====================================================================
// =====================================================================


void VectorProcessing::createDataSource(OGRwkbGeometryType GeometryType, std::string FilePath, std::string FileName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::string FullFilePath = FilePath + "/" + FileName;

	deleteFile(FullFilePath);

	OGRDataSource *DS = mp_Driver->CreateDataSource(FullFilePath.c_str());

	OGRLayer *Layer = DS->CreateLayer((getLayerNameFromFileName(FileName)).c_str(), mp_SRS, GeometryType);

	OGRDataSource::DestroyDataSource(DS);

}


// =====================================================================
// =====================================================================


void VectorProcessing::polygonizeRaster(std::pair<OGRFieldType,std::string> FieldInformation, const std::string &RasterFilePath, const std::string &RasterFileName, std::string FullFilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FullFilePath.empty())
	{
		FullFilePath = m_FullFilePath;
	}

	CPLErr ERR;

	createField(FieldInformation.first, (FieldInformation.second).c_str());

	RasterProcessing RP(RasterFilePath, RasterFileName, GA_ReadOnly);

	ERR = GDALPolygonize(RP.mp_Band, nullptr, mp_Layer, 0, nullptr, nullptr, nullptr);

	mp_Layer->SyncToDisk();

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		if (mp_Feature->GetFieldAsInteger((FieldInformation.second).c_str()) == -9999)
		{
			mp_Feature->SetField((FieldInformation.second).c_str(), 0);
			mp_Layer->SetFeature(mp_Feature);
		}
		OGRFeature::DestroyFeature(mp_Feature);
	}

}


// =====================================================================
// =====================================================================


std::pair <double, double> VectorProcessing::getCentroidCoordinates()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_Point = new OGRPoint;
	mp_Geometry->Centroid(mp_Point);

	std::pair <double, double> Coordinates;
	Coordinates = std::make_pair(mp_Point->getX(), mp_Point->getY());

	delete mp_Point;
	return Coordinates;

}


// =====================================================================
// =====================================================================


void VectorProcessing::fillFieldFromRaster(const std::string &FieldName, const std::string &RasterFilePath, const std::string &RasterFileName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	RasterProcessing RP(RasterFilePath, RasterFileName, GA_ReadOnly);

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		mp_Geometry = mp_Feature->GetGeometryRef()->clone();
		int Value = RP.extractFromRasterToPoint(getCentroidCoordinates());
		if (Value == -9999)
		{
			mp_Feature->SetField(FieldName.c_str(), 0);
			mp_Layer->SetFeature(mp_Feature);
		}
		else
		{
			mp_Feature->SetField(FieldName.c_str(), Value);
			mp_Layer->SetFeature(mp_Feature);
		}
		OGRFeature::DestroyFeature(mp_Feature);
		delete(mp_Geometry);
	}

}


// =====================================================================
// =====================================================================


void VectorProcessing::executeSQLRequest(const std::string &SQLExpression, std::string NewFileName, std::string NewFilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(NewFilePath.empty())
	{
		NewFilePath = m_FilePath;
	}

	std::string NewFileFullPath = NewFilePath + "/" + NewFileName;

	deleteFile(NewFileFullPath);

	std::string NewLayerName = getLayerNameFromFileName(NewFileName);

	OGRLayer *Layer = mp_DataSource->ExecuteSQL(SQLExpression.c_str(), nullptr, m_SQLDialect.c_str());

	mp_DataSource->CopyLayer(Layer, NewLayerName.c_str(), nullptr);

	mp_DataSource->ReleaseResultSet(Layer);

}


// =====================================================================
// =====================================================================


std::vector < std::vector <int>> VectorProcessing::transformToVector()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector < std::vector <int>> Table;

	Table.resize(mp_Layer->GetLayerDefn()->GetFieldCount());

	for(unsigned int i = 0 ; i < mp_Layer->GetLayerDefn()->GetFieldCount(); ++i)
	{
		Table[i].resize(mp_Layer->GetFeatureCount());
	}

	for (unsigned int i = 0; i < mp_Layer->GetLayerDefn()->GetFieldCount(); i++)
	{
		for (unsigned int j = 0; j < mp_Layer->GetFeatureCount(); j++)
		{
			mp_Feature = nullptr;
			mp_Feature = mp_Layer->GetFeature(j);
			Table[i][j] = mp_Feature->GetFieldAsInteger(i);
			OGRFeature::DestroyFeature(mp_Feature);
		  }
	  }

	return Table;
}


// =====================================================================
// =====================================================================


void VectorProcessing::deleteFeatures(std::vector<unsigned int> FIDList)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(!FIDList.empty())
	{
		mp_Layer->ResetReading();

		for (unsigned int i = 0; i < FIDList.size(); i++)
		{
			mp_Layer->DeleteFeature(FIDList[i]);
		}
	}
}


// =====================================================================
// =====================================================================


unsigned int VectorProcessing::getFeaturesCount(OGRLayer *Layer)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	unsigned int NFeatures = 0;

	if(!Layer)
	{
		NFeatures = mp_Layer->GetFeatureCount();
	}
	else
	{
		NFeatures = Layer->GetFeatureCount();
	}

	return NFeatures;
}


// =====================================================================
// =====================================================================


void VectorProcessing::copyDataSourceToDataSource(const std::string &NewFileName, const std::string &NewFilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	OGRDataSource *DS = mp_Driver->CopyDataSource(mp_DataSource, (NewFilePath + "/" + NewFileName).c_str(), nullptr);

	OGRDataSource::DestroyDataSource(DS);
}


// =====================================================================
// =====================================================================


void VectorProcessing::setDefaultFieldValue(const std::string &FieldName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_Layer->ResetReading();

	int FieldIndex = mp_Layer->GetLayerDefn()->GetFieldIndex(FieldName.c_str());

	OGRFieldType FieldType = mp_Layer->GetLayerDefn()->GetFieldDefn(FieldIndex)->GetType();

	if (FieldType == OFTInteger or FieldType == OFTReal)
	{
		mp_Layer->ResetReading();
		while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			mp_Feature->SetField(FieldName.c_str(), 0);
			mp_Layer->SetFeature(mp_Feature);
			OGRFeature::DestroyFeature(mp_Feature);
		}
	}
	else if(FieldType == OFTString)
	{
		while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			mp_Feature->SetField(FieldName.c_str(), "None");
			mp_Layer->SetFeature(mp_Feature);
			OGRFeature::DestroyFeature(mp_Feature);
		}
	}
}

