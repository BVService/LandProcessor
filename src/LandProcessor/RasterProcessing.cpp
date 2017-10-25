/*
 * RasterProcessing.cpp
 *
 *  Created on: 21 sept. 2017
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
#include <LandProcessor/Helpers.hpp>


RasterProcessing::RasterProcessing(const std::string &FilePath, const std::string &FileName, GDALAccess AccessType):
						m_FilePath(FilePath), m_FileName(FileName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	OGRRegisterAll();
	GDALAllRegister();

	std::setlocale(LC_NUMERIC, "C");

	if (!openfluid::tools::Filesystem::isDirectory(m_FilePath))
	{
		throw std::runtime_error("RasterProcessing::RasterProcessing(): input raster directory " + m_FilePath + " does not exists");
	}

	if (!openfluid::tools::Filesystem::isFile(m_FilePath + "/" + m_FileName))
	{
		throw std::runtime_error("RasterProcessing::RasterProcessing(): input raster file " + m_FileName + " does not exists in the "
				"input directory " + m_FilePath);
	}

	m_FullFilePath = m_FilePath + "/" + m_FileName;

	mp_Driver = GetGDALDriverManager()->GetDriverByName(m_DriverName.c_str());

	mp_Dataset = (GDALDataset *) GDALOpen( m_FullFilePath.c_str(), AccessType); //GA_Update or GA_ReadOnly

	if (mp_Dataset->GetDriver() != mp_Driver)
	{
		throw std::runtime_error("RasterProcessing::RasterProcessing(): raster driver " + ((std::string) mp_Dataset->GetDriver()->GetDescription()) + " is not supported");
	}

	if ((!mp_Dataset->GetProjectionRef()))
	{
		throw std::runtime_error("RasterProcessing::RasterProcessing(): " + m_FullFilePath + ": no spatial reference information is provided for the raster file");
	}

	m_ProjectionRef = mp_Dataset->GetProjectionRef();

	mp_SRS.SetFromUserInput(mp_Dataset->GetProjectionRef());

	mp_Band = mp_Dataset->GetRasterBand(1);

	GDALGetGeoTransform(mp_Dataset, m_GeoTransform);

}


// =====================================================================
// =====================================================================


RasterProcessing::~RasterProcessing()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(mp_Dataset)
	{
		GDALClose( mp_Dataset );
	}
}


// =====================================================================
// =====================================================================


void RasterProcessing::deleteFile(std::string FullFilePath)
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FullFilePath.empty())
	{
		FullFilePath = m_FullFilePath;
	}

	if(openfluid::tools::Filesystem::isFile(FullFilePath))
	{
		GDALDataset *DS = (GDALDataset *) GDALOpen( m_FullFilePath.c_str(), GA_ReadOnly); //GA_Update or GA_ReadOnly

		std::string DriverName = DS->GetDriver()->GetDescription();

		GDALClose( DS );

		if(DriverName != m_DriverName)
		{
			std::cout << "RasterProcessing::deleteFile(): raster file " << FullFilePath << " could not be deleted because raster driver " << DriverName << " "
					"is not supported" << std::endl;
		}
		else
		{
			mp_Driver->Delete(FullFilePath.c_str());
		}
	}
}


// =====================================================================
// =====================================================================


void RasterProcessing::checkEPSG()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_SRS.SetFromUserInput(mp_Dataset->GetProjectionRef());

	if (!mp_SRS.GetAttrValue("AUTHORITY",1))
	{
		throw std::runtime_error("RasterProcessing::RasterProcessing(): " + m_FullFilePath + ": no EPSG code provided for the raster data");
	}

	if ((std::string) (mp_SRS.GetAttrValue("AUTHORITY",1)) != m_EPSGCode.substr(m_EPSGCode.find(":")+1, m_EPSGCode.length()-1))
	{
		throw std::runtime_error("RasterProcessing::RasterProcessing(): " + m_FullFilePath + ": raster data EPSG code does not correspond to the default EPSG code");
	}

}


// =====================================================================
// =====================================================================


std::vector<double> RasterProcessing::getGeoTransformValues()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<double> GeoTransformValues;

	for (unsigned int i = 0; i < 6; i++)
	{
		GeoTransformValues.push_back(m_GeoTransform[i]);
	}

	return GeoTransformValues;
}


// =====================================================================
// =====================================================================


int RasterProcessing::extractFromRasterToPoint(std::pair <double, double> Coordinates)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	int Value;
	int *ScanValue = (int*) CPLMalloc(sizeof(int));

	GDALRasterIO(mp_Band, GF_Read,
			calculateOffset(Coordinates).first,
			calculateOffset(Coordinates).second,
			1, 1,
			ScanValue,
			1, 1, GDT_Int32,
			0, 0);

	Value = ScanValue[0];
	CPLFree(ScanValue);

	return Value;

}


//======================================================================
//======================================================================


std::pair <double, double> RasterProcessing::calculateOffset(std::pair <double, double> Coordinates)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::pair <double, double> Offset;

	Offset = std::make_pair(
      (Coordinates.first - m_GeoTransform[0])/m_GeoTransform[1],
      (Coordinates.second - m_GeoTransform[3])/m_GeoTransform[5]);

	return Offset;

}


//======================================================================
//======================================================================


std::pair <double, double> RasterProcessing::getMinMaxValues()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::pair <double, double> MinMaxValues;

	double MinValue, MaxValue;

	mp_Band->GetStatistics(1, 1, &MinValue, &MaxValue,nullptr,nullptr);

	MinMaxValues = std::make_pair(MinValue, MaxValue);

	return MinMaxValues;

}


//======================================================================
//======================================================================


std::pair <int, int> RasterProcessing::getRasterSize()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::pair <int, int> RasterSize;

	int XSize, YSize;

	XSize = mp_Dataset->GetRasterXSize();
	YSize = mp_Dataset->GetRasterYSize();

	RasterSize = std::make_pair(XSize, YSize);

	return RasterSize;

}


//======================================================================
//======================================================================


void RasterProcessing::createIDRaster(std::string FileName, std::string FilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FilePath.empty())
	{
		FilePath = m_FilePath;
	}

	std::string FullFilePath = FilePath + "/" + FileName;

	deleteFile(FullFilePath);

	double XSize = getRasterSize().first, YSize = getRasterSize().second;

	GDALDataset *IDRaster = mp_Driver->Create( FullFilePath.c_str(), XSize, YSize, 1, GDT_Int32, nullptr );
	GDALRasterBand *Band = IDRaster->GetRasterBand(1);

	int *Row = (int*) CPLMalloc(sizeof(int)*XSize);

	for (unsigned int i = 0; i < YSize; i++)
	{
		for (unsigned int j = 0; j < XSize; j++)
		{
			if (i == 0)
			{
				Row[j] = j+1;
			}
			else
			{
				Row[j] = XSize*i+(j+1);
			}
		}
		Band->RasterIO( GF_Write, 0, i, XSize, 1,
				Row, XSize, 1, GDT_Int32, 0, 0 );
	}

	CPLFree(Row);

	Band->FlushCache();
	Band->SetNoDataValue(-9999);
	IDRaster->SetGeoTransform(m_GeoTransform);
	IDRaster->SetProjection(m_ProjectionRef);

	GDALClose( IDRaster );

}


//======================================================================
//======================================================================


void RasterProcessing::createOutletRaster(std::string PlotsRasterName,
		std::string DrainageRasterName, std::string IDRasterName, std::string FileName, std::string FilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FilePath.empty())
	{
		FilePath = m_FilePath;
	}

	std::string FullFilePath = FilePath + "/" + FileName;

	deleteFile(FullFilePath);

	if (!openfluid::tools::Filesystem::isFile(FilePath + "/" + PlotsRasterName) or !openfluid::tools::Filesystem::isFile(FilePath + "/" + DrainageRasterName) or
			!openfluid::tools::Filesystem::isFile(FilePath + "/" + IDRasterName))
	{
		throw std::runtime_error("RasterProcessing::createOutletRaster(): operation could not be complete; one or several necessary files are missing or could not be opened");
	}

	double XSize = getRasterSize().first, YSize = getRasterSize().second;

	GDALDataset *Outlets = mp_Driver->Create( FullFilePath.c_str(), XSize, YSize, 1, GDT_Int32, nullptr );
	GDALRasterBand *Band = Outlets->GetRasterBand(1);

	GDALDataset *Plots = (GDALDataset *) GDALOpen((FilePath + "/" + PlotsRasterName).c_str(), GA_ReadOnly );
	GDALDataset *Drainage = (GDALDataset *) GDALOpen((FilePath + "/" + DrainageRasterName).c_str(), GA_ReadOnly );
	GDALDataset *IDRaster = (GDALDataset *) GDALOpen((FilePath + "/" + IDRasterName).c_str(), GA_ReadOnly );

	if (!Plots or !Drainage or !IDRaster)
	{
		if(Plots)
		{
			GDALClose( Plots );
		}
		if(Drainage)
		{
			GDALClose( Drainage );
		}
		if(IDRaster)
		{
			GDALClose( IDRaster );
		}
		throw std::runtime_error("RasterProcessing::createOutletRaster(): operation could not be complete because one or several "
				"necessary files could not be opened");
	}

	int *TargetRow, *PlotsUpperRow, *PlotsMiddleRow, *PlotsLowerRow, *DrainageRow, *IDRasterRow;

	TargetRow = (int*) CPLMalloc(sizeof(int)*XSize);
	PlotsUpperRow = (int*) CPLMalloc(sizeof(int)*XSize);
	PlotsMiddleRow = (int*) CPLMalloc(sizeof(int)*XSize);
	PlotsLowerRow = (int*) CPLMalloc(sizeof(int)*XSize);
	DrainageRow = (int*) CPLMalloc(sizeof(int)*XSize);
	IDRasterRow = (int*) CPLMalloc(sizeof(int)*XSize);

	for (unsigned int i = 0; i < YSize; i++)
	{
		for (unsigned int j = 0; j < XSize; j++)
		{
			if (i == 0 or j == 0 or i == YSize-1 or j == XSize-1)
			{
				TargetRow[j] = -9999;
			}
			else
			{
				Plots->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, PlotsMiddleRow, XSize, 1, GDT_Int32, 0, 0 );
				if (PlotsMiddleRow[j] == -9999)
				{
					TargetRow[j] = -9999;
				}
				else
				{
					Plots->GetRasterBand(1)->RasterIO( GF_Read, 0, i-1, XSize, 1, PlotsUpperRow, XSize, 1, GDT_Int32, 0, 0 );
					Plots->GetRasterBand(1)->RasterIO( GF_Read, 0, i+1, XSize, 1, PlotsLowerRow, XSize, 1, GDT_Int32, 0, 0 );
					Drainage->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, DrainageRow, XSize, 1, GDT_Int32, 0, 0 );
					IDRaster->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, IDRasterRow, XSize, 1, GDT_Int32, 0, 0 );
					if ((DrainageRow[j] == 1 and PlotsMiddleRow[j]!= PlotsUpperRow[j+1]) or
							(DrainageRow[j] == 2 and PlotsMiddleRow[j] != PlotsUpperRow[j]) or
							(DrainageRow[j] == 3 and PlotsMiddleRow[j] != PlotsUpperRow[j-1]) or
							(DrainageRow[j] == 4 and PlotsMiddleRow[j] != PlotsMiddleRow[j-1]) or
							(DrainageRow[j] == 5 and PlotsMiddleRow[j] != PlotsLowerRow[j-1]) or
							(DrainageRow[j] == 6 and PlotsMiddleRow[j] != PlotsLowerRow[j]) or
							(DrainageRow[j] == 7 and PlotsMiddleRow[j] != PlotsLowerRow[j+1]) or
							(DrainageRow[j] == 8 and PlotsMiddleRow[j] != PlotsMiddleRow[j+1]))
					{
						TargetRow[j] = IDRasterRow[j];
					}
					else
					{
						TargetRow[j] = -9999;
					}
				}
			}
		}
		Band->RasterIO( GF_Write, 0, i, XSize, 1, TargetRow, XSize, 1, GDT_Int32, 0, 0 );
	}

	CPLFree(TargetRow);
	CPLFree(PlotsUpperRow);
	CPLFree(PlotsMiddleRow);
	CPLFree(PlotsLowerRow);
	CPLFree(DrainageRow);
	CPLFree(IDRasterRow);

	Band->FlushCache();
	Band->SetNoDataValue(-9999);
	Outlets->SetGeoTransform(m_GeoTransform);
	Outlets->SetProjection(m_ProjectionRef);

	GDALClose( Outlets );
	GDALClose( Plots );
	GDALClose( Drainage );
	GDALClose( IDRaster );

}


//======================================================================
//======================================================================


void RasterProcessing::createReceiversRaster(std::string OutletsRasterName,
		std::string DrainageRasterName, std::string IDRasterName, std::string FileName, std::string FilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FilePath.empty())
	{
		FilePath = m_FilePath;
	}

	std::string FullFilePath = FilePath + "/" + FileName;

	deleteFile(FullFilePath);

	if (!openfluid::tools::Filesystem::isFile(FilePath + "/" + OutletsRasterName) or !openfluid::tools::Filesystem::isFile(FilePath + "/" + DrainageRasterName) or
			!openfluid::tools::Filesystem::isFile(FilePath + "/" + IDRasterName))
	{
		throw std::runtime_error("RasterProcessing::createReceiversRaster(): operation could not be complete because one or several "
				"necessary files are missing");
	}

	double XSize = getRasterSize().first, YSize = getRasterSize().second;

	GDALDataset *Receivers = mp_Driver->Create(FullFilePath.c_str(), XSize, YSize, 1, GDT_Int32, nullptr );
	GDALRasterBand *Band = Receivers->GetRasterBand(1);

	GDALDataset *Outlets = (GDALDataset *) GDALOpen((FilePath + "/" + OutletsRasterName).c_str(), GA_ReadOnly );
	GDALDataset *Drainage = (GDALDataset *) GDALOpen((FilePath + "/" + DrainageRasterName).c_str(), GA_ReadOnly );
	GDALDataset *IDRaster = (GDALDataset *) GDALOpen((FilePath + "/" + IDRasterName).c_str(), GA_ReadOnly );

	if (!Outlets or !Drainage or !IDRaster)
	{
		if(Outlets)
		{
			GDALClose( Outlets );
		}
		if(Drainage)
		{
			GDALClose( Drainage );
		}
		if(IDRaster)
		{
			GDALClose( IDRaster );
		}
		throw std::runtime_error("RasterProcessing::createReceiversRaster(): operation could not be complete because one or several "
				"necessary files could not be opened");
	}

	int *TargetRow, *OutletsUpperRow, *OutletsMiddleRow, *OutletsLowerRow, *DrainageUpperRow, *DrainageMiddleRow, *DrainageLowerRow, *IDRasterRow;

	TargetRow = (int*) CPLMalloc(sizeof(int)*XSize);
	OutletsUpperRow = (int*) CPLMalloc(sizeof(int)*XSize);
	OutletsMiddleRow = (int*) CPLMalloc(sizeof(int)*XSize);
	OutletsLowerRow = (int*) CPLMalloc(sizeof(int)*XSize);
	DrainageUpperRow = (int*) CPLMalloc(sizeof(int)*XSize);
	DrainageMiddleRow = (int*) CPLMalloc(sizeof(int)*XSize);
	DrainageLowerRow = (int*) CPLMalloc(sizeof(int)*XSize);
	IDRasterRow = (int*) CPLMalloc(sizeof(int)*XSize);

	for (unsigned int i = 0; i < YSize; i++)
	{
		for (unsigned int j = 0; j < XSize; j++)
		{
			  if (i == 0 or j == 0 or i == YSize-1 or j == XSize-1)
			  {
				  TargetRow[j] = -9999;
			  }
			  else
			  {
				  Outlets->GetRasterBand(1)->RasterIO( GF_Read, 0, i-1, XSize, 1, OutletsUpperRow, XSize, 1, GDT_Int32, 0, 0 );
				  Outlets->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, OutletsMiddleRow, XSize, 1, GDT_Int32, 0, 0 );
				  Outlets->GetRasterBand(1)->RasterIO( GF_Read, 0, i+1, XSize, 1, OutletsLowerRow, XSize, 1, GDT_Int32, 0, 0 );
				  Drainage->GetRasterBand(1)->RasterIO( GF_Read, 0, i-1, XSize, 1, DrainageUpperRow, XSize, 1, GDT_Int32, 0, 0 );
				  Drainage->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, DrainageMiddleRow, XSize, 1, GDT_Int32, 0, 0 );
				  Drainage->GetRasterBand(1)->RasterIO( GF_Read, 0, i+1, XSize, 1, DrainageLowerRow, XSize, 1, GDT_Int32, 0, 0 );
				  IDRaster->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, IDRasterRow, XSize, 1, GDT_Int32, 0, 0 );
				  if ((DrainageUpperRow[j-1] == 7 and OutletsUpperRow[j-1] != -9999) or
						  (DrainageUpperRow[j] == 6 and OutletsUpperRow[j] != -9999) or
						  (DrainageUpperRow[j+1] == 5 and OutletsUpperRow[j+1] != -9999) or
						  (DrainageMiddleRow[j+1] == 4 and OutletsMiddleRow[j+1] != -9999) or
						  (DrainageLowerRow[j+1] == 3 and OutletsLowerRow[j+1] != -9999) or
						  (DrainageLowerRow[j] == 2 and OutletsLowerRow[j] != -9999) or
						  (DrainageLowerRow[j-1] == 1 and OutletsLowerRow[j-1] != -9999) or
						  (DrainageMiddleRow[j-1] == 8 and OutletsMiddleRow[j-1] != -9999))
				  {
					  TargetRow[j] = IDRasterRow[j];
				  }
				  else
				  {
					  TargetRow[j] = -9999;
				  }
			  }
		  }
		  Band->RasterIO( GF_Write, 0, i, XSize, 1, TargetRow, XSize, 1, GDT_Int32, 0, 0 );
	  }

	  CPLFree(OutletsUpperRow);
	  CPLFree(OutletsMiddleRow);
	  CPLFree(OutletsLowerRow);
	  CPLFree(DrainageUpperRow);
	  CPLFree(DrainageMiddleRow);
	  CPLFree(DrainageLowerRow);
	  CPLFree(IDRasterRow);
	  CPLFree(TargetRow);

	  Band->FlushCache();
	  Band->SetNoDataValue(-9999);
	  Receivers->SetGeoTransform(m_GeoTransform);
	  Receivers->SetProjection(m_ProjectionRef);

	  GDALClose( Receivers );
	  GDALClose( Outlets );
	  GDALClose( Drainage );
	  GDALClose( IDRaster );

}


//======================================================================
//======================================================================


void RasterProcessing::createConnectingRaster(std::string OutletsRasterName,
		std::string DrainageRasterName, std::string IDRasterName, std::string FileName, std::string FilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FilePath.empty())
	{
		FilePath = m_FilePath;
	}

	std::string FullFilePath = FilePath + "/" + FileName;

	deleteFile(FullFilePath);

	if (!openfluid::tools::Filesystem::isFile(FilePath + "/" + OutletsRasterName) or !openfluid::tools::Filesystem::isFile(FilePath + "/" + DrainageRasterName) or
			!openfluid::tools::Filesystem::isFile(FilePath + "/" + IDRasterName))
	{
		throw std::runtime_error("RasterProcessing::createConnectingRaster(): operation could not be complete because one or several "
				"necessary files are missing");
	}

	double XSize = getRasterSize().first, YSize = getRasterSize().second;

	GDALDataset *Connecting = mp_Driver->Create( FullFilePath.c_str(), XSize, YSize, 1, GDT_Int32, nullptr );
	GDALRasterBand *Band = Connecting->GetRasterBand(1);

	GDALDataset *Outlets = (GDALDataset *) GDALOpen((FilePath + "/" + OutletsRasterName).c_str(), GA_ReadOnly );
	GDALDataset *Drainage = (GDALDataset *) GDALOpen((FilePath + "/" + DrainageRasterName).c_str(), GA_ReadOnly );
	GDALDataset *IDRaster = (GDALDataset *) GDALOpen((FilePath + "/" + IDRasterName).c_str(), GA_ReadOnly );

	if (!Outlets or !Drainage or !IDRaster)
	{
		if(Outlets)
		{
			GDALClose( Outlets );
		}
		if(Drainage)
		{
			GDALClose( Drainage );
		}
		if(IDRaster)
		{
			GDALClose( IDRaster );
		}
		throw std::runtime_error("RasterProcessing::createConnectingRaster(): operation could not be complete because one or several "
				"necessary files could not be opened");
	}

	int *TargetRow, *OutletsRow, *IDRasterUpperRow, *IDRasterMiddleRow, *IDRasterLowerRow, *DrainageRow;

	TargetRow = (int*) CPLMalloc(sizeof(int)*XSize);
	OutletsRow = (int*) CPLMalloc(sizeof(int)*XSize);
	DrainageRow = (int*) CPLMalloc(sizeof(int)*XSize);
	IDRasterUpperRow = (int*) CPLMalloc(sizeof(int)*XSize);
	IDRasterMiddleRow = (int*) CPLMalloc(sizeof(int)*XSize);
	IDRasterLowerRow = (int*) CPLMalloc(sizeof(int)*XSize);

	for (unsigned int i = 0; i < YSize; i++)
	{
		  for (unsigned int j = 0; j<XSize; j++)
		  {
			  if (i == 0 or j == 0 or i == YSize-1 or j == XSize-1)
			  {
				  TargetRow[j] = -9999;
			  }
			  else
			  {
				  Outlets->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, OutletsRow, XSize, 1, GDT_Int32, 0, 0 );
				  Drainage->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, DrainageRow, XSize, 1, GDT_Int32, 0, 0 );
				  IDRaster->GetRasterBand(1)->RasterIO( GF_Read, 0, i-1, XSize, 1, IDRasterUpperRow, XSize, 1, GDT_Int32, 0, 0 );
				  IDRaster->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, IDRasterMiddleRow, XSize, 1, GDT_Int32, 0, 0 );
				  IDRaster->GetRasterBand(1)->RasterIO( GF_Read, 0, i+1, XSize, 1, IDRasterLowerRow, XSize, 1, GDT_Int32, 0, 0 );
				  if (DrainageRow[j] == 1 and OutletsRow[j] != -9999)
				  {
					  TargetRow[j] = IDRasterUpperRow[j+1];
				  }
				  else if (DrainageRow[j] == 2 and OutletsRow[j] != -9999)
				  {
					  TargetRow[j] = IDRasterUpperRow[j];
				  }
				  else if (DrainageRow[j] == 3 and OutletsRow[j] != -9999)
				  {
					  TargetRow[j] = IDRasterUpperRow[j-1];
				  }
				  else if (DrainageRow[j] == 4 and OutletsRow[j] != -9999)
				  {
					  TargetRow[j] = IDRasterMiddleRow[j-1];
				  }
				  else if (DrainageRow[j] == 5 and OutletsRow[j] != -9999)
				  {
					  TargetRow[j] = IDRasterLowerRow[j-1];
				  }
				  else if (DrainageRow[j] == 6 and OutletsRow[j] != -9999)
				  {
					  TargetRow[j] = IDRasterLowerRow[j];
				  }
				  else if (DrainageRow[j] == 7 and OutletsRow[j] != -9999)
				  {
					  TargetRow[j] = IDRasterLowerRow[j+1];
				  }
				  else if (DrainageRow[j] == 8 and OutletsRow[j] != -9999)
				  {
					  TargetRow[j] = IDRasterMiddleRow[j+1];
				  }
				  else
				  {
					  TargetRow[j] = -9999;
				  }
			  }
		  }
		  Band->RasterIO( GF_Write, 0, i, XSize, 1, TargetRow, XSize, 1, GDT_Int32, 0, 0 );
	}

	CPLFree(OutletsRow);
	CPLFree(DrainageRow);
	CPLFree(IDRasterUpperRow);
	CPLFree(IDRasterMiddleRow);
	CPLFree(IDRasterLowerRow);
	CPLFree(TargetRow);

	Band->FlushCache();
	Band->SetNoDataValue(-9999);
	Connecting->SetGeoTransform(m_GeoTransform);
	Connecting->SetProjection(m_ProjectionRef);

	GDALClose( Connecting );
	GDALClose( Outlets );
	GDALClose( Drainage );
	GDALClose( IDRaster );

}


//======================================================================
//======================================================================


void RasterProcessing::createCatchmentsRaster(std::string OutletsRasterName, std::string DrainageRasterName, std::string PlotsRasterName, std::string FileName, std::string FilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(FilePath.empty())
	{
		FilePath = m_FilePath;
	}

	std::string FullFilePath = FilePath + "/" + FileName;

	deleteFile(FullFilePath);

	if (!openfluid::tools::Filesystem::isFile(FilePath + "/" + OutletsRasterName) or !openfluid::tools::Filesystem::isFile(FilePath + "/" + DrainageRasterName) or
			!openfluid::tools::Filesystem::isFile(FilePath + "/" + PlotsRasterName))
	{
		throw std::runtime_error("RasterProcessing::createCatchmentsRaster(): operation could not be complete because one or several "
				"necessary files are missing");
	}

	double XSize = getRasterSize().first, YSize = getRasterSize().second;

	GDALDataset *Outlets = (GDALDataset *) GDALOpen((FilePath + "/" + OutletsRasterName).c_str(), GA_ReadOnly );

	if (!Outlets)
	{
		throw std::runtime_error("RasterProcessing::createCatchmentsRaster(): operation could not be complete because one or several "
				"necessary files could not be opened");
	}

	GDALDataset *Catchments;

	Catchments = mp_Driver->CreateCopy( FullFilePath.c_str(), Outlets, FALSE,
	                                nullptr, nullptr, nullptr );
	GDALRasterBand *Band = Catchments->GetRasterBand(1);

	GDALClose( Outlets );

	GDALDataset *Plots = (GDALDataset *) GDALOpen((FilePath + "/" + PlotsRasterName).c_str(), GA_ReadOnly );
	GDALDataset *Drainage = (GDALDataset *) GDALOpen((FilePath + "/" + DrainageRasterName).c_str(), GA_ReadOnly );

	if (!Plots or !Drainage)
	{
		if(Drainage)
		{
			GDALClose( Drainage );
		}
		if(Plots)
		{
			GDALClose( Plots );
		}
		throw std::runtime_error("RasterProcessing::createCatchmentsRaster(): operation could not be complete because one or several "
				"necessary files could not be opened");
	}

	int *TargetRow, *CatchmentsUpperRow, *CatchmentsMiddleRow, *CatchmentsLowerRow,	*PlotsUpperRow, *PlotsMiddleRow, *PlotsLowerRow, *DrainageRow;

	TargetRow = (int*) CPLMalloc(sizeof(int) *XSize);
	CatchmentsUpperRow = (int*) CPLMalloc(sizeof(int) *XSize);
	CatchmentsMiddleRow = (int*) CPLMalloc(sizeof(int) *XSize);
	CatchmentsLowerRow = (int*) CPLMalloc(sizeof(int) *XSize);
	PlotsUpperRow = (int*) CPLMalloc(sizeof(int) *XSize);
	PlotsMiddleRow = (int*) CPLMalloc(sizeof(int) *XSize);
	PlotsLowerRow = (int*) CPLMalloc(sizeof(int) *XSize);
	DrainageRow = (int*) CPLMalloc(sizeof(int) *XSize);

	unsigned int PlotsCount = 0, CatchmentsCount = 0;

	for (unsigned int i = 0; i < YSize; i++)
	{
		Plots->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, PlotsMiddleRow, XSize,
				1, GDT_Int32, 0, 0 );
		for (unsigned int j = 0; j < XSize; j++)
		{
			if (PlotsMiddleRow[j] != -9999)
			{
				PlotsCount += 1;
			}
		}
	}

	for (unsigned int i = 0; i < YSize; i++)
	{
		Catchments->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, CatchmentsMiddleRow, XSize, 1, GDT_Int32, 0, 0 );
		for (unsigned int j = 0; j < XSize; j++)
		{
			if (CatchmentsMiddleRow[j] != -9999)
			{
				CatchmentsCount += 1;
			}
		}
	}

	while (CatchmentsCount != PlotsCount)
	{
		for (unsigned int i = 0; i < YSize; i++)
		{
			for (unsigned int j = 0; j < XSize; j++)
			{
				if (i == 0 or j == 0 or i == YSize-1 or j == XSize-1)
				{
					TargetRow[j] = -9999;
				}
				else
				{
					Catchments->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, CatchmentsMiddleRow, XSize, 1, GDT_Int32, 0, 0 );
					if (CatchmentsMiddleRow[j] == -9999)
					{
						Catchments->GetRasterBand(1)->RasterIO( GF_Read, 0, i-1, XSize, 1, CatchmentsUpperRow, XSize, 1, GDT_Int32, 0, 0 );
						Catchments->GetRasterBand(1)->RasterIO( GF_Read, 0, i+1, XSize, 1, CatchmentsLowerRow, XSize, 1, GDT_Int32, 0, 0 );
						Plots->GetRasterBand(1)->RasterIO( GF_Read, 0, i-1, XSize, 1, PlotsUpperRow, XSize, 1, GDT_Int32, 0, 0 );
						Plots->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, PlotsMiddleRow, XSize, 1, GDT_Int32, 0, 0 );
						Plots->GetRasterBand(1)->RasterIO( GF_Read, 0, i+1, XSize, 1, PlotsLowerRow, XSize, 1, GDT_Int32, 0, 0 );
						Drainage->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, DrainageRow, XSize, 1, GDT_Int32, 0, 0 );
						if (CatchmentsUpperRow[j+1] != -9999 and DrainageRow[j] == 1 and PlotsMiddleRow[j] == PlotsUpperRow[j+1])
						{
							TargetRow[j] = CatchmentsUpperRow[j+1];
						}
						else if (CatchmentsUpperRow[j] != -9999 and DrainageRow[j] == 2 and PlotsMiddleRow[j] == PlotsUpperRow[j])
						{
							TargetRow[j] = CatchmentsUpperRow[j];
						}
						else if (CatchmentsUpperRow[j-1] != -9999 and DrainageRow[j] == 3 and PlotsMiddleRow[j] == PlotsUpperRow[j-1])
						{
							TargetRow[j] = CatchmentsUpperRow[j-1];
						}
						else if (CatchmentsMiddleRow[j-1] != -9999 and DrainageRow[j] == 4 and PlotsMiddleRow[j] == PlotsMiddleRow[j-1])
						{
							TargetRow[j] = CatchmentsMiddleRow[j-1];
						}
						else if (CatchmentsLowerRow[j-1] != -9999 and DrainageRow[j] == 5 and PlotsMiddleRow[j] == PlotsLowerRow[j-1])
						{
							TargetRow[j] = CatchmentsLowerRow[j-1];
						}
						else if (CatchmentsLowerRow[j] != -9999 and DrainageRow[j] == 6 and PlotsMiddleRow[j] == PlotsLowerRow[j])
						{
							TargetRow[j] = CatchmentsLowerRow[j];
						}
						else if (CatchmentsLowerRow[j+1] != -9999 and DrainageRow[j] == 7 and PlotsMiddleRow[j] == PlotsLowerRow[j+1])
						{
							TargetRow[j] = CatchmentsLowerRow[j+1];
						}
						else if (CatchmentsMiddleRow[j+1] != -9999 and DrainageRow[j] == 8 and PlotsMiddleRow[j] == PlotsMiddleRow[j+1])
						{
							TargetRow[j] = CatchmentsMiddleRow[j+1];
						}
						else
						{
							TargetRow[j] = -9999;
						}
					}
					else
					{
						TargetRow[j] = CatchmentsMiddleRow[j];
					}
				}
			}
			Band->RasterIO( GF_Write, 0, i, XSize, 1, TargetRow, XSize, 1, GDT_Int32, 0, 0 );
		}

		PlotsCount = 0;

		for (unsigned int i = 0; i < YSize; i++)
		{
			Plots->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, PlotsMiddleRow, XSize, 1, GDT_Int32, 0, 0 );
			for (unsigned int j = 0; j < XSize; j++)
			{
				if (PlotsMiddleRow[j] != -9999)
				{
					PlotsCount += 1;
				}
			}
		}
		CatchmentsCount = 0;
		for (unsigned int i = 0; i < YSize; i++)
		{
			Catchments->GetRasterBand(1)->RasterIO( GF_Read, 0, i, XSize, 1, CatchmentsMiddleRow, XSize, 1, GDT_Int32, 0, 0 );
			for (unsigned int j = 0; j < XSize; j++)
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
	CPLFree(DrainageRow);
	CPLFree(TargetRow);

	Band->FlushCache();
	Band->SetNoDataValue(-9999);
	Catchments->SetGeoTransform(m_GeoTransform);
	Catchments->SetProjection(m_ProjectionRef);

	GDALClose( Catchments );
	GDALClose( Plots );
	GDALClose( Drainage );
}
