/*
 * LinearEntitiesProcessing.cpp
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
#include <LandProcessor/ArealEntitiesProcessing.hpp>
#include <LandProcessor/LinearEntitiesProcessing.hpp>
#include <LandProcessor/Helpers.hpp>


LinearEntitiesProcessing::LinearEntitiesProcessing(const std::string &FilePath, const std::string &FileName):
				VectorProcessing(FilePath, FileName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	OGRRegisterAll();
	GDALAllRegister();

	std::setlocale(LC_NUMERIC, "C");

}


// =====================================================================
// =====================================================================


LinearEntitiesProcessing::~LinearEntitiesProcessing()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);
}


// =====================================================================
// =====================================================================


void LinearEntitiesProcessing::setLinearEntitiesFaces(const std::string &BaseVectorFileName, std::string OutputVectorFilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(OutputVectorFilePath.empty())
	{
		OutputVectorFilePath = m_FilePath;
	}

	std::vector<std::string> FacesTable, AttributesNames = {"FaceA","FaceB"};;

	OGRDataSource *DS = OGRSFDriverRegistrar::Open( (OutputVectorFilePath + "/" + BaseVectorFileName).c_str(), TRUE );
	OGRLayer *Layer = DS->GetLayer(0);
	OGRFeature *FeatureA = nullptr, *FeatureB = nullptr;

	Layer->ResetReading();

	while ((FeatureA = Layer->GetNextFeature()) != nullptr)
	{
	    std::string IDA = FeatureA->GetFieldAsString("ID");
	    unsigned int IDPlotA = std::stoi(IDA.substr(0, IDA.find("N")));

	    if (IDPlotA != 0)
	    {
	    	OGRGeometry *GeometryA = FeatureA->GetGeometryRef()->clone();

	    	for(unsigned int i = 0; i < Layer->GetFeatureCount(); i++)
	    	{
	    		FeatureB = Layer->GetFeature(i);
	    		OGRGeometry *GeometryB = FeatureB->GetGeometryRef()->clone();
	    		std::string IDB = FeatureB->GetFieldAsString("ID");
	    		unsigned int IDPlotB = std::stoi(IDB.substr(0, IDB.find("N")));

	    		if(IDPlotA != IDPlotB)
	    		{
	    			std::string Face = IDB + "-" + IDA;

	    			if (std::find(FacesTable.begin(), FacesTable.end(), Face.c_str()) == FacesTable.end())
	    			{
	    				OGRGeometry *IntersectionGeometry = GeometryA->Intersection(GeometryB);

	    				if(!IntersectionGeometry->IsEmpty() and IntersectionGeometry->IsValid())
	    				{
	    					if(IntersectionGeometry->getGeometryType() == 2 or IntersectionGeometry->getGeometryType() == 5)
	    					{
	    						mp_Feature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());
	    						mp_Feature->SetField(AttributesNames[0].c_str(), IDA.c_str());
	    						mp_Feature->SetField(AttributesNames[1].c_str(), IDB.c_str());
	    						mp_Feature->SetGeometry(IntersectionGeometry);
	    						mp_Layer->CreateFeature(mp_Feature);
	    						OGRFeature::DestroyFeature(mp_Feature);
	    						Face = IDA + "-" + IDB;
	    						FacesTable.push_back(Face.c_str());
	    					}
	    					else if(IntersectionGeometry->getGeometryType() == 7)
	    					{
	    						for(unsigned int i = 0; i < ((OGRGeometryCollection *) IntersectionGeometry)->getNumGeometries(); i++)
	    						{
	    							if(((OGRGeometryCollection *) IntersectionGeometry)->getGeometryRef(i)->getGeometryType() == 2 or
	    									((OGRGeometryCollection *) IntersectionGeometry)->getGeometryRef(i)->getGeometryType() == 5)
	    							{
	    								mp_Feature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());
	    								mp_Feature->SetField(AttributesNames[0].c_str(), IDA.c_str());
	    								mp_Feature->SetField(AttributesNames[1].c_str(), IDB.c_str());
	    								mp_Feature->SetGeometry(((OGRGeometryCollection *) IntersectionGeometry)->getGeometryRef(i));
	    								mp_Layer->CreateFeature(mp_Feature);
	    								OGRFeature::DestroyFeature(mp_Feature);
	    								Face = IDA + "-" + IDB;
	    								FacesTable.push_back(Face.c_str());
	    							}
	    						}
	    					}
	    				}
	    				delete(IntersectionGeometry);
	    			}
	    		}
	    		OGRFeature::DestroyFeature(FeatureB);
	    		delete(GeometryB);
	    	}

	    	delete(GeometryA);
	    }
	    OGRFeature::DestroyFeature(FeatureA);
	}

	OGRDataSource::DestroyDataSource(DS);
}


// =====================================================================
// =====================================================================


void LinearEntitiesProcessing::setLinearEntitiesIDs(const std::string &ArealEntitiesVectorFileName, std::string OutputVectorFilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(OutputVectorFilePath.empty())
	{
		OutputVectorFilePath = m_FilePath;
	}

	ArealEntitiesProcessing AEP(OutputVectorFilePath, ArealEntitiesVectorFileName);

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		std::string FaceA = mp_Feature->GetFieldAsString("FaceA");
		std::string FaceB = mp_Feature->GetFieldAsString("FaceB");

		AEP.mp_Layer->ResetReading();

	    while ((AEP.mp_Feature = AEP.mp_Layer->GetNextFeature()) != nullptr)
	    {
	    	std::string ID = (std::string) AEP.mp_Feature->GetFieldAsString("ID");
	    	std::string ID2 = (std::string) AEP.mp_Feature->GetFieldAsString("ID2");

	    	if((ID == FaceA and ID2 == FaceB) or (ID == FaceB and ID2 == FaceA))
	    	{
	    		if(ID == FaceA and ID2 == FaceB)
	    		{
	    			std::string IDLNR = FaceA + "-" + FaceB;
	    			mp_Feature->SetField("ID", IDLNR.c_str());
	    			mp_Feature->SetField("IDTo", FaceB.c_str());
	    			mp_Layer->SetFeature(mp_Feature);
	    		}
	    		else
	    		{
	    			std::string IDLNR = FaceB + "-" + FaceA;
	    			mp_Feature->SetField("ID", IDLNR.c_str());
	    			mp_Feature->SetField("IDTo", FaceA.c_str());
	    			mp_Layer->SetFeature(mp_Feature);
	    		}
	    	}
	    	OGRFeature::DestroyFeature(AEP.mp_Feature);
	    }
    	OGRFeature::DestroyFeature(mp_Feature);
	}

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		std::string ID = mp_Feature->GetFieldAsString("ID");
		std::string ID2 = mp_Feature->GetFieldAsString("IDTo");
		if(ID.empty() and ID2.empty())
		{
			std::string FaceA = mp_Feature->GetFieldAsString("FaceA");
			std::string FaceB = mp_Feature->GetFieldAsString("FaceB");
			std::string IDLNR = FaceA + "-" + FaceB;
			std::string IDLNRTo = "None";
			mp_Feature->SetField("ID", IDLNR.c_str());
			mp_Feature->SetField("IDTo", IDLNRTo.c_str());
			mp_Layer->SetFeature(mp_Feature);
		}
    	OGRFeature::DestroyFeature(mp_Feature);
	}
}


// =====================================================================
// =====================================================================


void LinearEntitiesProcessing::createLinearEntitiesAttributes()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	createField(OFTReal, "Length");
}


// =====================================================================
// =====================================================================


void LinearEntitiesProcessing::setLinearEntitiesLength()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		mp_Geometry = mp_Feature->GetGeometryRef()->clone();
		if(mp_Geometry->getGeometryType() == 2)
		{
			mp_Feature->SetField("Length", ((OGRLineString *) mp_Geometry)->get_Length());
			mp_Layer->SetFeature(mp_Feature);
			OGRFeature::DestroyFeature(mp_Feature);
			delete(mp_Geometry);
		}
		else if (mp_Geometry->getGeometryType() == 5)
		{
			double Length = 0;
			OGRMultiLineString* LNRMultiLineString = (OGRMultiLineString*) mp_Geometry;
			for (unsigned int i = 0; i < LNRMultiLineString->getNumGeometries(); i++)
			{
				Length += ((OGRLineString*) LNRMultiLineString->getGeometryRef(i))->get_Length();
			}
			mp_Feature->SetField("Length", Length);
			mp_Layer->SetFeature(mp_Feature);
			OGRFeature::DestroyFeature(mp_Feature);
			delete(mp_Geometry);
		}
	}
}


// =====================================================================
// =====================================================================


void LinearEntitiesProcessing::creatLinearStructuresIntersectionVector(const std::string &IntersectionVectorFileName, std::string OutputVectorFilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(OutputVectorFilePath.empty())
	{
		OutputVectorFilePath = m_FilePath;
	}

	createDataSource(wkbLineString, OutputVectorFilePath, IntersectionVectorFileName);

	{
		VectorProcessing VP(OutputVectorFilePath, IntersectionVectorFileName);

		std::vector<std::pair<OGRFieldType, std::string>> FieldTypesAndNames;

		std::vector<std::string> FieldNames = {m_IDFieldName, "IDLS", "FIDPlot"};

		for (unsigned int i = 0; i < FieldNames.size(); i++)
		{
			FieldTypesAndNames.push_back(std::make_pair(OFTInteger, FieldNames[i]));
		}

		for (unsigned int i = 0; i < FieldTypesAndNames.size(); i++)
		{
			VP.createField(FieldTypesAndNames[i].first, FieldTypesAndNames[i].second);
		}
	}
}


// =====================================================================
// =====================================================================


void LinearEntitiesProcessing::setLinearStructuresIntersectionVector(const std::string &LinearStructureFileName, const std::string &ParcelVectorFileName, std::string OutputVectorFilePath) //Intersection
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(OutputVectorFilePath.empty())
	{
		OutputVectorFilePath = m_FilePath;
	}

	VectorProcessing Parcel(OutputVectorFilePath, ParcelVectorFileName);

	VectorProcessing LinearStructure(OutputVectorFilePath, LinearStructureFileName);

	OGRLayer *LinearStructureLayer = LinearStructure.getLayer();
	OGRLayer *ParcelLayer = Parcel.getLayer();

	OGRFeature *LinearStructureFeature = nullptr, *ParcelFeature = nullptr;

	while ((LinearStructureFeature = LinearStructureLayer->GetNextFeature()) != nullptr)
	{
		OGRGeometry *LinearStructureGeometry = LinearStructureFeature->GetGeometryRef()->clone();
		ParcelLayer->ResetReading();

		while ((ParcelFeature = ParcelLayer->GetNextFeature()) != nullptr)
		{
			OGRGeometry *ParcelGeometry = ParcelFeature->GetGeometryRef()->clone();
			OGRGeometry *NewGeometry = LinearStructureGeometry->Intersection(ParcelGeometry);

			if (!NewGeometry->IsEmpty())
			{
				if (NewGeometry->getGeometryType() == 2)
				{
					OGRFeature *NewFeature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());
					NewFeature->SetField(m_IDFieldName.c_str(), mp_Layer->GetFeatureCount()+1);
					NewFeature->SetField("IDLS", LinearStructureFeature->GetFieldAsInteger(m_IDFieldName.c_str()));
					NewFeature->SetField("FIDPlot", (int) ParcelFeature->GetFID());
					NewFeature->SetGeometry(NewGeometry);
					mp_Layer->CreateFeature(NewFeature);
					mp_Layer->SyncToDisk();
					OGRFeature::DestroyFeature(NewFeature);
				}
				else if (NewGeometry->getGeometryType() == 5)
				{
					for (unsigned int i = 0; i < ((OGRMultiLineString *) NewGeometry)->getNumGeometries(); i++)
					{
						OGRFeature *NewFeature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());
						NewFeature->SetField(m_IDFieldName.c_str(), mp_Layer->GetFeatureCount()+1);
						NewFeature->SetField("IDLS", LinearStructureFeature->GetFieldAsInteger(m_IDFieldName.c_str()));
						NewFeature->SetField("FIDPlot", (int) ParcelFeature->GetFID());
						NewFeature->SetGeometry(((OGRMultiLineString *) NewGeometry)->getGeometryRef(i));
						mp_Layer->CreateFeature(NewFeature);
						mp_Layer->SyncToDisk();
						OGRFeature::DestroyFeature(NewFeature);
					}
				}
				else if (NewGeometry->getGeometryType() == 7)
				{
					for (unsigned int i = 0; i < ((OGRGeometryCollection *) NewGeometry)->getNumGeometries(); i++)
					{
						if (((OGRGeometryCollection *) NewGeometry)->getGeometryRef(i)->getGeometryType() == 2)
						{
							OGRFeature *NewFeature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());
							NewFeature->SetField(m_IDFieldName.c_str(), mp_Layer->GetFeatureCount()+1);
							NewFeature->SetField("IDLS", LinearStructureFeature->GetFieldAsInteger(m_IDFieldName.c_str()));
							NewFeature->SetField("FIDPlot", (int) ParcelFeature->GetFID());
							NewFeature->SetGeometry(((OGRGeometryCollection *) NewGeometry)->getGeometryRef(i));
							mp_Layer->CreateFeature(NewFeature);
							mp_Layer->SyncToDisk();
							OGRFeature::DestroyFeature(NewFeature);
						}
						else if (((OGRGeometryCollection *) NewGeometry)->getGeometryRef(i)->getGeometryType() == 5)
						{
							for (unsigned int j = 0; j < ((OGRMultiLineString *) ((OGRGeometryCollection *) NewGeometry)->getGeometryRef(i))->getNumGeometries(); j++)
							{
								OGRFeature *NewFeature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());
								NewFeature->SetField(m_IDFieldName.c_str(), mp_Layer->GetFeatureCount()+1);
								NewFeature->SetField("IDLS", LinearStructureFeature->GetFieldAsInteger(m_IDFieldName.c_str()));
								NewFeature->SetField("FIDPlot", (int) ParcelFeature->GetFID());
								NewFeature->SetGeometry(((OGRMultiLineString *) ((OGRGeometryCollection *) NewGeometry)->getGeometryRef(i))->getGeometryRef(j));
								mp_Layer->CreateFeature(NewFeature);
								mp_Layer->SyncToDisk();
								OGRFeature::DestroyFeature(NewFeature);
							}
						}
					}
				}
				OGRFeature::DestroyFeature(ParcelFeature);
				delete(ParcelGeometry);
				delete(NewGeometry);
			}
			else
			{
				OGRFeature::DestroyFeature(ParcelFeature);
				delete(ParcelGeometry);
				delete(NewGeometry);
			}
		}
		OGRFeature::DestroyFeature(LinearStructureFeature);
		delete(LinearStructureGeometry);
	}
}


// =====================================================================
// =====================================================================


void LinearEntitiesProcessing::attributeLinearStructures(const std::string &ParcelVectorFileName, std::string OutputVectorFilePath)
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(OutputVectorFilePath.empty())
	{
		OutputVectorFilePath = m_FilePath;
	}

	ArealEntitiesProcessing AEP(OutputVectorFilePath, ParcelVectorFileName);

	std::vector<int> FIDToRemoveTable;

	FIDToRemoveTable.clear();

	mp_Layer->ResetReading();

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		mp_Geometry = mp_Feature->GetGeometryRef()->clone();
		mp_Geometry->segmentize(1);

		OGRLineString *IntersectionLineString = nullptr;
		IntersectionLineString = (OGRLineString *) mp_Geometry;

		AEP.mp_Feature = AEP.mp_Layer->GetFeature(mp_Feature->GetFieldAsInteger("FIDPlot"));
		AEP.mp_Geometry = AEP.mp_Feature->GetGeometryRef()->clone();

		AEP.mp_Geometry->segmentize(1);

		OGRPolygon *PlotsPolygon = nullptr;

		OGRLinearRing *PlotsRing = nullptr;

		PlotsPolygon = (OGRPolygon *) AEP.mp_Geometry;

		PlotsRing = PlotsPolygon->getExteriorRing();

		OGRPoint *sPoint, *ePoint, *RingPoint, *PointToAdd;

		sPoint = new OGRPoint;
		ePoint = new OGRPoint;

		int numPoints = IntersectionLineString->getNumPoints();

		IntersectionLineString->StartPoint(sPoint);
		IntersectionLineString->EndPoint(ePoint);

		std::vector <OGRPoint*> Points = {sPoint, ePoint};
		std::vector <int> PointsID = {0,0};

		PointsID.clear();
		double Distance = 99999999;

		for (unsigned int i = 0; i < Points.size(); i++)
		{
			for (unsigned int j = 0; j < PlotsRing->getNumPoints(); j++)
			{
				RingPoint = new OGRPoint;
				PlotsRing->getPoint(j, RingPoint);
				OGRGeometry *RingPointGeometry = nullptr;
				RingPointGeometry = (OGRGeometry *) RingPoint;
				if (((OGRGeometry *) Points[i])->Distance(RingPointGeometry) < Distance)
				{
					Distance = ((OGRGeometry *) Points[i])->Distance(RingPointGeometry);
					PointsID[i] = j;
				}
				delete RingPoint;
			}
			Distance = 99999999;
		}

		delete sPoint;
		delete ePoint;

		if(PointsID[0] == PointsID[1])
		{
			mp_Layer->DeleteFeature(mp_Feature->GetFID());
		}
		else
		{
			OGRLineString *FirstLineString = nullptr, *SecondLineString = nullptr;
			FirstLineString = new OGRLineString;
			SecondLineString = new OGRLineString;

			if (!PlotsRing->isClockwise())
			{
				PlotsRing->reverseWindingOrder();
			}

			if (PointsID[0] < PointsID[1])
			{
				for(unsigned int i = PointsID[0]; i <= PointsID[1]; i++)
				{
					PointToAdd = new OGRPoint;
					PlotsRing->getPoint(i, PointToAdd);
					FirstLineString->addPoint(PointToAdd->getX(), PointToAdd->getY());
					delete PointToAdd;
				}
				for(unsigned int i = PointsID[1]; i < PlotsRing->getNumPoints(); i++)
				{
					PointToAdd = new OGRPoint;
					PlotsRing->getPoint(i, PointToAdd);
					SecondLineString->addPoint(PointToAdd->getX(), PointToAdd->getY());
					delete PointToAdd;
				}
				for(unsigned int i = 0; i <= PointsID[0]; i++)
				{
					PointToAdd = new OGRPoint;
					PlotsRing->getPoint(i, PointToAdd);
					SecondLineString->addPoint(PointToAdd->getX(), PointToAdd->getY());
					delete PointToAdd;
				}
			}
			else if (PointsID[0] > PointsID[1])
			{
				for(unsigned int i = PointsID[0]; i < PlotsRing->getNumPoints(); i++)
				{
					PointToAdd = new OGRPoint;
					PlotsRing->getPoint(i, PointToAdd);
					FirstLineString->addPoint(PointToAdd->getX(), PointToAdd->getY());
					delete PointToAdd;
				}
				for(unsigned int i = 0; i <= PointsID[1]; i++)
				{
					PointToAdd = new OGRPoint;
					PlotsRing->getPoint(i, PointToAdd);
					FirstLineString->addPoint(PointToAdd->getX(), PointToAdd->getY());
					delete PointToAdd;
				}
				for(unsigned int i = PointsID[1]; i <= PointsID[0]; i++)
				{
					PointToAdd = new OGRPoint;
					PlotsRing->getPoint(i, PointToAdd);
					SecondLineString->addPoint(PointToAdd->getX(), PointToAdd->getY());
					delete PointToAdd;
				}
			}

			if (((OGRGeometry *) FirstLineString)->IsValid() and !((OGRGeometry *) FirstLineString)->IsEmpty() and
					((OGRGeometry *) FirstLineString)->IsValid() and !((OGRGeometry *) FirstLineString)->IsEmpty())
			{
				if(std::abs((IntersectionLineString->get_Length()) - FirstLineString->get_Length()) <= std::abs((IntersectionLineString->get_Length()) - SecondLineString->get_Length()))
				{
					mp_Feature->SetGeometry((OGRGeometry *) FirstLineString);
					mp_Layer->SetFeature(mp_Feature);
					mp_Layer->SyncToDisk();
				}
				else
				{
					mp_Feature->SetGeometry((OGRGeometry *) SecondLineString);
					mp_Layer->SetFeature(mp_Feature);
					mp_Layer->SyncToDisk();
				}
			}
			else
			{
				mp_Layer->DeleteFeature(mp_Feature->GetFID());
			}

			delete FirstLineString;
			delete SecondLineString;

			PointsID.clear();
			Points.clear();
		}

		OGRFeature::DestroyFeature(AEP.mp_Feature);
		OGRFeature::DestroyFeature(mp_Feature);
		delete(AEP.mp_Geometry);
		delete(mp_Geometry);
	}

	repack();
}


// =====================================================================
// =====================================================================


void LinearEntitiesProcessing::setAttributedLinearStructureLength(const std::string &LengthFieldName, const std::string &RelativeLengthFieldName, const std::string &IntersectionVectorFileName, std::string OutputVectorFilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(OutputVectorFilePath.empty())
	{
		OutputVectorFilePath = m_FilePath;
	}

	std::string LayerName = getLayerName();

	OGRDataSource *DS = OGRSFDriverRegistrar::Open( OutputVectorFilePath.c_str(), TRUE );

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		m_SQLRequest = "SELECT ST_Length(ST_Intersection(geometry,(SELECT ST_Buffer(geometry,0.05) FROM " + LayerName + " WHERE ROWID == " + std::to_string(mp_Feature->GetFID()) + "))) AS Length "
				"FROM " + getLayerNameFromFileName(IntersectionVectorFileName) + " "
				"WHERE ST_Length(ST_Intersection(geometry,(SELECT ST_Buffer(geometry,0.05) FROM " + LayerName + " WHERE ROWID == " + std::to_string(mp_Feature->GetFID()) + "))) > 0.1";

		OGRLayer *SQLLayer = DS->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
		if (SQLLayer->GetFeatureCount() > 0)
		{
			double Length = 0;
			for (unsigned int k = 0; k < SQLLayer->GetFeatureCount(); k++)
			{
				OGRFeature *SQLFeature = SQLLayer->GetFeature(k);
				Length += SQLFeature->GetFieldAsDouble("Length");
				OGRFeature::DestroyFeature(SQLFeature);
			}
			mp_Feature->SetField(LengthFieldName.c_str(), Length);
			mp_Layer->SetFeature(mp_Feature);
		}
		DS->ReleaseResultSet(SQLLayer);
		OGRFeature::DestroyFeature(mp_Feature);
	}

	mp_Layer->SyncToDisk();
	OGRDataSource::DestroyDataSource(DS);

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		double Length = mp_Feature->GetFieldAsDouble("Length");
		double LSLength = mp_Feature->GetFieldAsDouble(LengthFieldName.c_str());
		if (LSLength > Length)
		{
			mp_Feature->SetField(LengthFieldName.c_str(), Length);
			mp_Layer->SetFeature(mp_Feature);
		}
		OGRFeature::DestroyFeature(mp_Feature);
	}

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		double Length = mp_Feature->GetFieldAsDouble("Length");
		double LSLength = mp_Feature->GetFieldAsDouble(LengthFieldName.c_str());
		if (LSLength != Length and LSLength < 0.05*6)
		{
			mp_Feature->SetField(LengthFieldName.c_str(), 0);
			mp_Layer->SetFeature(mp_Feature);
		}
		OGRFeature::DestroyFeature(mp_Feature);
	}

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		if(mp_Feature->GetFieldAsDouble(LengthFieldName.c_str()) != 0)
		{
			double Length = mp_Feature->GetFieldAsDouble("Length");
			double LSLength = mp_Feature->GetFieldAsDouble(LengthFieldName.c_str());
			mp_Feature->SetField(RelativeLengthFieldName.c_str(), LSLength/Length);
			mp_Layer->SetFeature(mp_Feature);
		}
		mp_Layer->SyncToDisk();
		OGRFeature::DestroyFeature(mp_Feature);
	}

	mp_Layer->SyncToDisk();

}


// =====================================================================
// =====================================================================


void LinearEntitiesProcessing::setLinearEntitiesFlowDistance(const std::string &ArealEntityVectorFileName, std::string OutputVectorFilePath)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	if(OutputVectorFilePath.empty())
	{
		OutputVectorFilePath = m_FilePath;
	}

	OGRPoint *Centroid, *PointOnSurf;

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		std::string ID = mp_Feature->GetFieldAsString("ID");
		std::string IDTo = mp_Feature->GetFieldAsString("IDTo");

		mp_Geometry = mp_Feature->GetGeometryRef()->clone();

		if (IDTo.find("-") == std::string::npos)
		{
			ArealEntitiesProcessing AEP(OutputVectorFilePath, ArealEntityVectorFileName);

			AEP.mp_Layer->ResetReading();

			while((AEP.mp_Feature = AEP.mp_Layer->GetNextFeature()) != nullptr)
			{
				if((std::string)AEP.mp_Feature->GetFieldAsString("ID") == "SU#" + IDTo)
				{
					AEP.mp_Geometry = AEP.mp_Feature->GetGeometryRef()->clone();

					Centroid = new OGRPoint;

					AEP.mp_Geometry->Centroid(Centroid);

					geos::geom::Geometry *GeometryGEOS = (geos::geom::Geometry*) openfluid::landr::convertOGRGeometryToGEOS(AEP.mp_Geometry);

					geos::geom::Point *GeometryPointOnSurfGEOS = GeometryGEOS->getInteriorPoint();

					OGRGeometry *GeometryPointOnSurf = openfluid::landr::convertGEOSGeometryToOGR((GEOSGeom) ((geos::geom::Geometry *) GeometryPointOnSurfGEOS));
					PointOnSurf = (OGRPoint*) GeometryPointOnSurf->clone();

					delete GeometryGEOS;
					delete GeometryPointOnSurfGEOS;
					delete(GeometryPointOnSurf);

					if(Centroid->Within(AEP.mp_Geometry))
					{
						double Distance = mp_Geometry->Distance(Centroid);
						mp_Feature->SetField("FlowDist", Distance);
						mp_Layer->SetFeature(mp_Feature);
					}
					else
					{
						double Distance = mp_Geometry->Distance(PointOnSurf);
						mp_Feature->SetField("FlowDist", Distance);
						mp_Layer->SetFeature(mp_Feature);
					}
					delete Centroid;
					delete PointOnSurf;
					delete(AEP.mp_Geometry);
				}
				OGRFeature::DestroyFeature(AEP.mp_Feature);
			}
		}
		else
		{
			mp_Feature->SetField("FlowDist", 0);
			mp_Layer->SetFeature(mp_Feature);
		}
		OGRFeature::DestroyFeature(mp_Feature);
		delete(mp_Geometry);
	}
}
