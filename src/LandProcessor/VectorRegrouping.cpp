/*
 * VectorRegrouping.cpp
 *
 *  Created on: 4 oct. 2017
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

#include <LandProcessor/RasterProcessing.hpp>
#include <LandProcessor/VectorProcessing.hpp>
#include <LandProcessor/VectorRegrouping.hpp>
#include <LandProcessor/Helpers.hpp>


VectorRegrouping::VectorRegrouping(const std::string &FilePath, const std::string &FileName):
				VectorProcessing(FilePath, FileName)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	OGRRegisterAll();
	GDALAllRegister();

	std::setlocale(LC_NUMERIC, "C");

}


// =====================================================================
// =====================================================================


VectorRegrouping::~VectorRegrouping()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

}


// =====================================================================
// =====================================================================


std::vector<unsigned int> VectorRegrouping::getPlotIDWithLargeFIDList(unsigned int MinEntSize)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> PlotIDWithLargeFIDList;

	mp_Layer->ResetReading();

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		if (mp_Feature->GetGeometryRef()->getGeometryType() == 3)
		{
			if ((float) ((OGRPolygon *) mp_Feature->GetGeometryRef())->get_Area() > MinEntSize and
					std::find(PlotIDWithLargeFIDList.begin(), PlotIDWithLargeFIDList.end(), mp_Feature->GetFieldAsInteger("IDPlot")) == PlotIDWithLargeFIDList.end())
			{
				PlotIDWithLargeFIDList.push_back(mp_Feature->GetFieldAsInteger("IDPlot"));
			}
		}
		if (mp_Feature->GetGeometryRef()->getGeometryType() == 6)
		{
			if ((float) ((OGRMultiPolygon *) mp_Feature->GetGeometryRef())->get_Area() > MinEntSize and
					std::find(PlotIDWithLargeFIDList.begin(), PlotIDWithLargeFIDList.end(), mp_Feature->GetFieldAsInteger("IDPlot")) == PlotIDWithLargeFIDList.end())
			{
				PlotIDWithLargeFIDList.push_back(mp_Feature->GetFieldAsInteger("IDPlot"));
			}
		}
		OGRFeature::DestroyFeature(mp_Feature);
	}

	return PlotIDWithLargeFIDList;

}


//======================================================================
//======================================================================


std::vector<unsigned int> VectorRegrouping::getWithAllSmallFIDList(unsigned int MinEntSize)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> PlotIDWithAllSmallFIDList, PlotIDWithLargeFIDList = getPlotIDWithLargeFIDList(MinEntSize);

	mp_Layer->ResetReading();

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		if (mp_Feature->GetGeometryRef()->getGeometryType() == 3)
		{
			if ((float) ((OGRPolygon *) mp_Feature->GetGeometryRef())->get_Area() <= MinEntSize and
					std::find(PlotIDWithLargeFIDList.begin(), PlotIDWithLargeFIDList.end(), mp_Feature->GetFieldAsInteger("IDPlot")) == PlotIDWithLargeFIDList.end() and
					std::find(PlotIDWithAllSmallFIDList.begin(), PlotIDWithAllSmallFIDList.end(), mp_Feature->GetFieldAsInteger("IDPlot")) == PlotIDWithAllSmallFIDList.end())
			{
				PlotIDWithAllSmallFIDList.push_back(mp_Feature->GetFieldAsInteger("IDPlot"));
			}
		}
		if (mp_Feature->GetGeometryRef()->getGeometryType() == 6)
		{
			if ((float) ((OGRMultiPolygon *) mp_Feature->GetGeometryRef())->get_Area() <= MinEntSize and
					std::find(PlotIDWithLargeFIDList.begin(), PlotIDWithLargeFIDList.end(), mp_Feature->GetFieldAsInteger("IDPlot")) == PlotIDWithLargeFIDList.end() and
					std::find(PlotIDWithAllSmallFIDList.begin(), PlotIDWithAllSmallFIDList.end(), mp_Feature->GetFieldAsInteger("IDPlot")) == PlotIDWithAllSmallFIDList.end())
			{
				PlotIDWithAllSmallFIDList.push_back(mp_Feature->GetFieldAsInteger("IDPlot"));
			}
		}
		OGRFeature::DestroyFeature(mp_Feature);
	}

	return PlotIDWithAllSmallFIDList;

}


//======================================================================
//======================================================================


std::vector<unsigned int> VectorRegrouping::getSmallFIDList(unsigned int MinEntSize)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> SmallFIDList, PlotIDWithLargeFIDList = getPlotIDWithLargeFIDList(MinEntSize);

	mp_Layer->ResetReading();

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		if (mp_Feature->GetGeometryRef()->getGeometryType() == 3)
		{
			if ((float) ((OGRPolygon *) mp_Feature->GetGeometryRef())->get_Area() <= MinEntSize and
					std::find(PlotIDWithLargeFIDList.begin(), PlotIDWithLargeFIDList.end(), mp_Feature->GetFieldAsInteger("IDPlot")) != PlotIDWithLargeFIDList.end())
			{
				SmallFIDList.push_back(mp_Feature->GetFID());
			}
		}
		else if (mp_Feature->GetGeometryRef()->getGeometryType() == 6)
		{
			if ((float) ((OGRMultiPolygon *) mp_Feature->GetGeometryRef())->get_Area() <= MinEntSize and
					std::find(PlotIDWithLargeFIDList.begin(), PlotIDWithLargeFIDList.end(), mp_Feature->GetFieldAsInteger("IDPlot")) != PlotIDWithLargeFIDList.end())
			{
				SmallFIDList.push_back(mp_Feature->GetFID());
			}
		}
		OGRFeature::DestroyFeature(mp_Feature);
	}

	return SmallFIDList;
}


//======================================================================
//======================================================================


void VectorRegrouping::regroupSmallFeaturesThatConstituteSinglePlot(unsigned int MinEntSize)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> PlotIDList = getWithAllSmallFIDList(MinEntSize);

	if(!PlotIDList.empty())
	{
		std::string IDNew, ID2New, ID3New;
		int IDPlot2New, IDPlot3New, IDFID = -1;

		for (unsigned int i = 0; i < PlotIDList.size(); i++)
		{
			OGRDataSource *DataSource = OGRSFDriverRegistrar::Open( (m_FilePath + "/" + m_FileName).c_str(), TRUE );
			m_SQLRequest = "SELECT ROWID as RID, IDPlot, ID, IDPlot2, ID2, IDPlot3, ID3 FROM " + m_LayerName + " "
					"WHERE IDPlot == " + std::to_string(PlotIDList[i]) + " "
					"AND IDPlot != IDPlot3 ORDER BY ST_Area(geometry) DESC";
			OGRLayer *SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
			if (SQLLayer->GetFeatureCount() != 0)
			{
				OGRFeature *SQLFeature = SQLLayer->GetFeature(0);
				IDFID = SQLFeature->GetFieldAsInteger("RID");
				IDNew = SQLFeature->GetFieldAsString("ID");
				IDPlot2New = SQLFeature->GetFieldAsInteger("IDPlot2");
				ID2New = SQLFeature->GetFieldAsString("ID2");
				IDPlot3New = SQLFeature->GetFieldAsInteger("IDPlot3");
				ID3New = SQLFeature->GetFieldAsString("ID3");
				OGRFeature::DestroyFeature(SQLFeature);
				DataSource->ReleaseResultSet(SQLLayer);
			}
			else
			{
				DataSource->ReleaseResultSet(SQLLayer);
			}

			OGRDataSource::DestroyDataSource(DataSource);

			if (IDFID != -1)
			{
				std::vector <unsigned int> FIDList;
				mp_Layer->ResetReading();

				while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
				{
					if (mp_Feature->GetFieldAsInteger("IDPlot") == PlotIDList[i])
					{
						FIDList.push_back(mp_Feature->GetFID());
						OGRFeature::DestroyFeature(mp_Feature);
					}
					else
					{
						OGRFeature::DestroyFeature(mp_Feature);
					}
				}

				mp_DataSource->SyncToDisk();

				DataSource = OGRSFDriverRegistrar::Open( (m_FullFilePath).c_str(), TRUE );
				m_SQLRequest = "SELECT ST_Union(geometry) FROM " + m_LayerName + " WHERE IDPlot == " + std::to_string(PlotIDList[i]);
				SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
				OGRFeature *SQLFeature = SQLLayer->GetFeature(0);
				OGRGeometry *NewGeometry = SQLFeature->GetGeometryRef()->clone();
				DataSource->ReleaseResultSet(SQLLayer);
				OGRFeature::DestroyFeature(SQLFeature);
				OGRDataSource::DestroyDataSource(DataSource);

				OGRFeature *NewFeature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());
				NewFeature->SetField("IDPlot", (int) PlotIDList[i]);
				NewFeature->SetField("ID", IDNew.c_str());
				NewFeature->SetField("IDPlot2", IDPlot2New);
				NewFeature->SetField("ID2", ID2New.c_str());
				NewFeature->SetField("IDPlot3", IDPlot3New);
				NewFeature->SetField("ID3", ID3New.c_str());
				NewFeature->SetGeometry(NewGeometry);
				mp_Layer->CreateFeature(NewFeature);
				OGRFeature::DestroyFeature(NewFeature);
				delete(NewGeometry);

				deleteFeatures(FIDList);

				repack();

				mp_Layer->ResetReading();

				while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
				{
					if (mp_Feature->GetFieldAsInteger("IDPlot2") == PlotIDList[i])
					{
						mp_Feature->SetField("ID2", IDNew.c_str());
						mp_Feature->SetField("IDPlot3", IDPlot2New);
						mp_Feature->SetField("ID3", ID2New.c_str());
						mp_Layer->SetFeature(mp_Feature);
						OGRFeature::DestroyFeature(mp_Feature);
					}
					else
					{
						OGRFeature::DestroyFeature(mp_Feature);
					}
				}

				mp_Layer->ResetReading();

				while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
				{
					if (mp_Feature->GetFieldAsInteger("IDPlot3") == PlotIDList[i])
					{
						mp_Feature->SetField("ID3", IDNew.c_str());
						mp_Layer->SetFeature(mp_Feature);
						OGRFeature::DestroyFeature(mp_Feature);
					}
					else
					{
						OGRFeature::DestroyFeature(mp_Feature);
					}
				}
			}
		}
		mp_DataSource->SyncToDisk();
	}
}


//======================================================================
//======================================================================


void VectorRegrouping::regroupSmallFeaturesWithLargerFeatures(unsigned int MinEntSize)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> SmallFIDList = getSmallFIDList(MinEntSize);

	if (SmallFIDList.size() != 0)
	{
		std::vector <unsigned int> FIDToEraseFromFIDList, FIDToDeleteList, NewPlotIDList, NewPlotID2List;
		std::vector<std::string> OldIDList, NewIDList, NewID2List;

		FIDToDeleteList = SmallFIDList;

		while(SmallFIDList.size() != 0)
		{
			unsigned int Beginning = SmallFIDList.size();
			for (unsigned int i = 0; i < SmallFIDList.size(); i++)
			{
				mp_Feature = mp_Layer->GetFeature(SmallFIDList[i]);
				std::string ID = mp_Feature->GetFieldAsString("ID");
				mp_Geometry = mp_Feature->GetGeometryRef()->clone();
				int IDPlot = mp_Feature->GetFieldAsInteger("IDPlot");

				OGRDataSource *DataSource = OGRSFDriverRegistrar::Open((m_FullFilePath).c_str(), TRUE );
				m_SQLRequest = "SELECT ROWID as RID, * FROM " + m_LayerName + " WHERE "
						"IDPlot == " + std::to_string(IDPlot) + " AND ID3 != " + m_QMark + ID + m_QMark + " AND "
						"ST_Area(geometry) > " + std::to_string(MinEntSize) + " AND "
						"ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + m_LayerName + " WHERE ROWID == " + std::to_string(SmallFIDList[i])+ "))) > 0 "
						"ORDER BY ST_Length(ST_Intersection(geometry, (SELECT geometry FROM " + m_LayerName + " WHERE ROWID == " + std::to_string(SmallFIDList[i])+ "))) DESC, ST_Area(geometry) DESC";
				OGRLayer *SQLLayer = DataSource->ExecuteSQL(m_SQLRequest.c_str(), nullptr, m_SQLDialect.c_str());
				if (SQLLayer->GetFeatureCount() != 0)
				{
					OGRFeature *SQLFeature = SQLLayer->GetFeature(0);
					unsigned int FID = SQLFeature->GetFieldAsInteger("RID");
					OldIDList.push_back(mp_Feature->GetFieldAsString("ID"));
					NewIDList.push_back(SQLFeature->GetFieldAsString("ID"));
					NewID2List.push_back(SQLFeature->GetFieldAsString("ID2"));
					NewPlotIDList.push_back(SQLFeature->GetFieldAsInteger("IDPlot"));
					NewPlotID2List.push_back(SQLFeature->GetFieldAsInteger("IDPlot2"));
					OGRFeature *AggregatingFeature = mp_Layer->GetFeature(FID);
					OGRGeometry *AggregatingGeometry = AggregatingFeature->GetGeometryRef()->clone();
					OGRGeometry *InterGeometry = AggregatingGeometry->Union(mp_Geometry);
					AggregatingFeature->SetGeometry(InterGeometry);
					mp_Layer->SetFeature(AggregatingFeature);
					mp_Layer->SyncToDisk();
					FIDToEraseFromFIDList.push_back(SmallFIDList[i]);
					OGRFeature::DestroyFeature(SQLFeature);
					DataSource->ReleaseResultSet(SQLLayer);
					OGRDataSource::DestroyDataSource(DataSource);
					OGRFeature::DestroyFeature(AggregatingFeature);
					OGRFeature::DestroyFeature(mp_Feature);
					delete(mp_Geometry);
					delete(AggregatingGeometry);
					delete(InterGeometry);
				}
				else
				{
					DataSource->ReleaseResultSet(SQLLayer);
					OGRDataSource::DestroyDataSource(DataSource);
					OGRFeature::DestroyFeature(mp_Feature);
					delete(mp_Geometry);
				}
			}
			for (unsigned int k = 0; k < FIDToEraseFromFIDList.size(); k++)
			{
				if (std::find(SmallFIDList.begin(), SmallFIDList.end(), FIDToEraseFromFIDList[k]) != SmallFIDList.end())
				{
					unsigned int pos = std::find(SmallFIDList.begin(), SmallFIDList.end(), FIDToEraseFromFIDList[k]) - SmallFIDList.begin();
					SmallFIDList.erase(SmallFIDList.begin()+pos);
				}
			}
			FIDToEraseFromFIDList.clear();
			unsigned int End = SmallFIDList.size();
			if (Beginning == End)
			{
				break;
			}
		}

		for (unsigned int i = 0; i < FIDToDeleteList.size(); i++)
		{
			if (std::find(SmallFIDList.begin(), SmallFIDList.end(), FIDToDeleteList[i]) == SmallFIDList.end())
			{
				mp_Layer->DeleteFeature(FIDToDeleteList[i]);
			}
		}

		repack();

		if(!OldIDList.empty())
		{
			updateTextFieldValues("ID2", "ID3", OldIDList, NewID2List);
			updateIntegerFieldValues("ID2", "IDPlot3", OldIDList, NewPlotID2List);
			updateTextFieldValues("ID2", "ID2", OldIDList, NewIDList);

			updateTextFieldValues("ID3", "ID3", OldIDList, NewIDList);
		}

		mp_DataSource->SyncToDisk();
	}
}


//======================================================================
//======================================================================


std::vector<unsigned int> VectorRegrouping::findFeaturesThatShareAttributesAndLimits(unsigned int FID)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> FIDList;

	OGRDataSource *DataSource = OGRSFDriverRegistrar::Open( (m_FullFilePath).c_str(), TRUE );
	OGRLayer *Layer = DataSource->GetLayer(0);

	OGRFeature *FeatureInQuestion = Layer->GetFeature(FID), *Feature = nullptr;

	unsigned int IDPlot = FeatureInQuestion->GetFieldAsInteger("IDPlot");
	std::string ID = FeatureInQuestion->GetFieldAsString("ID");
	std::string ID2 = FeatureInQuestion->GetFieldAsString("ID2");

	if (ID != "0N1")
	{
		while((Feature = Layer->GetNextFeature()) != nullptr)
		{
			if(Feature->GetFieldAsInteger("IDPlot") == IDPlot and (std::string) Feature->GetFieldAsString("ID2") == ID2 and Feature != FeatureInQuestion)
			{
				OGRGeometry *GeometryInQuestion = FeatureInQuestion->GetGeometryRef()->clone(), *Geometry = Feature->GetGeometryRef()->clone(), *InterGeometry = nullptr;
				InterGeometry = GeometryInQuestion->Intersection(Geometry);
				if(InterGeometry->getGeometryType() == wkbLineString or InterGeometry->getGeometryType() == wkbMultiLineString)
				{
					FIDList.push_back(Feature->GetFID());
				}
				delete(GeometryInQuestion);
				delete(Geometry);
				delete(InterGeometry);
			}
			OGRFeature::DestroyFeature(Feature);
		}
	}

	OGRFeature::DestroyFeature(FeatureInQuestion);
	OGRDataSource::DestroyDataSource(DataSource);

	return FIDList;
}


//======================================================================
//======================================================================


void VectorRegrouping::regroupFeaturesThatShareAttributesAndLimits()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	unsigned int N = 1;

	while(N != 0)
	{
		std::vector<std::string> NewIDList, OldIDList;
		std::vector<unsigned int> FIDToRemoveList, NewPlotIDList, ShearingFIDList;

		mp_Layer->ResetReading();

		while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			unsigned int FID = mp_Feature->GetFID();
			unsigned int IDPlot = mp_Feature->GetFieldAsInteger("IDPlot");
			std::string ID = mp_Feature->GetFieldAsString("ID");
			std::string ID2 = mp_Feature->GetFieldAsString("ID2");

			if (ID != "0N1" and std::find(FIDToRemoveList.begin(), FIDToRemoveList.end(), FID) == FIDToRemoveList.end()
					and std::find(NewPlotIDList.begin(), NewPlotIDList.end(), IDPlot) == NewPlotIDList.end())
			{
				ShearingFIDList.clear();
				ShearingFIDList = findFeaturesThatShareAttributesAndLimits(FID);
				if (!ShearingFIDList.empty())
				{
					for (unsigned int i = 0; i < ShearingFIDList.size(); i++)
					{
						OGRFeature *Feature = mp_Layer->GetFeature(ShearingFIDList[i]);
						OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();

						if (std::find(FIDToRemoveList.begin(), FIDToRemoveList.end(), ShearingFIDList[i]) == FIDToRemoveList.end())
						{
							FIDToRemoveList.push_back(Feature->GetFID());
							OldIDList.push_back((std::string) Feature->GetFieldAsString("ID"));
							NewIDList.push_back(ID.c_str());
							NewPlotIDList.push_back(IDPlot);
							mp_Geometry = mp_Feature->GetGeometryRef()->clone();
							OGRGeometry *InterGeometry = nullptr;
							InterGeometry = mp_Geometry->Union(Geometry);
							mp_Feature->SetGeometry(InterGeometry);
							mp_Layer->SetFeature(mp_Feature);
							OGRFeature::DestroyFeature(Feature);
							delete(Geometry);
							delete(InterGeometry);
							delete(mp_Geometry);
						}
						else
						{
							OGRFeature::DestroyFeature(Feature);
							delete(Geometry);
						}
					}
				}
				OGRFeature::DestroyFeature(mp_Feature);
			}
			else
			{
				OGRFeature::DestroyFeature(mp_Feature);
			}
		}

		deleteFeatures(FIDToRemoveList);

		repack();

		mp_Layer->SyncToDisk();

		if(!OldIDList.empty())
		{
			updateIntegerFieldValues("ID2", "IDPlot2", OldIDList, NewPlotIDList);
			updateTextFieldValues("ID2", "ID2", OldIDList, NewIDList);

			updateIntegerFieldValues("ID3", "IDPlot3", OldIDList, NewPlotIDList);
			updateTextFieldValues("ID3", "ID3", OldIDList, NewIDList);
		}

		N = 0;

		mp_Layer->ResetReading();

		while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			if(!findFeaturesThatShareAttributesAndLimits(mp_Feature->GetFID()).empty())
			{
				N++;
			}
			OGRFeature::DestroyFeature(mp_Feature);
		}
	}
}


//======================================================================
//======================================================================


int VectorRegrouping::findAggregatingFeature(unsigned int FID)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	int AggregatingFID = -1;
	unsigned int NOfNeighborFeaturesWithinSamePlot = 0;

	OGRDataSource *DataSource = OGRSFDriverRegistrar::Open( (m_FullFilePath).c_str(), TRUE );
	OGRLayer *Layer = DataSource->GetLayer(0);
	OGRFeature *Feature = nullptr;

	std::vector<std::string> ID2List;

	Layer->ResetReading();

	while((Feature = Layer->GetNextFeature()) != nullptr)
	{
		if(std::find(ID2List.begin(), ID2List.end(), (std::string)(Feature->GetFieldAsString("ID2"))) == ID2List.end())
		{
			ID2List.push_back((std::string)(Feature->GetFieldAsString("ID2")));
		}
	    OGRFeature::DestroyFeature(Feature);
	}

	Layer->ResetReading();

	OGRFeature *FeatureInQuestion = Layer->GetFeature(FID);

	if(std::find(ID2List.begin(), ID2List.end(), (std::string)(FeatureInQuestion->GetFieldAsString("ID"))) == ID2List.end())
	{
		unsigned int IDPlot = FeatureInQuestion->GetFieldAsInteger("IDPlot");
		double PlotArea = 0;

		Layer->ResetReading();

		while((Feature = Layer->GetNextFeature()) != nullptr)
		{
			if(Feature->GetFieldAsInteger("IDPlot") == IDPlot)
			{
				OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();
				if(Geometry->getGeometryType() == wkbPolygon)
				{
					PlotArea+=((OGRPolygon *) Geometry)->get_Area();
				}
				else if(Geometry->getGeometryType() == wkbMultiPolygon)
				{
					for(unsigned int i = 0; i < ((OGRMultiPolygon *) Geometry)->getNumGeometries(); i++)
					{
						PlotArea+=((OGRPolygon *) (((OGRMultiPolygon *) Geometry)->getGeometryRef(i)))->get_Area();
					}
				}
				delete(Geometry);
			}
			OGRFeature::DestroyFeature(Feature);
		}

		Layer->ResetReading();

		while((Feature = Layer->GetNextFeature()) != nullptr)
		{
			if(Feature->GetFID() != FID)
			{
				OGRGeometry *GeometryInQuestion = FeatureInQuestion->GetGeometryRef()->clone(), *Geometry = Feature->GetGeometryRef()->clone();
				OGRGeometry *IntersectionGeometry = GeometryInQuestion->Intersection(Geometry);
				double GeometryInQuestionArea = 0;

				if(GeometryInQuestion->getGeometryType() == wkbPolygon)
				{
					GeometryInQuestionArea+=((OGRPolygon *) GeometryInQuestion)->get_Area();
				}
				else if(GeometryInQuestion->getGeometryType() == wkbMultiPolygon)
				{
					for(unsigned int i = 0; i < ((OGRMultiPolygon *) GeometryInQuestion)->getNumGeometries(); i++)
					{
						GeometryInQuestionArea+=((OGRPolygon *) (((OGRMultiPolygon *) GeometryInQuestion)->getGeometryRef(i)))->get_Area();
					}
				}

				if(GeometryInQuestion->Touches(Geometry) and Feature->GetFieldAsInteger("IDPlot") == IDPlot and (IntersectionGeometry->getGeometryType() == wkbLineString or
						IntersectionGeometry->getGeometryType() == wkbMultiLineString) and GeometryInQuestionArea < (PlotArea/6))
				{
					NOfNeighborFeaturesWithinSamePlot++;
					AggregatingFID = Feature->GetFID();
				}
				delete(Geometry);
				delete(GeometryInQuestion);
				delete(IntersectionGeometry);
			}
			OGRFeature::DestroyFeature(Feature);
		}
	}

	OGRFeature::DestroyFeature(FeatureInQuestion);

	OGRDataSource::DestroyDataSource(DataSource);

	if(NOfNeighborFeaturesWithinSamePlot == 1 and AggregatingFID != -1)
	{
		return AggregatingFID;
	}
	else
	{
		return -1;
	}
}


//======================================================================
//======================================================================


void VectorRegrouping::regroupLockedFeatures()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	unsigned int N = 1;

	while(N != 0)
	{
		std::vector<unsigned int> FIDToRemoveList;

		mp_Layer->ResetReading();

		while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			unsigned int FID = mp_Feature->GetFID();

			if ((std::string) (mp_Feature->GetFieldAsString("ID")) != "0N1" and (std::string) (mp_Feature->GetFieldAsString("ID2")) != "0N1" and
					std::find(FIDToRemoveList.begin(), FIDToRemoveList.end(), FID) == FIDToRemoveList.end())
			{
				int AggregatingFID = findAggregatingFeature(FID);

				if (AggregatingFID != -1 and std::find(FIDToRemoveList.begin(), FIDToRemoveList.end(), AggregatingFID) == FIDToRemoveList.end())
				{
					OGRFeature *Feature = mp_Layer->GetFeature(AggregatingFID);
					mp_Geometry = mp_Feature->GetGeometryRef()->clone();
					OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();
					OGRGeometry *InterGeometry = Geometry->Union(mp_Geometry);
					Feature->SetGeometry(InterGeometry);
					mp_Layer->SetFeature(Feature);
					FIDToRemoveList.push_back(FID);
					OGRFeature::DestroyFeature(Feature);
					delete(mp_Geometry);
					delete(Geometry);
					delete(InterGeometry);
				}
				OGRFeature::DestroyFeature(mp_Feature);
			}
			else
			{
				OGRFeature::DestroyFeature(mp_Feature);
			}
		}

		deleteFeatures(FIDToRemoveList);

		repack();

		N = 0;

		mp_Layer->ResetReading();

		while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			if((std::string) (mp_Feature->GetFieldAsString("ID")) != "0N1" and (std::string) (mp_Feature->GetFieldAsString("ID2")) != "0N1")
			{
				if(findAggregatingFeature(mp_Feature->GetFID()) != -1)
				{
					N++;
				}
			}
			OGRFeature::DestroyFeature(mp_Feature);
		}
	}
}


//======================================================================
//======================================================================


void VectorRegrouping::createOuterPolygon()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_Layer->ResetReading();

	OGRFeature *NewFeature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());

	for (unsigned int i = 0; i < NewFeature->GetFieldCount(); i++)
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

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		unsigned int IDField = mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str());
		if (!IDField)
		{
			mp_Geometry = mp_Feature->GetGeometryRef()->clone();
			NewFeature->SetGeometry(mp_Geometry);
			mp_Layer->DeleteFeature(mp_Feature->GetFID());
			mp_Layer->CreateFeature(NewFeature);
			OGRFeature::DestroyFeature(mp_Feature);
			delete(mp_Geometry);
			break;
		}
		else
		{
			OGRFeature::DestroyFeature(mp_Feature);
		}
	}

	OGRFeature::DestroyFeature(NewFeature);

	repack();

	mp_Layer->SyncToDisk();
	mp_DataSource->SyncToDisk();

	NewFeature = mp_Layer->GetFeature(mp_Layer->GetFeatureCount()-1);

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		unsigned int IDField = mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str());
		if (!IDField and mp_Feature->GetFID() != mp_Layer->GetFeatureCount()-1)
		{
			OGRGeometry *NewGeometry = NewFeature->GetGeometryRef()->clone();
			OGRGeometry *UnionGeometry = mp_Feature->GetGeometryRef()->clone();
			OGRGeometry *TempGeometry = NewGeometry->Union(UnionGeometry);
			NewFeature->SetGeometry(TempGeometry);
			mp_Layer->DeleteFeature(mp_Feature->GetFID());
			mp_Layer->SetFeature(NewFeature);
			OGRFeature::DestroyFeature(mp_Feature);
			delete(NewGeometry);
			delete(UnionGeometry);
			delete(TempGeometry);
		}
		else
		{
			OGRFeature::DestroyFeature(mp_Feature);
		}
	}

	OGRFeature::DestroyFeature(NewFeature);

	repack();
}


//======================================================================
//======================================================================


std::vector<unsigned int> VectorRegrouping::getPlotIDList()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> PlotIDList;

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		if (std::find(PlotIDList.begin(), PlotIDList.end(), mp_Feature->GetFieldAsInteger("IDPlot")) == PlotIDList.end())
		{
			PlotIDList.push_back(mp_Feature->GetFieldAsInteger("IDPlot"));
			OGRFeature::DestroyFeature(mp_Feature);
	    }
		else
		{
			OGRFeature::DestroyFeature(mp_Feature);
		}
	}

	return PlotIDList;
}


//======================================================================
//======================================================================


std::vector<unsigned int> VectorRegrouping::getMissingPlotIDList()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> MissingPlotIDList, PlotIDList;

	PlotIDList = getPlotIDList();

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		if (std::find(PlotIDList.begin(), PlotIDList.end(), mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str())) == PlotIDList.end() and
	        std::find(MissingPlotIDList.begin(), MissingPlotIDList.end(), mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str())) == MissingPlotIDList.end())
	    {
			MissingPlotIDList.push_back(mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()));
			OGRFeature::DestroyFeature(mp_Feature);
	    }
		else
		{
			OGRFeature::DestroyFeature(mp_Feature);
		}
	}

	return MissingPlotIDList;
}


//======================================================================
//======================================================================


std::vector<unsigned int> VectorRegrouping::getSliverFIDList()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> SliverFIDList, MissingPlotIDList;

	MissingPlotIDList = getMissingPlotIDList();

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		if (mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) != mp_Feature->GetFieldAsInteger("IDPlot") and
	        std::find(MissingPlotIDList.begin(), MissingPlotIDList.end(), mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str())) == MissingPlotIDList.end())
	    {
			SliverFIDList.push_back(mp_Feature->GetFID());
			OGRFeature::DestroyFeature(mp_Feature);
	    }
		else
		{
			OGRFeature::DestroyFeature(mp_Feature);
		}
	}

	return SliverFIDList;
}


//======================================================================
//======================================================================


void VectorRegrouping::regroupSliverFeaturesWithSingleNonSliverNeighbor()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> FIDToRemoveList, SliverFIDList = getSliverFIDList();

	for (unsigned int i = 0; i < SliverFIDList.size(); i++)
	{
    	OGRFeature *SliverFeature = mp_Layer->GetFeature(SliverFIDList[i]);
    	OGRGeometry *SliverGeometry = SliverFeature->GetGeometryRef()->clone();

    	int SliverFID = SliverFeature->GetFID(), UnionFID, N = 0;

	    mp_Layer->ResetReading();
    	mp_Feature = nullptr;

	    while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	    {
	    	mp_Geometry = mp_Feature->GetGeometryRef()->clone();

	    	OGRGeometry *IntersectionGeometry = mp_Geometry->Intersection(SliverGeometry);

	    	if (mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) == mp_Feature->GetFieldAsInteger("IDPlot") and
	    			mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) == SliverFeature->GetFieldAsInteger(m_IDFieldName.c_str()) and
					((IntersectionGeometry)->getGeometryType() == 2 or
							(IntersectionGeometry)->getGeometryType() == 5))
	    	{
	    		N++;
	    		UnionFID = mp_Feature->GetFID();
	    		OGRFeature::DestroyFeature(mp_Feature);
	    		delete(mp_Geometry);
	    		delete(IntersectionGeometry);
	    	}
	    	else
	    	{
	    		OGRFeature::DestroyFeature(mp_Feature);
	    		delete(mp_Geometry);
	    		delete(IntersectionGeometry);
	    	}
	    }
	    if(N == 1)
	    {
	    	mp_Layer->ResetReading();

	    	mp_Feature = mp_Layer->GetFeature(UnionFID);

	    	mp_Geometry = mp_Feature->GetGeometryRef()->clone();

	    	OGRGeometry *InterGeometry = nullptr;

	    	InterGeometry = mp_Geometry->Union(SliverGeometry);
	    	mp_Feature->SetGeometry(InterGeometry);
	    	mp_Layer->SetFeature(mp_Feature);
	    	FIDToRemoveList.push_back(SliverFIDList[i]);
	    	OGRFeature::DestroyFeature(mp_Feature);
	    	OGRFeature::DestroyFeature(SliverFeature);

	    	delete(mp_Geometry);
	    	delete(SliverGeometry);
	    	delete(InterGeometry);

	    }
	    else
	    {
	    	OGRFeature::DestroyFeature(mp_Feature);
	    	OGRFeature::DestroyFeature(SliverFeature);
	    	delete(SliverGeometry);
	    }
	}

    mp_Layer->ResetReading();

	deleteFeatures(FIDToRemoveList);

    repack();
}


//======================================================================
//======================================================================


void VectorRegrouping::regroupSliverFeaturesBasedOnGeometryAndAttributes()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> FIDToRemoveList, SliverFIDList = getSliverFIDList();

	for (unsigned int i = 0; i < SliverFIDList.size(); i++)
	{
		OGRFeature *SliverFeature = mp_Layer->GetFeature(SliverFIDList[i]);
		OGRGeometry *SliverGeometry = SliverFeature->GetGeometryRef()->clone();

		OGRGeometry *TempGeometry = nullptr;

		int SliverFID = SliverFeature->GetFID();
		std::string SliverID = SliverFeature->GetFieldAsString("ID");
		std::string SliverID2 = SliverFeature->GetFieldAsString("ID2");

		int UnionFIDN1 = -1, UnionFIDN2 = -1, UnionFIDN3 = -1;

		mp_Layer->ResetReading();

		while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			std::string UnionID2 = mp_Feature->GetFieldAsString("ID2");
			mp_Geometry = mp_Feature->GetGeometryRef()->clone();

			OGRGeometry *IntersectionGeometry = mp_Geometry->Intersection(SliverGeometry);

			if (mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) == mp_Feature->GetFieldAsInteger("IDPlot") and
					mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) == SliverFeature->GetFieldAsInteger(m_IDFieldName.c_str()) and
					(IntersectionGeometry->getGeometryType() == 2 or
							IntersectionGeometry->getGeometryType() == 5) and
							SliverID == UnionID2)
			{
				UnionFIDN1 = mp_Feature->GetFID();
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
				delete(IntersectionGeometry);
				break;
			}
			else
			{
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
				delete(IntersectionGeometry);
			}
		}

		mp_Layer->ResetReading();

		while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			std::string UnionID = mp_Feature->GetFieldAsString("ID");
			mp_Geometry = mp_Feature->GetGeometryRef()->clone();

			OGRGeometry *IntersectionGeometry = mp_Geometry->Intersection(SliverGeometry);

			if (mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) == mp_Feature->GetFieldAsInteger("IDPlot") and
					mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) == SliverFeature->GetFieldAsInteger(m_IDFieldName.c_str()) and
					(IntersectionGeometry->getGeometryType() == 2 or
							IntersectionGeometry->getGeometryType() == 5) and
							SliverID2 == UnionID)
			{
				UnionFIDN2 = mp_Feature->GetFID();
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
				delete(IntersectionGeometry);
				break;
			}
			else
			{
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
				delete(IntersectionGeometry);
			}
		}

		mp_Layer->ResetReading();

		while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			std::string UnionID2 = mp_Feature->GetFieldAsString("ID2");
			mp_Geometry = mp_Feature->GetGeometryRef()->clone();

			OGRGeometry *IntersectionGeometry = mp_Geometry->Intersection(SliverGeometry);

			if (mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) == mp_Feature->GetFieldAsInteger("IDPlot") and
					mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) == SliverFeature->GetFieldAsInteger(m_IDFieldName.c_str()) and
					(IntersectionGeometry->getGeometryType() == 2 or
							IntersectionGeometry->getGeometryType() == 5) and
							SliverID2 == UnionID2)
			{
				UnionFIDN3 = mp_Feature->GetFID();
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
				delete(IntersectionGeometry);
				break;
			}
			else
			{
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
				delete(IntersectionGeometry);
			}
		}

		mp_Layer->ResetReading();

		if(UnionFIDN1 != -1)
		{
			mp_Feature = mp_Layer->GetFeature(UnionFIDN1);
			mp_Geometry = mp_Feature->GetGeometryRef()->clone();
			TempGeometry = mp_Geometry->Union(SliverGeometry);
			mp_Feature->SetGeometry(TempGeometry);
			mp_Layer->SetFeature(mp_Feature);
			FIDToRemoveList.push_back(SliverFIDList[i]);
			OGRFeature::DestroyFeature(mp_Feature);
			OGRFeature::DestroyFeature(SliverFeature);
			delete(mp_Geometry);
			delete(SliverGeometry);
			delete(TempGeometry);
		}
		else
		{
			if(UnionFIDN2 != -1)
			{
				mp_Feature = mp_Layer->GetFeature(UnionFIDN2);
				mp_Geometry = mp_Feature->GetGeometryRef()->clone();
				TempGeometry = mp_Geometry->Union(SliverGeometry);
				mp_Feature->SetGeometry(TempGeometry);
				mp_Layer->SetFeature(mp_Feature);
				FIDToRemoveList.push_back(SliverFIDList[i]);
				OGRFeature::DestroyFeature(mp_Feature);
				OGRFeature::DestroyFeature(SliverFeature);
				delete(mp_Geometry);
				delete(SliverGeometry);
				delete(TempGeometry);
			}
			else
			{
				if(UnionFIDN3 != -1)
				{
					mp_Feature = mp_Layer->GetFeature(UnionFIDN3);
					mp_Geometry = mp_Feature->GetGeometryRef()->clone();
					TempGeometry = mp_Geometry->Union(SliverGeometry);
					mp_Feature->SetGeometry(TempGeometry);
					mp_Layer->SetFeature(mp_Feature);
					FIDToRemoveList.push_back(SliverFIDList[i]);
					OGRFeature::DestroyFeature(mp_Feature);
					OGRFeature::DestroyFeature(SliverFeature);
					delete(mp_Geometry);
					delete(SliverGeometry);
					delete(TempGeometry);
				}
				else
				{
					OGRFeature::DestroyFeature(SliverFeature);
					delete(SliverGeometry);
					delete(TempGeometry);
				}
			}
		}
	}

	mp_Layer->ResetReading();

	deleteFeatures(FIDToRemoveList);

	repack();

}


//======================================================================
//======================================================================


void VectorRegrouping::regroupSliverFeaturesBasedOnGeometry()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> FIDToRemoveList, SliverFIDList = getSliverFIDList();

    for (unsigned int i = 0; i < SliverFIDList.size(); i++)
    {
    	OGRFeature *SliverFeature = mp_Layer->GetFeature(SliverFIDList[i]);
    	OGRGeometry *SliverGeometry = SliverFeature->GetGeometryRef()->clone(), *TempGeometry = nullptr;;

    	int SliverFID = SliverFeature->GetFID(), UnionFID = -1, N = 0;
    	double Length = 0;

    	mp_Layer->ResetReading();

    	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
    	{
    		mp_Geometry = mp_Feature->GetGeometryRef()->clone();
    		OGRGeometry *IntersectionGeometry = mp_Geometry->Intersection(SliverGeometry);

    		if (mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) == mp_Feature->GetFieldAsInteger("IDPlot") and
    				mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) == SliverFeature->GetFieldAsInteger(m_IDFieldName.c_str()) and
					(IntersectionGeometry->getGeometryType() == 2 or
							IntersectionGeometry->getGeometryType() == 5))
    		{
    			double IntersectionLength = 0;
    			if(IntersectionGeometry->getGeometryType() == 2)
    			{
    				OGRLineString *LineString = ((OGRLineString *) IntersectionGeometry);
    				IntersectionLength = LineString->get_Length();
    			}
    			else if(IntersectionGeometry->getGeometryType() == 5)
    			{
    				OGRMultiLineString *MultiLineString = ((OGRMultiLineString *) IntersectionGeometry);
    				for(unsigned int i = 0; i < MultiLineString->getNumGeometries(); i++)
    				{
    					IntersectionLength += ((OGRLineString *) MultiLineString->getGeometryRef(i))->get_Length();
    				}
    			}
    			if(IntersectionLength > Length)
    			{
    				Length = IntersectionLength;
    				UnionFID = mp_Feature->GetFID();
    				OGRFeature::DestroyFeature(mp_Feature);
    				delete(mp_Geometry);
    				delete(IntersectionGeometry);
    			}
    			else
    			{
    				OGRFeature::DestroyFeature(mp_Feature);
    				delete(mp_Geometry);
    				delete(IntersectionGeometry);
    			}
    		}
    		else
    		{
    			OGRFeature::DestroyFeature(mp_Feature);
    			delete(mp_Geometry);
    			delete(IntersectionGeometry);
    		}
    	}

    	if(UnionFID != -1)
    	{
    		mp_Feature = mp_Layer->GetFeature(UnionFID);
    		mp_Geometry = mp_Feature->GetGeometryRef()->clone();
    		TempGeometry = mp_Geometry->Union(SliverGeometry);
    		mp_Feature->SetGeometry(TempGeometry);
    		mp_Layer->SetFeature(mp_Feature);
    		FIDToRemoveList.push_back(SliverFIDList[i]);
    		OGRFeature::DestroyFeature(mp_Feature);
    		OGRFeature::DestroyFeature(SliverFeature);
    		delete(mp_Geometry);
    		delete(SliverGeometry);
    		delete(TempGeometry);
    	}
    	else
    	{
    		OGRFeature::DestroyFeature(SliverFeature);
			delete(SliverGeometry);
			delete(TempGeometry);
    	}
    }

	mp_Layer->ResetReading();

	deleteFeatures(FIDToRemoveList);

	repack();
}


//======================================================================
//======================================================================


void VectorRegrouping::regroupSliverFeaturesBasedOnAttributes()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> FIDToRemoveList, SliverFIDList = getSliverFIDList();

	for (unsigned int i = 0; i < SliverFIDList.size(); i++)
	{
    	OGRFeature *SliverFeature = mp_Layer->GetFeature(SliverFIDList[i]);
    	OGRGeometry *SliverGeometry = SliverFeature->GetGeometryRef()->clone(), *TempGeometry = nullptr;;

		int SliverFID = SliverFeature->GetFID(), UnionFID = -1, N = 0;

		mp_Layer->ResetReading();

		while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			mp_Geometry = mp_Feature->GetGeometryRef()->clone();

			if (mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) == mp_Feature->GetFieldAsInteger("IDPlot") and
					mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()) == SliverFeature->GetFieldAsInteger(m_IDFieldName.c_str()) and
					mp_Geometry->Touches(SliverGeometry))
			{
				UnionFID = mp_Feature->GetFID();
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
				break;
			}
			else
			{
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
			}
		}

		if(UnionFID != -1)
		{
			mp_Feature = mp_Layer->GetFeature(UnionFID);
			mp_Geometry = mp_Feature->GetGeometryRef()->clone();
			TempGeometry = mp_Geometry->Union(SliverGeometry);
			mp_Feature->SetGeometry(TempGeometry);
			mp_Layer->SetFeature(mp_Feature);
			FIDToRemoveList.push_back(SliverFIDList[i]);
			OGRFeature::DestroyFeature(mp_Feature);
			OGRFeature::DestroyFeature(SliverFeature);
			delete(mp_Geometry);
			delete(SliverGeometry);
			delete(TempGeometry);
		}
		else
		{
			OGRFeature::DestroyFeature(SliverFeature);
			delete(SliverGeometry);
			delete(TempGeometry);
		}
	}

	mp_Layer->ResetReading();

	if (FIDToRemoveList.size() != 0)
	{
		for (unsigned int i = 0; i < FIDToRemoveList.size(); i++)
		{
			mp_Layer->DeleteFeature(FIDToRemoveList[i]);
		}
	}

	repack();
}


//======================================================================
//======================================================================


void VectorRegrouping::regroupFeaturesWithMissingPlotID()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> MissingPlotIDList = getMissingPlotIDList();

	for (unsigned int i = 0; i < MissingPlotIDList.size(); i++)
	{
		OGRFeature *NewFeature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());

		unsigned int IDPlot2, IDPlot3;
		double Area = 0;
		std::string IDNew = std::to_string(MissingPlotIDList[i]) + "N1", ID2, ID3;

		mp_Layer->ResetReading();

		while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			unsigned int IDField = mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str());
			mp_Geometry = mp_Feature->GetGeometryRef()->clone();

			if(IDField == MissingPlotIDList[i] and ((OGRPolygon *) mp_Geometry)->get_Area() > Area)
			{
				IDPlot2 = mp_Feature->GetFieldAsInteger("IDPlot2");
				ID2 = mp_Feature->GetFieldAsString("ID2");
				IDPlot3 = mp_Feature->GetFieldAsInteger("IDPlot3");
				ID3 = mp_Feature->GetFieldAsString("ID3");
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
			}
			else
			{
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
			}
		}

		NewFeature->SetField(m_IDFieldName.c_str(), (int) MissingPlotIDList[i]);
		NewFeature->SetField("IDPlot", (int) MissingPlotIDList[i]);
		NewFeature->SetField("ID", IDNew.c_str());
		NewFeature->SetField("IDPlot2", (int) IDPlot2);
		NewFeature->SetField("ID2", ID2.c_str());
		NewFeature->SetField("IDPlot3", (int) IDPlot3);
		NewFeature->SetField("ID3", ID3.c_str());

		mp_Layer->ResetReading();

		while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			unsigned int IDField = mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str());
			if (IDField == MissingPlotIDList[i])
			{
				mp_Geometry = mp_Feature->GetGeometryRef()->clone();
				NewFeature->SetGeometry(mp_Geometry);
				mp_Layer->DeleteFeature(mp_Feature->GetFID());
				mp_Layer->CreateFeature(NewFeature);
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
				break;
			}
			else
			{
				OGRFeature::DestroyFeature(mp_Feature);
			}
		}

		OGRFeature::DestroyFeature(NewFeature);

		repack();

		NewFeature = mp_Layer->GetFeature(mp_Layer->GetFeatureCount()-1);

		mp_Layer->ResetReading();

		while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			unsigned int IDField = mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str());
			if (IDField == MissingPlotIDList[i] and mp_Feature->GetFID() != mp_Layer->GetFeatureCount()-1)
			{
				OGRGeometry *NewGeometry = NewFeature->GetGeometryRef()->clone();
				mp_Geometry = mp_Feature->GetGeometryRef()->clone();
				OGRGeometry *TempGeometry = NewGeometry->Union(mp_Geometry);
				NewFeature->SetGeometry(TempGeometry);
				mp_Layer->DeleteFeature(mp_Feature->GetFID());
				mp_Layer->SetFeature(NewFeature);
				OGRFeature::DestroyFeature(mp_Feature);
				delete(NewGeometry);
				delete(mp_Geometry);
				delete(TempGeometry);
			}
			else
			{
				OGRFeature::DestroyFeature(mp_Feature);
			}
		}

		OGRFeature::DestroyFeature(NewFeature);

		repack();
	}
}


//======================================================================
//======================================================================


void VectorRegrouping::correctOuterPolygon()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector <unsigned int> FIDToRemoveList;

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		int FID = mp_Feature->GetFID();
		std::string ID = mp_Feature->GetFieldAsString("ID");
	    mp_Geometry = mp_Feature->GetGeometryRef();
	    if (ID == "0N1" and mp_Geometry->getGeometryType() == 6)
	    {
	    	OGRGeometry *NewGeometry = mp_Geometry->clone();
	    	for (unsigned int i = 0; i < ((OGRMultiPolygon *) NewGeometry)->getNumGeometries(); i++)
	    	{
	    		OGRFeature *NewFeature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());
	    		NewFeature->SetFrom(mp_Feature);
	    		NewFeature->SetGeometry((OGRGeometry *)(((OGRMultiPolygon*) NewGeometry)->getGeometryRef(i)));
	    		mp_Layer->CreateFeature(NewFeature);
	    		OGRFeature::DestroyFeature(NewFeature);
	    	}
	    	FIDToRemoveList.push_back(FID);
	    	delete(NewGeometry);
	    }
	    OGRFeature::DestroyFeature(mp_Feature);
	}

	mp_Layer->ResetReading();

	if (FIDToRemoveList.size() != 0)
	{
		for (unsigned int i = 0; i < FIDToRemoveList.size(); i++)
	    {
			mp_Layer->DeleteFeature(FIDToRemoveList[i]);
	    }
	}

	repack();

	double Area = 0;
	int ZeroFID = -1;

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		mp_Geometry = mp_Feature->GetGeometryRef()->clone();
	    std::string ID = mp_Feature->GetFieldAsString("ID");

	    if (ID == "0N1" and ((OGRPolygon *) mp_Geometry)->get_Area() > Area)
	    {
	    	Area = ((OGRPolygon *) mp_Geometry)->get_Area();
	    	ZeroFID = mp_Feature->GetFID();
	    }

	    OGRFeature::DestroyFeature(mp_Feature);
	    delete(mp_Geometry);
	}

	if (ZeroFID != -1)
	{
		mp_Layer->ResetReading();

		while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			mp_Geometry = mp_Feature->GetGeometryRef()->clone();
			std::string ID = mp_Feature->GetFieldAsString("ID");
			if (ID == "0N1" and mp_Feature->GetFID() != ZeroFID)
			{
				OGRDataSource *DataSource = OGRSFDriverRegistrar::Open((m_FullFilePath).c_str(), TRUE );
				OGRLayer *Layer = DataSource->GetLayer(0);
				OGRFeature *Feature = nullptr;
				Layer->ResetReading();

				double Length = 0;
				int FID = -1;

				while ((Feature = Layer->GetNextFeature()) != nullptr)
				{
					OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();
					OGRGeometry *IntersectionGeometry = mp_Geometry->Intersection(Geometry);
					if(IntersectionGeometry->getGeometryType() == 2 or IntersectionGeometry->getGeometryType() == 5)
					{
						double IntersectionLength = 0;
						if(IntersectionGeometry->getGeometryType() == 2)
						{
							OGRLineString *LineString = ((OGRLineString *) IntersectionGeometry);
							IntersectionLength = LineString->get_Length();
						}
						else if(IntersectionGeometry->getGeometryType() == 5)
						{
							OGRMultiLineString *MultiLineString = ((OGRMultiLineString *) IntersectionGeometry);
							for(unsigned int i = 0; i < MultiLineString->getNumGeometries(); i++)
							{
								IntersectionLength += ((OGRLineString *) MultiLineString->getGeometryRef(i))->get_Length();
							}
						}
						if(IntersectionLength > Length)
						{
							Length = IntersectionLength;
							FID = Feature->GetFID();
						}
					}
					OGRFeature::DestroyFeature(Feature);
					delete(Geometry);
					delete(IntersectionGeometry);
				}

				Layer->ResetReading();

				if(FID != -1)
				{
					Feature = Layer->GetFeature(FID);
					mp_Feature->SetField(m_IDFieldName.c_str(), Feature->GetFieldAsInteger(m_IDFieldName.c_str()));
					mp_Feature->SetField("IDPlot", Feature->GetFieldAsInteger("IDPlot"));
					mp_Feature->SetField("ID", Feature->GetFieldAsString("ID"));
					mp_Feature->SetField("IDPlot2", Feature->GetFieldAsInteger("IDPlot2"));
					mp_Feature->SetField("ID2", Feature->GetFieldAsString("ID2"));
					mp_Feature->SetField("IDPlot3", Feature->GetFieldAsInteger("IDPlot3"));
					mp_Feature->SetField("ID3", Feature->GetFieldAsString("ID3"));
					mp_Layer->SetFeature(mp_Feature);
					OGRFeature::DestroyFeature(Feature);
				}

				OGRDataSource::DestroyDataSource(DataSource);
			}
			OGRFeature::DestroyFeature(mp_Feature);
			delete(mp_Geometry);
		}
	}
}


//======================================================================
//======================================================================


void VectorRegrouping::correctMultiparts()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> FIDToRemoveList;
	std::vector< std::vector<int>> PlotIDMaxList(2,std::vector<int >(0));

	mp_Layer->ResetReading();

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		unsigned int IDPlot = mp_Feature->GetFieldAsInteger("IDPlot");
	    std::string ID = mp_Feature->GetFieldAsString("ID");
	    unsigned int IDN = std::stoi(ID.substr(ID.find("N")+1, ID.length()-1));
	    if (std::find(PlotIDMaxList[0].begin(), PlotIDMaxList[0].end(), IDPlot) == PlotIDMaxList[0].end())
	    {
	    	PlotIDMaxList[0].push_back(IDPlot);
	    	PlotIDMaxList[1].push_back(IDN);
	    }
	    else
	    {
	    	unsigned int pos = std::find(PlotIDMaxList[0].begin(), PlotIDMaxList[0].end(), IDPlot) - PlotIDMaxList[0].begin();
	    	if (PlotIDMaxList[1][pos] < IDN)
	    	{
	    		PlotIDMaxList[1][pos] = IDN;
	    	}
	    }
	    OGRFeature::DestroyFeature(mp_Feature);
	}

	mp_Layer->ResetReading();

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		mp_Geometry = mp_Feature->GetGeometryRef()->clone();
		OGRFeature *NewFeature = OGRFeature::CreateFeature(mp_Layer->GetLayerDefn());
	    if (mp_Geometry->getGeometryType() == 6)
	    {
	    	std::vector <unsigned int> NumberOfGeometriesList;
	    	for (unsigned int i = 0; i < ((OGRMultiPolygon*) mp_Geometry)->getNumGeometries(); i++)
	    	{
	    		if(std::find(NumberOfGeometriesList.begin(), NumberOfGeometriesList.end(), i) == NumberOfGeometriesList.end())
	    		{
	    			NumberOfGeometriesList.push_back(i);
	    			NewFeature->SetGeometry(((OGRGeometry *) ((OGRMultiPolygon*) mp_Geometry)->getGeometryRef(i)));
	    			for(unsigned int j = 0; j < NumberOfGeometriesList.size(); j++)
	    			{
	    				for(unsigned int k = 0; k < ((OGRMultiPolygon*) mp_Geometry)->getNumGeometries(); k++)
	    				{
	    					OGRGeometry *NewGeometry = NewFeature->GetGeometryRef()->clone();
	    					OGRGeometry *Geometry = ((OGRGeometry *) ((OGRMultiPolygon*) mp_Geometry)->getGeometryRef(k))->clone();
	    					if(k != NumberOfGeometriesList[j] and Geometry->Distance(NewGeometry) < 1 and std::find(NumberOfGeometriesList.begin(), NumberOfGeometriesList.end(), k) == NumberOfGeometriesList.end())
	    					{
	    						OGRGeometry *MultipolygonGeometry = OGRGeometryFactory::createGeometry(wkbMultiPolygon);
	    						OGRPolygon *GeometryPolygon = ((OGRPolygon *) Geometry);
	   							((OGRMultiPolygon *) MultipolygonGeometry)->addGeometry(GeometryPolygon);
	    						if(NewGeometry->getGeometryType() == 3)
	    						{
	    							OGRPolygon *NewGeometryPolygon = ((OGRPolygon *) NewGeometry);
	    							((OGRMultiPolygon *) MultipolygonGeometry)->addGeometry(NewGeometryPolygon);
	    						}
	    						else
	    						{
	    							OGRMultiPolygon *NewGeometryMultipolygon = ((OGRMultiPolygon *) NewGeometry);
	    							for(unsigned int l = 0; l < NewGeometryMultipolygon->getNumGeometries(); l++)
	    							{
	    								OGRPolygon *NewGeometryPoly = ((OGRPolygon *) NewGeometryMultipolygon->getGeometryRef(l));
	    								((OGRMultiPolygon *) MultipolygonGeometry)->addGeometry(NewGeometryPoly);
	    							}
	    						}
	    						NewFeature->SetGeometry(MultipolygonGeometry);
	    						NumberOfGeometriesList.push_back(k);
	    						delete(MultipolygonGeometry);
	    					}
	    					delete(NewGeometry);
	    					delete(Geometry);
	    				}
	    			}
	    			if(((OGRMultiPolygon*) NewFeature->GetGeometryRef())->getNumGeometries() == ((OGRMultiPolygon*) mp_Geometry)->getNumGeometries())
	    			{
	    				NumberOfGeometriesList.clear();
	    				break;
	    			}
	    			else
	    			{
	    				if(i == 0)
	    				{
	    					NewFeature->SetField(m_IDFieldName.c_str(), mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()));
	    					NewFeature->SetField("IDPlot", mp_Feature->GetFieldAsInteger("IDPlot"));
	    					NewFeature->SetField("ID",  mp_Feature->GetFieldAsString("ID"));
	    					NewFeature->SetField("IDPlot2", mp_Feature->GetFieldAsInteger("IDPlot2"));
	    					NewFeature->SetField("ID2", mp_Feature->GetFieldAsString("ID2"));
	    					NewFeature->SetField("IDPlot3", mp_Feature->GetFieldAsInteger("IDPlot3"));
	    					NewFeature->SetField("ID3", mp_Feature->GetFieldAsString("ID3"));
	    					mp_Layer->CreateFeature(NewFeature);
	    					if(std::find(FIDToRemoveList.begin(), FIDToRemoveList.end(), mp_Feature->GetFID()) == FIDToRemoveList.end())
	    					{
	    						FIDToRemoveList.push_back(mp_Feature->GetFID());
	    					}
	    				}
	    				else
	    				{
	    					unsigned int pos = std::find(PlotIDMaxList[0].begin(), PlotIDMaxList[0].end(), mp_Feature->GetFieldAsInteger("IDPlot")) - PlotIDMaxList[0].begin();
	    					unsigned int IDMax = PlotIDMaxList[1][pos];
	    					PlotIDMaxList[1][pos] = IDMax+1;
	    					std::string IDNew = std::to_string(mp_Feature->GetFieldAsInteger("IDPlot"))+"N"+std::to_string(IDMax+1);
	    					NewFeature->SetField(m_IDFieldName.c_str(), mp_Feature->GetFieldAsInteger(m_IDFieldName.c_str()));
	    					NewFeature->SetField("IDPlot", mp_Feature->GetFieldAsInteger("IDPlot"));
	    					NewFeature->SetField("ID", IDNew.c_str());
	    					NewFeature->SetField("IDPlot2", mp_Feature->GetFieldAsInteger("IDPlot2"));
	    					NewFeature->SetField("ID2", mp_Feature->GetFieldAsString("ID2"));
	    					NewFeature->SetField("IDPlot3", mp_Feature->GetFieldAsInteger("IDPlot3"));
	    					NewFeature->SetField("ID3", mp_Feature->GetFieldAsString("ID3"));
	    					mp_Layer->CreateFeature(NewFeature);
	    					if(std::find(FIDToRemoveList.begin(), FIDToRemoveList.end(), mp_Feature->GetFID()) == FIDToRemoveList.end())
	    					{
	    						FIDToRemoveList.push_back(mp_Feature->GetFID());
	    					}
	    				}
	    			}
	    		}
	    	}
	    }
		OGRFeature::DestroyFeature(NewFeature);
	    OGRFeature::DestroyFeature(mp_Feature);
	    delete(mp_Geometry);
	}

	deleteFeatures(FIDToRemoveList);

	repack();
}


//======================================================================
//======================================================================


std::vector<unsigned int> VectorRegrouping::getConnectedFIDList()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	OGRDataSource *DataSource  = OGRSFDriverRegistrar::Open((m_FullFilePath).c_str(), TRUE );
	OGRLayer *Layer = DataSource->GetLayer(0);

	std::vector<unsigned int> ConnectedFIDList;

	mp_Layer->ResetReading();

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		std::string ID2 = mp_Feature->GetFieldAsString("ID2");
	    mp_Geometry = mp_Feature->GetGeometryRef()->clone();
	    if ((std::string) mp_Feature->GetFieldAsString("ID") != "0N1")
	    {
	    	OGRFeature *Feature = nullptr;
	    	Layer->ResetReading();
	    	while((Feature = Layer->GetNextFeature()) != nullptr)
	    	{
	    		OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();
	    		std::string ID = Feature->GetFieldAsString("ID");
	    		if(ID == ID2 and mp_Geometry->Touches(Geometry))
	    		{
	    			ConnectedFIDList.push_back(mp_Feature->GetFID());
	    		}
	    		OGRFeature::DestroyFeature(Feature);
	    		delete(Geometry);
	    	}
	    }
	    OGRFeature::DestroyFeature(mp_Feature);
	    delete(mp_Geometry);
	}

	OGRDataSource::DestroyDataSource(DataSource);

	return ConnectedFIDList;

}


//======================================================================
//======================================================================


std::vector<unsigned int> VectorRegrouping::getFIDWithMissingConnectionList()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> FIDWithMissingConnectionList, ConnectedFIDList = getConnectedFIDList();

	mp_Layer->ResetReading();

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
	    if (std::find(ConnectedFIDList.begin(), ConnectedFIDList.end(), mp_Feature->GetFID()) == ConnectedFIDList.end() and
	    		(std::string) mp_Feature->GetFieldAsString("ID") != "0N1")
	    {
	    	FIDWithMissingConnectionList.push_back(mp_Feature->GetFID());
	    }
	    OGRFeature::DestroyFeature(mp_Feature);
	}

	return FIDWithMissingConnectionList;
}


//======================================================================
//======================================================================


void VectorRegrouping::setNewConnectionBasedOnGeometry()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> FIDWithMissingConnectionList = getFIDWithMissingConnectionList();

	for (unsigned int i = 0; i < FIDWithMissingConnectionList.size(); i++)
	{
		OGRFeature *FeatureInQuestion = mp_Layer->GetFeature(FIDWithMissingConnectionList[i]);
		OGRGeometry *GeometryInQuestion = FeatureInQuestion->GetGeometryRef()->clone();
		unsigned int IDPlot = FeatureInQuestion->GetFieldAsInteger("IDPlot");
		std::string ID = FeatureInQuestion->GetFieldAsString("ID");
		std::string ID2 = FeatureInQuestion->GetFieldAsString("ID2");

		mp_Layer->ResetReading();

		int ReceivingFeatureFID = -1;

		while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			if ((std::string) mp_Feature->GetFieldAsString("ID") == ID2)
			{
				ReceivingFeatureFID = mp_Feature->GetFID();
			}
			OGRFeature::DestroyFeature(mp_Feature);
		}

		if(ReceivingFeatureFID != -1)
		{
			mp_Layer->ResetReading();

			OGRFeature *ReceivingFeature = mp_Layer->GetFeature(ReceivingFeatureFID);
			OGRGeometry *ReceivingGeometry = ReceivingFeature->GetGeometryRef()->clone();
			mp_Layer->ResetReading();

			double LengthIntersection = 0;
			int NewReceivingFeatureFID = -1;

			while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
			{
				mp_Geometry = mp_Feature->GetGeometryRef()->clone();

				if (std::find(FIDWithMissingConnectionList.begin(), FIDWithMissingConnectionList.end(), mp_Feature->GetFID()) == FIDWithMissingConnectionList.end() and
						mp_Geometry->Touches(GeometryInQuestion) and mp_Feature->GetFieldAsInteger("IDPlot") != IDPlot  and
						(std::string) mp_Feature->GetFieldAsString("ID2") == ID2)
				{

					OGRGeometry *IntersectionGeometry = (mp_Geometry->Intersection(GeometryInQuestion));
					if(IntersectionGeometry->getGeometryType() == 2 and
							((OGRLineString *) IntersectionGeometry)->get_Length() > LengthIntersection)
					{
						LengthIntersection = ((OGRLineString *) IntersectionGeometry)->get_Length();
						NewReceivingFeatureFID = mp_Feature->GetFID();
					}
					else if(IntersectionGeometry->getGeometryType() == 5 and
							((OGRMultiLineString *) IntersectionGeometry)->get_Length() > LengthIntersection)
					{
						LengthIntersection = ((OGRMultiLineString *) IntersectionGeometry)->get_Length();
						NewReceivingFeatureFID = mp_Feature->GetFID();
					}
					else if(IntersectionGeometry->getGeometryType() == 7)
					{
						for(unsigned int i = 0; i < ((OGRGeometryCollection *) IntersectionGeometry)->getNumGeometries(); i++)
						{
							if(((OGRGeometryCollection *) IntersectionGeometry)->getGeometryRef(i)->getGeometryType() == 2 and
									((OGRLineString *) ((OGRGeometryCollection *) IntersectionGeometry)->getGeometryRef(i))->get_Length() > LengthIntersection)
							{
								LengthIntersection = ((OGRLineString *) ((OGRGeometryCollection *) IntersectionGeometry)->getGeometryRef(i))->get_Length();
								NewReceivingFeatureFID = mp_Feature->GetFID();
							}
						}
					}
					OGRFeature::DestroyFeature(mp_Feature);
					delete(mp_Geometry);
					delete(IntersectionGeometry);
				}
				else
				{
					OGRFeature::DestroyFeature(mp_Feature);
					delete(mp_Geometry);
				}
			}

			mp_Layer->ResetReading();

			if(NewReceivingFeatureFID != -1)
			{
				OGRFeature *NewReceivingFeature = mp_Layer->GetFeature(NewReceivingFeatureFID);
				FeatureInQuestion->SetField("IDPlot2", NewReceivingFeature->GetFieldAsInteger("IDPlot"));
				FeatureInQuestion->SetField("ID2", NewReceivingFeature->GetFieldAsString("ID"));
				FeatureInQuestion->SetField("IDPlot3", NewReceivingFeature->GetFieldAsInteger("IDPlot2"));
				FeatureInQuestion->SetField("ID3", NewReceivingFeature->GetFieldAsString("ID2"));
				mp_Layer->SetFeature(FeatureInQuestion);
				mp_Layer->ResetReading();

				while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
				{
					if((std::string) mp_Feature->GetFieldAsString("ID2") == ID)
					{
						mp_Feature->SetField("IDPlot3", NewReceivingFeature->GetFieldAsInteger("IDPlot"));
						mp_Feature->SetField("ID3", NewReceivingFeature->GetFieldAsString("ID"));
						mp_Layer->SetFeature(mp_Feature);
					}
					OGRFeature::DestroyFeature(mp_Feature);
				}

				OGRFeature::DestroyFeature(NewReceivingFeature);
			}

			OGRFeature::DestroyFeature(ReceivingFeature);
			delete(ReceivingGeometry);
		}
		OGRFeature::DestroyFeature(FeatureInQuestion);
		delete(GeometryInQuestion);
	}

	mp_DataSource->SyncToDisk();
}


//======================================================================
//======================================================================


void VectorRegrouping::setNewConnectionBasedOnDistance()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> FIDWithMissingConnectionList = getFIDWithMissingConnectionList();

	for (unsigned int i = 0; i < FIDWithMissingConnectionList.size(); i++)
	{
		OGRFeature *FeatureInQuestion = mp_Layer->GetFeature(FIDWithMissingConnectionList[i]);
		OGRGeometry *GeometryInQuestion = FeatureInQuestion->GetGeometryRef()->clone();

		unsigned int IDPlot = FeatureInQuestion->GetFieldAsInteger("IDPlot");
		std::string ID = FeatureInQuestion->GetFieldAsString("ID");
		std::string ID2 = FeatureInQuestion->GetFieldAsString("ID2");

		mp_Layer->ResetReading();

		int ReceivingFeatureFID = -1;

		while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
		{
			if ((std::string) mp_Feature->GetFieldAsString("ID") == ID2)
			{
				ReceivingFeatureFID = mp_Feature->GetFID();
			}
			OGRFeature::DestroyFeature(mp_Feature);
		}

		if(ReceivingFeatureFID != -1)
		{
			mp_Layer->ResetReading();

			OGRFeature *ReceivingFeature = mp_Layer->GetFeature(ReceivingFeatureFID);
			OGRGeometry *ReceivingGeometry = ReceivingFeature->GetGeometryRef()->clone();

			mp_Layer->ResetReading();

			double Distance = 9999999;
			int NewReceivingFeatureFID = -1;

			while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
			{
				mp_Geometry = mp_Feature->GetGeometryRef()->clone();
				if (std::find(FIDWithMissingConnectionList.begin(), FIDWithMissingConnectionList.end(), mp_Feature->GetFID()) == FIDWithMissingConnectionList.end() and
						mp_Geometry->Touches(GeometryInQuestion) and mp_Feature->GetFieldAsInteger("IDPlot") != IDPlot  and
						(std::string) mp_Feature->GetFieldAsString("ID2") != ID and (std::string) mp_Feature->GetFieldAsString("ID3") != ID)
				{
					if(mp_Geometry->Distance(ReceivingGeometry) < Distance)
					{
						Distance = mp_Geometry->Distance(ReceivingGeometry);
						NewReceivingFeatureFID = mp_Feature->GetFID();
					}
				}
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
			}

			mp_Layer->ResetReading();

			if(NewReceivingFeatureFID != -1)
			{
				OGRFeature *NewReceivingFeature = mp_Layer->GetFeature(NewReceivingFeatureFID);
				FeatureInQuestion->SetField("IDPlot2", NewReceivingFeature->GetFieldAsInteger("IDPlot"));
				FeatureInQuestion->SetField("ID2", NewReceivingFeature->GetFieldAsString("ID"));
				FeatureInQuestion->SetField("IDPlot3", NewReceivingFeature->GetFieldAsInteger("IDPlot2"));
				FeatureInQuestion->SetField("ID3", NewReceivingFeature->GetFieldAsString("ID2"));
				mp_Layer->SetFeature(FeatureInQuestion);
				mp_Layer->ResetReading();

				while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
				{
					if((std::string) mp_Feature->GetFieldAsString("ID2") == ID)
					{
						mp_Feature->SetField("IDPlot3", NewReceivingFeature->GetFieldAsInteger("IDPlot"));
						mp_Feature->SetField("ID3", NewReceivingFeature->GetFieldAsString("ID"));
						mp_Layer->SetFeature(mp_Feature);
					}
					OGRFeature::DestroyFeature(mp_Feature);
				}

				OGRFeature::DestroyFeature(NewReceivingFeature);
			}

			OGRFeature::DestroyFeature(FeatureInQuestion);
			OGRFeature::DestroyFeature(ReceivingFeature);
			delete(GeometryInQuestion);
			delete(ReceivingGeometry);
		}
		else
		{
			OGRFeature::DestroyFeature(FeatureInQuestion);
			delete(GeometryInQuestion);
		}
	}

	mp_DataSource->SyncToDisk();
}


//======================================================================
//======================================================================


void VectorRegrouping::dissolveInnerFeatures()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<unsigned int> FIDWithMissingConnectionList = getFIDWithMissingConnectionList(), FIDToRemoveList;

	for (unsigned int i = 0; i < FIDWithMissingConnectionList.size(); i++)
    {
		OGRFeature *FeatureInQuestion = mp_Layer->GetFeature(FIDWithMissingConnectionList[i]);
		OGRGeometry *GeometryInQuestion = FeatureInQuestion->GetGeometryRef()->clone();

		unsigned int IDPlot = FeatureInQuestion->GetFieldAsInteger("IDPlot");
    	std::string ID = FeatureInQuestion->GetFieldAsString("ID");
    	std::string ID2 = FeatureInQuestion->GetFieldAsString("ID2");

    	mp_Layer->ResetReading();

    	double LengthIntersection = 0;
    	int FeatureFID = -1;

    	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
    	{
    		mp_Geometry = mp_Feature->GetGeometryRef()->clone();
    		if (std::find(FIDWithMissingConnectionList.begin(), FIDWithMissingConnectionList.end(), mp_Feature->GetFID()) == FIDWithMissingConnectionList.end() and
    				mp_Geometry->Touches(GeometryInQuestion) and mp_Feature->GetFieldAsInteger("IDPlot") == IDPlot  and
					(std::string) mp_Feature->GetFieldAsString("ID2") != ID and (std::string) mp_Feature->GetFieldAsString("ID3") != ID)
    		{
    			OGRGeometry *IntersectionGeometry = mp_Geometry->Intersection(GeometryInQuestion);
    			if(IntersectionGeometry->getGeometryType() == 2 and ((OGRLineString *) IntersectionGeometry)->get_Length() > LengthIntersection)
    			{
    				LengthIntersection = ((OGRLineString *) IntersectionGeometry)->get_Length();
    				FeatureFID = mp_Feature->GetFID();
    			}
    			else if(IntersectionGeometry->getGeometryType() == 5 and
    					((OGRMultiLineString *) IntersectionGeometry)->get_Length() > LengthIntersection)
    			{
    				LengthIntersection = ((OGRLineString *) IntersectionGeometry)->get_Length();
    				FeatureFID = mp_Feature->GetFID();
    			}
    			else if(IntersectionGeometry->getGeometryType() == 7)
    			{
    				for(unsigned int i = 0; i < ((OGRGeometryCollection *) IntersectionGeometry)->getNumGeometries(); i++)
    				{
    					if(((OGRGeometryCollection *) IntersectionGeometry)->getGeometryRef(i)->getGeometryType() == 2 and
    							((OGRLineString *)((OGRGeometryCollection *) IntersectionGeometry)->getGeometryRef(i))->get_Length() > LengthIntersection)
    					{
    						LengthIntersection = ((OGRLineString *) IntersectionGeometry)->get_Length();
    						FeatureFID = mp_Feature->GetFID();
    					}
    				}
    			}
    			delete(IntersectionGeometry);
    		}
    		OGRFeature::DestroyFeature(mp_Feature);
    		delete(mp_Geometry);
    	}

    	mp_Layer->ResetReading();

    	if(FeatureFID != -1)
    	{
    		OGRFeature *Feature = mp_Layer->GetFeature(FeatureFID);
    		OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();

    		OGRGeometry *TempGeometry = Geometry->Union(GeometryInQuestion);

    		Feature->SetGeometry(TempGeometry);
    		mp_Layer->SetFeature(Feature);

    		std::string IDNew = Feature->GetFieldAsString("ID");
    		unsigned int IDPlot2New = Feature->GetFieldAsInteger("IDPlot2");
    		std::string ID2New = Feature->GetFieldAsString("ID2");
    		unsigned int IDPlot3New = Feature->GetFieldAsInteger("IDPlot3");
    		std::string ID3New = Feature->GetFieldAsString("ID3");

    		mp_Layer->ResetReading();

    		while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
    		{
    			if((std::string) mp_Feature->GetFieldAsString("ID2") == ID)
    			{
    				mp_Feature->SetField("ID2", IDNew.c_str());
    				mp_Feature->SetField("IDPlot3", (int) IDPlot2New);
    				mp_Feature->SetField("ID3", ID2New.c_str());
    				mp_Layer->SetFeature(mp_Feature);
    			}
    			OGRFeature::DestroyFeature(mp_Feature);
    		}

    		mp_Layer->ResetReading();

    		while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
    		{
    			if((std::string) mp_Feature->GetFieldAsString("ID3") == ID)
    			{
    				mp_Feature->SetField("ID3", IDNew.c_str());
    				mp_Layer->SetFeature(mp_Feature);
    			}
    			OGRFeature::DestroyFeature(mp_Feature);
    		}

    		OGRFeature::DestroyFeature(Feature);
    		delete(Geometry);
    		delete(TempGeometry);
    		FIDToRemoveList.push_back(FIDWithMissingConnectionList[i]);
    	}

    	OGRFeature::DestroyFeature(FeatureInQuestion);
    	delete(GeometryInQuestion);
    }

	deleteFeatures(FIDToRemoveList);

	repack();
}


//======================================================================
//======================================================================


void VectorRegrouping::regroupFeaturesBasedOnCommonReceivingFeature()
{

	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<int> FIDToMergeList;

	std::vector<std::string> OldIDList, NewIDList;

	OGRDataSource *DataSource  = OGRSFDriverRegistrar::Open((m_FullFilePath).c_str(), TRUE );

	mp_Layer->ResetReading();

	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		mp_Geometry = mp_Feature->GetGeometryRef()->clone();

		unsigned int FID = mp_Feature->GetFID();
		unsigned int IDPlot = mp_Feature->GetFieldAsInteger("IDPlot");
		std::string ID = mp_Feature->GetFieldAsString("ID");
		std::string ID2 = mp_Feature->GetFieldAsString("ID2");

		if (ID != "0N1" and std::find(FIDToMergeList.begin(), FIDToMergeList.end(), FID) == FIDToMergeList.end())
		{
			int ReceivingFeatureFID = -1;

			OGRLayer *Layer = DataSource->GetLayer(0);
			OGRFeature *Feature = nullptr;
			Layer->ResetReading();

			while((Feature = Layer->GetNextFeature()) != nullptr)
			{
				if((std::string) Feature->GetFieldAsString("ID") == ID2)
				{
					ReceivingFeatureFID = Feature->GetFID();
				}
				OGRFeature::DestroyFeature(Feature);
			}

			if(ReceivingFeatureFID != -1)
			{
				Layer = DataSource->GetLayer(0);
				Layer->ResetReading();

				OGRFeature *ReceivingFeature = Layer->GetFeature(ReceivingFeatureFID);
				OGRGeometry *ReceivingGeometry = ReceivingFeature->GetGeometryRef()->clone();

				OGRFeature::DestroyFeature(ReceivingFeature);

				Layer = DataSource->GetLayer(0);
				Layer->ResetReading();

				while((Feature = Layer->GetNextFeature()) != nullptr)
				{
					OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();
					if(Feature->GetFID() != FID and Feature->GetFieldAsInteger("IDPlot") == IDPlot and (std::string) Feature->GetFieldAsString("ID2") == ID2 and
							mp_Geometry->Touches(Geometry) and
							std::find(FIDToMergeList.begin(), FIDToMergeList.end(), Feature->GetFID()) == FIDToMergeList.end())
					{
						OGRGeometry *IntersectionGeometry1 = mp_Geometry->Intersection(ReceivingGeometry);
						OGRGeometry *IntersectionGeometry2 = Geometry->Intersection(ReceivingGeometry);
						if(IntersectionGeometry1->Touches(IntersectionGeometry2))
						{
							std::string IDNew = Feature->GetFieldAsString("ID");
							FIDToMergeList.push_back(FID);
							FIDToMergeList.push_back(Feature->GetFID());
							OGRGeometry *NewGeometry = mp_Geometry->Union(Geometry);
							Feature->SetGeometry(NewGeometry);
							mp_Layer->SetFeature(Feature);
							mp_Layer->DeleteFeature(FID);
							OldIDList.push_back(ID);
							NewIDList.push_back(IDNew);
							OGRFeature::DestroyFeature(Feature);
							delete(Geometry);
							delete(NewGeometry);
							delete(IntersectionGeometry1);
							delete(IntersectionGeometry2);
							break;
						}
						else
						{
							OGRFeature::DestroyFeature(Feature);
							delete(Geometry);
							delete(IntersectionGeometry1);
							delete(IntersectionGeometry2);
						}
					}
					else
					{
						OGRFeature::DestroyFeature(Feature);
						delete(Geometry);
					}
				}

				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
				delete(ReceivingGeometry);
			}
			else
			{
				OGRFeature::DestroyFeature(mp_Feature);
				delete(mp_Geometry);
			}
		}
		else
		{
			OGRFeature::DestroyFeature(mp_Feature);
			delete(mp_Geometry);
		}
	}

	OGRDataSource::DestroyDataSource(DataSource);

	repack();

	if(!OldIDList.empty())
	{
		updateTextFieldValues("ID2", "ID2", OldIDList, NewIDList);

		updateTextFieldValues("ID3", "ID3", OldIDList, NewIDList);
	}

	mp_DataSource->SyncToDisk();
}


//======================================================================
//======================================================================


bool VectorRegrouping::areThereFeaturesWithCommonReceivingFeature()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

   	unsigned int N = 0;

   	mp_Layer->ResetReading();

   	OGRDataSource *DataSource  = OGRSFDriverRegistrar::Open((m_FullFilePath).c_str(), TRUE );

   	while((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
   	{
   		mp_Geometry = mp_Feature->GetGeometryRef()->clone();

   		unsigned int FID = mp_Feature->GetFID();
   		unsigned int IDPlot = mp_Feature->GetFieldAsInteger("IDPlot");
   		std::string ID = mp_Feature->GetFieldAsString("ID");
   		std::string ID2 = mp_Feature->GetFieldAsString("ID2");

   		if (ID != "0N1")
   		{
   			int ReceivingFeatureFID = -1;

   			OGRLayer *Layer = DataSource->GetLayer(0);
   			OGRFeature *Feature = nullptr;
   			Layer->ResetReading();

   			while((Feature = Layer->GetNextFeature()) != nullptr)
   			{
   				if((std::string) Feature->GetFieldAsString("ID") == ID2)
   				{
   					ReceivingFeatureFID = Feature->GetFID();
   				}
   				OGRFeature::DestroyFeature(Feature);
   			}

   			if(ReceivingFeatureFID != -1)
   			{
   				Layer = DataSource->GetLayer(0);
   				Layer->ResetReading();

   				OGRFeature *ReceivingFeature = Layer->GetFeature(ReceivingFeatureFID);
   				OGRGeometry *ReceivingGeometry = ReceivingFeature->GetGeometryRef()->clone();

   				OGRFeature::DestroyFeature(ReceivingFeature);

   				Layer = DataSource->GetLayer(0);
   				Layer->ResetReading();

   				while((Feature = Layer->GetNextFeature()) != nullptr)
   				{
   					OGRGeometry *Geometry = Feature->GetGeometryRef()->clone();
   					if(Feature->GetFID() != FID and Feature->GetFieldAsInteger("IDPlot") == IDPlot and
   							(std::string) Feature->GetFieldAsString("ID2") == ID2 and
							mp_Geometry->Touches(Geometry))
   					{
   						OGRGeometry *IntersectionGeometry1 = mp_Geometry->Intersection(ReceivingGeometry);
   						OGRGeometry *IntersectionGeometry2 = Geometry->Intersection(ReceivingGeometry);
   						if(IntersectionGeometry1->Touches(IntersectionGeometry2))
   						{
   							N++;
   							OGRFeature::DestroyFeature(Feature);
   							delete(Geometry);
   							delete(IntersectionGeometry1);
   							delete(IntersectionGeometry2);
   							break;
   						}
   						else
   						{
   							OGRFeature::DestroyFeature(Feature);
   							delete(Geometry);
   							delete(IntersectionGeometry1);
   							delete(IntersectionGeometry2);
   						}
   					}
   					else
   					{
   						OGRFeature::DestroyFeature(Feature);
   						delete(Geometry);
   					}
   				}

   				OGRFeature::DestroyFeature(mp_Feature);
   				delete(mp_Geometry);
   				delete(ReceivingGeometry);
   			}
   			else
   			{
   				OGRFeature::DestroyFeature(mp_Feature);
   				delete(mp_Geometry);
   			}
   		}
   		else
   		{
   			OGRFeature::DestroyFeature(mp_Feature);
   			delete(mp_Geometry);
   		}
   	}

   	OGRDataSource::DestroyDataSource(DataSource);

   	if(N != 0)
   	{
   		return true;
   	}
   	else
   	{
   		return false;
   	}

}


//======================================================================
//======================================================================


void VectorRegrouping::setConsecutiveID()
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	std::vector<std::string> OldIDList, NewIDList;

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		unsigned int FID = mp_Feature->GetFID();
	    std::string ID = mp_Feature->GetFieldAsString("ID");
	    unsigned int IDPlot = std::stoi(ID.substr(0, ID.find("N")));
	    if (FID == 0)
	    {
	    	if (ID.substr(ID.find("N")+1, ID.length()-1) != "1")
	    	{
	    		OldIDList.push_back(ID.c_str());
	    		std::string IDNew = std::to_string(IDPlot) + "N1";
	    		NewIDList.push_back(IDNew);
	    		mp_Feature->SetField("ID", IDNew.c_str());
	    		mp_Layer->SetFeature(mp_Feature);
	    	}
	    }
	    else
	    {
	    	OGRFeature *FeatureAbove = mp_Layer->GetFeature(FID-1);
	    	std::string IDAbove = FeatureAbove->GetFieldAsString("ID");
	    	unsigned int IDPlotAbove = std::stoi(IDAbove.substr(0, IDAbove.find("N")));
	    	if (IDPlot == IDPlotAbove)
	    	{
	    		if(std::stoi(ID.substr(ID.find("N")+1, ID.length()-1)) != std::stoi(IDAbove.substr(IDAbove.find("N")+1, IDAbove.length()-1))+1)
	    		{
	    			OldIDList.push_back(ID.c_str());
	    			std::string IDNew = std::to_string(IDPlot) + "N" + std::to_string(std::stoi(IDAbove.substr(IDAbove.find("N")+1, IDAbove.length()-1))+1);
	    			NewIDList.push_back(IDNew.c_str());
	    			mp_Feature->SetField("ID", IDNew.c_str());
	    			mp_Layer->SetFeature(mp_Feature);
	    		}
	    	}
	    	else
	    	{
	    		if (ID.substr(ID.find("N")+1, ID.length()-1) != "1")
	    		{
	    			OldIDList.push_back(ID.c_str());
	    			std::string IDNew = std::to_string(IDPlot) + "N1";
	    			NewIDList.push_back(IDNew);
	    			mp_Feature->SetField("ID", IDNew.c_str());
	    			mp_Layer->SetFeature(mp_Feature);
	    		}
	    	}
	    	OGRFeature::DestroyFeature(FeatureAbove);
	    }
	    OGRFeature::DestroyFeature(mp_Feature);
	}

	mp_Layer->ResetReading();

	if(!OldIDList.empty())
	{
		updateTextFieldValues("ID2", "ID2", OldIDList, NewIDList);

		updateTextFieldValues("ID3", "ID3", OldIDList, NewIDList);
	}

	mp_DataSource->SyncToDisk();
}


//======================================================================
//======================================================================


void VectorRegrouping::updateTextFieldValues(std::string FieldNameToCheck, std::string FieldNameToUpdate, std::vector<std::string> OldIDList, std::vector<std::string> NewIDList)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		std::string Value = mp_Feature->GetFieldAsString(FieldNameToCheck.c_str());
		if (std::find(OldIDList.begin(), OldIDList.end(), Value.c_str()) != OldIDList.end())
		{
			unsigned int pos = std::find(OldIDList.begin(), OldIDList.end(), Value.c_str()) - OldIDList.begin();
			mp_Feature->SetField(FieldNameToUpdate.c_str(), NewIDList[pos].c_str());
			mp_Layer->SetFeature(mp_Feature);
		}
		OGRFeature::DestroyFeature(mp_Feature);
	}
}


//======================================================================
//======================================================================


void VectorRegrouping::updateIntegerFieldValues(std::string FieldNameToCheck, std::string FieldNameToUpdate, std::vector<std::string> OldIDList, std::vector<unsigned int> NewIDList)
{
	VERBOSE_MESSAGE(1,"Entering " << __PRETTY_FUNCTION__);

	mp_Layer->ResetReading();

	while ((mp_Feature = mp_Layer->GetNextFeature()) != nullptr)
	{
		std::string Value = mp_Feature->GetFieldAsString(FieldNameToCheck.c_str());
		if (std::find(OldIDList.begin(), OldIDList.end(), Value) != OldIDList.end())
		{
			unsigned int pos = std::find(OldIDList.begin(), OldIDList.end(), Value.c_str()) - OldIDList.begin();
			mp_Feature->SetField(FieldNameToUpdate.c_str(), (int) NewIDList[pos]);
			mp_Layer->SetFeature(mp_Feature);
		}
		OGRFeature::DestroyFeature(mp_Feature);
	}
}


//======================================================================
//======================================================================


