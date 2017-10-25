/*
 * VectorProcessing.hpp
 *
 *  Created on: 26 juil. 2017
 *      Author: zadonina
 */

#ifndef SRC_LANDPROCESSOR_VECTORPROCESSING_HPP_
#define SRC_LANDPROCESSOR_VECTORPROCESSING_HPP_


#include <string>
#include <vector>
#include <map>
#include <gdal/ogrsf_frmts.h>
#include <gdal/gdal_priv.h>
#include <gdal/cpl_conv.h>

class VectorProcessing
{
	friend class RasterProcessing;

  protected:

	std::string m_FilePath;

	std::string m_FileName;

	std::string m_FullFilePath;

	std::string m_LayerName;

	std::vector <unsigned int> m_FIDlist;

    OGRDataSource *mp_DataSource = nullptr;

    OGRSpatialReference *mp_SRS = nullptr;

    OGRLayer *mp_Layer = nullptr;

    OGRFeature *mp_Feature = nullptr;

    OGRGeometry *mp_Geometry = nullptr;

    OGRPoint *mp_Point = nullptr;

    OGRSFDriver *mp_Driver = nullptr;

	std::vector <int> m_GeometryTypes = {};

    const std::string m_DriverName = "ESRI Shapefile";

    const std::string m_EPSGCode = "EPSG:2154";

    std::string m_IDFieldName = "OFLD_ID";

    std::string m_SQLDialect = "SQLITE";

    std::string m_SQLRequest;

    std::string m_QMark = "'";

    std::vector<unsigned int> findWithinPolygons();

    std::vector<unsigned int> findOverlaps(unsigned int FID);

    std::vector<unsigned int> findDuplicates(unsigned int FID);

    void correctWithinPolygons();

    void correctOverlaps();

    void correctMultyparts();

    void correctHoles();

    void correctDuplicates();

    void correctEmpty();

  public:

	VectorProcessing(const std::string &FilePath, const std::string &FileName);

    virtual ~VectorProcessing();

    OGRLayer* getLayer();

    bool doesDataExist(const std::string &FileName, std::string FilePath = "");

    bool canDataBeOpened(const std::string &FileName, std::string FilePath = "");

    bool isDriverSupported(const std::string &FileName, std::string FilePath = "");

    bool isSRSInfoProvided(const std::string &FileName, std::string FilePath = "");

    bool isOfDeclaredGeometryType(std::vector <int> GeometryTypes, const std::string &FileName, std::string FilePath = "");

    bool areGeometriesConform(const std::string &FileName, std::string FilePath = "");

    void checkEPSG();

    void checkGeometryType(std::vector<int> GeometryTypes);

    void checkGeometries();

    void copyFile(std::string NewFilePath, std::string NewFileName = "");

    void deleteFile(std::string FullFilePath = "");

    void deleteField(std::string FieldName = "");

    void deleteUnnecessaryFields(std::vector<std::string> FieldNamesToKeep);

    bool doesFieldExist(std::string FieldName = "");

    void createField(OGRFieldType FieldType = OFTInteger, std::string FieldName = "");

    void createIDField();

    void setDefaultFieldValue(const std::string &FieldName);

    void fillFieldFromRaster(const std::string &FieldName, const std::string &RasterFilePath, const std::string &RasterFileName);

    void correctGeometry();

    std::vector<std::pair<double,double>> getLayerExtent();

    std::string getLayerName(OGRLayer *Layer = nullptr);

    std::string getLayerNameFromFileName(std::string FileName = "");

    std::string getLayerNameFromFilePath(std::string FullFilePath = "");

    void repack(const std::string &FullFilePath = "");

    void createDataSource(OGRwkbGeometryType GeometryType, std::string FilePath, std::string FileName);

    void copyDataSourceToDataSource(const std::string &NewFileName, const std::string &NewFilePath);

    void polygonizeRaster(std::pair<OGRFieldType,std::string> FieldInformation, const std::string &RasterFilePath, const std::string &RasterFileName, std::string FullFilePath = "");

    std::pair <double, double> getCentroidCoordinates();

    void executeSQLRequest(const std::string &SQLExpression, std::string NewFileName, std::string NewFilePath = "");

    std::vector < std::vector <int>> transformToVector();

    void deleteFeatures(std::vector<unsigned int> FIDList);

    unsigned int getFeaturesCount(OGRLayer *Layer = nullptr);


};


#endif /* SRC_LANDPROCESSOR_VECTORPROCESSING_HPP_ */
