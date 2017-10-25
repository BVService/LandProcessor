/*
 * RasterProcessing.hpp
 *
 *  Created on: 21 sept. 2017
 *      Author: zadonina
 */

#ifndef SRC_LANDPROCESSOR_RASTERPROCESSING_HPP_
#define SRC_LANDPROCESSOR_RASTERPROCESSING_HPP_

#include <string>
#include <vector>
#include <map>
#include <gdal/ogrsf_frmts.h>
#include <gdal/gdal_priv.h>
#include <gdal/cpl_conv.h>

class RasterProcessing
{
	friend class VectorProcessing;
	friend class ArealEntitiesProcessing;

  protected:

	std::string m_FilePath;

	std::string m_FileName;

	std::string m_FullFilePath = m_FilePath + "/" + m_FileName;

    const std::string m_EPSGCode = "EPSG:2154";

    const std::string m_DriverName = "GTiff";

    OGRSpatialReference mp_SRS;

    const char *m_ProjectionRef = nullptr;

    GDALDriver *mp_Driver = nullptr;

    GDALDataset *mp_Dataset = nullptr;

    GDALRasterBand *mp_Band = nullptr;

    double m_GeoTransform[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  public:

	RasterProcessing(const std::string &FilePath, const std::string &FileName, GDALAccess AccessType);

    virtual ~RasterProcessing();

    void deleteFile(std::string FullFilePath = "");

    std::vector<double> getGeoTransformValues();

    int extractFromRasterToPoint(std::pair <double, double> Coordinates);

    std::pair <double, double> calculateOffset(std::pair <double, double> Coordinates);

    std::pair <double, double> getMinMaxValues();

    std::pair <int, int> getRasterSize();

    void checkEPSG();

    void createIDRaster(std::string FileName, std::string FilePath = "");

    void createOutletRaster(std::string PlotsRasterName,
    		std::string DrainageRasterName, std::string IDRasterName, std::string FileName, std::string FilePath = "");

    void createReceiversRaster(std::string OutletsRasterName,
    		std::string DrainageRasterName, std::string IDRasterName, std::string FileName, std::string FilePath = "");

    void createConnectingRaster(std::string OutletsRasterName,
            		std::string DrainageRasterName, std::string IDRasterName, std::string FileName, std::string FilePath = "");

    void createCatchmentsRaster(std::string OutletsRasterName, std::string DrainageRasterName, std::string PlotsRasterName, std::string FileName, std::string FilePath = "");

};

#endif /* SRC_LANDPROCESSOR_RASTERPROCESSING_HPP_ */
