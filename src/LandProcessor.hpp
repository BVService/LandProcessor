/*
 * LandProcessor.hpp
 *
 *  Created on: 8 sept. 2016
 *      Author: ezadonina
 */

#ifndef __LANDPROCESSOR_HPP__
#define __LANDPROCESSOR_HPP__


#include <string>
#include <vector>
#include <gdal/ogrsf_frmts.h>
#include <gdal/gdal_priv.h>
#include <gdal/cpl_conv.h>


class LandProcessor
{
  private:

    std::string m_InputPath;

    std::string m_OutputPath;

    const std::string m_VectorDir = "vector";

    const std::string m_RasterDir = "raster";

    const std::string m_InputPlotsFile = "plots.shp";

    const std::string m_InputDitchesFile = "fosses.shp";

    const std::string m_InputHedgesFile = "haies.shp";

    const std::string m_InputGrassBandFile = "bandesenherbees.shp";

    const std::string m_InputRivesrFile = "coursdeau.shp";

    const std::string m_InputThalwegsFile = "talweg.shp";

    const std::string m_InputBenchesFile = "talus.shp";

    const std::string m_InputDEMFile = "dem.tif";

    const std::string m_OutputDEMFile = "dem.tif";

    const std::string m_OutputPlotsRasterFile = "plots.tif";

    const std::string m_OutputSlopeRasterFile = "slope.tif";

    const std::string m_OutputDrainageRasterFile = "drainage.tif";

    const std::string m_OutputIDRasterFile = "rasterID.tif";

    const std::string m_OutputOutletsRasterFile = "outlets.tif";

    const std::string m_OutputReceiversRasterFile = "receivers.tif";

    const std::string m_OutputDownslopeRasterFile = "downslope.tif";

    const std::string m_OutputCatchmentsRasterFile = "catchments.tif";

    const std::string m_OutputPlotsVectorFile = "plots.shp";

    const std::string m_OutputPlotsLimitsVectorFile = "plotslimits.shp";

    const std::string m_OutputCatchmentsVectorFile = "catchments.shp";

    const std::string m_OutputOutletsVectorFile = "outlets.shp";

    const std::string m_OutputReceiversVectorFile = "receivers.shp";

    const std::string m_OutputCatchmentsGroupedVectorFile = "catchmentsgrouped.shp";

    const std::string m_OutputEntitiesVectorFile = "entities.shp";

    const std::string m_OutputPlotsAndEntitiesUnionVectorFile = "unionplotsentities.shp";

    const std::string m_OutputEntitiesGroupedVectorFile = "entitiesgrouped.shp";

    const std::string m_OutputSurfaceEntitiesVectorFile = "SRF.shp";

    const std::string m_OutputLinearEntitiesVectorFile = "LNR.shp";

    const std::string m_OutputSUVectorFile = "SU.shp";

    const std::string m_OutputRSVectorFile = "RS.shp";

    const std::string m_OutputLIVectorFile = "LI.shp";

    const std::string m_VectorDriverName = "ESRI Shapefile";

    const std::string m_RasterDriverName = "GTiff";

    OGRSFDriver *mp_VectorDriver;

    OGRSpatialReference *mp_SRS;

    OGRPoint *mp_Point;

    GDALDriver *mp_RasterDriver;

    double *mp_GeoTransformVal;

    bool m_IsReady = true;

    unsigned int m_MinEntSize = 250;

    std::string m_SnapDistance = "1.0e-08";

    std::string m_IDFieldName = "OFLD_ID";

    std::string m_LandUseFieldName = "LandUse";

    const std::string m_EPSGCode = "EPSG:2154";

    std::string m_SQLDialect = "SQLITE";

    std::string m_SQLRequest;

    std::string m_QMark = "'";


  public:

	LandProcessor(const std::string& InputPath, const std::string& OutputPath);

	virtual ~LandProcessor();

    bool isReady() const
    { return m_IsReady; }

    std::string getOutputVectorPath(const std::string& Filename = "") const;

    std::string getInputVectorPath(const std::string& Filename = "") const;

    std::string getOutputRasterPath(const std::string& Filename = "") const;

    std::string getInputRasterPath(const std::string& Filename = "") const;

	void preprocessVectorData();

	void preprocessRasterData();

	void createSRFandLNR();

	void setSRFParameters();

	void setLNRParameters();

	void createSU();

	void createRS();

	void createLI();

	void setSUParameters();

	void setRSParameters();

	void setLIParameters();

	void extractPlotsLimits();

	void attributeLinearStructures();

	int extractFromRasterToPoint(GDALDataset *Dataset, unsigned int RasterBandIndex);

	void createField(OGRLayer *LayerName, std::string FieldName, OGRFieldType FieldType);

	void getCentroidPoint(OGRGeometry *Geometry);

	std::pair <double, double> getCoordinatesOfPoint();

	std::pair <double, double> calculateOffset(double XCoord, double YCoord);

	GDALRasterBand* getRasterBand(GDALDataset *Dataset, unsigned int RastreBandIndex);

	void getGeoTransform(GDALDataset *Dataset);

	void checkLinearStructuresVectorData();

	void checkPolygonVectorData();

};


#endif /* __LANDPROCESSOR_HPP__ */
