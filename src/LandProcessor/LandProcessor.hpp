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
  @file LandProcessor.hpp
  @author Ekaterina Zadonina <ekaterina.zadonina@inra.fr>
  @author Jean-Christophe FABRE <jean-christophe.fabre@inra.fr>
*/


#ifndef __LANDPROCESSOR_HPP__
#define __LANDPROCESSOR_HPP__


#include <string>
#include <vector>
#include <map>
#include <gdal/ogrsf_frmts.h>
#include <gdal/gdal_priv.h>
#include <gdal/cpl_conv.h>


class LandProcessor
{
  private:

    std::string m_InputPath;

    std::string m_OutputPath;

    std::string m_ReleasePath;

    const std::string m_VectorDir = "vector";

    const std::string m_RasterDir = "raster";


    // input files

    const std::string m_InputPlotsVectorFile = "plots.shp";

    const std::string m_InputDitchesVectorFile = "fosses.shp";

    const std::string m_InputHedgesVectorFile = "haies.shp";

    const std::string m_InputGrassBandVectorFile = "bandesenherbees.shp";

    const std::string m_InputRiversVectorFile = "coursdeau.shp";

    const std::string m_InputThalwegsVectorFile = "talweg.shp";

    const std::string m_InputBenchesVectorFile = "talus.shp";


    const std::string m_InputDEMFile = "dem.tif";


    // Output files

    const std::string m_OutputPlotsVectorFile = "plots.shp";

    const std::string m_OutputPlotsLimitsVectorFile = "plotslimits.shp";

    const std::string m_OutputPlotsRasterVectorFile = "plotsrastervector.shp";

    const std::string m_OutputCatchmentsVectorFile = "catchments.shp";

    const std::string m_OutputOutletsVectorFile = "outlets.shp";

    const std::string m_OutputReceiversVectorFile = "receivers.shp";

    const std::string m_OutputCatchmentsGroupedVectorFile = "catchmentsgrouped.shp";

    const std::string m_OutputEntitiesVectorFile = "entities.shp";

    const std::string m_OutputPlotsAndEntitiesUnionVectorFile = "unionplotsentities.shp";

    const std::string m_OutputEntitiesGroupedVectorFile = "entitiesgrouped.shp";

    const std::string m_OutputLinearStructureVectorFile = "linearstructure.shp";

    const std::string m_OutputIntersectionVectorFile = "intersection.shp";


    const std::string m_OutputDEMFile = "dem.tif";

    const std::string m_OutputPlotsRasterFile = "plots.tif";

    const std::string m_OutputSlopeRasterFile = "slope.tif";

    const std::string m_OutputDrainageRasterFile = "drainage.tif";

    const std::string m_OutputIDRasterFile = "rasterID.tif";

    const std::string m_OutputOutletsRasterFile = "outlets.tif";

    const std::string m_OutputReceiversRasterFile = "receivers.tif";

    const std::string m_OutputDownslopeRasterFile = "downslope.tif";

    const std::string m_OutputCatchmentsRasterFile = "catchments.tif";


    const std::string m_OutputSurfaceEntitiesVectorFile = "SRF.shp";

    const std::string m_OutputLinearEntitiesVectorFile = "LNR.shp";

    const std::string m_OutputSUVectorFile = "SU.shp";

    const std::string m_OutputRSVectorFile = "RS.shp";

    const std::string m_OutputLIVectorFile = "LI.shp";


    // Files to release

    std::vector<std::string> m_VectorFilesToRelease = {};

    std::vector<std::string> m_RasterFilesToRelease = {};


    const std::string m_VectorDriverName = "ESRI Shapefile";

    const std::string m_RasterDriverName = "GTiff";

    OGRSFDriver *mp_VectorDriver = nullptr;

    OGRSpatialReference *mp_SRS = nullptr;

    OGRPoint *mp_Point = nullptr;

    GDALDriver *mp_RasterDriver = nullptr;

    double m_GeoTransformVal[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    bool m_IsReady = true;

    unsigned int m_MinEntSize = 250; // minimal entity size (normally in square meters)

    std::string m_SnapDistance = "1.0e-08"; // snap distance in meters

    std::string m_IDFieldName = "OFLD_ID";

    std::string m_LandUseFieldName = "LandUse";

    const std::string m_EPSGCode = "EPSG:2154";

    std::string m_SQLDialect = "SQLITE";

    std::string m_SQLRequest;

    std::string m_QMark = "'";

    std::string getInputVectorPath(const std::string& Filename = "") const;

    std::string getInputRasterPath(const std::string& Filename = "") const;

    std::string getOutputVectorPath(const std::string& Filename = "") const;

    std::string getOutputRasterPath(const std::string& Filename = "") const;

    std::string getReleaseVectorPath(const std::string& Filename = "") const;

    std::string getReleaseRasterPath(const std::string& Filename = "") const;

    static std::string getLayerNameFromFilename(const std::string& Filename);

    void releaseVectorFile(const std::string& Filename);

    void releaseRasterFile(const std::string& Filename) const;


  public:

    LandProcessor(const std::string& InputPath,
                  const std::string& OutputPath,
                  const std::string& ReleasePath);

    virtual ~LandProcessor();

    bool isReady() const
    { return m_IsReady; }

    /**
      Preprocessing of vector data
    */
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

    void releaseFiles();

    int extractFromRasterToPoint(GDALDataset *Dataset, unsigned int RasterBandIndex);

    void createField(OGRLayer *LayerName, std::string FieldName, OGRFieldType FieldType);

    void getCentroidPoint(OGRGeometry *Geometry);

    std::pair <double, double> getCoordinatesOfPoint();

    std::pair <double, double> calculateOffset(double XCoord, double YCoord);

    GDALRasterBand* getRasterBand(GDALDataset *Dataset, unsigned int RastreBandIndex);

    void getGeoTransform(GDALDataset *Dataset);

    void checkLinearVectorData();

    void checkPolygonVectorData();

};


#endif /* __LANDPROCESSOR_HPP__ */