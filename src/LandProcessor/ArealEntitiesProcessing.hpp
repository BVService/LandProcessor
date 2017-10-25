/*
 * ArealEntitiesProcessing.hpp
 *
 *  Created on: 6 oct. 2017
 *      Author: zadonina
 */

#ifndef SRC_LANDPROCESSOR_AREALENTITIESPROCESSING_HPP_
#define SRC_LANDPROCESSOR_AREALENTITIESPROCESSING_HPP_

#include <string>
#include <vector>
#include <map>
#include <gdal/ogrsf_frmts.h>
#include <gdal/gdal_priv.h>
#include <gdal/cpl_conv.h>
#include <LandProcessor/VectorProcessing.hpp>
#include <LandProcessor/LinearEntitiesProcessing.hpp>
#include <LandProcessor/RasterProcessing.hpp>

class ArealEntitiesProcessing: public VectorProcessing
{
	friend class LinearEntitiesProcessing;
	friend class RasterProcessing;

  private:

  public:

	ArealEntitiesProcessing(const std::string &FilePath, const std::string &FileName);

    virtual ~ArealEntitiesProcessing();

    void setArealEntitiesIDs();

    void createArealEntitiesAttributes();

    void setArealEntitiesSurface();

    void setArealEntitiesLanduse(const std::string &LandUseFieldName, const std::string &ParcelVectorFileName, std::string OutputVectorFilePath = "");

    void setArealEntitiesSlope(const std::string &SlopeRasterFileName, const std::string &SlopeRasterFilePath);

    void setArealEntitiesFlowDistance(const std::string &LinearEntityVectorFileName, std::string OutputVectorFilePath = "");

};


#endif /* SRC_LANDPROCESSOR_AREALENTITIESPROCESSING_HPP_ */
