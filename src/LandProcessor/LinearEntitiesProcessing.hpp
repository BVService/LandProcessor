/*
 * LinearEntitiesProcessing.hpp
 *
 *  Created on: 6 oct. 2017
 *      Author: zadonina
 */

#ifndef SRC_LANDPROCESSOR_LINEARENTITIESPROCESSING_HPP_
#define SRC_LANDPROCESSOR_LINEARENTITIESPROCESSING_HPP_

#include <string>
#include <vector>
#include <map>
#include <gdal/ogrsf_frmts.h>
#include <gdal/gdal_priv.h>
#include <gdal/cpl_conv.h>
#include <LandProcessor/VectorProcessing.hpp>

class LinearEntitiesProcessing: public VectorProcessing
{
	friend class ArealEntitiesProcessing;

  private:

  public:

	LinearEntitiesProcessing(const std::string &FilePath, const std::string &FileName);

    virtual ~LinearEntitiesProcessing();

    void setLinearEntitiesFaces(const std::string &BaseVectorFileName, std::string OutputVectorFilePath = "");

    void setLinearEntitiesIDs(const std::string &ArealEntitiesVectorFileName, std::string OutputVectorFilePath = "");

    void createLinearEntitiesAttributes();

    void setLinearEntitiesLength();

    void creatLinearStructuresIntersectionVector(const std::string &IntersectionVectorFileName, std::string OutputVectorFilePath = "");

    void setLinearStructuresIntersectionVector(const std::string &LinearStructureFileName, const std::string &ParcelVectorFileName, std::string OutputVectorFilePath = "");

    void attributeLinearStructures(const std::string &ParcelVectorFileName, std::string OutputVectorFilePath = "");

    void setAttributedLinearStructureLength(const std::string &LengthFieldName, const std::string &RelativeLengthFieldName, const std::string &IntersectionVectorFileName, std::string OutputVectorFilePath = "");

    void setLinearEntitiesFlowDistance(const std::string &ArealEntityVectorFileName, std::string OutputVectorFilePath = "");
};


#endif /* SRC_LANDPROCESSOR_LINEARENTITIESPROCESSING_HPP_ */
