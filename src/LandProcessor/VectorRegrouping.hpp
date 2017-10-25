/*
 * VectorRegrouping.hpp
 *
 *  Created on: 4 oct. 2017
 *      Author: zadonina
 */

#ifndef SRC_VECTORREGROUPING_HPP_
#define SRC_VECTORREGROUPING_HPP_

#include <string>
#include <vector>
#include <map>
#include <gdal/ogrsf_frmts.h>
#include <gdal/gdal_priv.h>
#include <gdal/cpl_conv.h>
#include <LandProcessor/VectorProcessing.hpp>

class VectorRegrouping: public VectorProcessing
{
  private:

    int findAggregatingFeature(unsigned int FID);

    std::vector<unsigned int> findFeaturesThatShareAttributesAndLimits(unsigned int FID);

  public:

	VectorRegrouping(const std::string &FilePath, const std::string &FileName);

    virtual ~VectorRegrouping();

    std::vector<unsigned int> getPlotIDWithLargeFIDList(unsigned int MinEntSize);

    std::vector<unsigned int> getWithAllSmallFIDList(unsigned int MinEntSize);

    std::vector<unsigned int> getSmallFIDList(unsigned int MinEntSize);

    void regroupSmallFeaturesThatConstituteSinglePlot(unsigned int MinEntSize);

    void regroupSmallFeaturesWithLargerFeatures(unsigned int MinEntSize);

    void regroupFeaturesThatShareAttributesAndLimits();

    void regroupLockedFeatures();

    void createOuterPolygon();

    std::vector<unsigned int> getPlotIDList();

    std::vector<unsigned int> getMissingPlotIDList();

    std::vector<unsigned int> getSliverFIDList();

    void regroupSliverFeaturesWithSingleNonSliverNeighbor();

    void regroupSliverFeaturesBasedOnGeometryAndAttributes();

    void regroupSliverFeaturesBasedOnGeometry();

    void regroupSliverFeaturesBasedOnAttributes();

    void regroupFeaturesWithMissingPlotID();

    void correctOuterPolygon();

    void correctMultiparts();

    std::vector<unsigned int> getConnectedFIDList();

    std::vector<unsigned int> getFIDWithMissingConnectionList();

    void setNewConnectionBasedOnGeometry();

    void setNewConnectionBasedOnDistance();

    void dissolveInnerFeatures();

    void regroupFeaturesBasedOnCommonReceivingFeature();

    bool areThereFeaturesWithCommonReceivingFeature();

    void setConsecutiveID();

    void updateTextFieldValues(std::string FieldNameToCheck, std::string FieldNameToUpdate, std::vector<std::string> OldIDList, std::vector<std::string> NewIDList);

    void updateIntegerFieldValues(std::string FieldNameToCheck, std::string FieldNameToUpdate, std::vector<std::string> OldIDList, std::vector<unsigned int> NewIDList);

};



#endif /* SRC_VECTORREGROUPING_HPP_ */
