/*
 * test0.cpp
 *
 *  Created on: 8 sept. 2016
 *      Author: ezadonina
 */


#include <iostream>

#include <openfluid/base/Environment.hpp>
#include <openfluid/tools/FileHelpers.hpp>

#include "LandProcessor.hpp"

#include "tests-config.hpp"



int main(int argc, char *argv[])
{
  openfluid::base::Environment::init();


  openfluid::tools::emptyDirectoryRecursively(TESTS_RESULTS_PATH+"/DardaillonSmall_0");

  LandProcessor LP(TESTS_DATASETS_PATH+"/DardaillonSmall",
		           TESTS_RESULTS_PATH+"/DardaillonSmall_0");

  if (!LP.isReady())
  {
    std::cout << "LandProcessor is not ready" << std::endl;
    return -1;
  }

  try
  {
	  LP.preprocessVectorData();
	  LP.preprocessRasterData();
    LP.createSRFandLNR();
	  LP.setSRFParameters();
	  LP.setLNRParameters();
    LP.extractPlotsLimits();
    LP.attributeLinearStructures();
    LP.createSU();
    LP.createRS();
    LP.createLI();
    LP.setSUParameters();
    LP.setRSParameters();
    LP.setLIParameters();
  }
  catch (std::exception &E)
  {
    std::cout << E.what() << '\n';
    return -1;
  }

}

