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
  @file test0.cpp
  @author Ekaterina Zadonina <ekaterina.zadonina@inra.fr>
  @author Jean-Christophe FABRE <jean-christophe.fabre@inra.fr>
*/


#include <iostream>

#include <openfluid/base/Environment.hpp>
#include <openfluid/tools/FileHelpers.hpp>
#include <openfluid/tools/Filesystem.hpp>

#include <LandProcessor/LandProcessor.hpp>

#include "tests-config.hpp"
#include "TestHelpers.hpp"


int main(int argc, char *argv[])
{
  openfluid::base::Environment::init();


  openfluid::tools::emptyDirectoryRecursively(TESTS_EXECS_PATH+"/DardaillonSmall_0/output");
  openfluid::tools::emptyDirectoryRecursively(TESTS_EXECS_PATH+"/DardaillonSmall_0/release");


  try
  {
    LandProcessor LP(TESTS_DATASETS_PATH+"/DardaillonSmall",
                     TESTS_EXECS_PATH+"/DardaillonSmall_0/output",
                     TESTS_EXECS_PATH+"/DardaillonSmall_0/release");

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
    LP.releaseFiles();
  }
  catch (std::exception &E)
  {
    std::cout << E.what() << '\n';
    return -1;
  }


  std::vector<std::string> ShapefilesToCompare = {"SRF","LNR","SU","RS","LI"};
  std::string WorkTmpPath = openfluid::tools::Filesystem::makeUniqueSubdirectory(openfluid::base::Environment::getTempDir(),
                                                                                 "full-shpcompare");

  for (auto& ShpFile : ShapefilesToCompare)
  {

    if (CompareShapefiles(TESTS_REFERENCES_PATH+"/DardaillonSmall/"+ShpFile+".shp",
                          TESTS_EXECS_PATH+"/DardaillonSmall_0/release/vector/"+ShpFile+".shp",
                          WorkTmpPath))
    {
      std::cout << "OK: "+ShpFile + " shapefile is similar to reference file" << std::endl;
    }
    else
    {
      std::cout << "Error: "+ShpFile + " shapefile is different from reference file" << std::endl;
      return -1;
    }
  }

}
