/*
 * test0.cpp
 *
 *  Created on: 8 sept. 2016
 *      Author: fabrejc
 */


#include <iostream>

#include <openfluid/base/Environment.hpp>

#include "LandProcessor.hpp"


int main(int argc, char *argv[])
{

  openfluid::base::Environment::init();

  std::cout << openfluid::base::Environment::getTempDir() << std::endl;

  LandProcessor LP(argv[1],argv[2]);

  std::cout << LP.getInputVectorPath() << std::endl;
  std::cout << LP.getInputRasterPath() << std::endl;

  if (!LP.isReady())
    return -1;

  LP.preprocessRasterData();

}



