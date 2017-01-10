/*
 * Helpers.hpp
 *
 *  Created on: 6 janv. 2017
 *      Author: fabrejc
 */


#ifndef __HELPERS_HPP__
#define __HELPERS_HPP__


#include <iostream>
#include <sstream>

#ifndef VERBOSE_LEVEL
#  define VERBOSE_LEVEL 0
#endif


#define _STREAMTOSTRING(_stream) ((static_cast<std::ostringstream&>(std::ostringstream().flush() << _stream)).str())

#define VERBOSE_MESSAGE(level,msg) \
  if ((level) <= VERBOSE_LEVEL) std::cout << std::string(2*(level),' ') << (_STREAMTOSTRING(msg)) << std::endl;


#endif /* __HELPERS_HPP__ */
