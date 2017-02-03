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
  @file Helpers.hpp
  @author Jean-Christophe FABRE <jean-christophe.fabre@inra.fr>
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
