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
  @file TestHelpers.hpp
  @author Jean-Christophe FABRE <jean-christophe.fabre@inra.fr>
*/


#ifndef __TESTHELPERS_HPP__
#define __TESTHELPERS_HPP__


#include<cstdlib>

#include <openfluid/tools/Filesystem.hpp>


bool CompareShapefiles(const std::string& ReferenceFileName, const std::string& CandidateFileName,
                       const std::string& TmpWorkPath)
{
  std::string ConvertedRefFile = TmpWorkPath + "/ref-" + openfluid::tools::Filesystem::basename(ReferenceFileName)+".csv";
  std::string ConvertedCandFile = TmpWorkPath + "/cand-" + openfluid::tools::Filesystem::basename(CandidateFileName)+".csv";

  openfluid::tools::Filesystem::removeFile(ConvertedRefFile);
  openfluid::tools::Filesystem::removeFile(ConvertedCandFile);


  // -t_srs \"EPSG:4326\"
  std::string RefConversionCmd = "ogr2ogr -f \"CSV\" "+ConvertedRefFile+" "+ReferenceFileName+" -lco \"GEOMETRY=AS_WKT\" -lco \"LINEFORMAT=CRLF\" -lco \"SEPARATOR=SEMICOLON\"";
  if (std::system(RefConversionCmd.c_str()) != 0)
    return false;

  std::string CandConversionCmd = "ogr2ogr -f \"CSV\" "+ConvertedCandFile+" "+CandidateFileName+" -lco \"GEOMETRY=AS_WKT\" -lco \"LINEFORMAT=CRLF\" -lco \"SEPARATOR=SEMICOLON\"";
  if (std::system(CandConversionCmd.c_str()) != 0)
    return false;


  std::string DiffCmd = "diff "+ConvertedRefFile+" "+ConvertedCandFile;

  return (std::system(DiffCmd.c_str()) == 0);
}


#endif /* __SHAPEFILESTOOLS_HPP__ */
