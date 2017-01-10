
ADD_DEFINITIONS(-DVERBOSE_LEVEL=3)

ADD_LIBRARY(bvservice-landprocessor SHARED LandProcessor.cpp LandProcessor.hpp
                                           Helpers.cpp Helpers.hpp)