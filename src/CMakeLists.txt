
ADD_DEFINITIONS(-DVERBOSE_LEVEL=3)

SET(CPP_FILES LandProcessor.cpp Helpers.cpp)
SET(HPP_FILES LandProcessor.hpp Helpers.hpp)

ADD_LIBRARY(bvservice-landprocessor SHARED ${CPP_FILES} ${HPP_FILES})


INSTALL(FILES ${HPP_FILES}
        DESTINATION include)
                                           
INSTALL(TARGETS bvservice-landprocessor
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)