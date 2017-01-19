 
SET(TESTS_DATASETS_PATH "${CMAKE_CURRENT_LIST_DIR}/datasets")
SET(TESTS_REFERENCES_PATH "${CMAKE_CURRENT_LIST_DIR}/references")
SET(TESTS_EXECS_PATH "${CMAKE_BINARY_DIR}/tests-runs")



FILE(MAKE_DIRECTORY "${TESTS_RESULTS_PATH}")

CONFIGURE_FILE(tests-config.hpp.in tests-config.hpp @ONLY)

INCLUDE_DIRECTORIES(../src ${CMAKE_CURRENT_BINARY_DIR})



ADD_EXECUTABLE(test0 test0.cpp)
TARGET_LINK_LIBRARIES(test0 bvservice-landprocessor 
                            ${OpenFLUID_LIBRARIES} ${GDAL_LIBRARIES} ${QT_QTCORE_LIBRARIES} ${Boost_LIBRARIES})
ADD_TEST(test0 test0)


ADD_EXECUTABLE(LandscapeSegmentation LandscapeSegmentation.cpp)
TARGET_LINK_LIBRARIES(LandscapeSegmentation bvservice-landprocessor 
                      ${OpenFLUID_LIBRARIES} ${GDAL_LIBRARIES} ${QT_QTCORE_LIBRARIES} ${Boost_LIBRARIES})
ADD_TEST(LandscapeSegmentation LandscapeSegmentation)


ADD_EXECUTABLE(ScenarioTesting ScenarioTesting.cpp)
TARGET_LINK_LIBRARIES(ScenarioTesting bvservice-landprocessor 
                      ${OpenFLUID_LIBRARIES} ${GDAL_LIBRARIES} ${QT_QTCORE_LIBRARIES} ${Boost_LIBRARIES})
ADD_TEST(ScenarioTesting ScenarioTesting)
SET_PROPERTY(TEST ScenarioTesting APPEND PROPERTY DEPENDS LandscapeSegmentation)
