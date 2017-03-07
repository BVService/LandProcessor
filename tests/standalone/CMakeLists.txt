
SET(TESTS_DATASETS_PATH "${CMAKE_CURRENT_LIST_DIR}/datasets")
SET(TESTS_REFERENCES_PATH "${CMAKE_CURRENT_LIST_DIR}/references")
SET(TESTS_EXECS_PATH "${CMAKE_BINARY_DIR}/tests-standalone")



FILE(MAKE_DIRECTORY "${TESTS_EXECS_PATH}")

CONFIGURE_FILE(tests-config.hpp.in tests-config.hpp @ONLY)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR})



ADD_EXECUTABLE(standalone-Full Full.cpp)
TARGET_LINK_LIBRARIES(standalone-Full bvservice-landprocessor
                                      ${OpenFLUID_LIBRARIES}
                                      ${GDAL_LIBRARIES}
                                      ${QT_QTCORE_LIBRARIES})
ADD_TEST(NAME standalone-Full COMMAND standalone-Full)


ADD_EXECUTABLE(standalone-Preparation Preparation.cpp)
TARGET_LINK_LIBRARIES(standalone-Preparation bvservice-landprocessor
                                            ${OpenFLUID_LIBRARIES}
                                            ${GDAL_LIBRARIES}
                                            ${QT_QTCORE_LIBRARIES})
ADD_TEST(NAME standalone-Preparation COMMAND standalone-Preparation)


ADD_EXECUTABLE(standalone-Scenario Scenario.cpp)
TARGET_LINK_LIBRARIES(standalone-Scenario bvservice-landprocessor
                                          ${OpenFLUID_LIBRARIES}
                                          ${GDAL_LIBRARIES}
                                          ${QT_QTCORE_LIBRARIES})
ADD_TEST(NAME standalone-Scenario COMMAND standalone-Scenario)
SET_PROPERTY(TEST standalone-Scenario APPEND PROPERTY DEPENDS standalone-Preparation)
