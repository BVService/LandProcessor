

SET(OPENFLUID_RUN_OPTS "--simulators-paths=${CMAKE_BINARY_DIR}/src/openfluid-simulators/landprocessor.preparation-stage:${CMAKE_BINARY_DIR}/src/openfluid-simulators/landprocessor.scenario-stage")

ADD_TEST(NAME openfluid-Full
         COMMAND "${OpenFLUID_CMD_PROGRAM}" run "${CMAKE_CURRENT_SOURCE_DIR}/datasets/Full" "${CMAKE_BINARY_DIR}/tests-openfluid/Full"
                                            ${OPENFLUID_RUN_OPTS})
                                            
ADD_TEST(NAME openfluid-FullWithParam
         COMMAND "${OpenFLUID_CMD_PROGRAM}" run "${CMAKE_CURRENT_SOURCE_DIR}/datasets/FullWithParam" "${CMAKE_BINARY_DIR}/tests-openfluid/FullWithParam"
                                            ${OPENFLUID_RUN_OPTS})
