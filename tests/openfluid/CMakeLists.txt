

SET(OPENFLUID_RUN_OPTS "--simulators-paths=${CMAKE_BINARY_DIR}/src/openfluid-simulators/landprocessor.preprocess-stage:${CMAKE_BINARY_DIR}/src/openfluid-simulators/landprocessor.scenario-stage")

ADD_TEST(NAME openfluid-Full
         COMMAND "${OpenFLUID_CMD_PROGRAM}" run "${CMAKE_CURRENT_SOURCE_DIR}/datasets/Full" "${CMAKE_BINARY_DIR}/tests-openfluid/Full"
                                            ${OPENFLUID_RUN_OPTS})
