

ADD_SUBDIRECTORY(simulators)


SET(OPENFLUID_RUN_OPTS "--simulators-paths=${CMAKE_CURRENT_BINARY_DIR}/simulators/tests.preprocess-stage:${CMAKE_CURRENT_BINARY_DIR}/simulators/tests.scenario-stage")

ADD_TEST(NAME openfluid-Full 
         COMMAND "${OpenFLUID_CMD_PROGRAM}" run "${CMAKE_CURRENT_SOURCE_DIR}/datasets/Full" "${CMAKE_BINARY_DIR}/tests-openfluid/Full"
                                            ${OPENFLUID_RUN_OPTS})