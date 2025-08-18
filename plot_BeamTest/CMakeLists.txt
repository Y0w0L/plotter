# plot_histogram/CMakeLists.txt

SET(LIB_SOURCES
    plot_BeamTest.cpp)

ADD_LIBRARY(plot_BeamTest_lib STATIC ${LIB_SOURCES})

TARGET_INCLUDE_DIRECTORIES(plot_BeamTest_lib
                        PUBLIC
                            ${PROJECT_SOURCE_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}
                            ${ROOT_INCLUDE_DIRS})