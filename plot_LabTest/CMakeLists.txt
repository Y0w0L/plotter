# plot_inPixel/CMakeLists.txt

SET(LIB_SOURCES
    plot_LabTest.cpp)

ADD_LIBRARY(plot_LabTest_lib STATIC ${LIB_SOURCES})

TARGET_INCLUDE_DIRECTORIES(plot_LabTest_lib
                        PUBLIC
                            ${PROJECT_SOURCE_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}
                            ${ROOT_INCLUDE_DIRS})