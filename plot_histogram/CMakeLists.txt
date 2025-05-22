# plot_histogram/CMakeLists.txt

SET(LIB_SOURCES
    plot_histogram.cpp)

ADD_LIBRARY(plot_histogram_lib STATIC ${LIB_SOURCES})

TARGET_INCLUDE_DIRECTORIES(plot_histogram_lib
                        PUBLIC
                            ${PROJECT_SOURCE_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}
                            ${ROOT_INCLUDE_DIRS})