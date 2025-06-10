# plot_inPixel/CMakeLists.txt

SET(LIB_SOURCES
    plot_inPixel.cpp)

ADD_LIBRARY(plot_inPixel_lib STATIC ${LIB_SOURCES})

TARGET_INCLUDE_DIRECTORIES(plot_inPixel_lib
                        PUBLIC
                            ${PROJECT_SOURCE_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}
                            ${ROOT_INCLUDE_DIRS})