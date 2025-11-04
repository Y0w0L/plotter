# plot_histogram/CMakeLists.txt

SET(LIB_SOURCES
    plot_Check.cpp)

ADD_LIBRARY(plot_Check_lib STATIC ${LIB_SOURCES})

TARGET_INCLUDE_DIRECTORIES(plot_Check_lib
                        PUBLIC
                            ${PROJECT_SOURCE_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}
                            ${ROOT_INCLUDE_DIRS})