# plot_MobilityModel/CMakeLists.txt

SET(LIB_SOURCES
    plot_MobilityModel.cpp)

ADD_LIBRARY(plot_MobilityModel_lib STATIC ${LIB_SOURCES})

TARGET_INCLUDE_DIRECTORIES(plot_MobilityModel_lib
                        PUBLIC
                            ${PROJECT_SOURCE_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}
                            ${ROOT_INCLUDE_DIRS})