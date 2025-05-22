# plot_ComparigVoltage/CMakeLists.txt

SET(LIB_SOURCES
    plot_ComparingVoltage.cpp)

ADD_LIBRARY(plot_ComparingVoltage_lib STATIC ${LIB_SOURCES})

TARGET_INCLUDE_DIRECTORIES(plot_ComparingVoltage_lib
                        PUBLIC
                            ${PROJECT_SOURCE_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}
                            ${ROOT_INCLUDE_DIRS})