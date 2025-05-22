# plot_simulation/CMakeLists.txt

SET(LIB_SOURCES
    plot_simulation.cpp)

ADD_LIBRARY(plot_simulation_lib STATIC ${LIB_SOURCES})

TARGET_INCLUDE_DIRECTORIES(plot_simulation_lib
                        PUBLIC
                            ${PROJECT_SOURCE_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}
                            ${ROOT_INCLUDE_DIRS})