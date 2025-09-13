# plot_histogram/CMakeLists.txt

SET(LIB_SOURCES
    plot_AllpixData.cpp)

ADD_LIBRARY(plot_AllpixData_lib STATIC ${LIB_SOURCES})

TARGET_INCLUDE_DIRECTORIES(plot_AllpixData_lib
                        PUBLIC
                            ${PROJECT_SOURCE_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}
                            ${ROOT_INCLUDE_DIRS}
                            /home/towa/package/allpix/allpix_dev
)

# TARGET_INCLUDE_LIBRARIES(plot_AllpixData_lib
#                         PUBLIC
#                             Allpix::Objects
#                             ROOT::Core)