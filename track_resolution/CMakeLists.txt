# track_resolution/CMakeLists.txt

SET(LIB_SOURCES
    track_resolution.cpp)

ADD_LIBRARY(track_resolution_lib STATIC ${LIB_SOURCES})

TARGET_INCLUDE_DIRECTORIES(track_resolution_lib
                        PUBLIC
                            ${PROJECT_SOURCE_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}
                            ${ROOT_INCLUDE_DIRS})