# Messenger/CMakeLists.txt

SET(LIB_SOURCES
    Messenger.cpp)

ADD_LIBRARY(Messenger_lib STATIC ${LIB_SOURCES})

TARGET_INCLUDE_DIRECTORIES(Messenger_lib
                        PUBLIC
                            ${PROJECT_SOURCE_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}
                            ${ROOT_INCLUDE_DIRS})