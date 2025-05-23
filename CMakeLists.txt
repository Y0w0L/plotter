CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)

PROJECT(MYPROJECT VERSION 1.0 LANGUAGES CXX)

## Find ROOT
FIND_PACKAGE(ROOT REQUIRED COMPONENTS)
IF(NOT ROOT_FOUND)
    MESSAGE(FATAL_ERROR "Could not find ROOT, check the ROOT path or install\n"
                        "$ source YOUR_ROOT_DIR/bin/thisroot.sh")
ENDIF()

# Check ROOT version
IF(NOT ${ROOT_VERSION} VERSION_GREATER "6.0")
    MESSAGE(FATAL_ERROR "ROOT version below 6.0 are not supported")
ENDIF()

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF)

# Set global build flags
SET(SOMMON_CXX_FLAGS
    -pedantic
    -Wall
    -Wextra
    -Wcast-align
    -Wcast-qual
    -Wconversion
    -Wuseless-cast
    -Wold-style-cast
    -Wsign-conversion
    -Wctor-dtor-privacy
    -Wnoexcept
    -Woverloaded-virtual
    -Wstrict-null-sentinel
    -Wzero-as-null-pointer-constant
    -Wdisabled-optimization
    -Wformat=2
    -Wformat-security
    -Winit-self
    -Wlogical-op
    -Wmissing-declaretions
    -Wmissing-include-dirs
    -Wredundant-decls
    -Wsign-promo
    -Wswitch-default
    -Wundef
    -Wshadow
    -Wdeprecated
    -Wheader-hygiene
    -fdiagnostics-color=auto)

LIST(APPEND CMAKE_CXX_FLAGS "${COMMON_CXX_FLAGS}")

SET(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -pg")
SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -s -DNDEBUG -march=native")
SET(CMAKE_CXX_FLAGS_RELWITAHDEBINFO "-g3 -Og -pg")
SET(CMAKE_CXX_FLAGS_MINISIZEREL "-Os -s -DNDEBUG -march=native")

# Add Subdirectory
ADD_SUBDIRECTORY(plot_histogram)
ADD_SUBDIRECTORY(plot_simulation)
ADD_SUBDIRECTORY(track_resolution)
ADD_SUBDIRECTORY(plot_MobilityModel)
ADD_SUBDIRECTORY(plot_ComparingVoltage)
ADD_SUBDIRECTORY(Messenger)

ADD_EXECUTABLE(ce65_model
    ce65.cpp)

SET_TARGET_PROPERTIES(ce65_model PROPERTIES
                    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")

# Add root include directory to target
TARGET_INCLUDE_DIRECTORIES(ce65_model 
                    PUBLIC
                        ${PROJECT_SOURCE_DIR}
                        ${ROOT_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(ce65_model
                    PUBLIC
                        plot_histogram_lib
                        plot_simulation_lib
                        track_resolution_lib
                        plot_MobilityModel_lib
                        plot_ComparingVoltage_lib
                        Messenger_lib
                        ${ROOT_LIBRARIES})

MESSAGE(STATUS "Compiler Flags for ce65_model: $<TARGET_PROPERTY:ce65_model,COMPILE_OPTIONS>")