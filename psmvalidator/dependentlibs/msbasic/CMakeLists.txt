cmake_minimum_required(VERSION 3.2)
project(msbasicproj)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wreturn-type  -fPIC -fpermissive -Werror=return-type ")

set(SOURCE_FILES
        BasePSMFeature.cpp
        CDebugMode.cpp
        CFragScore.cpp
        classifier.cpp
        commonfun.cpp
        ConcretePSMFeatures.cpp
        CPeakPair.cpp
        CPeakPairsImporter.cpp
        FragmentationScore.cpp
        gnuplot_functions.cpp
        liblinearWrapper.cpp
        MGFReader.cpp
        PeakInfo.cpp
        PeakPairScore.cpp
        peak_annotation.cpp
        PSMScoreCalculator.cpp
        CFlow.cpp
        )
include_directories(../../../External/spdlog-1.x/include)
link_directories(${CMAKE_BINARY_DIR}/Release)


find_package(Boost REQUIRED COMPONENTS iostreams system filesystem program_options timer)
message("MSBasic: boost found libraries" ${Boost_LIBRARIES})
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS} ../../../External/liblinear-2.11 ./)
    add_library(msbasic SHARED ${SOURCE_FILES})
    target_link_libraries(msbasic ${Boost_LIBRARIES} pthread
            msms
            spectrastlib
            linear)
    set_target_properties(msbasic PROPERTIES POSITION_INDEPENDENT_CODE ON)

    add_library(msbasic_static STATIC ${SOURCE_FILES})
    target_link_libraries(msbasic_static ${Boost_LIBRARIES} pthread msms_static
            spectrastlib_static
            linear)
    set_target_properties(msbasic_static PROPERTIES POSITION_INDEPENDENT_CODE ON)

endif ()

