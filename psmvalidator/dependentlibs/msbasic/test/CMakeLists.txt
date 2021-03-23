cmake_minimum_required(VERSION 3.2)
project(testMsBasic)

# tried many different method. does not working...

set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS}  --std=c++11 -fpermissive -Werror=return-type ")

include_directories(../../../../External/gtest-1.7.0/include)
include_directories(../../../../External/spdlog-1.x/include)
include_directories(../../../../psmvalidator/dependentlibs/msbasic/)

set(majorPath ${CMAKE_CURRENT_SOURCE_DIR}/../../../../psmvalidator)
MESSAGE(${majorPath})

add_executable(testCFlow ../testCFlow.cpp ../CFlow.cpp
        ../BasePSMFeature.cpp
        ../CDebugMode.cpp
        ../CFragScore.cpp
        ../classifier.cpp
        ../commonfun.cpp
        ../ConcretePSMFeatures.cpp
        ../CPeakPair.cpp
        ../CPeakPairsImporter.cpp
        ../FragmentationScore.cpp
        ../gnuplot_functions.cpp
        ../liblinearWrapper.cpp
        ../MGFReader.cpp
        ../PeakInfo.cpp
        ../PeakPairScore.cpp
        ../peak_annotation.cpp
        ../PSMScoreCalculator.cpp
        ../CFlow.cpp
        ${majorPath}/../librarymsms/CThreadsPool.cpp

        ${majorPath}/../librarymsms/SpectraST_ramp.cpp # SpectraST_cramp requires
        ${majorPath}/../librarymsms/SpectraST_cramp.cpp # proteomicsdatatypes requires

        ${majorPath}/../librarymsms/ICGtInfoUpdate.cpp # proteomicdatatypes requires
        ${majorPath}/../librarymsms/ProteomicsDataTypes.cpp #CFlow Requires

        ${majorPath}/../librarymsms/PeakList.cpp # XMLFileParser requires
        ${majorPath}/../librarymsms/XMLFileParser.cpp # CFlow requires

        ${majorPath}/../librarymsms/DatabaseManager.cpp # Util requires
        ${majorPath}/../librarymsms/Util.cpp
        ${majorPath}/../librarymsms/CPSMAnnotation.cpp
        )
message(${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(testCFlow
        pthread
        sqlite3
        spectrastlib_static
        ${CMAKE_CURRENT_SOURCE_DIR}/../../../../Release/Release/liblinear_static.a
        blas
        z
        dl
        gsl gslcblas
        gtest gtest_main spdlog
        boost_system
        boost_filesystem
        boost_iostreams
        boost_program_options
        boost_timer
        boost_chrono
        )
target_include_directories(testCFlow PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../../../Release/Release)
target_link_directories(testCFlow PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../../../Release/Release)
