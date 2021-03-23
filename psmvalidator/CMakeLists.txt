cmake_minimum_required(VERSION 3.2)
project(psmvalidator)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wreturn-type  -fopenmp -fPIC -fpermissive -Werror=return-type ")
# message(PSMValidator "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS} AND HE BOOST DIR IS  ${Boost_INCLUDE_DIR}")


set(SOURCE_FILES
        dependentlibs/msbasic/BasePSMFeature.cpp
        dependentlibs/msbasic/CDebugMode.cpp
        dependentlibs/msbasic/CFlow.cpp
        dependentlibs/msbasic/CFragScore.cpp
        dependentlibs/msbasic/classifier.cpp # 5
        dependentlibs/msbasic/commonfun.cpp
        dependentlibs/msbasic/ConcretePSMFeatures.cpp
        dependentlibs/msbasic/CPeakPair.cpp
        dependentlibs/msbasic/CPeakPairsImporter.cpp
        dependentlibs/msbasic/FragmentationScore.cpp # 10
        dependentlibs/msbasic/gnuplot_functions.cpp
        dependentlibs/msbasic/liblinearWrapper.cpp
        dependentlibs/msbasic/MGFReader.cpp
        dependentlibs/msbasic/peak_annotation.cpp
        dependentlibs/msbasic/PeakInfo.cpp # 15
        dependentlibs/msbasic/PeakPairScore.cpp
        dependentlibs/msbasic/PSMScoreCalculator.cpp

        ../librarymsms/CThreadsPool.cpp

        ../librarymsms/SpectraST_ramp.cpp # SpectraST_cramp requires
            ../librarymsms/SpectraST_cramp.cpp # proteomicsdatatypes requires

        ../librarymsms/ICGtInfoUpdate.cpp # proteomicdatatypes requires
            ../librarymsms/ProteomicsDataTypes.cpp #CFlow Requires

        ../librarymsms/PeakList.cpp # XMLFileParser requires
            ../librarymsms/XMLFileParser.cpp # CFlow requires

        ../librarymsms/DatabaseManager.cpp # Util requires
            ../librarymsms/Util.cpp
        ../librarymsms/CPSMAnnotation.cpp
        )

set(Boost_USE_STATIC_LIBS ON)

message("Boost version is:  ${Boost_LIB_VERSION} and the libraries are  ${Boost_LIBRARIES}")

add_executable(psmvalidator ${SOURCE_FILES} psmvalidator_main.cpp )
add_dependencies(psmvalidator boost)
include_directories(../External/spdlog-1.x/include )

target_link_libraries(psmvalidator
        ${Boost_LIBRARIES}
        pthread
        spectrastlib_static
        linear_static
        sqlite3
        )

target_link_directories(psmvalidator PRIVATE ${CMAKE_BINARY_DIR}/Release)

add_dependencies(psmvalidator ranger comet   linear_static  spdlog boost)

add_subdirectory(dependentlibs/msbasic/test)

install(TARGETS psmvalidator ranger comet  DESTINATION bin)
install(FILES psmvalidator.conf DESTINATION param)
install(PROGRAMS train_small_demo_mgf.bash validate_psm.bash DESTINATION scripts)
install(FILES model/nist_human_hcd_selected_features.txtmtry_8_ntree_900_trN_20000.forest model/nist_human_hcd_selected_features.txt_2.000000_frag.model DESTINATION model)

install(FILES
        model/interact-01625b_GD2-TUM_first_pool_12_01_01-3xHCD-1h-R1.ipro.pep.xml
        model/01625b_GD2-TUM_first_pool_12_01_01-3xHCD-1h-R1.mzML
        model/TUM_first_pool_12.pep.list
        DESTINATION predict)
install(FILES
        train/comet16low.param
        train/human_hcd_selected_new_small.mgf

       # human database 2020
       train/uniprot-human-2020-12.fasta
       train/uniprot-human-2020-12TD_only_decoy.fasta

       # yeast database
       train/uniprot_yeast_reviewed_6721_Nov152016TD_only_decoy.fasta
       train/uniprot_yeast_reviewed_6721_Nov152016.fasta

       train/featurelist.txt
       train/features.txt
        DESTINATION train)

include(InstallRequiredSystemLibraries)
include(CPack)

