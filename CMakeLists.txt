cmake_minimum_required(VERSION 3.2)
project(PSMValidator)


# after I update my project in to a newer version, I got some problem
# I find the solution from the following pages:
# http://stackoverflow.com/questions/6594796/how-do-i-make-cmake-output-into-a-bin-dir
# I just put everything in one folder, which makes it easy to find a file.
# Maybe later i will change it in the future
#set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)

set( boost_INSTALL ${CMAKE_BINARY_DIR}/External/boost1.65.1/build )
set( boost_INCLUDE_DIR ${boost_INSTALL}/include )
set( boost_LIB_DIR ${boost_INSTALL}/lib )
set(BOOST_ROOT ${boost_INSTALL})
# Note it is important that we put system after filesystem !!!
#find_package(Boost REQUIRED COMPONENTS system filesystem iostreams program_options timer chrono)
set(Boost_LIBRARIES
        ${boost_LIB_DIR}/libboost_system.a
        ${boost_LIB_DIR}/libboost_filesystem.a
        ${boost_LIB_DIR}/libboost_iostreams.a
        ${boost_LIB_DIR}/libboost_program_options.a
        ${boost_LIB_DIR}/libboost_timer.a
        ${boost_LIB_DIR}/libboost_chrono.a 
        )

add_subdirectory(External/boost1.65.1)
ExternalProject_Get_Property(boost install_dir)
include_directories(${install_dir}/include)

set(Boost_NO_SYSTEM_PATHS TRUE)
set(Boost_NO_BOOST_CMAKE TRUE)

add_subdirectory(External/comet_source_2016012)
add_subdirectory(External/gtest-1.7.0)
add_subdirectory(External/liblinear-2.11)
add_subdirectory(External/ranger)
add_subdirectory(External/spdlog-1.x)
add_subdirectory(External/SpectraST)
file(TOUCH_NOCREATE psmvalidator_main.cpp)
add_subdirectory(psmvalidator)
add_dependencies(psmvalidator boost ranger train predict)
add_dependencies(testCFlow boost)
#add_subdirectory(librarymsms)
include_directories(External/spdlog-1.x/include)

set(CMAKE_CXX_FLAGS " -Wreturn-type ${CMAKE_CXX_FLAGS} -std=c++11  ")

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
include(CPack)
