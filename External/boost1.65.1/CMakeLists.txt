
include( ExternalProject )

set( boost_URL "https://udomain.dl.sourceforge.net/project/boost/boost/1.65.1/boost_1_65_1.tar.gz" )
set( boost_SHA1 "c066ac5c2f42fa2b870362c3c931ef73ffc6f24f" )
find_package(Python3 COMPONENTS Interpreter Development)
message( "--- the root of boost is  ${BOOST_ROOT} The python root is ${Python3_INCLUDE_DIRS} the install dir is : ${install_dir}/include")




ExternalProject_Add( boost
        PREFIX boost
        URL ${boost_URL}
        URL_HASH SHA1=${boost_SHA1}
        BUILD_IN_SOURCE 1
        DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/External/boost1.65.1/tmp
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/External/boost1.65.1/src
        CONFIGURE_COMMAND
        ./bootstrap.sh
        --with-libraries=filesystem
        --with-libraries=system
        --with-libraries=iostreams
        --with-libraries=program_options
        --with-libraries=timer
        --with-libraries=chrono
        --prefix=<INSTALL_DIR>
       --with-python=${Python3_INCLUDE_DIRS}
        BUILD_COMMAND
#        &&
        ./b2 install
        #        link=static variant=release threading=multi runtime-link=static
        -j 30
#        echo done
        INSTALL_COMMAND ""
        INSTALL_DIR ${boost_INSTALL}
        )


