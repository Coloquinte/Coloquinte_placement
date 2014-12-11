# -*- explicit-buffer-name: "CMakeLists.txt<Coloquinte/src>" -*-

 project(Coloquinte)

 set(CMAKE_C_FLAGS_DEBUG     " -Wall -O2 -g"       CACHE STRING "C Compiler Debug options."   FORCE)
 set(CMAKE_C_FLAGS_RELEASE   " -Wall -O3 -fopenmp" CACHE STRING "C Compiler Release options." FORCE)
 set(CMAKE_CXX_FLAGS_DEBUG   " -Wall -Og -g"   CACHE STRING "C++ Compiler Debug options."   FORCE)
 set(CMAKE_CXX_FLAGS_RELEASE " -Wall -O3 -fopenmp" CACHE STRING "C++ Compiler Release options." FORCE)

 cmake_minimum_required(VERSION 2.4.0)

# list(INSERT CMAKE_MODULE_PATH 0 "${DESTDIR}/$ENV{COLOQUINTE_TOP}/share/cmake/Modules/")

 add_definitions(-std=c++11)

#   include_directories ( 
#                         ${Coloquinte_SOURCE_DIR}/src
#                       )
                   set ( includes          Coloquinte/circuit.hxx
                                           Coloquinte/circuit_helper.hxx
                                           Coloquinte/common.hxx
                                           Coloquinte/netlist.hxx
                                           Coloquinte/solvers.hxx
                                           Coloquinte/rough_legalizers.hxx
                                           Coloquinte/transportation.hxx
                                           Coloquinte/detailed.hxx
                       )	           
                   set ( cpps              circuit.cxx
                                           rough_legalizers.cxx
                                           solvers.cxx
                                           transportation.cxx
                                           orientation.cxx
                                           detailed.cxx
                                           legalizer.cxx
                       )
                   set ( coloquintecpps    main.cxx )
					   
            add_library ( coloquinte        ${cpps} )
            add_executable ( coloquinte.bin    ${coloquintecpps})
            target_link_libraries ( coloquinte.bin    coloquinte )
