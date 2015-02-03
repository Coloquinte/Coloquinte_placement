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
                   set ( includes          src/coloquinte/circuit.hxx
                                           src/coloquinte/circuit_helper.hxx
                                           src/coloquinte/common.hxx
                                           src/coloquinte/netlist.hxx
                                           src/coloquinte/solvers.hxx
                                           src/coloquinte/rough_legalizers.hxx
                                           src/coloquinte/detailed.hxx
                                           src/coloquinte/topologies.hxx
                                           src/coloquinte/optimization_subproblems.hxx
                       )	           
                   set ( cpps              src/circuit.cxx
                                           src/rough_legalizers.cxx
                                           src/solvers.cxx
                                           src/optimization_subproblems.cxx
                                           src/orientation.cxx
                                           src/detailed.cxx
                                           src/topologies.cxx
                                           src/lookup_table.cxx
                                           src/legalizer.cxx
                       )
                   set ( coloquintecpps    src/main.cxx )
					   
            add_library ( coloquinte        ${cpps} )
            add_executable ( coloquinte.bin    ${coloquintecpps})
            target_link_libraries ( coloquinte.bin    coloquinte )
