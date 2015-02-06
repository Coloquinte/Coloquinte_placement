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
                   set ( includes          coloquinte/circuit.hxx
                                           coloquinte/circuit_helper.hxx
                                           coloquinte/common.hxx
                                           coloquinte/netlist.hxx
                                           coloquinte/solvers.hxx
                                           coloquinte/rough_legalizers.hxx
                                           coloquinte/legalizer.hxx
                                           coloquinte/detailed.hxx
                                           coloquinte/topologies.hxx
                                           coloquinte/optimization_subproblems.hxx
                       )	           
                   set ( cpps              circuit.cxx
                                           rough_legalizers.cxx
                                           solvers.cxx
                                           optimization_subproblems.cxx
                                           orientation.cxx
                                           detailed.cxx
                                           topologies.cxx
                                           lookup_table.cxx
                                           legalizer.cxx
                       )
                   set ( coloquintecpps    main.cxx )
					   
            add_library ( coloquinte        ${cpps} )
            add_executable ( coloquinte.bin    ${coloquintecpps})
            target_link_libraries ( coloquinte.bin    coloquinte )

               install( TARGETS coloquinte  DESTINATION lib${LIB_SUFFIX} )
               install( FILES ${includes}   DESTINATION include/Coloquinte ) 
