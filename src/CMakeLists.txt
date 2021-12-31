# -*- explicit-buffer-name: "CMakeLists.txt<Coloquinte/src>" -*-

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
                        checkers.cxx
                        rough_legalizers.cxx
                        solvers.cxx
                        optimization_subproblems.cxx
                        orientation.cxx
                        detailed.cxx
                        cell_swapping.cxx
                        MCF_opt.cxx
                        row_opt.cxx
                        topologies.cxx
                        lookup_table.cxx
                        legalizer.cxx
    )
set ( coloquintecpps    main.cxx )

# Make sure this points to where lemon is installed!
include_directories($ENV{HOME}/usr/lemon/include)

add_library ( coloquinte        ${cpps} )
add_executable ( coloquinte.bin    ${coloquintecpps})
target_link_libraries ( coloquinte.bin    coloquinte )
set_property(TARGET coloquinte PROPERTY CXX_STANDARD 11)

install( TARGETS coloquinte  DESTINATION lib${LIB_SUFFIX} )
install( FILES ${includes}   DESTINATION include/coloquinte ) 
