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
