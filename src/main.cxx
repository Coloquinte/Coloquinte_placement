
/*
 * This file is an example of how to interface Coloquinte in your own tool
 *
 * Create the objects netlist, placement_t and the box<int_t> representing the placement surface just as in input_stdin
 * Then run some optimization schedule as the one written below (you'll need to give the standard cell height)
 * Convert back your result
 */

#include "coloquinte/circuit.hxx"
#include "coloquinte/legalizer.hxx"

#include <iostream>
#include <vector>
#include <ctime>

using namespace coloquinte::gp;
using namespace coloquinte;

// Input of the circuit
void input_stdin(netlist & circuit, placement_t & pl, box<int_t> & surface){
    // Will need the lower left (integer) positions of the cells' and boolean orientations, and the placement surface
    // The orientations represent mirroring of the cells (centered)
    // The pins are given relative to the lower-left corner in the standard (true, true) orientation
    // The netlist object is constructed from objects temporary_cell, temporary_net and temporary_pin

    int_t x_min, x_max, y_min, y_max;

    std::cin >> x_min >> x_max >> y_min >> y_max;

    index_t cell_cnt, net_cnt;
    std::cin >> cell_cnt;

    std::vector<temporary_cell> cells(cell_cnt);
    std::vector<point<int_t> > positions(cell_cnt);
    for(index_t i=0; i<cell_cnt; ++i){
        index_t cell_ind;
        int_t x_s, y_s;
        std::cin >> cell_ind >> x_s >> y_s;
        assert(cell_ind == i);
        cells[i].size = point<int_t>(x_s, y_s);
        cells[i].list_index = cell_ind;
        cells[i].area = static_cast<capacity_t>(x_s) * static_cast<capacity_t>(y_s);
    }
    for(index_t i=0; i<cell_cnt; ++i){
        index_t cell_ind;
        int_t x, y, fixed;
        std::cin >> cell_ind >> x >> y >> fixed;
        assert(cell_ind == i);
        positions[i] = point<int_t>(x, y);
        if(fixed != 0){
            cells[i].attributes = 0;
        }
        else{
            cells[i].attributes = XMovable|YMovable|XFlippable|YFlippable;
        }
    }
    std::cin >> net_cnt;
    std::vector<temporary_net> nets(net_cnt);
    std::vector<temporary_pin> pins;
    for(index_t i=0; i<net_cnt; ++i){
        index_t pin_cnt;
        std::cin >> pin_cnt;
        nets[i] = temporary_net(i, 1);
        for(index_t j=0; j<pin_cnt; ++j){
            index_t cell_ind;
            float_t x, y;
            std::cin >> cell_ind >> x >> y;
            // The pin's information is only given in the temporary_pin structure: nets and cells create pointers to the pins during the construction
            pins.push_back(temporary_pin(point<int_t>(point<float_t>(x,y) + 0.5f * point<float_t>(cells[cell_ind].size)), cell_ind, i));
        }
    }
    std::vector<point<bool> > orientations(cell_cnt, point<bool>(true, true));

    surface = box<int_t>(x_min, x_max, y_min, y_max);
    pl.positions_ = positions;
    pl.orientations_ = orientations;
    circuit = netlist(cells, nets, pins);
}

void output_stdout(netlist const & circuit, placement_t const & pl, box<int_t> surface){
    std::cout << surface.x_min_ << " " << surface.x_max_ << std::endl;
    std::cout << surface.y_min_ << " " << surface.y_max_ << std::endl;
    std::cout << std::endl;
    for(index_t i=0; i<circuit.cell_cnt(); ++i){
        std::cout << pl.positions_[i].x_ << " " << pl.positions_[i].y_ << std::endl;
    }
}

void output_report(netlist const & circuit, placement_t const & LB_pl, placement_t const & UB_pl){
    std::cout << "HPWL: " << get_HPWL_wirelength(circuit, UB_pl);
    std::cout << "\tRSMT: " << get_RSMT_wirelength(circuit, UB_pl);
    //std::cout << "\tMST: " << get_MST_wirelength(circuit, UB_pl);
    std::cout << "\tTime: " << time(NULL) << "\tLinear D: " << get_mean_linear_disruption(circuit, LB_pl, UB_pl) << "\tQuad D: " << get_mean_quadratic_disruption(circuit, LB_pl, UB_pl) << std::endl;
    //std::cout << "HPWL: " << get_HPWL_wirelength(circuit, UB_pl) << "\tTime: " << time(NULL) << "\tLinear D: " << get_mean_linear_disruption(circuit, LB_pl, UB_pl) << "\tQuad D: " << get_mean_quadratic_disruption(circuit, LB_pl, UB_pl) << std::endl;
}
void output_report(netlist const & circuit, placement_t const & LB_pl){
    std::cout << "HPWL: " << get_HPWL_wirelength(circuit, LB_pl);
    std::cout << "\tRSMT: " << get_RSMT_wirelength(circuit, LB_pl);
    //std::cout << "\tMST: " << get_MST_wirelength(circuit, LB_pl);
    std::cout << "\tTime: " << time(NULL) << std::endl;
    //std::cout << "HPWL: " << get_HPWL_wirelength(circuit, LB_pl) << "\tTime: " << time(NULL) << std::endl;
}


int main(){
    box<int_t> surface;
    netlist circuit;
    placement_t LB_pl, UB_pl;
    input_stdin(circuit, LB_pl, surface);
    UB_pl = LB_pl;
    circuit.selfcheck();
    LB_pl.selfcheck();

    std::cout << "The initial wirelength is " << get_HPWL_wirelength(circuit, LB_pl) << " at " << time(NULL) << std::endl;
    
    auto first_legalizer = get_rough_legalizer(circuit, LB_pl, surface);
    first_legalizer.selfcheck();
    get_rough_legalization(circuit, UB_pl, first_legalizer);
    UB_pl.selfcheck();
    std::cout << "The simply legalized wirelength is " << get_HPWL_wirelength(circuit, UB_pl) << " at " << time(NULL) << " with linear disruption " << get_mean_linear_disruption(circuit, LB_pl, UB_pl) << " and quadratic disruption " << get_mean_quadratic_disruption(circuit, LB_pl, UB_pl) << std::endl;
    LB_pl = UB_pl;

    // Early topology-independent solution
    auto solv = get_star_linear_system(circuit, LB_pl, 1.0, 0, 10000)
            + get_pulling_forces(circuit, UB_pl, 1000000.0); // Big distance: doesn't pull strongly, but avoids problems if there are no fixed pins
    std::cout << "Star optimization at time " << time(NULL) << std::endl;
    solve_linear_system(circuit, LB_pl, solv, 200); // number of iterations=200
    output_report(circuit, LB_pl);

    float_t pulling_force = 0.01;
    for(int i=0; i<50; ++i, pulling_force += 0.03){
        // Create a legalizer and bipartition it until we have sufficient precision (~2 to 10 standard cell widths)
        auto legalizer = get_rough_legalizer(circuit, LB_pl, surface);
        for(int quad_part =0; 10u * (1u << (2*quad_part)) < circuit.cell_cnt(); quad_part++){ // Here, approximately 10 cells in each region
            legalizer.x_bipartition();
            legalizer.y_bipartition();
            legalizer.redo_diagonal_bipartitions();
            legalizer.redo_line_partitions();
            legalizer.redo_diagonal_bipartitions();
            legalizer.redo_line_partitions();
            legalizer.redo_diagonal_bipartitions();
            legalizer.selfcheck();
        }
        UB_pl = LB_pl; // Keep the orientation between LB and UB

        get_rough_legalization(circuit, UB_pl, legalizer);
        std::cout << "Roughly legalized" << std::endl;
        output_report(circuit, LB_pl, UB_pl);

        auto LEG = dp::legalize(circuit, UB_pl, surface, 12);
        dp::get_result(circuit, LEG, UB_pl);
        std::cout << "Legalized" << std::endl;
        output_report(circuit, LB_pl, UB_pl);

        // Get the system to optimize (tolerance, maximum and minimum pin counts) and the pulling forces (threshold distance)
        auto solv = get_HPWLF_linear_system(circuit, LB_pl, 0.01, 2, 100000)
            + get_linear_pulling_forces(circuit, UB_pl, LB_pl, pulling_force, 40.0);
        std::cout << "Got the linear system at time " << time(NULL) << std::endl;
        solve_linear_system(circuit, LB_pl, solv, 400); // number of iterations
        output_report(circuit, LB_pl);

        // Optimize orientation sometimes
        if(i%5 == 0){
            optimize_exact_orientations(circuit, LB_pl);
            std::cout << "Oriented" << std::endl;
            output_report(circuit, LB_pl);
        }
    }

    std::cout << "Now let's detailed place" << std::endl; 
    for(index_t i=0; i<2; ++i){
        optimize_exact_orientations(circuit, UB_pl);
        std::cout << "Oriented" << std::endl;
        output_report(circuit, LB_pl, UB_pl);

        auto LEG = dp::legalize(circuit, UB_pl, surface, 12);
        dp::get_result(circuit, LEG, UB_pl);
        std::cout << "Legalized" << std::endl;
        output_report(circuit, LB_pl, UB_pl);

        dp::swaps_global_HPWL(circuit, LEG, 3, 4);
        dp::get_result(circuit, LEG, UB_pl);
        std::cout << "Global swaps" << std::endl;
        output_report(circuit, LB_pl, UB_pl);

        //dp::OSRP_convex_HPWL(circuit, LEG);
        dp::OSRP_noncvx_HPWL(circuit, LEG);
        dp::get_result(circuit, LEG, UB_pl);
        std::cout << "Ordered row optimization" << std::endl;
        output_report(circuit, LB_pl, UB_pl);

        //dp::swaps_row_convex_HPWL(circuit, LEG, 4);
        dp::swaps_row_noncvx_HPWL(circuit, LEG, 4);
        dp::get_result(circuit, LEG, UB_pl);
        std::cout << "Local swaps" << std::endl;
        output_report(circuit, LB_pl, UB_pl);
    }
    return 0;
}

