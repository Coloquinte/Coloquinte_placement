
#include "Coloquinte/circuit.hxx"
#include "Coloquinte/legalizer.hxx"

#include <iostream>
#include <vector>
#include <ctime>

using namespace coloquinte::gp;
using namespace coloquinte;

// Input of the circuit
void input_stdin(netlist & circuit, placement_t & pl, box<int_t> & surface){
    // Will need the surface to legalize, the positions of the cells' centers and orientations between -1.0 and 1.0
    // The netlist object is constructed from objects temporary_cell, temporary_net and temporary_pin
    // The pin's information is only given in the temporary_pin structure: nets and cells don't have pointers to the pins during the construction

    int_t x_min, x_max, y_min, y_max;

    std::cin >> x_min >> x_max >> y_min >> y_max;

    index_t cell_cnt, net_cnt;
    std::cin >> cell_cnt;

    std::vector<temporary_cell> cells(cell_cnt);
    std::vector<point<float_t> > positions(cell_cnt);
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
        positions[i] = point<float_t>(x, y) + static_cast<float_t>(0.5) * static_cast<point<float_t> >(cells[i].size);
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
        nets[i] = temporary_net(i, 1.0);
        for(index_t j=0; j<pin_cnt; ++j){
            index_t cell_ind;
            float_t x, y;
            std::cin >> cell_ind >> x >> y;
            pins.push_back(temporary_pin(point<float_t>(x,y), cell_ind, i));
        }
    }
    std::vector<point<float_t> > orientations(cell_cnt, point<float_t>(1.0, 1.0));

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


    std::cout << "The initial wirelength is " << get_HPWL_wirelength(circuit, LB_pl) << " at " << time(NULL) << std::endl;
    
    auto first_legalizer = get_rough_legalizer(circuit, LB_pl, surface);
    get_result(circuit, UB_pl, first_legalizer);
    std::cout << "The simply legalized wirelength is " << get_HPWL_wirelength(circuit, UB_pl) << " at " << time(NULL) << " with linear disruption " << get_mean_linear_disruption(circuit, LB_pl, UB_pl) << " and quadratic disruption " << get_mean_quadratic_disruption(circuit, LB_pl, UB_pl) << std::endl;
    LB_pl = UB_pl;

    // No orientation at the beginning
    zero_orientations(circuit, LB_pl);

    // Early topology-independent solution
    auto solv = get_star_linear_system(circuit, LB_pl, 1.0, 0, 10000);
    std::cout << "Star optimization at time " << time(NULL) << std::endl;
    get_result(circuit, LB_pl, solv, 200); // number of iterations=200
    output_report(circuit, LB_pl);
    
    for(int i=0; i<10; ++i){
        //auto solv = get_clique_linear_system(circuit, LB_pl, 1.0, 0, 10) + get_HPWLF_linear_system(circuit, LB_pl, 1.0, 10, 100000);
        auto solv = get_HPWLF_linear_system(circuit, LB_pl, 1.0, 2, 100000);
        std::cout << "Got the linear system at time " << time(NULL) << std::endl;
        get_result(circuit, LB_pl, solv, 300); // number of iterations
        output_report(circuit, LB_pl);
    }

    float_t pulling_force = 0.03;

    for(int i=0; i<50; ++i, pulling_force += 0.03){
        // Create a legalizer and bipartition it until we have sufficient precision (~2 to 10 standard cell widths)
        auto legalizer = get_rough_legalizer(circuit, LB_pl, surface);
        for(int quad_part =0; quad_part < 8; quad_part++){
            legalizer.x_bipartition();
            legalizer.y_bipartition();
            legalizer.redo_line_partitions();
            legalizer.redo_diagonal_bipartitions();
            legalizer.redo_line_partitions();
            legalizer.redo_diagonal_bipartitions();
            legalizer.selfcheck();
        }
        if(i < 10){
            spread_orientations(circuit, LB_pl);
        }
        UB_pl = LB_pl; // Keep the orientation between LB and UB

        get_result(circuit, UB_pl, legalizer);
        std::cout << "Roughly legalized" << std::endl;
        output_report(circuit, LB_pl, UB_pl);

        if(i >= 30){
            auto LEG = dp::legalize(circuit, UB_pl, surface, 12);
            dp::get_result(circuit, LEG, UB_pl);
            std::cout << "Legalized" << std::endl;
            output_report(circuit, LB_pl, UB_pl);
        }

        // Get the system to optimize (tolerance, maximum and minimum pin counts) and the pulling forces (threshold distance)
        auto solv = get_HPWLF_linear_system(circuit, LB_pl, 0.01, 2, 100000) + get_linear_pulling_forces(circuit, UB_pl, LB_pl, pulling_force, 40.0);
        std::cout << "Got the linear system at time " << time(NULL) << std::endl;
        get_result(circuit, LB_pl, solv, 400); // number of iterations
        output_report(circuit, LB_pl);

        // Optimize orientation sometimes
        if(i>=10 and i%5 == 0){
            optimize_exact_orientations(circuit, LB_pl);
            std::cout << "Oriented" << std::endl;
            output_report(circuit, LB_pl);
        }
    }

    std::cout << "Now let's detailed place" << std::endl; 
    for(index_t i=0; i<10; ++i){
        optimize_exact_orientations(circuit, UB_pl);
        std::cout << "Oriented" << std::endl;
        output_report(circuit, LB_pl, UB_pl);

        auto LEG = dp::legalize(circuit, UB_pl, surface, 12);
        dp::get_result(circuit, LEG, UB_pl);
        std::cout << "Legalized" << std::endl;
        output_report(circuit, LB_pl, UB_pl);

        dp::swaps_global(circuit, LEG, 3, 4);
        dp::get_result(circuit, LEG, UB_pl);
        std::cout << "Global swaps" << std::endl;
        output_report(circuit, LB_pl, UB_pl);

        dp::OSRP_convex(circuit, LEG);
        dp::get_result(circuit, LEG, UB_pl);
        std::cout << "Ordered row optimization" << std::endl;
        output_report(circuit, LB_pl, UB_pl);

        dp::swaps_row(circuit, LEG, 4);
        dp::get_result(circuit, LEG, UB_pl);
        std::cout << "Local swaps" << std::endl;
        output_report(circuit, LB_pl, UB_pl);
    }
    return 0;
}

