#include "Coloquinte/circuit_helper.hxx"

#include <set>

namespace coloquinte{
namespace gp{

void add_force(pin_1D const p1, pin_1D const p2, linear_system & L, float_t force){
    if(p1.movable && p2.movable){
        L.add_force(
            force,
            p1.cell_ind, p2.cell_ind,
            p1.offs,     p2.offs
        );
    }
    else if(p1.movable){
        L.add_fixed_force(
            force,
            p1.cell_ind,
            p2.pos,
            p1.offs
        );
    }
    else if(p2.movable){
        L.add_fixed_force(
            force,
            p2.cell_ind,
            p1.pos,
            p2.offs
        );
    }
}

void add_force(pin_1D const p1, pin_1D const p2, linear_system & L, float_t tol, float_t scale){
    add_force(p1, p2, L, scale/std::max(tol, std::abs(p2.pos-p1.pos)));
}

point<linear_system> empty_linear_systems(netlist const & circuit, placement_t const & pl){
    point<linear_system> ret = point<linear_system>(linear_system(circuit.cell_cnt()), linear_system(circuit.cell_cnt()));

    for(index_t i=0; i<circuit.cell_cnt(); ++i){
        if( (XMovable & circuit.get_cell(i).attributes) == 0){
            ret.x_.add_triplet(i, i, 1.0);
            ret.x_.add_doublet(i, pl.positions_[i].x_);
        }
        if( (YMovable & circuit.get_cell(i).attributes) == 0){
            ret.y_.add_triplet(i, i, 1.0);
            ret.y_.add_doublet(i, pl.positions_[i].y_);
        }
    }

    return ret;
}

namespace{ // Anonymous namespace for helper functions

void get_HPWLF(std::vector<pin_1D> const & pins, linear_system & L, float_t tol){
    if(pins.size() >= 2){
        auto min_elt = std::min_element(pins.begin(), pins.end()), max_elt = std::max_element(pins.begin(), pins.end());

        for(auto it = pins.begin(); it != pins.end(); ++it){
            // Just comparing the iterator is poorer due to redundancies in the benchmarks!
            if(it != min_elt){
                add_force(*it, *min_elt, L, tol, 1.0/(pins.size()-1));
                if(it != max_elt){ // Hopefully only one connexion between the min and max pins
                    add_force(*it, *max_elt, L, tol, 1.0/(pins.size()-1));
                }
            }
        }
    }
}

void get_HPWLR(std::vector<pin_1D> const & pins, linear_system & L, float_t tol){
    std::vector<pin_1D> sorted_pins = pins;
    std::sort(sorted_pins.begin(), sorted_pins.end());
    // Pins are connected to the pin two places away
    for(index_t i=0; i+2<sorted_pins.size(); ++i){
        add_force(sorted_pins[i], sorted_pins[i+2], L, tol, 0.5);
    }
    // The extreme pins are connected with their direct neighbour too
    if(sorted_pins.size() > 1){
        add_force(sorted_pins[0], sorted_pins[1], L, tol, 0.5);
        add_force(sorted_pins[sorted_pins.size()-1], sorted_pins[sorted_pins.size()-2], L, tol, 0.5);
    }
}

void get_star(std::vector<pin_1D> const & pins, linear_system & L, float_t tol, index_t star_index){
    // The net is empty, but we still populate the diagonal to avoid divide by zeros
    if(pins.size() < 2){
        L.add_triplet(star_index, star_index, 1.0);
        return;
    }

    for(pin_1D p : pins){
        pin_1D star_pin = pin_1D(star_index, std::numeric_limits<float_t>::quiet_NaN(), 0.0, true);
        add_force(p, star_pin, L, 1.0/pins.size());
    }
}

void get_clique(std::vector<pin_1D> const & pins, linear_system & L, float_t tol){
    // Pins are connected to the pin two places away
    for(index_t i=0; i+1<pins.size(); ++i){
        for(index_t j=i+1; j<pins.size(); ++j){
            add_force(pins[i], pins[j], L, tol, 1.0/(pins.size()-1));
        }
    }
}

inline void northeast_octant_neighbours(std::vector<pin_2D> pins, std::vector<std::pair<index_t, index_t> > & edges){
    struct indexed_pt : point<float_t>{
        index_t index;
        indexed_pt(point<float_t> pt, index_t pos) : point<float_t>(pt), index(pos) {}
    };

    std::vector<indexed_pt> point_list;
    for(index_t i=0; i<pins.size(); ++i){
        point_list.push_back(indexed_pt(pins[i].pos, i));
    }

    std::sort(point_list.begin(), point_list.end(),
                [](indexed_pt const a, indexed_pt const b){ return a.x_ + a.y_ < b.x_ + b.y_; }
              );

    // Decreasing order of x and y; multiset not necessary because no two elements have same coordinate
    std::set<indexed_pt, std::function<bool (indexed_pt const, indexed_pt const)> >
                      active_upper_octant([](indexed_pt const a, indexed_pt const b)->bool{return a.x_ > b.x_;}),
                      active_lower_octant([](indexed_pt const a, indexed_pt const b)->bool{return a.y_ > b.y_;});

    for(indexed_pt const current : point_list){
        { // North to north-east region
            auto first_it = active_upper_octant.lower_bound(current); // Largest x with x <= current.x
            auto it = first_it;
            for(; it != active_upper_octant.end() && it->x_ - it->y_ >= current.x_ - current.y_; ++it){
                edges.push_back(std::pair<index_t, index_t>(current.index, it->index));
            }
            if(first_it != active_upper_octant.end()){ active_upper_octant.erase(first_it, it); }
            active_upper_octant.insert(it, current); // Hint to insert the element since it is the correct position
        } // End region
        { // North-east to east region
            auto first_it = active_lower_octant.lower_bound(current); // Largest y with y <= current.y
            auto it = first_it;
            for(; it != active_lower_octant.end() && it->y_ - it->x_ >= current.y_ - current.x_; ++it){
                edges.push_back(std::pair<index_t, index_t>(current.index, it->index));
            }
            if(first_it != active_lower_octant.end()){ active_lower_octant.erase(first_it, it); }
            active_lower_octant.insert(it, current); // Hint to insert the element since it is the correct position
        } // End region
    }
}

// Gets the nearest octant neighbour for each point in the south-east quadrant
inline void southeast_octant_neighbours(std::vector<pin_2D> pins, std::vector<std::pair<index_t, index_t> > & edges){
    for(auto & pin : pins){
        pin.pos.y_ = - pin.pos.y_;
    }
    northeast_octant_neighbours(pins, edges);
}

std::vector<std::pair<index_t, index_t> > get_spanning_tree(std::vector<pin_2D> const & pins){
    typedef std::pair<index_t, index_t> edge_t;

	std::vector<edge_t> edges;
    
    if(pins.size() <= 2){
        if(pins.size() == 2){
            edges.push_back(edge_t(0, 1));
        }
        if(pins.size() == 3){
            auto dists = std::array<float_t, 3>({dist(pins[1], pins[2]), dist(pins[1], pins[2]), dist(pins[0], pins[1])});
            index_t mx = std::max_element(dists.begin(), dists.end()) - dists.begin();
            for(index_t i=0; i<3; ++i){
                if(i != mx)
                    edges.push_back(edge_t((i+1) % 3, (i+2) % 3));
            }
        }
        return edges;
    }
    
    northeast_octant_neighbours(pins, edges);
    southeast_octant_neighbours(pins, edges);

	std::vector<edge_t> returned_edges;

    class union_find{
    	index_t* connex_representants;
    	index_t  sz;
    
    	public:
    	void merge(index_t a, index_t b){
    		connex_representants[find(a)] = b;
    	}
    
    	index_t find(index_t ind){
    		if(connex_representants[ind] != ind){
    			connex_representants[ind] = find(connex_representants[ind]);
    		}
    		return connex_representants[ind];
    	}
    
    	union_find(index_t s){
    		sz = s;
    		connex_representants = new index_t[size()];
    		for(index_t i=0; i<size(); ++i){
    			connex_representants[i] = i;
    		}
    	}
    	
    	~union_find(){
    		delete connex_representants;
    	}
    
        bool is_connex(){
            bool connex = true;
            for(index_t i=0; i+1<size(); ++i){
                connex = connex && (find(i) == find(i+1));
            }
            return connex;
        }
    	
    	index_t size() const { return sz; }
    };

    auto edge_length = [&](edge_t E){
        point<float_t> p1 = pins[E.first].pos,
                       p2 = pins[E.second].pos;
        return std::abs(p1.x_ - p2.x_) + std::abs(p1.y_ - p2.y_);
    };
	// Perform Kruskal to get the tree
	std::sort(edges.begin(), edges.end(), [&](edge_t a, edge_t b){ return edge_length(a) < edge_length(b); });

	union_find merger(pins.size());

	for(index_t i=0; i<edges.size() && returned_edges.size()+1 < pins.size(); ++i){
		edge_t E = edges[i];
		if(merger.find(E.first) != merger.find(E.second)){
			merger.merge(E.first, E.second);
            assert(merger.find(E.first) == merger.find(E.second));
			returned_edges.push_back(E);
		}
	}
	assert(returned_edges.size() + 1 == pins.size());
    assert(merger.is_connex());
	return returned_edges;
}

} // End anonymous namespace

point<linear_system> get_HPWLF_linear_system (netlist const & circuit, placement_t const & pl, float_t tol, index_t min_s, index_t max_s){
    point<linear_system> L = empty_linear_systems(circuit, pl);
    for(index_t i=0; i<circuit.net_cnt(); ++i){
        // Has the net the right pin count?
        index_t pin_cnt = circuit.get_net(i).pin_cnt;
        if(pin_cnt < min_s or pin_cnt >= max_s) continue;

        auto pins = get_pins_1D(circuit, pl, i);
        get_HPWLF(pins.x_, L.x_, tol);
        get_HPWLF(pins.y_, L.y_, tol);
    }
    return L;
}

point<linear_system> get_HPWLR_linear_system (netlist const & circuit, placement_t const & pl, float_t tol, index_t min_s, index_t max_s){
    point<linear_system> L = empty_linear_systems(circuit, pl);
    for(index_t i=0; i<circuit.net_cnt(); ++i){
        // Has the net the right pin count?
        index_t pin_cnt = circuit.get_net(i).pin_cnt;
        if(pin_cnt < min_s or pin_cnt >= max_s) continue;

        auto pins = get_pins_1D(circuit, pl, i);
        get_HPWLR(pins.x_, L.x_, tol);
        get_HPWLR(pins.y_, L.y_, tol);
    }
    return L;
}

point<linear_system> get_star_linear_system  (netlist const & circuit, placement_t const & pl, float_t tol, index_t min_s, index_t max_s){
    point<linear_system> L = empty_linear_systems(circuit, pl);
    L.x_.add_variables(circuit.net_cnt());
    L.y_.add_variables(circuit.net_cnt());
    for(index_t i=0; i<circuit.net_cnt(); ++i){
        // Has the net the right pin count?
        index_t pin_cnt = circuit.get_net(i).pin_cnt;
        if(pin_cnt < min_s or pin_cnt >= max_s){
            // Put a one in the intermediate variable in order to avoid non-invertible matrices
            L.x_.add_triplet(i+circuit.cell_cnt(), i+circuit.cell_cnt(), 1.0);
            L.y_.add_triplet(i+circuit.cell_cnt(), i+circuit.cell_cnt(), 1.0);
            continue;
        }

        auto pins = get_pins_1D(circuit, pl, i);
        // Provide the index of the star's central pin in the linear system
        get_star(pins.x_, L.x_, tol, i+circuit.cell_cnt());
        get_star(pins.y_, L.y_, tol, i+circuit.cell_cnt());
    }
    return L;
}

point<linear_system> get_clique_linear_system (netlist const & circuit, placement_t const & pl, float_t tol, index_t min_s, index_t max_s){
    point<linear_system> L = empty_linear_systems(circuit, pl);
    for(index_t i=0; i<circuit.net_cnt(); ++i){
        // Has the net the right pin count?
        index_t pin_cnt = circuit.get_net(i).pin_cnt;
        if(pin_cnt < min_s or pin_cnt >= max_s) continue;

        auto pins = get_pins_1D(circuit, pl, i);
        get_clique(pins.x_, L.x_, tol);
        get_clique(pins.y_, L.y_, tol);
    }
    return L;
}

point<linear_system> get_MST_linear_system(netlist const & circuit, placement_t const & pl, float_t tol, index_t min_s, index_t max_s){
    point<linear_system> L = empty_linear_systems(circuit, pl);
    for(index_t i=0; i<circuit.net_cnt(); ++i){
        // Has the net the right pin count?
        index_t pin_cnt = circuit.get_net(i).pin_cnt;
        if(pin_cnt < min_s or pin_cnt >= max_s or pin_cnt <= 1) continue;
            
        auto pins = get_pins_2D(circuit, pl, i);
        auto edges = get_spanning_tree(pins);
        for(auto E : edges){
            add_force(pins[E.first].x(), pins[E.second].x(), L.x_, tol, 1.0);
            add_force(pins[E.first].y(), pins[E.second].y(), L.y_, tol, 1.0);
        }
    }
    return L;
}

float_t get_HPWL_wirelength(netlist const & circuit, placement_t const & pl){
    float_t sum = 0.0;
    for(index_t i=0; i<circuit.net_cnt(); ++i){
        if(circuit.get_net(i).pin_cnt <= 1) continue;

        auto pins = get_pins_1D(circuit, pl, i);
        auto minmaxX = std::minmax_element(pins.x_.begin(), pins.x_.end()), minmaxY = std::minmax_element(pins.y_.begin(), pins.y_.end());
        sum += ((minmaxX.second->pos - minmaxX.first->pos) + (minmaxY.second->pos - minmaxY.first->pos));
    }
    return sum;
}

// The true wirelength with minimum spanning trees, except for very small nets (<= 3) where we have HPWL == true WL
float_t get_MST_wirelength(netlist const & circuit, placement_t const & pl){
    float_t sum = 0.0;
    for(index_t i=0; i<circuit.net_cnt(); ++i){
        if(circuit.get_net(i).pin_cnt <= 1) continue;

        if(circuit.get_net(i).pin_cnt <= 3){
            auto pins = get_pins_1D(circuit, pl, i);
            auto minmaxX = std::minmax_element(pins.x_.begin(), pins.x_.end()), minmaxY = std::minmax_element(pins.y_.begin(), pins.y_.end());
            sum += ((minmaxX.second->pos - minmaxX.first->pos) + (minmaxY.second->pos - minmaxY.first->pos));
        }
        else{
            auto pins = get_pins_2D(circuit, pl, i);
            auto edges = get_spanning_tree(pins);
            for(auto E : edges){
                sum += std::abs(pins[E.first].pos.x_ - pins[E.second].pos.x_);
                sum += std::abs(pins[E.first].pos.y_ - pins[E.second].pos.y_);
            }
        }
    }
    return sum;
}

void get_result(netlist const & circuit, placement_t & pl, point<linear_system> & L, float_t tol){
    std::vector<float_t> x_sol, y_sol;
    std::vector<float_t> x_guess(pl.cell_cnt()), y_guess(pl.cell_cnt());
    
    assert(L.x_.internal_size() == x_guess.size());
    assert(L.y_.internal_size() == y_guess.size());

    for(index_t i=0; i<pl.cell_cnt(); ++i){
        x_guess[i] = pl.positions_[i].x_;
        y_guess[i] = pl.positions_[i].y_;
    }
    #pragma omp parallel sections num_threads(2)
    {
    #pragma omp section
    x_sol = L.x_.solve_CG(x_guess, tol);
    #pragma omp section
    y_sol = L.y_.solve_CG(y_guess, tol);
    }
    for(index_t i=0; i<pl.cell_cnt(); ++i){
        if( (circuit.get_cell(i).attributes & XMovable) != 0){
            pl.positions_[i].x_ = x_sol[i];
        }
        if( (circuit.get_cell(i).attributes & YMovable) != 0){
            pl.positions_[i].y_ = y_sol[i];
        }
    }
}

// Intended to be used by pulling forces to adapt the forces to the cell's areas
std::vector<float_t> get_area_scales(netlist const & circuit){
    std::vector<float_t> ret(circuit.cell_cnt());
    capacity_t int_tot_area = 0;
    for(index_t i=0; i<circuit.cell_cnt(); ++i){
        capacity_t A = circuit.get_cell(i).area;
        ret[i] = static_cast<float_t>(A);
        int_tot_area += A;
    }
    float_t average_area = static_cast<float_t>(int_tot_area) / circuit.cell_cnt();
    for(index_t i=0; i<circuit.cell_cnt(); ++i){
        ret[i] /= average_area;
    }
    return ret;
}

point<linear_system> get_pulling_forces (netlist const & circuit, placement_t const & pl, float_t typical_distance){
    point<linear_system> L = empty_linear_systems(circuit, pl);
    float_t typical_force = 1.0 / typical_distance;
    std::vector<float_t> scaling = get_area_scales(circuit);
    for(index_t i=0; i<pl.cell_cnt(); ++i){
        L.x_.add_anchor(
            typical_force * scaling[i],
            i, pl.positions_[i].x_
        );
        L.y_.add_anchor(
            typical_force * scaling[i],
            i, pl.positions_[i].y_
        );
    }
    
    return L;
}

point<linear_system> get_linear_pulling_forces (netlist const & circuit, placement_t const & UB_pl, placement_t const & LB_pl, float_t force, float_t min_distance){
    point<linear_system> L = empty_linear_systems(circuit, UB_pl);
    assert(LB_pl.cell_cnt() == UB_pl.cell_cnt());
    std::vector<float_t> scaling = get_area_scales(circuit);
    for(index_t i=0; i<LB_pl.cell_cnt(); ++i){
        L.x_.add_anchor(
            force * scaling[i] / (std::max(std::abs(UB_pl.positions_[i].x_ - LB_pl.positions_[i].x_), min_distance)),
            i, UB_pl.positions_[i].x_
        );
        L.y_.add_anchor(
            force * scaling[i] / (std::max(std::abs(UB_pl.positions_[i].y_ - LB_pl.positions_[i].y_), min_distance)),
            i, UB_pl.positions_[i].y_
        );
    }
    

    return L;
}

region_distribution get_rough_legalizer(netlist const & circuit, placement_t const & pl, box<int_t> surface){
    std::vector<region_distribution::movable_cell> movable_cells;
    std::vector<region_distribution::fixed_cell>   fixed_cells;

    for(index_t i=0; i<circuit.cell_cnt(); ++i){
        auto C = circuit.get_cell(i);
        if((C.attributes & (XMovable|YMovable)) != 0){
            movable_cells.push_back(region_distribution::movable_cell(C.area, pl.positions_[i], i));
        }
        else{
            fixed_cells.push_back(region_distribution::fixed_cell(C.size, pl.positions_[i]));
        }
    }

    return region_distribution(surface, movable_cells, fixed_cells);
}

void get_result(netlist const & circuit, placement_t & pl, region_distribution const & legalizer){
    auto exportation = legalizer.export_spread_positions_linear();
    for(auto const C : exportation){
        pl.positions_[C.index_in_placement_] = C.pos_;
    }
}

float_t get_mean_linear_disruption(netlist const & circuit, placement_t const & LB_pl, placement_t const & UB_pl){
    float_t tot_cost = 0.0;
    float_t tot_area = 0.0;
    for(index_t i=0; i<circuit.cell_cnt(); ++i){
        float_t area = static_cast<float_t>(circuit.get_cell(i).area);
        point<float_t> diff = LB_pl.positions_[i] - UB_pl.positions_[i];

        if( (circuit.get_cell(i).attributes & XMovable) == 0.0) assert(diff.x_ == 0.0);
        if( (circuit.get_cell(i).attributes & YMovable) == 0.0) assert(diff.y_ == 0.0);

        tot_cost += area * (std::abs(diff.x_) + std::abs(diff.y_));
        tot_area += area;
    }
    return tot_cost / tot_area;
}

float_t get_mean_quadratic_disruption(netlist const & circuit, placement_t const & LB_pl, placement_t const & UB_pl){
    float_t tot_cost = 0.0;
    float_t tot_area = 0.0;
    for(index_t i=0; i<circuit.cell_cnt(); ++i){
        float_t area = static_cast<float_t>(circuit.get_cell(i).area);
        point<float_t> diff = LB_pl.positions_[i] - UB_pl.positions_[i];

        if( (circuit.get_cell(i).attributes & XMovable) == 0.0) assert(diff.x_ == 0.0);
        if( (circuit.get_cell(i).attributes & YMovable) == 0.0) assert(diff.y_ == 0.0);

        float_t manhattan = (std::abs(diff.x_) + std::abs(diff.y_));
        tot_cost += area * manhattan * manhattan;
        tot_area += area;
    }
    return std::sqrt(tot_cost / tot_area);
}

} // namespace gp
} // namespace coloquinte


