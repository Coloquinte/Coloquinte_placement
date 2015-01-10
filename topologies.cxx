
#include "Coloquinte/topologies.hxx"

#include <array>
#include <algorithm>
#include <cassert>
#include <set>

namespace coloquinte{
namespace {

template<int n, int array_size>
float_t get_wirelength(std::vector<point<float_t> > const & pins, std::array<Hconnectivity<n>, array_size> const & lookups){
    std::array<point<float_t>, n> points;
    for(index_t i=0; i<n; ++i){
        points[i] = pins[i];
    }

    std::sort(points.begin(), points.end(), [](point<float_t> a , point<float_t> b){return a.x_ < b.x_; });
    float_t cost = std::numeric_limits<float_t>::infinity();
    for(auto const L : lookups){
        cost = std::min(cost, L.get_wirelength(points));
    }
    return cost;
}

template<int n, int array_size>
std::vector<std::pair<index_t, index_t> > get_topology(std::vector<point<float_t> > const & pins, std::array<Hconnectivity<n>, array_size> const & lookups){
    std::array<point<float_t>, n> points;
    for(index_t i=0; i<n; ++i){
        points[i] = pins[i];
    }

    std::sort(points.begin(), points.end(), [](point<float_t> a , point<float_t> b){return a.x_ < b.x_; });
    float_t cost = std::numeric_limits<float_t>::infinity();
    index_t ind  = std::numeric_limits<index_t>::max();
    for(index_t i=0; i<array_size; ++i){
        float_t this_cost = lookups[i].get_wirelength(points);
        if(this_cost < cost){
            cost = this_cost;
            ind = i;
        }
    }
    assert(ind < std::numeric_limits<index_t>::max());
    auto topo = lookups[ind].get_topology();
    return std::vector<std::pair<index_t, index_t> >(topo.begin(), topo.end());
}

inline void northeast_octant_neighbours(std::vector<point<float_t> > pins, std::vector<std::pair<index_t, index_t> > & edges){
    struct indexed_pt : point<float_t>{
        index_t index;
        indexed_pt(point<float_t> pt, index_t pos) : point<float_t>(pt), index(pos) {}
    };

    std::vector<indexed_pt> point_list;
    for(index_t i=0; i<pins.size(); ++i){
        point_list.push_back(indexed_pt(pins[i], i));
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
inline void southeast_octant_neighbours(std::vector<point<float_t> > pins, std::vector<std::pair<index_t, index_t> > & edges){
    for(auto & pin : pins){
        pin.y_ = - pin.y_;
    }
    northeast_octant_neighbours(pins, edges);
}

} // End anonymous namespace

std::vector<std::pair<index_t, index_t> > get_MST_topology(std::vector<point<float_t> > const & pins){
    typedef std::pair<index_t, index_t> edge_t;

	std::vector<edge_t> edges;
    
    if(pins.size() <= 2){
        if(pins.size() == 2){
            edges.push_back(edge_t(0, 1));
        }
        if(pins.size() == 3){
            auto D = [](point<float_t> a, point<float_t> b){ return std::abs(a.x_ - b.x_) + std::abs(a.y_ - b.y_); };
            auto dists = std::array<float_t, 3>({D(pins[1], pins[2]), D(pins[1], pins[2]), D(pins[0], pins[1])});
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
        point<float_t> p1 = pins[E.first],
                       p2 = pins[E.second];
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

float_t MST_length(std::vector<point<float_t> > const & pins){
    auto edges = get_MST_topology(pins);
    float_t sum = 0.0;
    for(auto E : edges){
        sum += std::abs(pins[E.first].x_ - pins[E.second].x_);
        sum += std::abs(pins[E.first].y_ - pins[E.second].y_);
    }
    return sum;
}

float_t RSMT_length(std::vector<point<float_t> > const & pins, index_t exactitude_limit){
    assert(exactitude_limit <= 10 and exactitude_limit >= 3);
    if(pins.size() <= exactitude_limit){
        if(pins.size() <= 3){
            if(pins.size() == 2){
                return std::abs(pins[0].x_ - pins[1].x_) + std::abs(pins[0].y_ - pins[1].y_);
            }
            else if(pins.size() == 3){
                auto minmaxX = std::minmax_element(pins.begin(), pins.end(), [](point<float_t> a, point<float_t> b){ return a.x_ < b.x_; }), 
                     minmaxY = std::minmax_element(pins.begin(), pins.end(), [](point<float_t> a, point<float_t> b){ return a.y_ < b.y_; });
                return (minmaxX.second->x_ - minmaxX.first->x_) + (minmaxY.second->y_ - minmaxY.first->y_);
                
            }
            else{
                return 0.0;
            }
        }
        switch(pins.size()){
            case 4:
                return get_wirelength<4, 2>(pins, steiner_lookup::topologies_4);
            case 5:
                return get_wirelength<5, 6>(pins, steiner_lookup::topologies_5);
            case 6:
                return get_wirelength<6, 23>(pins, steiner_lookup::topologies_6);
            case 7:
                return get_wirelength<7, 111>(pins, steiner_lookup::topologies_7);
            case 8:
                return get_wirelength<8, 642>(pins, steiner_lookup::topologies_8);
            case 9:
                return get_wirelength<9, 4334>(pins, steiner_lookup::topologies_9);
            case 10:
                return get_wirelength<10, 33510>(pins, steiner_lookup::topologies_10);
            default:
                abort();
        }
    }
    // Silent and dumb fallback to the spanning tree
    return MST_length(pins);
}

} // Namespace coloquinte

