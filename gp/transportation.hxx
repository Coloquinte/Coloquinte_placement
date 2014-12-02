
#include "common.hxx"

#include <queue>
#include <vector>
#include <cassert>

namespace coloquinte{
namespace gp{


class current_allocation{
    static const index_t null_ind = std::numeric_limits<index_t>::max();

    // Internal data structures

    // Priority queue element to determine the source to be used between regions
    struct movable_source{
        index_t source;
        float_t cost;
        bool operator<(movable_source const o) const{
               return cost > o.cost // Sorted by cost
            || (cost == o.cost && source < o.source); // And by index to limit the number of fractional elements between two regions
        }
        movable_source(index_t s, float_t c) : source(s), cost(c) {}
    };

    // Member data


    // The current state
    std::vector<std::vector<capacity_t>  > sr_allocations; // For each region, for each source, the capacity allocated by the region
    std::vector<std::vector<float_t> >     sr_costs;       // The costs from a region to a source
    std::vector<capacity_t>                r_capacities;   // The remaining capacities of the regions

    // Shortest path data
    std::vector<float_t>                   r_costs;        // The costs of allocating to a region
    std::vector<index_t>                   r_parents;      // The parents of the regions i.e. the regions where we push sources first (or null_ind)
    std::vector<index_t>                   r_sources;      // The sources involved in these edges
    std::vector<capacity_t>                arc_capacities; // The capacities of the edges to the parents, or of the region if no parent

    // Best edges data
    std::vector<std::vector<std::priority_queue<movable_source> > > best_interregions_costs; // What is the best source to move to go from region k1 to region k2?
    index_t dijkstra_cnt;


    // Helper functions

    // Number of regions
    index_t region_cnt() const{
        assert(sr_costs.size() == sr_allocations.size());
        return sr_costs.size();
    }

    // Update the edge between two regions
    void update_edge(index_t r1, index_t r2);
    // Add a source to all heaps of a region
    void add_source_to_heaps(index_t r, index_t source);
    // Initialize the heaps of a region
    void create_heaps(index_t reg);

    // Update the edge and returns if the path has been modified so that we need to rerun Dijkstra
    bool push_edge(index_t reg, capacity_t flow);
    // Run the shortest path algorithm to update the cost of each region
    void dijkstra_update();
    // Update a full path when pushing an element, updating the paths
    capacity_t push_path(index_t pushed_reg, capacity_t demanded);

    public:
    // Add a new source to the transportation problem; should be done in decreasing order of demand to keep low complexity
    void add_source(capacity_t demand, std::vector<float_t> const & costs);

    current_allocation(std::vector<capacity_t> caps)
        :
        sr_allocations(caps.size()),
        sr_costs(caps.size()),
        r_capacities(caps),
        r_costs(caps.size(), 0.0),
        r_parents(caps.size(), null_ind),
        r_sources(caps.size(), null_ind),
        arc_capacities(caps),
        best_interregions_costs(caps.size(), std::vector<std::priority_queue<movable_source> >(caps.size())),
        dijkstra_cnt(0)
        {
            assert(caps.size() > 0);
            dijkstra_update();
        }

    std::vector<std::vector<capacity_t> > get_allocations() const{ return sr_allocations; }
    index_t get_iterations_cnt() const { return dijkstra_cnt; }
};

} // Namespace gp
}

