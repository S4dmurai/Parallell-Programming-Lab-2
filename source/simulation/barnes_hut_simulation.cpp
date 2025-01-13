#include "simulation/barnes_hut_simulation.h"
#include "simulation/naive_parallel_simulation.h"
#include "physics/gravitation.h"
#include "physics/mechanics.h"

#include <cmath>

void BarnesHutSimulation::simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    for(int i = 0; i < num_epochs; i++){
        simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs);
    }
}

void BarnesHutSimulation::simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    Quadtree qt(universe, universe.get_bounding_box(), 0);
    qt.calculate_cumulative_masses();
    qt.calculate_center_of_mass();
    calculate_forces(universe, qt);
    NaiveParallelSimulation::calculate_velocities(universe);
    NaiveParallelSimulation::calculate_positions(universe);
    universe.current_simulation_epoch++;
    if (create_intermediate_plots) {
        if ((universe.current_simulation_epoch % plot_intermediate_epochs) == 0) {
            plotter.add_bodies_to_image(universe);
            plotter.write_and_clear();
        }
    }
    
}

void BarnesHutSimulation::get_relevant_nodes(Universe& universe, Quadtree& quadtree, std::vector<QuadtreeNode*>& relevant_nodes, Vector2d<double>& body_position, std::int32_t body_index, double threshold_theta){
    relevant_nodes.push_back(quadtree.root);

    for (int i = 0; i < relevant_nodes.size();) {
        //std::cout << "starti" << i << "\n";
        //Quadrant only containing K body is trivially irrelevant
        if (relevant_nodes[i]->body_identifier == body_index) {
            relevant_nodes.erase(relevant_nodes.begin() + i);
            continue;
        }
        //Quadrant containing body K and others needs to be subdivided
        if (relevant_nodes[i]->bounding_box.contains(body_position)) {
            relevant_nodes.insert(relevant_nodes.end(), relevant_nodes[i]->children.begin(), relevant_nodes[i]->children.end());
            relevant_nodes.erase(relevant_nodes.begin() + i);
        }
        //Processing of Quadrant not containing body K
        if (relevant_nodes[i]->body_identifier != -1) {
            i++;
            continue;
                        //Quadrant only containing single body (not K) is trivially relevant
        }else {
		    Vector2d<double> direction_vector = relevant_nodes[i]->calculate_node_center_of_mass() - body_position;
            double distance = sqrt(pow(direction_vector[0], 2) + pow(direction_vector[1], 2));
            double theta = relevant_nodes[i]->bounding_box.get_diagonal() / distance;
            //std::cout <<"i" << i << "\n" << relevant_nodes.size() << "\n";
            //d/r<THETA, quadrant is relevant
            if (theta < threshold_theta) {
                i++;
		    }
            //d/r>=THETA, Quadrant needs to be subdivided
		    else {

                //std::cout << "SPLIT childnum " << relevant_nodes[i]->children.size() << "\n";
                relevant_nodes.insert(relevant_nodes.end(), relevant_nodes[i]->children.begin(), relevant_nodes[i]->children.end());
                relevant_nodes.erase(relevant_nodes.begin() + i);
                
		    }
	    }
    }
}

void BarnesHutSimulation::calculate_forces(Universe& universe, Quadtree& quadtree){
    const double threshold_theta = 0.2;
#pragma omp parallel for
    for (int body_idx = 0; body_idx < universe.num_bodies; body_idx++){
        Vector2d<double> body_position = universe.positions[body_idx];
        std::vector<QuadtreeNode*> relevant_nodes;
        get_relevant_nodes(universe, quadtree, relevant_nodes, body_position, body_idx, threshold_theta);

        double body_mass = universe.weights[body_idx];

        Vector2d<double> applied_force_vector;

        for (int i = 0; i < relevant_nodes.size(); ++i) {
            if (!relevant_nodes[i]->center_of_mass_ready) {
                relevant_nodes[i]->calculate_node_center_of_mass();
            }
            // get distant body positions
            Vector2d<double> distant_body_position = relevant_nodes[i]->center_of_mass;
            // calculate vector between bodies to get the direction of the gravitational force
            Vector2d<double> direction_vector = distant_body_position - body_position;

            // calculate the distance between the bodies
            double distance = sqrt(pow(direction_vector[0], 2) + pow(direction_vector[1], 2));

            if (!relevant_nodes[i]->cumulative_mass_ready) {
                relevant_nodes[i]->calculate_node_cumulative_mass();
            }

            // calculate gravitational force between the bodies
            double force = gravitational_force(body_mass, relevant_nodes[i]->cumulative_mass, distance);

            // create the force vector
            Vector2d<double> force_vector = direction_vector * (force / distance);

            // sum forces applied to body
            applied_force_vector = applied_force_vector + force_vector;

        }
        // store applied force 
        universe.forces[body_idx] = applied_force_vector;
    }

}