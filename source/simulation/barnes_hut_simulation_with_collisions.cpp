#include "simulation/barnes_hut_simulation_with_collisions.h"
#include "simulation/barnes_hut_simulation.h"
#include "simulation/naive_parallel_simulation.h"
#include <omp.h>

void BarnesHutSimulationWithCollisions::simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    for(int i = 0; i < num_epochs; i++){
        simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs);
    }
}

void BarnesHutSimulationWithCollisions::simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    Quadtree qt(universe, universe.get_bounding_box(), 0);
    qt.calculate_cumulative_masses();
    qt.calculate_center_of_mass();
    calculate_forces(universe, qt);
    NaiveParallelSimulation::calculate_velocities(universe);
    NaiveParallelSimulation::calculate_positions(universe);
    find_collisions(universe);
    universe.current_simulation_epoch++;
    if (create_intermediate_plots) {
        if ((universe.current_simulation_epoch % plot_intermediate_epochs) == 0) {
            plotter.add_bodies_to_image(universe);
            plotter.write_and_clear();
        }
    }
}

void BarnesHutSimulationWithCollisions::find_collisions(Universe& universe){

    std::vector<std::pair<int, int>> collisions;
    for (int body_idx = 0; body_idx < universe.num_bodies; body_idx++) {

        Vector2d<double> body_position = universe.positions[body_idx];

        double body_mass = universe.weights[body_idx];

        Vector2d<double> applied_force_vector;

        for (int distant_body_idx = 0; distant_body_idx < universe.num_bodies; distant_body_idx++) {
            if (body_idx == distant_body_idx) {
                continue;
            }
            Vector2d<double> distant_body_position = universe.positions[distant_body_idx];

            Vector2d<double> direction_vector = distant_body_position - body_position;

            double distance = sqrt(pow(direction_vector[0], 2) + pow(direction_vector[1], 2));
            if (distance < 100000000000) {
	            if (universe.weights[body_idx] <= universe.weights[distant_body_idx]) {
                    collisions.emplace_back(distant_body_idx, body_idx);

	            } else {
                    collisions.emplace_back(body_idx, distant_body_idx);
	            }
            }
        }
    }

    std::ranges::sort(collisions, [universe](const std::pair<int, int> a, const std::pair<int, int>& b) {
	    return universe.weights[a.first] > universe.weights[b.first];
    });

    
    std::vector<int> del;
    for (int i = 0; i < collisions.size(); i++) {
        int index1 = collisions[i].first;
        int index2 = collisions[i].second;

        if (universe.weights[index1] == 0.0 || universe.weights[index2] == 0.0) {
            continue;
        }

        double mass1 = universe.weights[index1];
        double mass2 = universe.weights[index2];
        Vector2d<double> v1 = universe.velocities[index1];
        Vector2d<double> v2 = universe.velocities[index2];

        double res1 = (mass1 * v1.operator[](0) + mass2 * v2.operator[](0)/ mass1 + mass2);
        double res2 = (mass1 * v1.operator[](1) + mass2 * v2.operator[](1) / mass1 + mass2);
        universe.velocities[index1] = Vector2d<double>(res1, res2);
        universe.weights[index1] = mass1 + mass2;
        universe.weights[index2] = 0.0;
        del.emplace_back(index2);
    }

    std::ranges::sort(del, [](int a, int b) {return a > b;});

    for (int i = 0; i < del.size(); ++i) {
        int index = del[i];
    	universe.positions.erase(universe.positions.begin() + index);
        universe.forces.erase(universe.forces.begin() + index);
        universe.weights.erase(universe.weights.begin() + index);
        universe.velocities.erase(universe.velocities.begin() + index);
        universe.num_bodies--;
    }
}

void BarnesHutSimulationWithCollisions::find_collisions_parallel(Universe& universe){
    std::vector<std::pair<int, int>> collisions;
#pragma omp parallel for
    for (int body_idx = 0; body_idx < universe.num_bodies; body_idx++) {

        Vector2d<double> body_position = universe.positions[body_idx];

        double body_mass = universe.weights[body_idx];

        Vector2d<double> applied_force_vector;

        for (int distant_body_idx = 0; distant_body_idx < universe.num_bodies; distant_body_idx++) {
            if (body_idx == distant_body_idx) {
                continue;
            }
            Vector2d<double> distant_body_position = universe.positions[distant_body_idx];

            Vector2d<double> direction_vector = distant_body_position - body_position;

            double distance = sqrt(pow(direction_vector[0], 2) + pow(direction_vector[1], 2));
            if (distance < 100000000000) {
                if (universe.weights[body_idx] <= universe.weights[distant_body_idx]) {
                    collisions.emplace_back(distant_body_idx, body_idx);

                }
                else {
                    collisions.emplace_back(body_idx, distant_body_idx);
                }
            }
        }
    }

    std::ranges::sort(collisions, [universe](const std::pair<int, int> a, const std::pair<int, int>& b) {
        return universe.weights[a.first] > universe.weights[b.first];
        });


    std::vector<int> del;
    for (int i = 0; i < collisions.size(); i++) {
        int index1 = collisions[i].first;
        int index2 = collisions[i].second;

        if (universe.weights[index1] == 0.0 || universe.weights[index2] == 0.0) {
            continue;
        }

        double mass1 = universe.weights[index1];
        double mass2 = universe.weights[index2];
        Vector2d<double> v1 = universe.velocities[index1];
        Vector2d<double> v2 = universe.velocities[index2];

        double res1 = (mass1 * v1.operator[](0) + mass2 * v2.operator[](0) / mass1 + mass2);
        double res2 = (mass1 * v1.operator[](1) + mass2 * v2.operator[](1) / mass1 + mass2);
        universe.velocities[index1] = Vector2d<double>(res1, res2);
        universe.weights[index1] = mass1 + mass2;
        universe.weights[index2] = 0.0;
        del.emplace_back(index2);
    }

    std::ranges::sort(del, [](int a, int b) {return a > b;});

    for (int i = 0; i < del.size(); ++i) {
        int index = del[i];
        universe.positions.erase(universe.positions.begin() + index);
        universe.forces.erase(universe.forces.begin() + index);
        universe.weights.erase(universe.weights.begin() + index);
        universe.velocities.erase(universe.velocities.begin() + index);
        universe.num_bodies--;
    }
}