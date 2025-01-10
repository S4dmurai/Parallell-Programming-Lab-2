#include "quadtreeNode.h"

#include <iostream>


double QuadtreeNode::calculate_node_cumulative_mass(){
    if (cumulative_mass_ready) { return cumulative_mass;}//Activates in both leaves and nodes that have already been calculated
    else if (body_identifier == -1) {
        double result = 0.0;
        for (int i=0; i < children.size(); i++) {
            result += children[i]->calculate_node_cumulative_mass();
        }
        cumulative_mass = result; //Save mass for later use
        cumulative_mass_ready = true;
        return result;
    }
    else { std::runtime_error("Leaf  has no valid cumulative mass!"); } //Error that should NEVER happen, as leaf notes get constructed WITH mass
}

QuadtreeNode::QuadtreeNode(BoundingBox arg_bounding_box){
    bounding_box = arg_bounding_box;
    children = {};
    body_identifier = -1;
}

QuadtreeNode::~QuadtreeNode(){
    for (int i = 0; i < children.size(); i++) {
        delete children[i];
    }
}

Vector2d<double> QuadtreeNode::calculate_node_center_of_mass(){
    if (center_of_mass_ready) { return center_of_mass; }//Activates in both leaves and nodes that have already been calculated
    else if (body_identifier == -1) {
        Vector2d<double> result = Vector2d(0.0, 0.0);
        for (int i = 0; i < children.size(); i++) {
            result = result + (children[i]->calculate_node_center_of_mass() * children[i]->calculate_node_cumulative_mass());
        }
        result = result/calculate_node_cumulative_mass();
        center_of_mass = result; //Save mass for later use
        center_of_mass_ready = true;
        return result;
    }
    else { std::runtime_error("Leaf has no valid center of mass!"); } //Error that should NEVER happen, as leaf notes get constructed WITH center of mass
}