#include "quadtree.h"

#include "quadtreeNode.h"
#include <set>
#include <algorithm>
#include <stdexcept>
#include <omp.h>


Quadtree::Quadtree(Universe& universe, BoundingBox bounding_box, std::int8_t construct_mode) {
    root = new QuadtreeNode(bounding_box);
    std::vector<int32_t> body_indices = {};
    //filling body index with all indices in given universe
    for (int32_t i=0; i < universe.num_bodies; i++) { body_indices.push_back(i); }
    if (construct_mode == 0) {
       root->children=Quadtree::construct(universe, bounding_box, body_indices);
    }
    else if (construct_mode == 1) {
        root->children = Quadtree::construct_task(universe, bounding_box, body_indices);
    }
    else if (construct_mode == 2) {
        root->children = Quadtree::construct_task_with_cutoff(universe, bounding_box, body_indices);
    }
}

Quadtree::~Quadtree(){
    delete root;

}

void Quadtree::calculate_cumulative_masses(){
    root->calculate_node_cumulative_mass();
}

void Quadtree::calculate_center_of_mass(){
    root->calculate_node_center_of_mass();
}


std::vector<QuadtreeNode*> Quadtree::construct(Universe& universe, BoundingBox BB, std::vector<std::int32_t> body_indices){
    std::vector<QuadtreeNode*> resulting_children = {};
    for (int i = 0; i < 4; i++) {
        BoundingBox sub_BB = BB.get_quadrant(i);
        std::vector<std::int32_t> sub_body_indices = {};
        //Create subset of body indices that are in that quadrant for further use
        for (int j = 0; j < body_indices.size(); j++) {
            int32_t current_indice = body_indices[j];
            if (sub_BB.contains(universe.positions[current_indice])) {
                sub_body_indices.push_back(current_indice);
            }
        }
        //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
        if (!sub_body_indices.empty()) {
            QuadtreeNode* new_child = new QuadtreeNode(sub_BB);
            resulting_children.push_back(new_child);
            //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
            //Cumulative Mass and Center of Mass are the mass/position of the sole body
            if (sub_body_indices.size() == 1) {
                new_child->body_identifier = sub_body_indices[0]; 
                new_child->cumulative_mass = universe.weights[sub_body_indices[0]];
                new_child->cumulative_mass_ready = true;
                new_child->center_of_mass = universe.positions[sub_body_indices[0]];
                new_child->center_of_mass_ready = true;

            }
            else { new_child->children = construct(universe, sub_BB, sub_body_indices); }
        }
        
    }
    //Return the children of the Node that the Node corresponding to the current Bounding Box is supposed to get.
    return resulting_children;
}

std::vector<QuadtreeNode*> Quadtree::construct_task(Universe& universe, BoundingBox BB, std::vector<std::int32_t> body_indices){
    std::vector<QuadtreeNode*> resulting_children = {};

    
    for (int i = 0; i < 4; i++) {
        #pragma omp task shared(universe) default(private)
        {
        BoundingBox sub_BB = BB.get_quadrant(i);
        std::vector<std::int32_t> sub_body_indices = {};
        //Create subset of body indices that are in that quadrant for further use
        for (int j = 0; j < body_indices.size(); j++) {
            int32_t current_indice = body_indices[j];
            if (sub_BB.contains(universe.positions[current_indice])) {
                sub_body_indices.push_back(current_indice);
            }           
        }
        //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
        if (!sub_body_indices.empty()) {
            QuadtreeNode* new_child = new QuadtreeNode(sub_BB);
            resulting_children.push_back(new_child);
            //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
            //Cumulative Mass and Center of Mass are the mass/position of the sole body
            if (sub_body_indices.size() == 1) {
                new_child->body_identifier = sub_body_indices[0];
                new_child->cumulative_mass = universe.weights[sub_body_indices[0]];
                new_child->cumulative_mass_ready = true;
                new_child->center_of_mass = universe.positions[sub_body_indices[0]];
                new_child->center_of_mass_ready = true;

            }

            else {
                new_child->children = construct(universe, sub_BB, sub_body_indices);
                
            }
        }
    }
          //Return the children of the Node that the Node corresponding to the current Bounding Box is supposed to get.
    return resulting_children;
    }
}

std::vector<QuadtreeNode*> Quadtree::construct_task_with_cutoff(Universe& universe, BoundingBox& BB, std::vector<std::int32_t>& body_indices){
    return std::vector<QuadtreeNode*>();
}

std::vector<BoundingBox> Quadtree::get_bounding_boxes(QuadtreeNode* qtn){
    // traverse quadtree and collect bounding boxes
    std::vector<BoundingBox> result;
    // collect bounding boxes from children
    for(auto child: qtn->children){
        for(auto bb: get_bounding_boxes(child)){
            result.push_back(bb);
        }
    }
    result.push_back(qtn->bounding_box);
    return result;
}








