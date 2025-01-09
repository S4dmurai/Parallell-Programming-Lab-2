#include "quadtree.h"

#include "quadtreeNode.h"
#include <set>
#include <algorithm>
#include <stdexcept>
#include <omp.h>

Quadtree::Quadtree(Universe& universe, BoundingBox bounding_box, std::int8_t construct_mode) {
    QuadtreeNode root_node = QuadtreeNode(bounding_box);
    root = &root_node;
    std::vector<int32_t> body_indices = {};
    if (construct_mode == 0) {
        Quadtree::construct(universe, bounding_box, body_indices);
    }
    else if (construct_mode == 1) {
        Quadtree::construct_task(universe, bounding_box, body_indices);
    }
    else if (construct_mode == 2) {
        Quadtree::construct_task_with_cutoff(universe, bounding_box, body_indices);
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
    return std::vector<QuadtreeNode*>();
}

std::vector<QuadtreeNode*> Quadtree::construct_task(Universe& universe, BoundingBox BB, std::vector<std::int32_t> body_indices){
    return std::vector<QuadtreeNode*>();
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








