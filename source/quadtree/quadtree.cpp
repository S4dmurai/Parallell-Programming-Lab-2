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
#pragma omp parallel
        {
#pragma omp single
            {
                root->children = Quadtree::construct_task(universe, bounding_box, body_indices);
            }
        }
    }
    else if (construct_mode == 2) {
#pragma omp parallel
        {
#pragma omp single
            {
        root->children = Quadtree::construct_task_with_cutoff(universe, bounding_box, body_indices);
            }
        }
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
    //Init possible children as NONE/Nullpointer to later be merged into resulting_children
    QuadtreeNode* q0 = nullptr;
    QuadtreeNode* q1 = nullptr;
    QuadtreeNode* q2 = nullptr;
    QuadtreeNode* q3 = nullptr;
    BoundingBox q0_BB = BB.get_quadrant(0);
    std::vector<std::int32_t> q0_body_indices = {};
    //Create subset of body indices that are in that quadrant for further use
    for (int j = 0; j < body_indices.size(); j++) {
        int32_t current_indice = body_indices[j];
        if (q0_BB.contains(universe.positions[current_indice])) {
            q0_body_indices.push_back(current_indice);
        }
    }
    //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
    if (!q0_body_indices.empty()) {
        q0 = new QuadtreeNode(q0_BB);

        //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
        //Cumulative Mass and Center of Mass are the mass/position of the sole body
        if (q0_body_indices.size() == 1) {
            q0->body_identifier = q0_body_indices[0];
            q0->cumulative_mass = universe.weights[q0_body_indices[0]];
            q0->cumulative_mass_ready = true;
            q0->center_of_mass = universe.positions[q0_body_indices[0]];
            q0->center_of_mass_ready = true;

        }

        else {
            q0->children = construct(universe, q0_BB, q0_body_indices);

        }
    }
 
    BoundingBox q1_BB = BB.get_quadrant(1);
    std::vector<std::int32_t> q1_body_indices = {};
    //Create subset of body indices that are in that quadrant for further use
    for (int j = 0; j < body_indices.size(); j++) {
        int32_t current_indice = body_indices[j];
        if (q1_BB.contains(universe.positions[current_indice])) {
            q1_body_indices.push_back(current_indice);
        }
    }
    //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
    if (!q1_body_indices.empty()) {
        q1 = new QuadtreeNode(q1_BB);

        //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
        //Cumulative Mass and Center of Mass are the mass/position of the sole body
        if (q1_body_indices.size() == 1) {
            q1->body_identifier = q1_body_indices[0];
            q1->cumulative_mass = universe.weights[q1_body_indices[0]];
            q1->cumulative_mass_ready = true;
            q1->center_of_mass = universe.positions[q1_body_indices[0]];
            q1->center_of_mass_ready = true;

        }

        else {
            q1->children = construct(universe, q1_BB, q1_body_indices);

        }
    }
    BoundingBox q2_BB = BB.get_quadrant(2);
    std::vector<std::int32_t> q2_body_indices = {};
    //Create subset of body indices that are in that quadrant for further use
    for (int j = 0; j < body_indices.size(); j++) {
        int32_t current_indice = body_indices[j];
        if (q2_BB.contains(universe.positions[current_indice])) {
            q2_body_indices.push_back(current_indice);
        }
    }
    //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
    if (!q2_body_indices.empty()) {
        q2 = new QuadtreeNode(q2_BB);

        //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
        //Cumulative Mass and Center of Mass are the mass/position of the sole body
        if (q2_body_indices.size() == 1) {
            q2->body_identifier = q2_body_indices[0];
            q2->cumulative_mass = universe.weights[q2_body_indices[0]];
            q2->cumulative_mass_ready = true;
            q2->center_of_mass = universe.positions[q2_body_indices[0]];
            q2->center_of_mass_ready = true;

        }

        else {
            q2->children = construct(universe, q2_BB, q2_body_indices);

        }
    }

    BoundingBox q3_BB = BB.get_quadrant(3);
    std::vector<std::int32_t> q3_body_indices = {};
    //Create subset of body indices that are in that quadrant for further use
    for (int j = 0; j < body_indices.size(); j++) {
        int32_t current_indice = body_indices[j];
        if (q3_BB.contains(universe.positions[current_indice])) {
            q3_body_indices.push_back(current_indice);
        }
    }
    //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
    if (!q3_body_indices.empty()) {
        q3 = new QuadtreeNode(q3_BB);
        //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
        //Cumulative Mass and Center of Mass are the mass/position of the sole body
        if (q3_body_indices.size() == 1) {
            q3->body_identifier = q3_body_indices[0];
            q3->cumulative_mass = universe.weights[q3_body_indices[0]];
            q3->cumulative_mass_ready = true;
            q3->center_of_mass = universe.positions[q3_body_indices[0]];
            q3->center_of_mass_ready = true;

        }

        else {
            q3->children = construct(universe, q3_BB, q3_body_indices);

        }
    }

    //Return the children of the Node that the Node corresponding to the current Bounding Box is supposed to get by merging teh children (if there are any)
    std::vector<QuadtreeNode*> resulting_children = {};

    if (q0 != nullptr) { resulting_children.push_back(q0); }
    if (q1 != nullptr) { resulting_children.push_back(q1); }
    if (q2 != nullptr) { resulting_children.push_back(q2); }
    if (q3 != nullptr) { resulting_children.push_back(q3); }
    return resulting_children;
}

std::vector<QuadtreeNode*> Quadtree::construct_task(Universe& universe, BoundingBox BB, std::vector<std::int32_t> body_indices){
    //Init possible children as NONE/Nullpointer to later be merged into resulting_children
    QuadtreeNode* q0 = nullptr;
    QuadtreeNode* q1 = nullptr;
    QuadtreeNode* q2 = nullptr;
    QuadtreeNode* q3 = nullptr;
#pragma omp task shared(universe, q0) //calculating q0
            {
                BoundingBox q0_BB = BB.get_quadrant(0);
                std::vector<std::int32_t> q0_body_indices = {};
                //Create subset of body indices that are in that quadrant for further use

                for (int j = 0; j < body_indices.size(); j++) {
                    int32_t current_indice = body_indices[j];
                    if (q0_BB.contains(universe.positions[current_indice])) {
                        q0_body_indices.push_back(current_indice);
                    }
                }
                //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
                if (!q0_body_indices.empty()) {
                    q0 = new QuadtreeNode(q0_BB);

                    //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
                    //Cumulative Mass and Center of Mass are the mass/position of the sole body
                    if (q0_body_indices.size() == 1) {
                        q0->body_identifier = q0_body_indices[0];
                        q0->cumulative_mass = universe.weights[q0_body_indices[0]];
                        q0->cumulative_mass_ready = true;
                        q0->center_of_mass = universe.positions[q0_body_indices[0]];
                        q0->center_of_mass_ready = true;

                    }

                    else {
                        q0->children = construct_task(universe, q0_BB, q0_body_indices);

                    }
                }
            }
#pragma omp task shared(universe, q1) //calculating q1
            {
                BoundingBox q1_BB = BB.get_quadrant(1);
                std::vector<std::int32_t> q1_body_indices = {};
                //Create subset of body indices that are in that quadrant for further use

                for (int j = 0; j < body_indices.size(); j++) {
                    int32_t current_indice = body_indices[j];
                    if (q1_BB.contains(universe.positions[current_indice])) {
                        q1_body_indices.push_back(current_indice);
                    }
                }
                //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
                if (!q1_body_indices.empty()) {
                    q1 = new QuadtreeNode(q1_BB);

                    //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
                    //Cumulative Mass and Center of Mass are the mass/position of the sole body
                    if (q1_body_indices.size() == 1) {
                        q1->body_identifier = q1_body_indices[0];
                        q1->cumulative_mass = universe.weights[q1_body_indices[0]];
                        q1->cumulative_mass_ready = true;
                        q1->center_of_mass = universe.positions[q1_body_indices[0]];
                        q1->center_of_mass_ready = true;

                    }

                    else {
                        q1->children = construct_task(universe, q1_BB, q1_body_indices);

                    }
                }
            }
#pragma omp task shared(universe, q2) //calculating q2
            {
                BoundingBox q2_BB = BB.get_quadrant(2);
                std::vector<std::int32_t> q2_body_indices = {};
                //Create subset of body indices that are in that quadrant for further use

                for (int j = 0; j < body_indices.size(); j++) {
                    int32_t current_indice = body_indices[j];
                    if (q2_BB.contains(universe.positions[current_indice])) {
                        q2_body_indices.push_back(current_indice);
                    }
                }
                //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
                if (!q2_body_indices.empty()) {
                    q2 = new QuadtreeNode(q2_BB);

                    //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
                    //Cumulative Mass and Center of Mass are the mass/position of the sole body
                    if (q2_body_indices.size() == 1) {
                        q2->body_identifier = q2_body_indices[0];
                        q2->cumulative_mass = universe.weights[q2_body_indices[0]];
                        q2->cumulative_mass_ready = true;
                        q2->center_of_mass = universe.positions[q2_body_indices[0]];
                        q2->center_of_mass_ready = true;

                    }

                    else {
                        q2->children = construct_task(universe, q2_BB, q2_body_indices);

                    }
                }
            }
#pragma omp task shared(q3, universe) //calculating q3
            {
                BoundingBox q3_BB = BB.get_quadrant(3);
                std::vector<std::int32_t> q3_body_indices = {};
                //Create subset of body indices that are in that quadrant for further use
                for (int j = 0; j < body_indices.size(); j++) {
                    int32_t current_indice = body_indices[j];
                    Vector2d<double> test= universe.positions[current_indice];
                    if (q3_BB.contains(universe.positions[current_indice])) {
                        q3_body_indices.push_back(current_indice);
                    }
                }
                //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
                if (!q3_body_indices.empty()) {
                    q3 = new QuadtreeNode(q3_BB);
                    //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
                    //Cumulative Mass and Center of Mass are the mass/position of the sole body
                    if (q3_body_indices.size() == 1) {
                        q3->body_identifier = q3_body_indices[0];
                        q3->cumulative_mass = universe.weights[q3_body_indices[0]];
                        q3->cumulative_mass_ready = true;
                        q3->center_of_mass = universe.positions[q3_body_indices[0]];
                        q3->center_of_mass_ready = true;

                    }

                    else {
                        q3->children = construct_task(universe, q3_BB, q3_body_indices);

                    }
                }
            }

          //Return the children of the Node that the Node corresponding to the current Bounding Box is supposed to get by merging teh children (if there are any)
#pragma omp taskwait
    std::vector<QuadtreeNode*> resulting_children = {};

    if (q0 != nullptr) { resulting_children.push_back(q0); }
    if (q1 != nullptr) { resulting_children.push_back(q1); }
    if (q2 != nullptr) { resulting_children.push_back(q2); }
    if (q3 != nullptr) { resulting_children.push_back(q3); }
    return resulting_children;
}

std::vector<QuadtreeNode*> Quadtree::construct_task_with_cutoff(Universe& universe, BoundingBox BB, std::vector<std::int32_t> body_indices) {
    uint32_t const CUTOFF = 1000;
    //Init possible children as NONE/Nullpointer to later be merged into resulting_children
    QuadtreeNode* q0 = nullptr;
    QuadtreeNode* q1 = nullptr;
    QuadtreeNode* q2 = nullptr;
    QuadtreeNode* q3 = nullptr;
#pragma omp task shared(universe, q0) if (body_indices.size()>CUTOFF) //calculating q0
    {
        BoundingBox q0_BB = BB.get_quadrant(0);
        std::vector<std::int32_t> q0_body_indices = {};
        //Create subset of body indices that are in that quadrant for further use

        for (int j = 0; j < body_indices.size(); j++) {
            int32_t current_indice = body_indices[j];
            if (q0_BB.contains(universe.positions[current_indice])) {
                q0_body_indices.push_back(current_indice);
            }
        }
        //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
        if (!q0_body_indices.empty()) {
            q0 = new QuadtreeNode(q0_BB);

            //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
            //Cumulative Mass and Center of Mass are the mass/position of the sole body
            if (q0_body_indices.size() == 1) {
                q0->body_identifier = q0_body_indices[0];
                q0->cumulative_mass = universe.weights[q0_body_indices[0]];
                q0->cumulative_mass_ready = true;
                q0->center_of_mass = universe.positions[q0_body_indices[0]];
                q0->center_of_mass_ready = true;

            }

            else {
                q0->children = construct_task_with_cutoff(universe, q0_BB, q0_body_indices);

            }
        }
    }
#pragma omp task shared(universe, q1) if (body_indices.size()>CUTOFF) //calculating q1
    {
        BoundingBox q1_BB = BB.get_quadrant(1);
        std::vector<std::int32_t> q1_body_indices = {};
        //Create subset of body indices that are in that quadrant for further use

        for (int j = 0; j < body_indices.size(); j++) {
            int32_t current_indice = body_indices[j];
            if (q1_BB.contains(universe.positions[current_indice])) {
                q1_body_indices.push_back(current_indice);
            }
        }
        //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
        if (!q1_body_indices.empty()) {
            q1 = new QuadtreeNode(q1_BB);

            //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
            //Cumulative Mass and Center of Mass are the mass/position of the sole body
            if (q1_body_indices.size() == 1) {
                q1->body_identifier = q1_body_indices[0];
                q1->cumulative_mass = universe.weights[q1_body_indices[0]];
                q1->cumulative_mass_ready = true;
                q1->center_of_mass = universe.positions[q1_body_indices[0]];
                q1->center_of_mass_ready = true;

            }

            else {
                q1->children = construct_task_with_cutoff(universe, q1_BB, q1_body_indices);

            }
        }
    }
#pragma omp task shared(universe, q2) if (body_indices.size()>CUTOFF) //calculating q2
    {
        BoundingBox q2_BB = BB.get_quadrant(2);
        std::vector<std::int32_t> q2_body_indices = {};
        //Create subset of body indices that are in that quadrant for further use

        for (int j = 0; j < body_indices.size(); j++) {
            int32_t current_indice = body_indices[j];
            if (q2_BB.contains(universe.positions[current_indice])) {
                q2_body_indices.push_back(current_indice);
            }
        }
        //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
        if (!q2_body_indices.empty()) {
            q2 = new QuadtreeNode(q2_BB);

            //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
            //Cumulative Mass and Center of Mass are the mass/position of the sole body
            if (q2_body_indices.size() == 1) {
                q2->body_identifier = q2_body_indices[0];
                q2->cumulative_mass = universe.weights[q2_body_indices[0]];
                q2->cumulative_mass_ready = true;
                q2->center_of_mass = universe.positions[q2_body_indices[0]];
                q2->center_of_mass_ready = true;

            }

            else {
                q2->children = construct_task_with_cutoff(universe, q2_BB, q2_body_indices);

            }
        }
    }
#pragma omp task shared(universe, q3) if (body_indices.size()>CUTOFF) //calculating q3
    {
        BoundingBox q3_BB = BB.get_quadrant(3);
        std::vector<std::int32_t> q3_body_indices = {};
        //Create subset of body indices that are in that quadrant for further use
        for (int j = 0; j < body_indices.size(); j++) {
            int32_t current_indice = body_indices[j];
            if (q3_BB.contains(universe.positions[current_indice])) {
                q3_body_indices.push_back(current_indice);
            }
        }
        //If quadrant isn't empty, create a Quadtree node with that Quadtrant, then recursively use condtruct to create that node's children. If it's empty, ignore quadrant.
        if (!q3_body_indices.empty()) {
            q3 = new QuadtreeNode(q3_BB);
            //If only one body in child (thus a leaf node), update body identifier, center of mass and cumulative mass, else recursive construction.
            //Cumulative Mass and Center of Mass are the mass/position of the sole body
            if (q3_body_indices.size() == 1) {
                q3->body_identifier = q3_body_indices[0];
                q3->cumulative_mass = universe.weights[q3_body_indices[0]];
                q3->cumulative_mass_ready = true;
                q3->center_of_mass = universe.positions[q3_body_indices[0]];
                q3->center_of_mass_ready = true;

            }

            else {
                q3->children = construct_task_with_cutoff(universe, q3_BB, q3_body_indices);

            }
        }
    }

    //Return the children of the Node that the Node corresponding to the current Bounding Box is supposed to get by merging teh children (if there are any)
#pragma omp taskwait
    std::vector<QuadtreeNode*> resulting_children = {};

    if (q0 != nullptr) { resulting_children.push_back(q0); }
    if (q1 != nullptr) { resulting_children.push_back(q1); }
    if (q2 != nullptr) { resulting_children.push_back(q2); }
    if (q3 != nullptr) { resulting_children.push_back(q3); }
    return resulting_children;
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








