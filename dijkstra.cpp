/*
 * dijkstra.cpp
 *
 *  Created on: Oct 7, 2021
 *      Author: d-w-h
 */

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include "fib_heap.hpp"
#include "memory.hpp"
#include "user_types.hpp"

int map_index(int n, int index, int s) {
    int r;

    if(index >= s) { r = index - s; }
    else { r = n - s + index; }

    return r;
}

int map_inverse(int n, int index, int s) {
    int r;

    r = s + index;
    if(r > n - 1) {
        r = r - n;
    }

    return r;
}

void relax(node* u, node* v, int** w, FibHeap* H) {

    if(v->key > u->key + w[u->index][v->index]) {
        int weight = u->key + w[u->index][v->index];
        fib_heap_decrease_key(H, v, weight);
        v->key = weight;
    }
}

void populate_weight_and_ref(int size_graph,
                             std::vector< std::vector<int> >& edges,
                             int start_vertex,
                             FibHeap* H,
                             int** weight_mat,
                             node** node_refs) {

    //Create heap
    for(int i = 0; i < size_graph; ++i) {
        node_refs[i] = new node;
        node_refs[i]->key = inf;
        node_refs[i]->index = i;
        node_refs[i]->index_og = map_inverse(size_graph, i, start_vertex);
        if(i == 0) {
            node_refs[i]->key = 0;
        }
        fib_heap_insert(H, node_refs[i]);
    }

    //Add references to adjacent nodes and set weight matrix
    int num_edges = (int) edges.size();
    for(int i = 0; i < num_edges; ++i) {
        int start_index = edges[i][0] - 1;
        int end_index = edges[i][1] - 1;
        int weight = edges[i][2];

        int start = map_index(size_graph, start_index, start_vertex);
        int end = map_index(size_graph, end_index, start_vertex);

        node_refs[start]->adj_nodes.push_back(end);
        node_refs[end]->adj_nodes.push_back(start);

        weight_mat[start][end] = weight;
        weight_mat[end][start] = weight;
    }

    //Traverse edges again to pick minimum weight
    for(int i = 0; i < num_edges; ++i) {
        int start_index = edges[i][0] - 1;
        int end_index = edges[i][1] - 1;
        int weight = edges[i][2];

        int start = map_index(size_graph, start_index, start_vertex);
        int end = map_index(size_graph, end_index, start_vertex);

        bool is_greater = weight_mat[start][end] >= weight;
        if(is_greater) {
            weight_mat[start][end] = weight;
            weight_mat[end][start] = weight;
        }
    }
}

void dijkstra(FibHeap* H, int** w, node** node_refs) {

    //Perform Dijkstra's algorithm
    while(H->n > 0) {
        node* u = fib_heap_extract_min(H);

        int num_adj_nodes = (int) u->adj_nodes.size();
        for(int i = 0; i < num_adj_nodes; ++i) {
            int index_ref = u->adj_nodes[i];
            node* v = node_refs[index_ref];
            relax(u, v, w, H);
        }
    }
}

void reorder_results(int n, int s, node** node_refs, std::vector<int>& results) {

    for(int i = 0; i < n; ++i) {
        int index = map_index(n, i, s);
        if(node_refs[index]->index_og != s) {
            int key = node_refs[index]->key;
            if(key == inf) { key = -1; }
            results.push_back(key);
        }
    }
}

std::vector<int> shortest_reach(int n, std::vector< std::vector<int> >& edges, int s) {

    //Declarations
    FibHeap H;
    std::vector<int> results;

    //Map start index s to s - 1
    s = s - 1;

    //Initialize heap reference and weight mat
    node** node_refs = new node*[n];
    int** weight_mat = int2D(n);

    //Set weight mat and heap references and create heap
    populate_weight_and_ref(n, edges, s, &H, weight_mat, node_refs);

    //Perform Dijkstra's algorithm
    dijkstra(&H, weight_mat, node_refs);

    //Reorder results
    reorder_results(n, s, node_refs, results);

    //Deallocate memory
    free_int2D(weight_mat, n);
    free_node_refs(node_refs, n);

    return results;
}


