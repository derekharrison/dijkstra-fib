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

#include "fibheap.hpp"
#include "memory.hpp"
#include "support.hpp"
#include "usertypes.hpp"

void relax(Node * u, Node * v, int w, FibHeap * H) {

    if(v->key > u->key + w) {
        int weight = u->key + w;
        fibHeapDecreaseKey(H, v, weight);
    }
}

int mapIndex(int n, int index, int s) {
    int r;

    if(index >= s) { r = index - s; }
    else { r = n - s + index; }

    return r;
}

void setWeightMatAndRef(int sizeGraph,
                        std::vector<std::vector<int>> & edges,
                        int startVertex,
                        FibHeap * H,
                        Node ** NodeRefs) {


    //Initialize and construct heap
    for(int i = 0; i < sizeGraph; ++i) {
        NodeRefs[i] = ((Node *) calloc(1, sizeof(Node)));
        NodeRefs[i]->key = inf;
        NodeRefs[i]->index = i;
        if(i == 0) {
            NodeRefs[i]->key = 0;
        }
        fibHeapInsert(H, NodeRefs[i]);
    }

    //Set weight  matrix and adjacent Nodes
    int numEdges = (int) edges.size();
    for(int i = 0; i < numEdges; ++i) {
        int startIndex = edges[i][0] - 1;
        int endIndex = edges[i][1] - 1;
        int weight = edges[i][2];

        int start = mapIndex(sizeGraph, startIndex, startVertex);
        int end = mapIndex(sizeGraph, endIndex, startVertex);

        NodeRefs[start]->adjNodes.push_back(end);
        NodeRefs[start]->adjWeights.push_back(weight);
    }
}

void dijkstra(FibHeap * H, Node ** NodeRefs) {

    //Perform Dijkstra's algorithm
    while(H->n > 0) {
        Node * u = fibHeapExtractMin(H);

        for(int i = 0; i < u->adjNodes.size(); i++) {
            Node * v = NodeRefs[u->adjNodes[i]];
            int weight = u->adjWeights[i];
            relax(u, v, weight, H);
        }
    }
}

std::vector<int> reorderResults(FibHeap * H, int n, int s, Node ** NodeRefs) {

    std::vector<int> results;
    for(int i = 0; i < n; ++i) {
        if(i != s) {
            int index = mapIndex(n, i, s);
            if(NodeRefs[index]->key == inf) {
                results.push_back(-1);
            }
            else {
                results.push_back(NodeRefs[index]->key);
            }
        }
    }

    return results;
}

int max(int a, int b) {
    int res = 0;
    
    if(a > b) res = a;
    else res = b;
    
    return res;
}

std::vector<int> dijkstra(int n, std::vector<std::vector<int>> & edges, int s) {

    //Heap and Results
    FibHeap H;
    std::vector<int> results;

    //Map start vertex index
    s = s - 1;

    //Node references
    Node ** NodeRefs = getNodeRef(n);

    //Set weight matrix and create heap
    setWeightMatAndRef(n, edges, s, &H, NodeRefs);

    //Perform Dijkstra's algorithm
    dijkstra(&H, NodeRefs);

    //Reorder results
    results = reorderResults(&H, n, s, NodeRefs);

    //Deallocate memory
    freeNodeRef(NodeRefs, n);

    return results;
} 


