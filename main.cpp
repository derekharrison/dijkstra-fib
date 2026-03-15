/*
 * main.cpp
 *
 *  Created on: Sep 17, 2021
 *      Author: d-w-h
 */

#include <iostream>
#include <vector>

#include "dijkstra.hpp"

int main(int argc, char* argv[]) {

    //Declarations
    int s = 2; //Start vertex. The minimum index for vertices is 1
    int n = 2499; //Number of vertices
    int numEdges = 33125; //Number of edges

    //Create edges
    std::vector< std::vector<int> > edges;
    for(int i = 0; i < numEdges; ++i) {
        int startVertex = rand() % n + 1;
        int endVertex = rand() % n + 1;
        int weight = rand() % 200 + 1;

        std::vector<int> edge;
        edge.push_back(startVertex);
        edge.push_back(endVertex);
        edge.push_back(weight);
        edges.push_back(edge);
    }

    //Compute distances to Nodes from start vertex
    std::vector<int> results = dijkstra(n, edges, s);

    //Print results
    int size_results = (int) results.size();
    for(int i = 0; i < size_results; ++i) {
        std::cout << results[i] << " ";
    }
    
    return 0;
}
