/*
 * dijkstra.cpp
 *
 *  Created on: Oct 7, 2021
 *      Author: d-w-h
 */

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

    if(index >= s)  r = index - s;
    else r = n - s + index;

    return r;
}

void createFibonacciHeap(int sizeGraph, Node ** NodeRefs, FibHeap * H) {
    for(int i = 0; i < sizeGraph; ++i)
        fibHeapInsert(H, NodeRefs[i]);
}

void createNodeReferences(int sizeGraph, std::vector<std::vector<int>> & edges, int startVertex, Node ** NodeRefs) {
    for(int i = 0; i < sizeGraph; ++i) {
        NodeRefs[i] = ((Node *) calloc(1, sizeof(Node)));
        NodeRefs[i]->key = inf;
        NodeRefs[i]->index = i;
        if(i == 0)
            NodeRefs[i]->key = 0;
    }
    
    for(std::vector<int> edge : edges) {
        int startIndex = edge[0] - 1;
        int endIndex = edge[1] - 1;
        int weight = edge[2];

        int start = mapIndex(sizeGraph, startIndex, startVertex);
        int end = mapIndex(sizeGraph, endIndex, startVertex);

        NodeRefs[start]->adjNodes.push_back(end);
        NodeRefs[start]->adjWeights.push_back(weight);
    }
}

void dijkstra(FibHeap * H, Node ** NodeRefs) {
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
    FibHeap H;
    std::vector<int> results;

    s = s - 1;

    Node ** NodeRefs = getNodeRef(n);
    
    createNodeReferences(n, edges, s, NodeRefs);
    
    createFibonacciHeap(n, NodeRefs, &H);

    dijkstra(&H, NodeRefs);

    results = reorderResults(&H, n, s, NodeRefs);

    freeNodeRef(NodeRefs, n);

    return results;
} 


