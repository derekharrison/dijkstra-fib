/*
 * UserTypes.hpp
 *
 *  Created on: Oct 7, 2021
 *      Author: d-w-h
 */

#ifndef USER_TYPES_HPP_
#define USER_TYPES_HPP_

#include <vector>

const int inf = 3e+8;

struct FibHeapProperties {
    bool degreeIsEqualToNumChildren;
    int numNodes;
};

struct Node {
    Node * left;
    Node * right;
    Node * p;
    Node * child;

    std::vector<int> adjNodes;
    std::vector<int> adjWeights;
    
    int key;
    int degree;
    int index;
    int index_og;
    bool mark;
};

class FibHeap {
public:
    int n;
    Node * min;
    FibHeap() { min = NULL; n = 0; }
};

#endif /* USER_TYPES_HPP_ */
