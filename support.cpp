/*
 * support.cpp
 *
 *  Created on: Oct 7, 2021
 *      Author: d-w-h
 */

#include <iostream>
#include <vector>

#include "memory.hpp"
#include "usertypes.hpp"

void printRootList(Node* z) {
    Node* xt = z;
    if(xt != NULL) {
        do {
            std::cout << "xt->key: " << xt->key;
            std::cout << ", xt->degree: " << xt->degree << std::endl;
            xt = xt->right;
        } while(xt != z);
    }
}

void printChildList(Node* child) {
    Node* xt = child;
    if(xt != NULL) {
        do {
            std::cout << "xt->child->key: " << xt->key;
            std::cout << ", xt->child->degree: " << xt->degree << std::endl;
            if(xt->child != NULL) {
                std::cout << "xt->child->child->key: " << xt->child->key << std::endl;
            }
            xt = xt->right;
        } while(xt != child);
    }
}

void printList(Node* z) {
    Node* xt = z;
    if(xt != NULL) {
        do {
            std::cout << "xt->key: " << xt->key;
            std::cout << ", xt->degree: " << xt->degree << std::endl;
            if(xt->child != NULL) {
                printChildList(xt->child);
            }
            xt = xt->right;
        } while(xt != z);
    }
}

bool numbersChildrenMatch(Node* z, int& numNodes1) {
    bool numsMatch = true;
    int numNodes = 0;

    Node* xt = z->child;
    if(xt != NULL) {
        do {
            numNodes++;
            if(xt->child != NULL) {
                numsMatch = numbersChildrenMatch(xt, numNodes1);
                if(!numsMatch) { return false; }
            }
            xt = xt->right;
        } while(xt != z->child);

        numNodes1 = numNodes1 + numNodes;

        if(numNodes == z->degree) { numsMatch = true; }
        else { numsMatch = false; }
    }

    return numsMatch;
}

FibHeapProperties numbersMatch(Node* z) {
    bool numsMatch = true;
    int numNodes = 0;
    FibHeapProperties props = { numsMatch, numNodes };

    Node* xt = z;
    if(xt != NULL) {
        do {
            numNodes++;
            numsMatch = numbersChildrenMatch(xt, numNodes);
            props.degreeIsEqualToNumChildren = numsMatch;
            props.numNodes = numNodes;
            if(!numsMatch) { return props; }
            xt = xt->right;
        } while(xt != z);
    }

    props.degreeIsEqualToNumChildren = numsMatch;
    props.numNodes = numNodes;

    return props;
}

bool childrenAreFibHeap(Node * z) {
    bool isFibHeap = true;

    Node* xt = z->child;
    if(xt != NULL) {
        do {
            if(xt->p->key > xt->key) {
                return isFibHeap = false;
            }
            if(xt->child != NULL) {
                isFibHeap = childrenAreFibHeap(xt);
                if(!isFibHeap) { return false; }
            }
            xt = xt->right;
        } while(xt != z->child);
    }

    return isFibHeap;
}

bool isFibonacciHeap(Node * z) {
    bool isFibHeap = true;

    Node* xt = z;
    if(xt != NULL) {
        do {
            isFibHeap = childrenAreFibHeap(xt);
            if(!isFibHeap) { return false; }
            xt = xt->right;
        } while(xt != z);
    }

    return isFibHeap;
}

bool checkFibHeap(FibHeap * H) {
    //This is the general test for the Fibonacci heap.
    //The function returns true if the heap satisfies
    //the Fibonacci heap properties

    //Compute heap properties
    FibHeapProperties props = numbersMatch(H->min);
    bool isFibHeap = isFibonacciHeap(H->min);

    //Check if number of children equal Node degrees
    bool degreeIsEqualToNumChildren = props.degreeIsEqualToNumChildren;

    //Check if number of Nodes counted in heap equals H.n
    int numNodes = props.numNodes;
    bool numNodesMatch = (numNodes == H->n);

    //Check to see if heap is properly structured
    bool heapIsOK = numNodesMatch && degreeIsEqualToNumChildren && isFibHeap;

    return heapIsOK;
}
