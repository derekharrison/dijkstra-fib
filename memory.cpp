/*
 * memory.cpp
 *
 *  Created on: Oct 7, 2021
 *      Author: d-w-h
 */

#include "usertypes.hpp"

Node ** getNodeRef(int n) {
    
    return (Node **) malloc(sizeof(Node *) * n);
}

void freeNodeRef(Node ** vRef, int size) {
    for(int i = 0; i < size; ++i) {
        free(vRef[i]);
    }

    free(vRef);
}

