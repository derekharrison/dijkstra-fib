/*
 * fib_heap.cpp
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

#include "memory.hpp"
#include "user_types.hpp"

void fib_heap_insert(FibHeap* H, node* x) {
    x->degree = 0;
    x->p = NULL;
    x->child = NULL;
    x->mark = false;

    if(H->min == NULL) {
        x->left = x;
        x->right = x;
        H->min = x;
        H->n = 0;
    }
    else {
        x->left = H->min;
        x->right = H->min->right;
        H->min->right->left = x;
        H->min->right = x;
        if(x->key < H->min->key) {
            H->min = x;
        }
    }

    H->n = H->n + 1;
}

void make_child_of(FibHeap* H, node* y, node* x) {

    //Remove node from root list
    y->left->right = y->right;
    y->right->left = y->left;

    if(x->child == NULL) {
        x->child = y;
        y->p = x;
        y->left = y;
        y->right = y;
    }
    else {
        y->left = x->child;
        y->right = x->child->right;
        y->p = x;
        x->child->right->left = y;
        x->child->right = y;
    }

    //Set mark
    y->mark = false;

    x->degree = x->degree + 1;
}

void link(FibHeap* H, node** A, node* x, node* y) {

    int d = x->degree;

    //Make y child of x;
    make_child_of(H, y, x);

    A[d] = NULL;
    A[d+1] = x;

    if(y == H->min) {
        H->min = x;
    }
}

void link_dup_deg(FibHeap* H, node** A, node*& x, bool& there_is_dup) {
    int d = x->degree;
    //There is a node with the same degree and A[d] is not node x
    if(A[d] != NULL && A[d] != x) {
        there_is_dup = true;
        node* y = A[d];
        //Link x and y
        if(y->key > x->key) {
            //Make y child of x;
            link(H, A, x, y);
        }
        else {
            //Make x child of y;
            link(H, A, y, x);
            x = y;
        }
    }
    //There is no node with the same degree or A[d] is node x
    else {
        A[d] = x;
    }
}

void consolidate(FibHeap* H) {

    //Compute upper bound root list
    double golden = (1.0 + sqrt(5.0)) / 2.0;
    double f = log(H->n) / log(golden);
    int D = floor(f + 0.01) + 1;

    //Allocate memory for root list construction
    node** A = new node*[D + 1];
    for(int i = 0; i < D + 1; ++i) {
        A[i] = NULL;
    }

    //Ensure all root nodes have unique degrees
    node* x = H->min;
    if(x != NULL) {
        bool there_is_dup = true;
        while(there_is_dup) {
            there_is_dup = false;
            x = H->min;
            do {
                link_dup_deg(H, A, x, there_is_dup);
                x = x->right;
            }  while(x != H->min);
        }
    }

    //Reconstruct root list
    H->min = NULL;
    for(int i = 0; i < D + 1; ++i) {
        if(A[i] != NULL) {
            if(H->min == NULL) {
                A[i]->left = A[i];
                A[i]->right = A[i];
                A[i]->p = NULL;
                H->min = A[i];
            }
            else {
                A[i]->left = H->min;
                A[i]->right = H->min->right;
                H->min->right->left = A[i];
                H->min->right = A[i];
                A[i]->p = NULL;
                if(A[i]->key < H->min->key) {
                    H->min = A[i];
                }
            }
        }
    }

    //Free root list references
    delete [] A;
}



void nullify_children_parent_node(node* z) {
    node* xt = z->child;
    if(xt != NULL) {
        do {
            xt->p = NULL;
            xt = xt->right;
        } while(xt != z->child);
    }
}

node* fib_heap_extract_min(FibHeap* H) {

    node* z = H->min;

    if(z != NULL) {
        //Add each child of z to root list
        node* y = z->child;
        if(y != NULL) {
            //Set children's parent node to NULL
            nullify_children_parent_node(z);

            y->left->right = z->right;
            z->right->left = y->left;
            y->left = z;
            z->right = y;
            z->degree = 0;

            z->child = NULL;
        }

        //Remove z from root list
        z->left->right = z->right;
        z->right->left = z->left;

        if(z == z->right) {
            H->min = NULL;
        }
        else {
            H->min = z->right;
            consolidate(H);
        }

        H->n = H->n - 1;

    }

    return z;

}

void cut(FibHeap* H, node* x, node* y) {

    //If x is only child set child of parent to null
    if(x == x->right) {
        y->child = NULL;
        y->degree = 0;
    }
    else {
        y->child = x->right;
        y->degree = y->degree - 1;
    }

    //Remove x from child list of y and add x to root list of H
    x->left->right = x->right;
    x->right->left = x->left;

    x->right = H->min->right;
    x->left = H->min;

    H->min->right->left = x;
    H->min->right = x;

    x->p = NULL;
    x->mark = false;
}

void cascading_cut(FibHeap* H, node* y) {
    node* z = y->p;
    if(z != NULL) {
        if(y->mark == false) {
            y->mark = true;
        }
        else {
            cut(H, y, z);
            cascading_cut(H, z);
        }
    }
}

void fib_heap_decrease_key(FibHeap* H, node* x, int k) {
    if(k > x->key) {
        const char* s = "new key is greater than current key";
        std::cout << s << std::endl;
        throw s;
    }

    x->key = k;
    node* y = x->p;
    if(y != NULL && x->key < y->key) {
        cut(H, x, y);
        cascading_cut(H, y);
    }

    if(x->key < H->min->key) {
        H->min = x;
    }
}




