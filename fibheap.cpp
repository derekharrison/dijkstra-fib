/*
 * FibHeap.cpp
 *
 *  Created on: Oct 7, 2021
 *      Author: d-w-h
 */

#include <vector>

#include "memory.hpp"
#include "usertypes.hpp"

void fibHeapInsert(FibHeap * H, Node * x) {
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

void makeChildOf(FibHeap * H, Node * y, Node * x) {
    //Remove Node from root list
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

void linkDupDegree(FibHeap * H, Node ** A, Node *& x) {
    int d = x->degree;

    if(A[d] == x) //Don't link Nodes to themselves
        return;
        
    while(A[d] != NULL) {
        Node * y = A[d];
        //Link x and y
        if(x != y) {
            if(y->key > x->key) {
                //Make y child of x
                makeChildOf(H, y, x);
                
                if(y == H->min) {
                    H->min = x;
                }
            }
            else {
                //Make x child of y
                makeChildOf(H, x, y);
                
                //Reset root Node and root list tracker
                H->min = y;
                x = H->min;
            }
        }
        A[d] = NULL;
        d = d + 1;
    }
    A[d] = x;

}

void consolidate(FibHeap * H) {
    //Compute upper bound root list
    double golden = (1.0 + sqrt(5.0)) / 2.0;
    double f = log(H->n) / log(golden);
    int D = floor(f + 0.01) + 1;

    //Allocate memory for root list construction
    Node ** A = getNodeRef(D + 1);
    for(int i = 0; i < D + 1; ++i) {
        A[i] = NULL;
    }

    //Ensure all root Nodes have unique degrees
    Node * x = H->min;
    if(x != NULL) {
        do {
            linkDupDegree(H, A, x);
            x = x->right;
        } while(x != H->min);
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

    //Free root list reference
    free(A);
}

void nullifyChildrenParentNode(Node * z) {
    Node * xt = z->child;
    if(xt != NULL) {
        do {
            xt->p = NULL;
            xt = xt->right;
        } while(xt != z->child);
    }
}

Node * fibHeapExtractMin(FibHeap * H) {
    Node * z = H->min;

    if(z != NULL) {
        //Add each child of z to root list
        Node * y = z->child;
        if(y != NULL) {
            //Set children's parent Node to NULL
            nullifyChildrenParentNode(z);

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

void cut(FibHeap * H, Node * x, Node * y) {
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

void cascadingCut(FibHeap * H, Node * y) {
    Node * z = y->p;
    if(z != NULL) {
        if(y->mark == false) {
            y->mark = true;
        }
        else {
            cut(H, y, z);
            cascadingCut(H, z);
        }
    }
}

void fibHeapDecreaseKey(FibHeap * H, Node * x, int k) {
    if(k > x->key) {
        throw "key invalid";
    }

    x->key = k;
    Node * y = x->p;
    if(y != NULL && x->key < y->key) {
        cut(H, x, y);
        cascadingCut(H, y);
    }

    if(x->key < H->min->key) {
        H->min = x;
    }
}




