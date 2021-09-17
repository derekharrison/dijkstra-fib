/*
 * main.cpp
 *
 *  Created on: Sep 17, 2021
 *      Author: d-w-h
 */

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <stdlib.h>
#include <time.h>

typedef struct Node {
	Node* left;
	Node* right;
	Node* p;
	Node* child;

	int key;
	int degree;
	bool mark;
} node;

class FibHeap {
public:
	int n;
	int n_root;
	node* min = NULL;
};

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

	if(y == H->min) {
		H->min = x;
	}

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

	x->degree = x->degree + 1;
}

void fib_heap_link(FibHeap* H, node* y, node* x) {

	//Remove node from root list
	y->left->right = y->right;
	y->right->left = y->left;

	//Make y child of x
	make_child_of(H, y, x);

	//Set mark
	y->mark = false;
}

void consolidate(FibHeap* H) {

	double golden = (1.0 + sqrt(5.0)) / 2.0;
	double f = log(H->n) / log(golden);
	int D = floor(f + 0.01) + 1;

	node** A = new node*[D + 2];
	for(int i = 0; i < D + 2; ++i) {
		A[i] = NULL;
	}

	node* x = H->min;
	while(x->right != H->min) {
		int d = x->degree;
		while(A[d] != NULL && A[d] != x) {
			node* y = A[d];

			if(x->key > y->key) {
				node* dummy = x;
				x = y;
				y = dummy;
			}

			fib_heap_link(H, y, x);

			A[d] = NULL;
			d = d + 1;
		}
		A[d] = x;
		x = x->right;
	}

	if(x->right == H->min) {
		int d = x->degree;
		while(A[d] != NULL && A[d] != x) {
			node* y = A[d];

			if(x->key > y->key) {
				node* dummy = x;
				x = y;
				y = dummy;
			}

			fib_heap_link(H, y, x);

			A[d] = NULL;
			d = d + 1;
		}
		A[d] = x;
	}

	H->min = NULL;
	for(int i = 0; i < D + 2; ++i) {
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

}

int num_times_exec = 0;
int tot_num_nodes = 0;

void print_child_circle(node* child);
void print_circle(node* z);

void print_child_circle(node* child) {
	node* xt = child;
	if(xt != NULL && num_times_exec > 80) {
		if(xt->right != child) {
			while(xt->right != child) {
				std::cout << "xt->child->key: " << xt->key;
				std::cout << ", xt->child->degree: " << xt->degree << std::endl;
				if(xt->child != NULL) {
					std::cout << "xt->child->child->key: " << xt->child->key << std::endl;
//					print_circle(xt->child);
				}
				xt = xt->right;
			}
			if(xt->right == child) {
				std::cout << "xt->child->key: " << xt->key;
				std::cout << ", xt->child->degree: " << xt->degree << std::endl;
				if(xt->child != NULL) {
					if(xt->child != NULL) {
						std::cout << "xt->child->child->key: " << xt->child->key << std::endl;
//						print_circle(xt->child);
					}
				}
			}
		}
		else {
			std::cout << "X->CHILD == X->CHILD->RIGHT" << std::endl;
			std::cout << "xt->child->key: " << xt->key;
			std::cout << ", xt->child->degree: " << xt->degree << std::endl;
		}
	}
}

void print_circle(node* z) {
	node* xt = z;
	if(xt != NULL && num_times_exec > 80) {
		if(xt->right != z) {
			while(xt->right != z) {
				std::cout << "xt->key: " << xt->key;
				std::cout << ", xt->degree: " << xt->degree << std::endl;
				if(xt->child != NULL) {
					print_child_circle(xt->child);
				}
				xt = xt->right;
			}
			if(xt->right == z) {
				std::cout << "xt->key: " << xt->key;
				std::cout << ", xt->degree: " << xt->degree << std::endl;
				if(xt->child != NULL) {
					if(xt->child != NULL) {
						print_child_circle(xt->child);
					}
				}
			}
		}
		else {
			std::cout << "X == X->RIGHT" << std::endl;
			std::cout << "xt->key: " << xt->key;
			std::cout << ", xt->degree: " << xt->degree << std::endl;
			if(xt->child != NULL) {
				print_child_circle(xt->child);
			}
		}
	}
}

bool numbers_children_match(node* z) {
	bool nums_match = true;
	int num_of_nodes = 0;

	node* xt = z->child;
	if(xt != NULL) {
		while(xt->right != z->child) {
			num_of_nodes++;
			tot_num_nodes++;
			if(xt->child != NULL) {
				nums_match = numbers_children_match(xt);
				if(!nums_match) { return false; }
			}
			xt = xt->right;
		}
		if(xt->right == z->child) {
			num_of_nodes++;
			tot_num_nodes++;
			if(xt->child != NULL) {
				nums_match = numbers_children_match(xt);
				if(!nums_match) { return false; }
			}
		}

		if(num_of_nodes == z->degree) { nums_match = true; }
		else { nums_match = false; }
	}

	return nums_match;
}

bool numbers_match(node* z) {
	bool nums_match = true;

	node* xt = z;
	if(xt != NULL) {
		while(xt->right != z) {
			tot_num_nodes++;
			nums_match = numbers_children_match(xt);
			if(!nums_match) { return false; }
			xt = xt->right;
		}
		if(xt->right == z) {
			tot_num_nodes++;
			nums_match = numbers_children_match(xt);
			if(!nums_match) { return false; }
		}
	}

	return nums_match;
}

bool is_fib_heap_children(node* z) {
	bool is_fibheap = true;

	node* xt = z->child;
	if(xt != NULL) {
		while(xt->right != z->child) {
			if(xt->p->key > xt->key) {
				return is_fibheap = false;
			}
			if(xt->child != NULL) {
				is_fibheap = is_fib_heap_children(xt);
				if(!is_fibheap) { return false; }
			}
			xt = xt->right;
		}
		if(xt->right == z->child) {
			if(xt->p->key > xt->key) {
				return is_fibheap = false;
			}
			if(xt->child != NULL) {
				is_fibheap = is_fib_heap_children(xt);
				if(!is_fibheap) { return false; }
			}
		}
	}

	return is_fibheap;
}

bool is_fib_heap(node* z) {
	bool is_fibheap = true;

	node* xt = z;
	if(xt != NULL) {
		while(xt->right != z) {
			is_fibheap = is_fib_heap_children(xt);
			if(!is_fibheap) { return false; }
			xt = xt->right;
		}
		if(xt->right == z) {
			is_fibheap = is_fib_heap_children(xt);
			if(!is_fibheap) { return false; }
		}
	}

	return is_fibheap;
}

node* fib_heap_extract_min(FibHeap* H) {

	node* z = H->min;

	if(z != NULL) {

		//Add each child of z to root list
		node* y = z->child;
		if(y != NULL) {
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
	//Remove x from child list of y and add x to root list of H
	x->left->right = x->right;
	x->right->left = x->left;

	x->right = H->min->right;
	x->left = H->min;

	H->min->right->left = x;
	H->min->right = x;

	//If x is only child set child of parent to null
	if(x == x->right) {
		y->child = NULL;
	}

	y->degree = y->degree - 1;

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

int main(int argc, char* argv[]) {

	//Declarations
	FibHeap H;
	int num_nodes = 357;
	int num_calls = 357;

	//Create nodes
	srand (time(NULL));
	node** v = new node*[num_nodes];
	for(int i = 0; i < num_nodes; ++i) {
		v[i] = new node;
		v[i]->key = rand() % num_nodes;
		std::cout << "keys: " << v[i]->key << std::endl;
	}

	//Insert nodes
	for(int i = 0; i < num_nodes; ++i) {
		fib_heap_insert(&H, v[i]);
	}

	//Extract nodes
	node** z = new node*[num_nodes];
	for(int i = 0; i < num_calls; ++i) {
		z[i] = fib_heap_extract_min(&H);
	}

	bool is_fibheap = is_fib_heap(H.min);

	//Print nodes
	for(int i = 0; i < num_calls; ++i) {
		if(z[i] != NULL) {
		    std::cout << "z " << i << " key: " << z[i]->key << " , degree: " << z[i]->degree << std::endl;
		}
		else {
			std::cout << "z key " << i << " is null" << std::endl;
		}
	}

	std::cout << "is fib heap: " << is_fibheap << std::endl;

	return 0;
}





