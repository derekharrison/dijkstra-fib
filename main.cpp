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
#include <vector>

const int SETVAR = 314159;

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
#include <vector>

int node_counter = 0;

typedef struct Node {
	Node* left;
	Node* right;
	Node* p;
	Node* child;
	Node* pi;

	std::vector<int> adj_nodes;

	int key;
	int d;
	int degree;
	int index;
	int index_og;
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

void print_root_circle(node* z) {
	node* xt = z;
	if(xt != NULL) {
		if(xt->right != z) {
			while(xt->right != z) {
				std::cout << "xt->key: " << xt->key;
				std::cout << ", xt->degree: " << xt->degree << std::endl;
				xt = xt->right;
			}
			if(xt->right == z) {
				std::cout << "xt->key: " << xt->key;
				std::cout << ", xt->degree: " << xt->degree << std::endl;
			}
		}
		else {
			std::cout << "X == X->RIGHT" << std::endl;
			std::cout << "xt->key: " << xt->key;
			std::cout << ", xt->degree: " << xt->degree << std::endl;
		}
	}
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

void fib_heap_link(FibHeap* H, node* y, node* x) {

	//Make y child of x
	make_child_of(H, y, x);
}

bool numbers_match(node* z);

int counter = 0;

void consolidate(FibHeap* H) {

	counter++;
	double golden = (1.0 + sqrt(5.0)) / 2.0;
	double f = log(H->n) / log(golden);
	int D = floor(f + 0.01) + 1;

	node** A = new node*[D + 2];
	for(int i = 0; i < D + 2; ++i) {
		A[i] = NULL;
	}

	std::cout << "PRINTING ROOT CIRCLE BEFORE: " << counter << std::endl;
	print_root_circle(H->min);

	//Test code
	node_counter = 0;
	int nums_match = numbers_match(H->min);
	std::cout << "numbers match bro before: " << nums_match << std::endl;
	std::cout << "H->n == node_counter before: " << (node_counter == (H->n - 1)) << std::endl;
	//End test code

	node* x = H->min;
	bool there_is_dup = true;
	if(x != NULL) {
		if(x->right != H->min) {

			//Check for duplicate degree
			there_is_dup = false;
			while(x->right != H->min) {
				int d = x->degree;
				if(A[d] != NULL) {
					there_is_dup = true;
				}
				else {
					A[d] = x;
				}
				x = x->right;
			}

			if(x->right == H->min) {
				int d = x->degree;
				if(A[d] != NULL) {
					there_is_dup = true;
				}
				else {
					A[d] = x;
				}
			}

			//Ensure all root nodes have unique degrees
			while(there_is_dup) {
				for(int i = 0; i < D + 2; ++i) {
					A[i] = NULL;
				}

				there_is_dup = false;
				x = H->min;
				while(x->right != H->min) {
					int d = x->degree;
					if(A[d] != NULL && A[d] != x) {
						there_is_dup = true;
						node* y = A[d];
						if(y->key > x->key) {
							//Make y child of x;

							 make_child_of(H, y, x);

							 A[d] = NULL;
							 A[d+1] = x;

							if(y == H->min) {
								H->min = x;
							}
						}
						else {
							//Make x child of y;
							make_child_of(H, x, y);

							A[d] = NULL;
							A[d+1] = y;

							if(x == H->min) {
								H->min = y;
							}

							x = y;
						}
					}
					else {
						A[d] = x;
					}
					x = x->right;
				}

				if(x->right == H->min) {
					int d = x->degree;
					if(A[d] != NULL && A[d] != x) {
						there_is_dup = true;
						node* y = A[d];
						if(y->key > x->key) {

							//Make y child of x;

							make_child_of(H, y, x);

							A[d] = NULL;
							A[d+1] = x;

							if(y == H->min) {
								H->min = x;
							}
						}
						else {
							//Make x child of y;

							make_child_of(H, x, y);

							A[d] = NULL;
							A[d+1] = y;

							if(x == H->min) {
								H->min = y;
							}

							x = y;
						}
					}
					else {
						A[d] = x;
					}
				}
			}
		}
		else {
			int d = x->degree;
			A[d] = x;
		}
	}

	//Test code
	std::cout << "PRINTING ROOT CIRCLE" << std::endl;
	print_root_circle(H->min);

	node_counter = 0;
	nums_match = numbers_match(H->min);
	std::cout << "numbers match bro after: " << nums_match << std::endl;
	std::cout << "H->n == node_counter after: " << (node_counter == (H->n - 1)) << std::endl;
	//End test code

//	//Test code below, remove when done testing
//	for(int i = 0; i < D + 2; ++i) {
//		A[i] = NULL;
//	}
//
//	x = H->min;
//	if(x != NULL) {
//		if(x->right != H->min) {
//			while(x->right != H->min) {
//				int d = x->degree;
//				std::cout << "x->degree: " << d << std::endl;
//				x->p = NULL;
//				if(A[d] != NULL) {
//					std::cout << "THERE IS DUPLICATE" << std::endl;
//					std::cout << "duplicate = " << d << std::endl;
//
//					exit(314159);
//				}
//				A[d] = x;
//				x = x->right;
//			}
//
//			if(x->right == H->min) {
//				int d = x->degree;
//				std::cout << "x->degree: " << d << std::endl;
//				if(A[d] != NULL) {
//					std::cout << "THERE IS DUPLICATE" << std::endl;
//					std::cout << "duplicate = " << d << std::endl;
//
//					exit(314159);
//				}
//				x->p = NULL;
//				A[d] = x;
//			}
//		}
//		else {
//			int d = x->degree;
//			A[d] = x = x;
//		}
//	}
//	//End test code

	//Test code
	std::cout << "START HERE" << std::endl;
	for(int i = 0; i < D + 2; ++i) {
		if(A[i] != NULL) {
		    std::cout << "A[" << i << "].degree: " << A[i]->degree << std::endl;
		}
		else {
			std::cout << "A[" << i << "] is null" << std::endl;
		}
	}

	std::cout << "H->n == node_counter 1: " << (node_counter == (H->n - 1)) << std::endl;
	if((node_counter != (H->n - 1))) { exit(314); }
	std::cout << "node_counter 1: " << (node_counter) << std::endl;
	std::cout << "H->n - 1, 1: " << (H->n - 1) << std::endl;
//	std::cout << "numbers match 1: " << nums_match << std::endl;
	std::cout << "H->min->right 1: " << (H->min->right) << std::endl;
	std::cout << "H->min->right->right 1: " << (H->min->right->right) << std::endl;
	std::cout << "H->min 1: " << (H->min) << std::endl;
	std::cout << "H->min->left 1: " << (H->min->left) << std::endl;
	//End test code

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


	//Test code
	node_counter = 0;
	nums_match = numbers_match(H->min);
	std::cout << "H->n == node_counter 2: " << (node_counter == (H->n - 1)) << std::endl;
	std::cout << "node_counter 2: " << (node_counter) << std::endl;
	std::cout << "H->n - 1, 2: " << (H->n - 1) << std::endl;
	std::cout << "numbers match 2: " << nums_match << std::endl;
	std::cout << "H->min->right 2: " << (H->min->right) << std::endl;
	std::cout << "H->min->right->right 2: " << (H->min->right->right) << std::endl;
	std::cout << "H->min 2: " << (H->min) << std::endl;
	std::cout << "H->min->left 2: " << (H->min->left) << std::endl;
	//End test code
}

void print_child_circle(node* child) {
	node* xt = child;
	if(xt != NULL) {
		if(xt->right != child) {
			while(xt->right != child) {
				std::cout << "xt->child->key: " << xt->key;
				std::cout << ", xt->child->degree: " << xt->degree << std::endl;
				if(xt->child != NULL) {
					std::cout << "xt->child->child->key: " << xt->child->key << std::endl;
				}
				xt = xt->right;
			}
			if(xt->right == child) {
				std::cout << "xt->child->key: " << xt->key;
				std::cout << ", xt->child->degree: " << xt->degree << std::endl;
				if(xt->child != NULL) {
					if(xt->child != NULL) {
						std::cout << "xt->child->child->key: " << xt->child->key << std::endl;
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
	if(xt != NULL) {
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
			node_counter++;
			if(xt->child != NULL) {
				nums_match = numbers_children_match(xt);
				if(!nums_match) { return false; }
			}
			xt = xt->right;
		}
		if(xt->right == z->child) {
			num_of_nodes++;
			node_counter++;
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
			node_counter++;
			nums_match = numbers_children_match(xt);
			if(!nums_match) { return false; }
			xt = xt->right;
		}
		if(xt->right == z) {
			node_counter++;
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
		else {
			std::cout << "y == null" << std::endl;
		}

		//Set all root parents to zero. Remove code when done testing
		node* x_track = H->min;
		while(x_track->right != H->min) {
			x_track->p = NULL;
			x_track = x_track->right;
		}

		if(x_track->right == H->min) {
			x_track->p = NULL;
		}

		//Remove z from root list
		z->left->right = z->right;
		z->right->left = z->left;

		if(z == z->right) {
			H->min = NULL;
			std::cout << "z == z->right" << std::endl;
		}
		else {

			H->min = z->right;
			consolidate(H);
		}

		H->n = H->n - 1;


	}

	return z;

}

bool** bool2D(const int size) {
    bool** p = new bool*[size];

    for(int i = 0; i < size; ++i)
        p[i] = new bool[size];

    return p;
}

int** int2D(const int size) {
    int** p = new int*[size];

    for(int i = 0; i < size; ++i)
        p[i] = new int[size];

    return p;
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

void relax(node* u, node* v, int** w, FibHeap* H) {

    if(v->key > u->key + w[u->index][v->index]) {
        int weight = u->key + w[u->index][v->index];
        fib_heap_decrease_key(H, v, weight);
        v->pi = u;
        v->d = weight;
        v->key = weight;
    }
}

void dijkstra(FibHeap* H) {
	int n = H->n;
	int inf = 3e+8;

    //Initialize heap
    int num_nodes = n;
	node** v_ref = new node*[num_nodes];
	for(int i = 0; i < num_nodes; ++i) {
		v_ref[i] = new node;
        v_ref[i]->key = inf;
        v_ref[i]->pi = NULL;
        v_ref[i]->d = inf;
        v_ref[i]->index = i;
		if(i == 0) {
			v_ref[0]->key = 0;
			v_ref[0]->d = 0;
		}
	}

    int** weight_mat = int2D(n);

    //Populate adjacency and weight matrices
    for(int i = 0; i < n; ++i) {
    	for(int j = 0; j < n; ++j) {
    		if(i != j) {
				bool populate_elem = (rand() % n) > n/2;
				if(populate_elem) {
					weight_mat[i][j] = rand() % n + 1;
					v_ref[i]->adj_nodes.push_back(j);
				}
				else {
					weight_mat[i][j] = inf;
				}
    		}
    		else {
				weight_mat[i][j] = inf;
    		}
    	}
    }

    //Insert nodes into heap
	for(int i = 0; i < num_nodes; ++i) {
		fib_heap_insert(H, v_ref[i]);
	}

    //Perform Dijkstra's algorithm
    while(H->n > 0) {
    	node* u = fib_heap_extract_min(H);

		int num_adj_nodes = u->adj_nodes.size();
		for(int i = 0; i < num_adj_nodes; ++i) {
			int index_ref = u->adj_nodes[i];
			node* v = v_ref[index_ref];
			relax(u, v, weight_mat, H);
		}
    }
}

void set_index_map(int size_graph, int* index_map, int* index_map_inverse, int s) {
	//Make 0 elements point to zero
	index_map[0] = index_map_inverse[0] = 0;

	int index_track = 1;
    for(int i = s; i <= size_graph; ++i) {
        index_map[i] = index_track;
        index_map_inverse[index_track] = i;
        index_track++;
    }
    for(int i = 1; i <= s - 1; ++i) {
        index_map[i] = index_track;
        index_map_inverse[index_track] = i;
        index_track++;
    }
}

void populate_adj_and_weight_hr(int* index_map, int** adj_mat, int** weight_mat, int size_graph, std::vector<std::vector<int>> edges, int s) {

    int** elem_is_set = int2D(size_graph + 1);

    int num_edges = edges.size();
    for(int i = 0; i < num_edges; ++i) {
        int start = index_map[edges[i][0]];
        int end = index_map[edges[i][1]];
        int weight = edges[i][2];
        if(elem_is_set[start][end] != SETVAR) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
            elem_is_set[start][end] = elem_is_set[end][start] = SETVAR;
        }
        else if(elem_is_set[start][end] == SETVAR && weight_mat[start][end] >= weight) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
        }
        adj_mat[start][end] = adj_mat[end][start] = SETVAR;
    }
}

std::vector<int> shortest_reach(int n, std::vector<std::vector<int>> edges, int s) {
	//Declarations
	FibHeap H;
	std::vector<node> rs_S;
    const int inf = 3e+8;

    //Set index map
    int* index_map = new int[n+1];
    int* index_map_inverse = new int[n+1];
    set_index_map(n, index_map, index_map_inverse, s);

    int* index_map_end = new int[n+1];
    for(int i = 0; i < n; ++i) {
        index_map_end[i] = 0;
    }

    //Initialize heap
    int num_nodes = n;
	node** v_ref = new node*[num_nodes];
	for(int i = 0; i < num_nodes; ++i) {
		v_ref[i] = new node;
        v_ref[i]->key = inf;
        v_ref[i]->pi = NULL;
        v_ref[i]->d = inf;
        v_ref[i]->index = i + 1;
        v_ref[i]->index_og = index_map_inverse[i+1];
		if(i == 0) {
			v_ref[0]->key = 0;
			v_ref[0]->d = 0;
		}
		fib_heap_insert(&H, v_ref[i]);
	}

    //Add references to adjacent nodes
    int num_edges = edges.size();
    for(int i = 0; i < num_edges; ++i) {
        int start_index = edges[i][0];
        int end_index = edges[i][1];

        int start_index_reordered = index_map[start_index];
        int end_index_reordered = index_map[end_index];
        v_ref[start_index_reordered-1]->adj_nodes.push_back(end_index_reordered);
        v_ref[end_index_reordered-1]->adj_nodes.push_back(start_index_reordered);
    }

    //Initialize weight and adjacency matrices
    int** adj_mat = int2D(n+1);
    int** weight_mat = int2D(n+1);

    populate_adj_and_weight_hr(index_map, adj_mat, weight_mat, n, edges, s);

    //Perform Dijkstra's algorithm
    while(H.n > 0) {
    	node* u = fib_heap_extract_min(&H);

    	bool is_fibheap = is_fib_heap(H.min);
    	bool nums_match = numbers_match(H.min);



        int num_adj_nodes = u->adj_nodes.size();
        for(int i = 0; i < num_adj_nodes; ++i) {
        	int index_ref = u->adj_nodes[i];
        	node* v = v_ref[index_ref-1];
        	relax(u, v, weight_mat, &H);
        	std::cout << "stuck somewhere else" << std::endl;
        }

    	std::cout << "is fibheap: " << is_fibheap << std::endl;
    	std::cout << "numbers match: " << nums_match << std::endl;
    }

    //Reorder results
    std::vector<int> rs_S_reordered;

    return rs_S_reordered;
}

int main(int argc, char* argv[]) {

	//Declarations
	FibHeap H;
	H.n = 411;
//	int num_nodes = H.n;
//	int num_calls = num_nodes - 1;

	//Execute algorithm
	dijkstra(&H);

	std::cout << "done" << std::endl;

//	//Create nodes
//	srand (time(NULL));
//	node** v = new node*[num_nodes];
//	for(int i = 0; i < num_nodes; ++i) {
//		v[i] = new node;
//		v[i]->key = rand() % num_nodes;
//		std::cout << "keys: " << v[i]->key << std::endl;
//	}
//
//	//Insert nodes
//	for(int i = 0; i < num_nodes; ++i) {
//		fib_heap_insert(&H, v[i]);
//	}
//
//	//Extract nodes
//	node** z = new node*[num_nodes];
//	for(int i = 0; i < num_calls; ++i) {
//		z[i] = fib_heap_extract_min(&H);
//	}
//
//	if(H.min->child != NULL) {
//	    fib_heap_decrease_key(&H, H.min->child, -10);
//	}
//
//	//Print stuff
//	for(int i = 0; i < num_calls; ++i) {
//		if(z[i] != NULL) {
//		    std::cout << "z " << i << " key: " << z[i]->key << " , degree: " << z[i]->degree << std::endl;
//		}
//		else {
//			std::cout << "z key " << i << " is null" << std::endl;
//		}
//	}
//
//	//Test code
//	bool is_fibheap = is_fib_heap(H.min);
//	bool nums_match = numbers_match(H.min);
//	std::cout << "is fib heap: " << is_fibheap << std::endl;
//	std::cout << "numbers match: " << nums_match << std::endl;
//	//End test code

	return 0;
}





