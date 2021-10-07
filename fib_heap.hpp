/*
 * fib_heap.hpp
 *
 *  Created on: Oct 7, 2021
 *      Author: d-w-h
 */

#ifndef FIB_HEAP_HPP_
#define FIB_HEAP_HPP_

#include <vector>

#include "user_types.hpp"

std::vector<int> shortest_reach(int n, std::vector< std::vector<int> >& edges, int s);
void fib_heap_decrease_key(FibHeap* H, node* x, int k);
void fib_heap_insert(FibHeap* H, node* x);
node* fib_heap_extract_min(FibHeap* H);

#endif /* FIB_HEAP_HPP_ */
