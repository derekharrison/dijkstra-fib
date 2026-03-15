/*
 * FibHeap.hpp
 *
 *  Created on: Oct 7, 2021
 *      Author: d-w-h
 */

#ifndef FIBHEAP_HPP_
#define FIBHEAP_HPP_

#include <vector>

#include "usertypes.hpp"

void fibHeapDecreaseKey(FibHeap * H, Node * x, int k);
void fibHeapInsert(FibHeap * H, Node * x);
Node * fibHeapExtractMin(FibHeap * H);

#endif /* FIBHEAP_HPP_ */
