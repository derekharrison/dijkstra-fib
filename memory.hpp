/*
 * memory.hpp
 *
 *  Created on: Oct 7, 2021
 *      Author: d-w-h
 */

#ifndef MEMORY_HPP_
#define MEMORY_HPP_

#include "user_types.hpp"

bool** bool2D(const int size);
int** int2D(const int size);
void free_bool2D(bool** p, const int size);
void free_int2D(int** p, const int size);
void free_node_refs(node** p, int size);

#endif /* MEMORY_HPP_ */
