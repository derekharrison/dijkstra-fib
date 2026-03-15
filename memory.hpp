/*
 * memory.hpp
 *
 *  Created on: Oct 7, 2021
 *      Author: d-w-h
 */

#ifndef MEMORY_HPP_
#define MEMORY_HPP_

#include "usertypes.hpp"

Node ** getNodeRef(int n);
void freeNodeRef(Node ** vRef, int size);

#endif /* MEMORY_HPP_ */
