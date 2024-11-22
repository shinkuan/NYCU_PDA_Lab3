//############################################################################
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   Cell Class Implementation Header File
//   
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   File Name   : cell.h
//   Release Version : V1.0
//   Description : 
//      
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   Key Features:
//   
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   Author          : shinkuan
//   Creation Date   : 2024-10-31
//   Last Modified   : 2024-10-31
//   Compiler        : g++/clang++
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   Usage Example:
//   #include "cell.h"
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   License: 
//
//############################################################################

#ifndef CELL_H
#define CELL_H

#include <string>
#include <vector>
#include "common.h"

#define CELL_PROPERTY_NOT_FIXED 0
#define CELL_PROPERTY_FIXED 1
#define CELL_PROPERTY_BLOCKAGE 2


class Cell {
public:
    std::string name;
    Point<double> lower_left;
    Size<double> size;
    bool fixed;

    int row_idx;
    int row_height;

    Cell();
    ~Cell();
};


#endif // CELL_H