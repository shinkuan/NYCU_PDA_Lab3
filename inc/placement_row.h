//############################################################################
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   PlacementRow Class Implementation Header File
//   
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   File Name   : placement_row.h
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
//   #include "placement_row.h"
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   License: 
//
//############################################################################

#ifndef PLACEMENT_ROW_H
#define PLACEMENT_ROW_H

#include <string>
#include "common.h"
#include "cell.h"
#include "intervaltree.h"


struct PlacementRow {
    Point<double> lower_left;
    int num_sites;

    IntervalTree<double, Cell*> itree;
};



#endif // PLACEMENT_ROW_H
