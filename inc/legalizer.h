//############################################################################
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   Legalizer Class Implementation Header File
//   
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   File Name   : legalizer.h
//   Release Version : V1.0
//   Description : 
//      This header file defines the class `Legalizer` for legalizing the placement
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   Key Features:
//   1. `Legalizer` Class:
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
//   #include "legalizer.h"
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   License: 
//
//############################################################################

#ifndef LEGALIZER_H
#define LEGALIZER_H

#define WINDOW_SIZE_X 6
#define WINDOW_SIZE_Y 6

#include <vector>
#include <string>
#include <list>
#include <deque>
#include <unordered_set>
#include <unordered_map>
#include <random>
#include "common.h"
#include "cell.h"
#include "placement_row.h"


struct OptimizeStep {
    Cell* to;
    std::vector<Cell*> from;
};

struct CellLRPack {
    Cell* cell;
    bool left_packed;
    double left_x;
    bool right_packed;
    double right_x;
};

struct InsertionInterval {
    int row;
    double left_x;
    double right_x;
    CellLRPack* left_cell;
    CellLRPack* right_cell;
};

struct InsertionPoint {
    int row;
    double left_x;
    double right_x;
    double best_x;
    std::vector<InsertionInterval*> involved_intervals;
};

struct Endpoint {
    double x;
    bool is_left;
    bool cleared;
    InsertionInterval* ii;
    Endpoint* pair;
};

struct LocalRegionRow {
    int row_idx;
    double left_x;
    double right_x;
    std::vector<InsertionInterval> insertion_intervals;
    IntervalTree<double, Cell*> itree;
};


struct OptimizedResult {\
    OptimizeStep step;
    Point<double> lower_left;
    int success;
    Point<double> affected_lower_left;
    Point<double> affected_upper_right;
    InsertionPoint* best_ip;    
    std::vector<Cell*> lrecored;
    std::vector<Cell*> rrecored;
};


class Legalizer {
private:
    int window_size_x;
    int window_size_y;

    double alpha, beta;
    Point<double> die_lower_left;
    Point<double> die_upper_right;

    double row_start_y;
    double site_height;

    std::unordered_map<std::string, Cell*> cells;
    std::vector<PlacementRow*> prows;

    std::vector<std::vector<int>> bin_score;
    double bin_width;
    double bin_height;

    std::random_device rd;
    std::mt19937 gen;

public:
    Legalizer();
    ~Legalizer();
    std::vector<OptimizeStep> optimize_steps;

    void load_lg(std::string filename);
    void load_opt(std::string filename);
    void dump_loaded_row(std::string filename);
    void dump_optimize_steps(std::string filename);

    void legalize(std::string filename);
    void remove_interval(std::vector<LocalRegionRow>& local_region_rows, int row_idx, double start_x, bool mode, std::unordered_set<int>& removed_intervals);
    double find_lpack_cell_x(std::vector<LocalRegionRow> &local_region_rows, CellLRPack *to_find_clp, std::vector<std::vector<CellLRPack *>> &lrpacks, int row_idx, std::vector<int> &packed_idx);
    int find_lpack_cell(std::vector<LocalRegionRow> &local_region_rows, CellLRPack *to_find_clp, std::vector<std::vector<CellLRPack *>> &lrpacks, int row_idx, std::vector<int> &packed_idx);
    double find_rpack_cell_x(std::vector<LocalRegionRow> &local_region_rows, CellLRPack *to_find_clp, std::vector<std::vector<CellLRPack *>> &lrpacks, int row_idx, std::vector<int> &packed_idx);
    int find_rpack_cell(std::vector<LocalRegionRow> &local_region_rows, CellLRPack *to_find_clp, std::vector<std::vector<CellLRPack *>> &lrpacks, int row_idx, std::vector<int> &packed_idx);
    std::vector<std::vector<CellLRPack*>> find_lrpack(std::vector<LocalRegionRow>& local_region_rows);
    void enumerate_insertion_points(Endpoint& e, int row, std::vector<std::deque<InsertionInterval*>*>& involved_queues, size_t index, std::vector<InsertionInterval*>& current_intervals, double current_min_right_x, std::vector<InsertionPoint*>& result);
    double cost_calculator(InsertionPoint* ip, OptimizeStep step);
    void cell_pusher(Cell* cell, double target_x, std::vector<Cell*>& moved_cells, bool push_left);
    std::vector<Cell*> result_realizer(OptimizedResult result);
    OptimizedResult MLL(OptimizeStep step, bool dump_local_region);
    OptimizedResult _MLL(OptimizeStep step, bool dump_local_region, Point<double> local_region);
};


#endif // LEGALIZER_H