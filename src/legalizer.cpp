#include "legalizer.h"
#include <cfloat>
#include <cstdint>
#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <map>
#include <unordered_set>
#include <set>
#include <list>
#include <queue>
#include <deque>
#include <random>
#include <omp.h>


Legalizer::Legalizer() {
    window_size_x = WINDOW_SIZE_X;
    window_size_y = WINDOW_SIZE_Y;
    alpha = 0.0;
    beta = 0.0;
    die_lower_left = {0.0, 0.0};
    die_upper_right = {0.0, 0.0};
    gen = std::mt19937(rd());
}

Legalizer::~Legalizer() {
    for (auto& cell : cells) {
        delete cell.second;
    }
    for (auto& prow : prows) {
        delete prow;
    }
}

void Legalizer::load_lg(std::string filename) {
    std::cout << "Loading " << filename << std::endl;

    // Load the *.lg file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open block file " << filename << std::endl;
        return;
    }

    row_start_y = DBL_MAX;
    site_height = DBL_MAX;

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string cmd_str;
        iss >> cmd_str;
        if (cmd_str == "Alpha") {
            iss >> alpha;
        } else if (cmd_str == "Beta") {
            iss >> beta;
        } else if (cmd_str == "DieSize") {
            iss >> die_lower_left.x >> die_lower_left.y >> die_upper_right.x >> die_upper_right.y;
        } else if (cmd_str == "PlacementRows") {
            double startX, startY, siteWidth, siteHeight;
            int numSites;
            iss >> startX >> startY >> siteWidth >> siteHeight >> numSites;
            PlacementRow* prow = new PlacementRow;
            prow->lower_left = {startX, startY};
            this->site_height = siteHeight; // Assume all the rows have the same site height
            prow->num_sites = numSites;
            prows.push_back(prow);

            if (startY < row_start_y) {
                row_start_y = startY;
            }
        } else {
            double lowerLeftX, lowerLeftY, width, height;
            std::string fixed;
            iss >> lowerLeftX >> lowerLeftY >> width >> height >> fixed;
            Cell* cell = new Cell;
            cell->lower_left = {lowerLeftX, lowerLeftY};
            cell->size = {width, height};
            cell->fixed = fixed == "FIX";
            cell->name = cmd_str;
            cells[cmd_str] = cell;
        }
    }

    file.close();

    // std::cout << "Sorting placement rows" << std::endl;
    std::sort(prows.begin(), prows.end(), [](PlacementRow* a, PlacementRow* b) {
        return a->lower_left.y < b->lower_left.y;
    });

    // std::cout << "Building interval trees" << std::endl;
    for (auto& cell : cells) {
        int prow_idx = (cell.second->lower_left.y - row_start_y) / this->site_height;
        if (prow_idx < 0) {
            std::cerr << "Error: Cell " << cell.second->name << " is out of the placement rows" << std::endl;
            prow_idx = 0;
        } else if (prow_idx >= (int)(prows.size())) {
            std::cerr << "Error: Cell " << cell.second->name << " is out of the placement rows" << std::endl;
            prow_idx = prows.size() - 1;
        }
        int cell_height = (int)ceil(cell.second->size.height / this->site_height);
        bool inserted = false;
        for (int i = 0; i < cell_height; i++) {
            if (prow_idx + i >= (int)prows.size()) {
                break;
            }
            PlacementRow* prow = prows[prow_idx + i];
            if (cell.second->lower_left.x >= prow->lower_left.x && cell.second->lower_left.x < prow->lower_left.x + (double)prow->num_sites) {
                prow->itree.insert({cell.second->lower_left.x, cell.second->lower_left.x + cell.second->size.width, cell.second});
                inserted = true;
            }
            if (!inserted) {
                std::cerr << "Error: Cell " << cell.second->name << " is out of the placement rows" << std::endl;
            }
        }
        cell.second->row_idx = prow_idx;
        cell.second->row_height = cell_height;
    }
    std::cout << "Loaded " << cells.size() << " cells and " << prows.size() << " placement rows" << std::endl;
}

void Legalizer::dump_loaded_row(std::string filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open output file " << filename << std::endl;
        return;
    }

    file << "Alpha " << alpha << std::endl;
    file << "Beta " << beta << std::endl;
    file << "DieSize " << static_cast<int>(die_lower_left.x) << " " << static_cast<int>(die_lower_left.y) << " " << static_cast<int>(die_upper_right.x) << " " << static_cast<int>(die_upper_right.y) << std::endl;

    for (int i = 0; i < (int)prows.size(); i++) {
        PlacementRow* prow = prows[i];
        file << "PlacementRows " << static_cast<int>(prow->lower_left.x) << " " << static_cast<int>(prow->lower_left.y) << " " << static_cast<int>(1) << " " << static_cast<int>(this->site_height) << " " << prow->num_sites << std::endl;
        for (auto interval_cell : prow->itree.intervals()) {
            Cell* cell = interval_cell.value;
            if (cell->row_idx == i) {
                file << cell->name << " " << static_cast<int>(cell->lower_left.x) << " " << static_cast<int>(cell->lower_left.y) << " " << static_cast<int>(cell->size.width) << " " << static_cast<int>(cell->size.height) << " " << (cell->fixed ? "FIX" : "NOTFIX") << std::endl;
            }
        }
    }
}

void Legalizer::load_opt(std::string filename) {
    std::cout << "Loading " << filename << std::endl;

    // Load the *.opt file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open block file " << filename << std::endl;
        return;
    }

    optimize_steps.clear();
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string cmd_str;
        iss >> cmd_str;
        if (cmd_str == "Banking_Cell:") {
            OptimizeStep step;
            std::string from;
            iss >> from;
            while (from != "-->") {
                step.from.push_back(cells[from]);
                iss >> from;
            }
            std::string to;
            double x, y, w, h;
            iss >> to;
            iss >> x >> y >> w >> h;
            Cell* to_cell = new Cell;
            to_cell->name = to;
            to_cell->lower_left = {x, y};
            to_cell->size = {w, h};
            to_cell->fixed = false;

            int prow_idx = (to_cell->lower_left.y - row_start_y) / this->site_height;
            if (prow_idx < 0) {
                std::cerr << "Warning: Banked Cell " << to_cell->name << "'s target location is out of the placement rows" << std::endl;
                prow_idx = 0;
            } else if (prow_idx >= (int)(prows.size())) {
                std::cerr << "Warning: Banked Cell " << to_cell->name << "'s target location is out of the placement rows" << std::endl;
                prow_idx = prows.size() - 1;
            }
            int cell_height = (int)ceil(to_cell->size.height / this->site_height);
            to_cell->row_idx = prow_idx;
            to_cell->row_height = cell_height;

            cells[to] = to_cell;
            step.to = to_cell;
            optimize_steps.push_back(step);
        }
    }

    file.close();

    std::cout << "Loaded " << optimize_steps.size() << " optimize steps" << std::endl;
    std::cout << "Creating Bins" << std::endl;
    bin_width = optimize_steps[0].to->size.width * WINDOW_SIZE_X;
    bin_height = site_height * WINDOW_SIZE_Y;
    bin_score.resize((int)ceil((die_upper_right.x - die_lower_left.x) / bin_width), std::vector<int>((int)ceil((die_upper_right.y - die_lower_left.y) / bin_height), 0));
    for (int i = 0; i < (int)bin_score.size(); i++) {
        for (int j = 0; j < (int)bin_score[i].size(); j++) {
            bin_score[i][j] = 0;
        }
    }
}

void Legalizer::dump_optimize_steps(std::string filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open output file " << filename << std::endl;
        return;
    }
    
    for (auto& step : optimize_steps) {
        file << "Banking_Cell:";
        for (auto& from : step.from) {
            file << " " << from->name;
        }
        file << " --> " << step.to->name << " " << static_cast<int>(step.to->lower_left.x) << " " << static_cast<int>(step.to->lower_left.y) << " " << static_cast<int>(step.to->size.width) << " " << static_cast<int>(step.to->size.height) << std::endl;
    }
}

void Legalizer::legalize(std::string filename) {
    std::cout << "Legalizing" << std::endl;

    std::ofstream file(filename);
    bool dump_local_region = true;
    int cnt = 0;
    for (auto& step : optimize_steps) {
        std::cout << "Optimizing " << step.to->name << std::endl;
        dump_local_region = cnt == -1;
        OptimizedResult result = MLL(step, dump_local_region);
        std::vector<Cell*> moved_cells = result_realizer(result);
        delete result.best_ip;
        file << static_cast<int>(result.lower_left.x) << " " << static_cast<int>(result.lower_left.y) << std::endl;
        file << moved_cells.size() << std::endl;
        for (auto& cell : moved_cells) {
            file << cell->name << " " << static_cast<int>(cell->lower_left.x) << " " << static_cast<int>(cell->lower_left.y) << std::endl;
        }
        cnt++;
    }
    file.close();

    std::cout << "Done!" << std::endl;
}

void Legalizer::remove_interval(std::vector<LocalRegionRow>& local_region_rows, int row_idx, double start_x, bool mode, std::unordered_set<int>& removed_intervals) {
    // mode = true: remove left interval
    // mode = false: remove right interval
    if (row_idx < 0 || row_idx >= (int)local_region_rows.size()) {
        // std::cerr << "Error: Row index out of range" << std::endl;
        return;
    }
    struct Remove {
        bool mode;
        double x;
    };
    std::unordered_map<int, Remove> to_be_removed;
    if (mode) {
        for (auto& interval_to_be_removed : local_region_rows[row_idx].itree.findOverlappingIntervals({(double)local_region_rows[row_idx].left_x, start_x, nullptr}, false)) {
            Cell* cell = interval_to_be_removed.value;
            int cell_y_low_diff = cell->row_idx - local_region_rows[row_idx].row_idx;
            int cell_y_high_diff = cell->row_idx + cell->row_height - 1 - local_region_rows[row_idx].row_idx;
            for (int i = row_idx + cell_y_low_diff; i < row_idx; i++) {
                if (i < 0) {
                    continue;
                }
                if (local_region_rows[i].left_x < cell->lower_left.x + cell->size.width) {
                    if (to_be_removed.count(i) == 0) {
                        to_be_removed[i] = {true, cell->lower_left.x + cell->size.width};
                    } else {
                        if (to_be_removed[i].x < cell->lower_left.x + cell->size.width) {
                            to_be_removed[i].x = cell->lower_left.x + cell->size.width;
                        }
                    }
                }
            }
            local_region_rows[row_idx].itree.remove(interval_to_be_removed);
            for (int i = row_idx + 1; i <= row_idx + cell_y_high_diff; i++) {
                if (i >= (int)local_region_rows.size()) {
                    continue;
                }
                if (local_region_rows[i].left_x < cell->lower_left.x + cell->size.width) {
                    if (to_be_removed.count(i) == 0) {
                        to_be_removed[i] = {true, cell->lower_left.x + cell->size.width};
                    } else {
                        if (to_be_removed[i].x < cell->lower_left.x + cell->size.width) {
                            to_be_removed[i].x = cell->lower_left.x + cell->size.width;
                        }
                    }
                }
            }
        }
        if (local_region_rows[row_idx].left_x < start_x) {
            local_region_rows[row_idx].left_x = start_x;
        }
        for (auto& remove : to_be_removed) {
            remove_interval(local_region_rows, remove.first, remove.second.x, true, removed_intervals);
        }
    } else {
        for (auto& interval_to_be_removed : local_region_rows[row_idx].itree.findOverlappingIntervals({start_x, (double)local_region_rows[row_idx].right_x, nullptr}, false)) {
            Cell* cell = interval_to_be_removed.value;
            int cell_y_low_diff = cell->row_idx - local_region_rows[row_idx].row_idx;
            int cell_y_high_diff = cell->row_idx + cell->row_height - 1 - local_region_rows[row_idx].row_idx;
            for (int i = row_idx + cell_y_low_diff; i < row_idx; i++) {
                if (i < 0) {
                    continue;
                }
                if (local_region_rows[i].right_x > cell->lower_left.x) {
                    if (to_be_removed.count(i) == 0) {
                        to_be_removed[i] = {false, cell->lower_left.x};
                    } else {
                        if (to_be_removed[i].x > cell->lower_left.x) {
                            to_be_removed[i].x = cell->lower_left.x;
                        }
                    }
                }
            }
            local_region_rows[row_idx].itree.remove(interval_to_be_removed);
            for (int i = row_idx + 1; i <= row_idx + cell_y_high_diff; i++) {
                if (i >= (int)local_region_rows.size()) {
                    continue;
                }
                if (local_region_rows[i].right_x > cell->lower_left.x) {
                    if (to_be_removed.count(i) == 0) {
                        to_be_removed[i] = {false, cell->lower_left.x};
                    } else {
                        if (to_be_removed[i].x > cell->lower_left.x) {
                            to_be_removed[i].x = cell->lower_left.x;
                        }
                    }
                }
            }
        }
        if (local_region_rows[row_idx].right_x > start_x) {
            local_region_rows[row_idx].right_x = start_x;
        }
        for (auto& remove : to_be_removed) {
            remove_interval(local_region_rows, remove.first, remove.second.x, false, removed_intervals);
        }
    }
}

double Legalizer::find_lpack_cell_x(std::vector<LocalRegionRow>& local_region_rows, CellLRPack* to_find_clp, std::vector<std::vector<CellLRPack*>>& lrpacks, int row_idx, std::vector<int>& packed_idx) {
    while (lrpacks[row_idx][packed_idx[row_idx]] != to_find_clp) {
        find_lpack_cell(local_region_rows, lrpacks[row_idx][packed_idx[row_idx]], lrpacks, row_idx, packed_idx);
    }
    packed_idx[row_idx]++;
    if (packed_idx[row_idx]-1 == 0) {
        return local_region_rows[row_idx].left_x;
    } else {
        return lrpacks[row_idx][packed_idx[row_idx] - 2]->left_x + lrpacks[row_idx][packed_idx[row_idx] - 2]->cell->size.width;
    }
}

int Legalizer::find_lpack_cell(std::vector<LocalRegionRow>& local_region_rows, CellLRPack* to_find_clp, std::vector<std::vector<CellLRPack*>>& lrpacks, int row_idx, std::vector<int>& packed_idx) {
    if (to_find_clp->left_packed) {
        return 1;
    }
    double max_x = DBL_MIN;
    for (int i = row_idx; i < row_idx + to_find_clp->cell->row_height; i++) {
        max_x = std::max(max_x, find_lpack_cell_x(local_region_rows, to_find_clp, lrpacks, i, packed_idx));
    }
    to_find_clp->left_packed = true;
    to_find_clp->left_x = max_x;
    return 0;
}

double Legalizer::find_rpack_cell_x(std::vector<LocalRegionRow>& local_region_rows, CellLRPack* to_find_clp, std::vector<std::vector<CellLRPack*>>& lrpacks, int row_idx, std::vector<int>& packed_idx) {
    while (lrpacks[row_idx][packed_idx[row_idx]] != to_find_clp) {
        find_rpack_cell(local_region_rows, lrpacks[row_idx][packed_idx[row_idx]], lrpacks, row_idx, packed_idx);
    }
    packed_idx[row_idx]--;
    if (packed_idx[row_idx]+1 == (int)lrpacks[row_idx].size()-1) {
        return local_region_rows[row_idx].right_x - to_find_clp->cell->size.width;
    } else {
        return lrpacks[row_idx][packed_idx[row_idx]+2]->right_x - to_find_clp->cell->size.width;
    }
}

int Legalizer::find_rpack_cell(std::vector<LocalRegionRow>& local_region_rows, CellLRPack* to_find_clp, std::vector<std::vector<CellLRPack*>>& lrpacks, int row_idx, std::vector<int>& packed_idx) {
    if (to_find_clp->right_packed) {
        return 1;
    }
    double min_x = DBL_MAX;
    for (int i = row_idx; i < row_idx + to_find_clp->cell->row_height; i++) {
        min_x = std::min(min_x, find_rpack_cell_x(local_region_rows, to_find_clp, lrpacks, i, packed_idx));
    }
    to_find_clp->right_packed = true;
    to_find_clp->right_x = min_x;
    return 0;
}

std::vector<std::vector<CellLRPack*>> Legalizer::find_lrpack(std::vector<LocalRegionRow>& local_region_rows) {
    std::vector<std::vector<CellLRPack*>> lrpacks(local_region_rows.size());
    std::list<CellLRPack*> temp_lrpacks;
    std::vector<int> packed_idx(local_region_rows.size(), -1);
    for (int i = 0; i < (int)local_region_rows.size(); i++) {
        LocalRegionRow& lrow = local_region_rows[i];
        for (auto& interval_cell : lrow.itree.intervals()) {
            Cell* cell = interval_cell.value;
            if ((cell->row_idx) == (lrow.row_idx)) {
                CellLRPack* clp = new CellLRPack;
                clp->cell = cell;
                clp->left_packed = false;
                clp->left_x = DBL_MAX;
                clp->right_packed = false;
                clp->right_x = DBL_MIN;
                lrpacks[i].push_back(clp);
                for (int j = i+1; j < i + cell->row_height; j++) {
                    temp_lrpacks.push_back(clp);
                }
            } else {
                for (auto it = temp_lrpacks.begin(); it != temp_lrpacks.end(); ++it) {
                    CellLRPack* clp = *it;
                    if (clp->cell == cell) {
                        lrpacks[i].push_back(clp);
                        temp_lrpacks.erase(it);
                        break;
                    }
                }
            }
        }
    }

    std::vector<int> packed_idx_left(local_region_rows.size(), 0);
    for (int i = 0; i < (int)local_region_rows.size(); i++) {
        for (auto& clp : lrpacks[i]) {
            find_lpack_cell(local_region_rows, clp, lrpacks, i, packed_idx_left);
        }
    }

    std::vector<int> packed_idx_right(local_region_rows.size());
    for (int i = 0; i < (int)local_region_rows.size(); i++) {
        packed_idx_right[i] = (int)lrpacks[i].size() - 1;
    }
    for (int i = 0; i < (int)local_region_rows.size(); i++) {
        for (auto it = lrpacks[i].rbegin(); it != lrpacks[i].rend(); ++it) {
            CellLRPack* clp = *it;
            find_rpack_cell(local_region_rows, clp, lrpacks, i, packed_idx_right);
        }
    }

    return lrpacks;
}

void Legalizer::enumerate_insertion_points(
    Endpoint& e, 
    int row,
    std::vector<std::deque<InsertionInterval*>*>& involved_queues, 
    size_t index, 
    std::vector<InsertionInterval*>& current_intervals, 
    double current_min_right_x,
    std::vector<InsertionPoint*>& result
) {
    if (index == involved_queues.size()) {
        if (e.x > current_min_right_x) {
            return;
        }
        InsertionPoint* ip = new InsertionPoint;
        ip->row = row;
        ip->left_x = e.x;
        ip->right_x = current_min_right_x;
        ip->involved_intervals = current_intervals;
        result.push_back(ip);
        return;
    }

    for (InsertionInterval* ii : *involved_queues[index]) {
        current_intervals.push_back(ii);
        double new_min_right_x = std::min(current_min_right_x, ii->right_x);
        enumerate_insertion_points(e, row, involved_queues, index+1, current_intervals, new_min_right_x, result);
        current_intervals.pop_back();
    }
}

double Legalizer::cost_calculator(InsertionPoint* ip, OptimizeStep step) {
    double cost = 0.0;

    struct CellPoint {
        double x;
        bool is_left;
    };
    std::vector<CellPoint> cell_points;
    std::vector<Cell*> lrecored;
    std::vector<Cell*> rrecored;

    for (auto& ii : ip->involved_intervals) {
        if (ii->left_cell != nullptr) {
            Cell* lcell = ii->left_cell->cell;
            if (std::find(lrecored.begin(), lrecored.end(), lcell) == lrecored.end()) {
                lrecored.push_back(lcell);
                cell_points.push_back({(lcell->lower_left.x + lcell->size.width), true});
            }
        }
        if (ii->right_cell != nullptr) {
            Cell* rcell = ii->right_cell->cell;
            if (std::find(rrecored.begin(), rrecored.end(), rcell) == rrecored.end()) {
                rrecored.push_back(rcell);
                cell_points.push_back({(rcell->lower_left.x - step.to->size.width), false});
            }
        }
    }

    std::sort(cell_points.begin(), cell_points.end(), [](CellPoint a, CellPoint b) {
        return a.x < b.x;
    });

    double best_x = ip->left_x;
    int x_cost = 0;
    int best_x_cost = INT_MAX;
    for (auto& cp : cell_points) {
        if (cp.is_left) x_cost--;
        else            x_cost++;
        if (ip->left_x <= cp.x && cp.x <= ip->right_x) {
            if (x_cost < best_x_cost) {
                best_x = cp.x;
                best_x_cost = x_cost;
            }
        }
    }
    ip->best_x = best_x;

    for (auto& lcell : lrecored) {
        if (lcell == nullptr) continue;
        if (lcell->lower_left.x + lcell->size.width > best_x) {
            cost += alpha + beta * (lcell->lower_left.x + lcell->size.width - best_x);
        }
    }
    for (auto& rcell : rrecored) {
        if (rcell == nullptr) continue;
        if (rcell->lower_left.x - step.to->size.width < best_x) {
            cost += alpha + beta * (best_x - rcell->lower_left.x + step.to->size.width);
        }
    }
    cost += beta * (abs(best_x - step.to->lower_left.x) + abs(ip->row - step.to->row_idx) * this->site_height);
    return cost;
}

void Legalizer::cell_pusher(Cell* cell, double target_x, std::vector<Cell*>& moved_cells, bool push_left) {
    if (push_left) {
        if (cell->lower_left.x <= target_x) {
            return;
        }
        for (int i = 0; i < cell->row_height; i++) {
            Interval<double, Cell*>* pred = prows[cell->row_idx + i]->itree.findPredecessor({cell->lower_left.x, cell->lower_left.x + cell->size.width, cell});
            std::vector<Interval<double, Cell*>> intervals = prows[cell->row_idx + i]->itree.intervals();
            if (pred == nullptr) {
                prows[cell->row_idx + i]->itree.remove({cell->lower_left.x, cell->lower_left.x + cell->size.width, cell});
                prows[cell->row_idx + i]->itree.insert({target_x, target_x + cell->size.width, cell});
                continue;
            }
            cell_pusher(pred->value, target_x-pred->value->size.width, moved_cells, push_left);
            prows[cell->row_idx + i]->itree.remove({cell->lower_left.x, cell->lower_left.x + cell->size.width, cell});
            prows[cell->row_idx + i]->itree.insert({target_x, target_x + cell->size.width, cell});
        }
        cell->lower_left.x = target_x;
        // std::cout << "Pushing Left " << cell->name << " to " << target_x << std::endl;
        if (std::find(moved_cells.begin(), moved_cells.end(), cell) == moved_cells.end()) {
            moved_cells.push_back(cell);
        }
    } else {
        if (cell->lower_left.x >= target_x) {
            return;
        }
        for (int i = 0; i < cell->row_height; i++) {
            Interval<double, Cell*>* succ = prows[cell->row_idx + i]->itree.findSuccessor({cell->lower_left.x, cell->lower_left.x + cell->size.width, cell});
            std::vector<Interval<double, Cell*>> intervals = prows[cell->row_idx + i]->itree.intervals();
            if (succ == nullptr) {
                prows[cell->row_idx + i]->itree.remove({cell->lower_left.x, cell->lower_left.x + cell->size.width, cell});
                prows[cell->row_idx + i]->itree.insert({target_x, target_x + cell->size.width, cell});
                continue;
            }
            cell_pusher(succ->value, target_x + cell->size.width, moved_cells, push_left);
            prows[cell->row_idx + i]->itree.remove({cell->lower_left.x, cell->lower_left.x + cell->size.width, cell});
            prows[cell->row_idx + i]->itree.insert({target_x, target_x + cell->size.width, cell});
        }
        cell->lower_left.x = target_x;
        // std::cout << "Pushing Right " << cell->name << " to " << target_x << std::endl;
        if (std::find(moved_cells.begin(), moved_cells.end(), cell) == moved_cells.end()) {
            moved_cells.push_back(cell);
        }
    }
}

std::vector<Cell*> Legalizer::result_realizer(OptimizedResult result) {
    std::vector<Cell*> moved_cells;

    // Left
    for (auto& lcell : result.lrecored) {
        cell_pusher(lcell, result.best_ip->best_x - lcell->size.width, moved_cells, true);
    }

    // Right
    for (auto& rcell : result.rrecored) {
        cell_pusher(rcell, result.best_ip->best_x + result.step.to->size.width, moved_cells, false);
    }

    result.step.to->lower_left.x = result.best_ip->best_x;
    result.step.to->lower_left.y = result.best_ip->row * this->site_height + row_start_y;
    result.step.to->row_idx = result.best_ip->row;
    result.step.to->fixed = false;
    for (int r = result.step.to->row_idx; r < result.step.to->row_idx + result.step.to->row_height; r++) {
        if (r < 0 || r >= (int)prows.size()) continue;
        PlacementRow* prow = prows[r];
        prow->itree.insert({(result.best_ip->best_x), result.best_ip->best_x + result.step.to->size.width, result.step.to});
    }
    return moved_cells;
}

OptimizedResult Legalizer::MLL(OptimizeStep step, bool dump_local_region = false) {
    // Remove Cells
    for (auto& from : step.from) {
        for (int i = 0; i < from->row_height; i++) {
            prows[from->row_idx + i]->itree.remove({from->lower_left.x, from->lower_left.x + from->size.width, from});
        }
        cells.erase(from->name);
        delete from;
    }

    Point<double> original_lower_left = step.to->lower_left;

    window_size_x = WINDOW_SIZE_X;
    window_size_y = WINDOW_SIZE_Y;
    OptimizedResult result;
    for (int i = 0; i < 7; i++) {
        // std::cout << "Trying " << step.to->lower_left.x << " " << step.to->lower_left.y << std::endl;
        result = _MLL(step, dump_local_region, original_lower_left);
        if (result.success > 0) {
            return result;
        }
        if (result.success == -1) {
            i--;
        } else {
            window_size_x += WINDOW_SIZE_X;
            window_size_y += WINDOW_SIZE_Y;
        }
    }
    window_size_x = WINDOW_SIZE_X;
    window_size_y = WINDOW_SIZE_Y;
    // std::cout << "Local Failed, moving xy" << std::endl;

    int try_per_step = 100;
    double sigma_x = WINDOW_SIZE_X * step.to->size.width;
    double sigma_y = WINDOW_SIZE_Y * this->site_height;
    double max_sigma_x = (die_upper_right.x - die_lower_left.x);
    double max_sigma_y = (die_upper_right.y - die_lower_left.y);
    while (result.success <= 0) {
        bool found = false;

        double local_sigma_x = sigma_x;
        double local_sigma_y = sigma_y;

        #pragma omp parallel
        {
            std::mt19937 thread_gen;
            thread_gen.seed(gen() + omp_get_thread_num());

            std::normal_distribution<double> dist_x(original_lower_left.x, local_sigma_x);
            std::normal_distribution<double> dist_y(original_lower_left.y, local_sigma_y);

            #pragma omp for
            for (int i = 0; i < try_per_step; i++) {
                if (found) continue;
                
                Point<double> local_region_lower_left = {dist_x(thread_gen), dist_y(thread_gen)};
                int local_region_row_idx = (int)((local_region_lower_left.y - row_start_y) / this->site_height);

                if (
                    local_region_lower_left.y < row_start_y ||
                    local_region_lower_left.y > row_start_y + (double)prows.size() * this->site_height ||
                    local_region_lower_left.x < prows[local_region_row_idx]->lower_left.x ||
                    local_region_lower_left.x > prows[local_region_row_idx]->lower_left.x + (double)prows[local_region_row_idx]->num_sites
                ) {
                    continue;
                }

                OptimizedResult local_result = _MLL(step, dump_local_region, local_region_lower_left);

                if (local_result.success > 0) {
                    #pragma omp critical
                    {
                        if (!found || local_result.success > result.success) {
                            result = local_result;
                            found = true;
                        }
                    }
                    #pragma omp cancel for
                }

                #pragma omp cancellation point for
            }
        } 

        if (found) {
            if (window_size_x > WINDOW_SIZE_X) {
                window_size_x -= WINDOW_SIZE_X;
            }
            if (window_size_y > WINDOW_SIZE_Y) {
                window_size_y -= WINDOW_SIZE_Y;
            }
            return result;
        }

        // std::cout << "Failed with sigma_x = " << sigma_x << " sigma_y = " << sigma_y << std::endl;
        sigma_x *= 1.3;
        sigma_y *= 1.3;
        if (sigma_x > max_sigma_x || sigma_y > max_sigma_y) {
            window_size_x += WINDOW_SIZE_X;
            window_size_y += WINDOW_SIZE_Y;
            sigma_x = max_sigma_x;
            sigma_y = max_sigma_y;
        }
    }

    return result;
}

// Multi-row Local Legalization algorithm (MLL)
OptimizedResult Legalizer::_MLL(OptimizeStep step, bool dump_local_region, Point<double> local_region) {
    OptimizedResult result;
    result.step = step;
    result.lower_left = {0.0, 0.0};
    result.success = 0;
    result.affected_lower_left = {0.0, 0.0};
    result.affected_upper_right = {0.0, 0.0};
    result.best_ip = nullptr;

    // Extract the local region
    int local_region_row = (int)((local_region.y - row_start_y) / this->site_height);
    int local_region_start_row = local_region_row - window_size_y;
    int local_region_end_row = local_region_row + step.to->row_height + window_size_y;
    if (local_region_start_row < 0) {
        local_region_start_row = 0;
        local_region_end_row += -local_region_start_row;
    }
    if (local_region_end_row > (int)prows.size()) {
        local_region_start_row -= local_region_end_row - (int)prows.size();
        local_region_end_row = (int)prows.size();
    }
    if (local_region_start_row < 0) {
        local_region_start_row = 0;
    }

    double local_region_x_low = local_region.x - window_size_x * step.to->size.width;
    double local_region_x_high = local_region.x + step.to->size.width + window_size_x * step.to->size.width;
    if (local_region_x_low < prows[local_region_row]->lower_left.x) {
        local_region_x_low = prows[local_region_row]->lower_left.x;
        local_region_x_high = local_region_x_low + step.to->size.width + window_size_x * step.to->size.width;
    }
    if (local_region_x_high > prows[local_region_row]->lower_left.x + (double)prows[local_region_row]->num_sites) {
        local_region_x_high = prows[local_region_row]->lower_left.x + (double)prows[local_region_row]->num_sites;
        local_region_x_low = local_region_x_high - step.to->size.width - window_size_x * step.to->size.width;
    }
    if (local_region_x_low < prows[local_region_row]->lower_left.x) {
        local_region_x_low = prows[local_region_row]->lower_left.x;
    }

    Interval<double, Cell*> local_region_x = {
        local_region_x_low,
        local_region_x_high,
        nullptr
    };

    struct RemoveInterval {
        bool remove_left;
        double left_to_be_removed;
        bool remove_right;
        double right_to_be_removed;
    };

    int local_region_height = local_region_end_row - local_region_start_row;
    std::vector<RemoveInterval> remove_intervals(local_region_height);
    std::vector<LocalRegionRow> local_region_rows(local_region_height);
    for (int ridx = local_region_start_row; ridx < local_region_end_row; ridx++) {
        PlacementRow* prow = prows[ridx];
        LocalRegionRow lrow;
        lrow.row_idx = ridx;
        lrow.left_x = (local_region_x.low < prow->lower_left.x ? prow->lower_left.x : local_region_x.low);
        lrow.right_x = (local_region_x.high > prow->lower_left.x + (double)prow->num_sites ? prow->lower_left.x + (double)prow->num_sites : local_region_x.high);

        // Interval tree
        RemoveInterval remove_interval = {false, DBL_MIN, false, DBL_MAX};
        for (auto& interval_cell : prow->itree.findOverlappingIntervals(local_region_x, false)) {
            Cell* cell = interval_cell.value;
            if ((cell->row_height + cell->row_idx >= local_region_end_row) ||
                (cell->row_idx < local_region_start_row)) {
                if (cell->lower_left.x < local_region.x) {
                    remove_interval.remove_left = true;
                    double tbr = cell->lower_left.x + cell->size.width;
                    if (tbr > remove_interval.left_to_be_removed) {
                        remove_interval.left_to_be_removed = tbr;
                    }
                } else {
                    remove_interval.remove_right = true;
                    double tbr = cell->lower_left.x;
                    if (tbr < remove_interval.right_to_be_removed) {
                        remove_interval.right_to_be_removed = tbr;
                    }
                }
                continue;
            }
            if (cell->lower_left.x < lrow.left_x) {
                remove_interval.remove_left = true;
                double tbr = cell->lower_left.x + cell->size.width;
                if (tbr > remove_interval.left_to_be_removed) {
                    remove_interval.left_to_be_removed = tbr;
                }
                continue;
            }
            if (cell->lower_left.x + cell->size.width > lrow.right_x) {
                remove_interval.remove_right = true;
                double tbr = cell->lower_left.x;
                if (tbr < remove_interval.right_to_be_removed) {
                    remove_interval.right_to_be_removed = tbr;
                }
                continue;
            }
            if (cell->fixed) {
                if (cell->lower_left.x < local_region.x) {
                    remove_interval.remove_left = true;
                    double tbr = cell->lower_left.x + cell->size.width;
                    if (tbr > remove_interval.left_to_be_removed) {
                        remove_interval.left_to_be_removed = tbr;
                    }
                } else {
                    remove_interval.remove_right = true;
                    double tbr = cell->lower_left.x;
                    if (tbr < remove_interval.right_to_be_removed) {
                        remove_interval.right_to_be_removed = tbr;
                    }
                }
                continue;
            }
            // if cell in step.from
            // if (std::find(step.from.begin(), step.from.end(), cell) != step.from.end()) {
            //     continue;
            // }
            lrow.itree.insert({cell->lower_left.x, cell->lower_left.x + cell->size.width, cell});
        }
        remove_intervals[ridx - local_region_start_row] = remove_interval;
        local_region_rows[ridx - local_region_start_row] = lrow;
    }

    // Remove the intervals
    for (int ridx = local_region_start_row; ridx < local_region_end_row; ridx++) {
        if (remove_intervals[ridx - local_region_start_row].remove_left) {
            if (remove_intervals[ridx - local_region_start_row].left_to_be_removed > local_region_rows[ridx - local_region_start_row].left_x) {
                std::unordered_set<int> removed_intervals = {ridx - local_region_start_row};
                remove_interval(local_region_rows, ridx - local_region_start_row, remove_intervals[ridx - local_region_start_row].left_to_be_removed, true, removed_intervals);
            }
        }
        if (remove_intervals[ridx - local_region_start_row].remove_right) {
            if (remove_intervals[ridx - local_region_start_row].right_to_be_removed < local_region_rows[ridx - local_region_start_row].right_x) {
                std::unordered_set<int> removed_intervals = {ridx - local_region_start_row};
                remove_interval(local_region_rows, ridx - local_region_start_row, remove_intervals[ridx - local_region_start_row].right_to_be_removed, false, removed_intervals);
            }
        }
    }

    // Dump the local region
    if (dump_local_region) {
        std::ofstream file("local_region.txt");
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open output file local_region.txt" << std::endl;
            return result;
        }


        file << "RowStartY " << row_start_y << std::endl;
        file << "RowHeight " << this->site_height << std::endl;
        for (auto& lrow : local_region_rows) {
            file << "Row " << lrow.row_idx << " " << lrow.left_x << " " << lrow.right_x << std::endl;
            for (auto& interval_cell : lrow.itree.intervals()) {
                Cell* cell = interval_cell.value;
                if ((int)(cell->lower_left.y - row_start_y) == (int)(lrow.row_idx * this->site_height)) {
                    file << cell->name << " " << static_cast<int>(cell->lower_left.x) << " " << static_cast<int>(cell->lower_left.y) << " " << static_cast<int>(cell->size.width) << " " << static_cast<int>(cell->size.height) << std::endl;
                }
            }
        }
    }

    // Find the LR Pack
    std::vector<std::vector<CellLRPack*>> lrpacks = find_lrpack(local_region_rows);

    // Dump the LR Pack
    if (dump_local_region) {
        std::ofstream file("lpack.txt");
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open output file lrpack.txt" << std::endl;
            return result;
        }

        file << "RowStartY " << row_start_y << std::endl;
        file << "RowHeight " << this->site_height << std::endl;
        for (int i = 0; i < (int)local_region_rows.size(); i++) {
            file << "Row " << local_region_rows[i].row_idx << " " << local_region_rows[i].left_x << " " << local_region_rows[i].right_x << std::endl;
            for (auto& clp : lrpacks[i]) {
                file << clp->cell->name << " " << clp->left_x << " " << static_cast<int>(clp->cell->lower_left.y) << " " << static_cast<int>(clp->cell->size.width) << " " << static_cast<int>(clp->cell->size.height) << std::endl;
            }
        }
        file.close();

        file.open("rpack.txt");
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open output file rpack.txt" << std::endl;
            return result;
        }

        file << "RowStartY " << row_start_y << std::endl;
        file << "RowHeight " << this->site_height << std::endl;
        for (int i = 0; i < (int)local_region_rows.size(); i++) {
            file << "Row " << local_region_rows[i].row_idx << " " << local_region_rows[i].left_x << " " << local_region_rows[i].right_x << std::endl;
            for (auto& clp : lrpacks[i]) {
                file << clp->cell->name << " " << clp->right_x << " " << static_cast<int>(clp->cell->lower_left.y) << " " << static_cast<int>(clp->cell->size.width) << " " << static_cast<int>(clp->cell->size.height) << std::endl;
            }
        }
        file.close();
    }

    // Get insertion intervals
    int continous_intervals = 0;
    int max_continous_intervals = 0;
    std::vector<InsertionInterval*> insertion_intervals;
    for (int i = 0; i < (int)local_region_rows.size(); i++) {
        double x_r, x_l;
        bool continous_added = false;
        if (lrpacks[i].size() == 0) {
            x_r = local_region_rows[i].right_x - step.to->size.width;
            x_l = local_region_rows[i].left_x;
            if (x_r - x_l >= 0) {
                InsertionInterval* ii = new InsertionInterval;
                ii->row = i;
                ii->left_x = x_l;
                ii->right_x = x_r;
                ii->left_cell = nullptr;
                ii->right_cell = nullptr;
                insertion_intervals.push_back(ii);
                continous_intervals++;
            } else {
                continous_intervals = 0;
            }
            continue;
        }
        x_r = lrpacks[i][0]->right_x - step.to->size.width;
        x_l = local_region_rows[i].left_x;
        if (x_r - x_l >= 0) {
            InsertionInterval* ii = new InsertionInterval;
            ii->row = i;
            ii->left_x = x_l;
            ii->right_x = x_r;
            ii->left_cell = nullptr;
            ii->right_cell = lrpacks[i][0];
            insertion_intervals.push_back(ii);
            if (!continous_added) {
                continous_intervals++;
                continous_added = true;
            }
        }
        for (int j = 1; j < (int)lrpacks[i].size(); j++) {
            x_r = lrpacks[i][j]->right_x - step.to->size.width;
            x_l = lrpacks[i][j-1]->left_x + lrpacks[i][j-1]->cell->size.width;
            if (x_r - x_l >= 0) {
                InsertionInterval* ii = new InsertionInterval;
                ii->row = i;
                ii->left_x = x_l;
                ii->right_x = x_r;
                ii->left_cell = lrpacks[i][j-1];
                ii->right_cell = lrpacks[i][j];
                insertion_intervals.push_back(ii);
                if (!continous_added) {
                    continous_intervals++;
                    continous_added = true;
                }
            }
        }
        x_r = local_region_rows[i].right_x - step.to->size.width;
        x_l = lrpacks[i][lrpacks[i].size()-1]->left_x + lrpacks[i][lrpacks[i].size()-1]->cell->size.width;
        if (x_r - x_l >= 0) {
            InsertionInterval* ii = new InsertionInterval;
            ii->row = i;
            ii->left_x = x_l;
            ii->right_x = x_r;
            ii->left_cell = lrpacks[i][lrpacks[i].size()-1];
            ii->right_cell = nullptr;
            insertion_intervals.push_back(ii);
            if (!continous_added) {
                continous_intervals++;
                continous_added = true;
            }
        }
        if (!continous_added) {
            continous_intervals = 0;
        }
        if (continous_intervals > max_continous_intervals) {
            max_continous_intervals = continous_intervals;
        }
    }

    // Early exit
    if (max_continous_intervals < step.to->row_height) {
        for (int i = (int)local_region_rows.size()-1; i >= 0; i--) {
            for (auto& clp : lrpacks[i]) {
                if (clp->cell->row_idx == local_region_rows[i].row_idx) {
                    delete clp;
                }
            }
        }

        for (auto& ii : insertion_intervals) {
            delete ii;
        }
        return result;
    }

    // Sort the endpoints
    std::vector<Endpoint> endpoints;
    for (auto& ii : insertion_intervals) {
        Endpoint e;
        e.x = ii->left_x;
        e.is_left = true;
        e.ii = ii;
        endpoints.push_back(e);
        Endpoint er;
        er.x = ii->right_x;
        er.is_left = false;
        er.ii = ii;
        endpoints.push_back(er);
    }
    std::sort(endpoints.begin(), endpoints.end(), [](const Endpoint& a, const Endpoint& b) {
        return a.x < b.x;
    });

    // Create Queue
    std::map<std::pair<int, int>, std::deque<InsertionInterval*>> queues;
    for (int s = 0; s < local_region_height; ++s) {
        for (int r = 0; r < local_region_height; ++r) {
            if (s != r && std::abs(r - s) < step.to->row_height) {
                std::deque<InsertionInterval*> q;
                queues[std::make_pair(s, r)] = q;
            }
        }
    }

    // Insertion Point Enumeration
    // Psuedo code:
    // for each endpoint e in sorted endpoints:
    //     if e is a left endpoint:
    //         generate the insertion points involving e.ii
    //         if left cell is multi-row:
    //             let S be the set of rows that the left cell occupies
    //             for each segment s in S:
    //                clear the queue of (e.ii.r, s)
    //         push e.ii to all queues of (r, s) where s = e.ii.r
    //     else:
    //         pop e.ii from all queues of (r, s) where s = e.ii.r
    std::vector<InsertionPoint*> insertion_points;
    for (auto& e : endpoints) {
        if (e.is_left) {
            for (int i = e.ii->row - step.to->row_height + 1; i < e.ii->row + 1; i++) {
                if (i < 0 || i + step.to->row_height > local_region_height) continue;
                std::vector<std::deque<InsertionInterval*>*> involved_queues;
                std::deque<InsertionInterval*>* q = new std::deque<InsertionInterval*>;
                for (int j = 0; j < step.to->row_height; j++) {
                    if (i + j == e.ii->row) {
                        q->push_back(e.ii);
                        involved_queues.push_back(q);
                        continue;
                    }
                    involved_queues.push_back(&queues[std::make_pair(e.ii->row, i + j)]);
                }
                int index = 0;
                std::vector<InsertionInterval*> current_intervals;
                enumerate_insertion_points(e, i+local_region_start_row, involved_queues, index, current_intervals, DBL_MAX, insertion_points);
                delete q;
                if (insertion_points.size() > 100000) {
                    // std::cerr << "Error: Too many insertion points" << std::endl;
                    result.success = -1;
                    window_size_x -= 1;
                    window_size_y -= 1;
                    for (int i = (int)local_region_rows.size()-1; i >= 0; i--) {
                        for (auto& clp : lrpacks[i]) {
                            if (clp->cell->row_idx == local_region_rows[i].row_idx) {
                                delete clp;
                            }
                        }
                    }

                    for (auto& ii : insertion_intervals) {
                        delete ii;
                    }

                    for (auto& ip : insertion_points) {
                        delete ip;
                    }
                    return result;
                }
            }
            if (e.ii->left_cell != nullptr)
            if (e.ii->left_cell->cell->row_height > 1) {
                for (int s = e.ii->left_cell->cell->row_idx - local_region_start_row; s < e.ii->left_cell->cell->row_idx - local_region_start_row + e.ii->left_cell->cell->row_height; s++) {
                    if (s == e.ii->row) continue;
                    // queues[std::make_pair(e.ii->row, s)].clear();
                    while (true) {
                        if (queues[std::make_pair(e.ii->row, s)].empty()) break;
                        if (queues[std::make_pair(e.ii->row, s)].front()->left_x >= e.ii->left_x) break;
                        queues[std::make_pair(e.ii->row, s)].pop_front();
                    }
                }
            }
            for (int r = e.ii->row - step.to->row_height + 1; r < e.ii->row + step.to->row_height; r++) {
                if (r < 0 || r >= local_region_height) continue;
                if (r == e.ii->row) continue;
                queues[std::make_pair(r, e.ii->row)].push_back(e.ii);
            }
        } else {
            for (int r = e.ii->row - step.to->row_height + 1; r < e.ii->row + step.to->row_height; r++) {
                if (r < 0 || r >= local_region_height) continue;
                if (r == e.ii->row) continue;
                // queues[std::make_pair(r, e.ii->row)].remove(e.ii);
                if (queues[std::make_pair(r, e.ii->row)].empty()) {
                    continue;
                }
                if (queues[std::make_pair(r, e.ii->row)].front() != e.ii) {
                    continue;
                }
                queues[std::make_pair(r, e.ii->row)].pop_front();
            }
        }
    }

    if (insertion_points.size() == 0) {
        for (int i = (int)local_region_rows.size()-1; i >= 0; i--) {
            for (auto& clp : lrpacks[i]) {
                if (clp->cell->row_idx == local_region_rows[i].row_idx) {
                    delete clp;
                }
            }
        }

        for (auto& ii : insertion_intervals) {
            delete ii;
        }

        for (auto& ip : insertion_points) {
            delete ip;
        }
        return result;
    }

    double best_cost = DBL_MAX;
    InsertionPoint* best_ip = nullptr;
    for (auto& ip : insertion_points) {
        double cost = cost_calculator(ip, step);
        if (cost < best_cost) {
            best_cost = cost;
            best_ip = ip;
        }
    }

    result.lower_left = {(best_ip->best_x), best_ip->row * this->site_height + row_start_y};
    result.success = 1;
    result.affected_lower_left = {local_region_x.low, (double)local_region_rows[0].row_idx};
    result.affected_upper_right = {local_region_x.high, (double)local_region_rows[local_region_rows.size()-1].row_idx};
    result.best_ip = best_ip;
    for (auto& ii : best_ip->involved_intervals) {
        if (ii->left_cell != nullptr) {
            Cell* lcell = ii->left_cell->cell;
            if (std::find(result.lrecored.begin(), result.lrecored.end(), lcell) == result.lrecored.end()) {
                result.lrecored.push_back(lcell);
            }
        }
        if (ii->right_cell != nullptr) {
            Cell* rcell = ii->right_cell->cell;
            if (std::find(result.rrecored.begin(), result.rrecored.end(), rcell) == result.rrecored.end()) {
                result.rrecored.push_back(rcell);
            }
        }
    }

    if (dump_local_region) {
        dump_loaded_row("placing_row.txt");
    }

    for (int i = (int)local_region_rows.size()-1; i >= 0; i--) {
        for (auto& clp : lrpacks[i]) {
            if (clp->cell->row_idx == local_region_rows[i].row_idx) {
                delete clp;
            }
        }
    }

    for (auto& ii : insertion_intervals) {
        delete ii;
    }

    for (auto& ip : insertion_points) {
        if (ip != best_ip)
            delete ip;
    }

    return result;

    // // Placement Realization
    // std::vector<Cell*> lrecored;
    // std::vector<Cell*> rrecored;
    // for (auto& ii : best_ip->involved_intervals) {
    //     if (ii->left_cell != nullptr) {
    //         Cell* lcell = ii->left_cell->cell;
    //         if (std::find(lrecored.begin(), lrecored.end(), lcell) == lrecored.end()) {
    //             lrecored.push_back(lcell);
    //         }
    //     }
    //     if (ii->right_cell != nullptr) {
    //         Cell* rcell = ii->right_cell->cell;
    //         if (std::find(rrecored.begin(), rrecored.end(), rcell) == rrecored.end()) {
    //             rrecored.push_back(rcell);
    //         }
    //     }
    // }

    // // Left
    // for (auto& lcell : lrecored) {
    //     cell_pusher(lcell, best_ip->best_x - lcell->size.width, result.moved_cells, true);
    // }

    // // Right
    // for (auto& rcell : rrecored) {
    //     cell_pusher(rcell, best_ip->best_x + step.to->size.width, result.moved_cells, false);
    // }

    // step.to->lower_left.x = best_ip->best_x;
    // step.to->lower_left.y = best_ip->row * this->site_height + row_start_y;
    // step.to->row_idx = best_ip->row;
    // step.to->fixed = false;
    // for (int r = step.to->row_idx; r < step.to->row_idx + step.to->row_height; r++) {
    //     if (r < 0 || r >= (int)prows.size()) continue;
    //     PlacementRow* prow = prows[r];
    //     prow->itree.insert({(best_ip->best_x), best_ip->best_x + step.to->size.width, step.to});
    // }

    // result.lower_left = {(best_ip->best_x), best_ip->row * this->site_height + row_start_y};
    // result.success = 1;

    // if (dump_local_region) {
    //     dump_loaded_row("placing_row.txt");
    // }

    // for (int i = (int)local_region_rows.size()-1; i >= 0; i--) {
    //     for (auto& clp : lrpacks[i]) {
    //         if (clp->cell->row_idx == local_region_rows[i].row_idx) {
    //             delete clp;
    //         }
    //     }
    // }

    // for (auto& ii : insertion_intervals) {
    //     delete ii;
    // }

    // for (auto& ip : insertion_points) {
    //     delete ip;
    // }


    // return result;
}