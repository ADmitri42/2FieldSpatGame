//
// Created by Dima on 26.01.2020.
//

#ifndef DYNAMIC_FRACTALS_CLUSTERING_H
#define DYNAMIC_FRACTALS_CLUSTERING_H

#include <vector>
#include <sys/types.h>
#include "games.h"

typedef struct LBF {
    std::vector<int> labeled_field;
    std::vector<int> cluster_sizes;
    LBF(int field_size = 1): labeled_field(field_size, 0), cluster_sizes(1, 0) { }
} LabeledField;


std::vector<int> n_m_distribution(MeanGame &game);
LabeledField** clustering(const std::vector<int>& field, int N, int M);

std::vector<int> py_n_m_distribution(MeanGame *game);
#endif //DYNAMIC_FRACTALS_CLUSTERING_H
