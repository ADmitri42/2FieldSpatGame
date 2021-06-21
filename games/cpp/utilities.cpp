#include <map>
#include <assert.h>
#include "utilities.h"
#include <iostream>
#include <algorithm>
#include <unordered_set>

/******************************
 ******************************
 ****                      ****
 ****   N-M distribution   ****
 ****                      ****
 ******************************
 ******************************/

std::vector<int> n_m_distribution(MeanGame &game){
    size_t L = game.size();
    std::vector<int> field = game.get_field();
    std::vector<double> score(L*L, 0);
    game.calculate_scores(score);
    std::vector<int> nmdistr(9*9, 0);
    int m, n, is, x1, x2, x3, y1, y2, y3;
    double n_sc, m_sc;

    for(size_t x = 0; x < L; ++x){
        for (size_t y = 0; y < L; ++y) {
            n = m = -1;
            n_sc = m_sc = -22;
            for (int i = -1; i < 2; ++i) {
                for (int j = -1; j < 2; ++j) {
                    is = ((L+y+j)%L)*L+(L + x+i)%L;
                    x1 = (L+x+i-1)%L;
                    x2 = (L+x+i)%L;
                    x3 = (L+x+i+1)%L;

                    y1 = (L+y+j-1)%L;
                    y2 = (L+y+j)%L;
                    y3 = (L+y+j+1)%L;

                    if((field[is] == 0)&&(score[is] > m_sc)){
                        m = field[y1*L + x1] + field[y1*L + x2] + field[y1*L + x3]
                                + field[y2*L + x1] +  field[y2*L + x3]
                                + field[y3*L + x1] + field[y3*L + x2] + field[y3*L + x3];

                        m_sc = score[is];
                    }

                    if((field[is] == 1)&&(score[is] > n_sc)){
                        n = field[y1*L + x1] + field[y1*L + x2] + field[y1*L + x3]
                            + field[y2*L + x1] +  field[y2*L + x3]
                            + field[y3*L + x1] + field[y3*L + x2] + field[y3*L + x3];

                        n_sc = score[is];
                    }
                }
            }

            if((m >= 0)&&(n >= 0)){
                nmdistr[m*9 + n] += 1;
            }
//            std::cout << n << " " << m << std::endl;
        }
    }
    return nmdistr;
}


std::vector<int> py_n_m_distribution(MeanGame* game){
    MeanGame& g = *game;
    return n_m_distribution(g);
}


/************************
 ************************
 ****                ****
 ****   CLUSTERING   ****
 ****                ****
 ************************
 ************************/

/*
 * Search for proper cluster label, since neighbor_cl can be not optimal(optimal == smallest os possible)
 */
int classify(int neighbor_cl, std::vector<int> &N){
    int t, r, c = 0;
    r = neighbor_cl;
    t = r;
    t = -N[t];
    while(t >= 0){
        r = t;
        t = -N[t];
        c++;
    }
    if(c > 2){
        N[neighbor_cl] = -r;
    }
    return r;
}

/*
 * Return proper label base on labels of top and left neighbours(0 if not proper class)
 */
int assign_label(int left, int top, int right, int bottom, std::vector<int> &N){
    int label;
    switch (!!left + !!top + !!right + !!bottom){
        case 0:
            N.push_back(1);
            label = N.size() - 1;
            break;
        case 1:
            label = std::max(std::max(top, bottom), std::max(left, right));
            label = classify(label, N);
            N[label] += 1;
            break;
        default:
            std::unordered_set<int> labels = {left, top, right, bottom};
            std::unordered_set<int> right_labels;
            for(auto t: labels){
                if(t > 0){
                    right_labels.insert(classify(t, N));
                }
            }
            label = *std::min_element(right_labels.begin(), right_labels.end());
            int size = 0;
            for(auto t: right_labels){
                size += N[t];
                N[t] = -label;
            }
            N[label] = size + 1;
            break;
    }
    return label;
}


LabeledField* fix_labels(LabeledField* lbf){
    std::vector<int> w_cluster_sizes; // For every label on the field it contains size of the cluster
    std::map<int, int> fixed_label;

    w_cluster_sizes.assign(lbf->cluster_sizes.begin(), lbf->cluster_sizes.end());
    lbf->cluster_sizes.clear();
    lbf->cluster_sizes.push_back(0);

    for(size_t i = 0, counter = 1; i < w_cluster_sizes.size(); ++i){
        if(w_cluster_sizes[i] > 0){
            fixed_label[i] = lbf->cluster_sizes.size();
            lbf->cluster_sizes.push_back(w_cluster_sizes[i]);
            counter++;
        }
    }

    for(size_t i = 0; i < lbf->labeled_field.size(); ++i){
        if(lbf->labeled_field[i] > 0){
            lbf->labeled_field[i] = fixed_label[classify(lbf->labeled_field[i],
                                                         w_cluster_sizes)];
        }
    }
    return lbf;
}


/*
 * Return field with labels
 * N - horizontal size of the field
 * M - vertical size of the field
 * field - vector of size N*M
 */
LabeledField **clustering(const std::vector<int> &field, int N, int M) {
    assert((N * M == field.size()));

    // Variables
    int ileft, itop, iright, ibottom;
    std::vector<int>   w_cluster_sizes(1, 0); // For every label on the field it contains size of the cluster
    std::map<int, int> fixed_label;

    LabeledField **lbf = new LabeledField *[2];
    lbf[0] = new LabeledField(field.size());
    lbf[1] = new LabeledField(field.size());

    for (size_t i = 0; i < field.size(); ++i) {
        ileft   = (N + i - 1) % N + (i / N) * N;
        itop    = (N * M + i - N) % (N * M);
        iright  = (N + i + 1) % N + (i / N) * N;
        ibottom = (N * M + i + N) % (N * M);
        lbf[field[i]]->labeled_field[i] = assign_label(lbf[field[i]]->labeled_field[ileft],
                                                       lbf[field[i]]->labeled_field[itop],
                                                       lbf[field[i]]->labeled_field[iright],
                                                       lbf[field[i]]->labeled_field[ibottom],
                                                       lbf[field[i]]->cluster_sizes);
    }



    // Create map object where current label placed against
    lbf[0] = fix_labels(lbf[0]);
    lbf[1] = fix_labels(lbf[1]);
    return lbf;
}