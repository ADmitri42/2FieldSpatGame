#include "games.h"
#include <cmath>
#include <iostream>
#include <numeric>


void NovakMayGame::calculate_scores(std::vector<double> &scores){
    scores.assign(L*L, 0);
    double density = densities.back();

    //Payoffs
    for (size_t k = 0; k < L*L; k++) {
        int y = k / L; // Row
        int x = k % L; // Col

        for (int i = -1; i <= 1; i++) //Row
        {
            for (int j = -1; j <= 1; j++) //Col
            {
                size_t memberIndex = (x + i + L) % L + L * ((y + j + L) % L);
                scores[k] += field[memberIndex];
            }
        }

        if (field[k] == 0)
        {
            scores[k] = scores[k] * b;
        }
    }
}

void MeanGame::calculate_scores(std::vector<double> &scores){
    scores.assign(L*L, 0);
    double density = densities.back();

    //Payoffs
    for (size_t k = 0; k < L*L; k++) {
        int y = k / L; // Row
        int x = k % L; // Col

        for (int i = -1; i <= 1; i++) //Row
        {
            for (int j = -1; j <= 1; j++) //Col
            {
                size_t memberIndex = (x + i + L) % L + L * ((y + j + L) % L);
                if((i == 0)&&(j == 0)){
                    scores[k] += density;
                } else {
                    scores[k] += field[memberIndex];// == 0 ? 1 : 0;
                }
            }
        }

        if (field[k] == 0)
        {
            scores[k] = scores[k] * b;
        }
    }
}

/*
 * Triangular lattice
 */

void NovakMayTriangularGame::calculate_scores(std::vector<double> &scores){
    scores.assign(L*L, 0);
//    double density = densities.back();

    //Payoffs
    for (size_t k = 0; k < L*L; k++) {
        int y = k / L; // Row
        int x = k % L; // Col

        for (int i = -1; i <= 1; i++) //Row
        {
            for (int j = -1; j <= 1; j++) //Col
            {
                if(i*j == 1){
                    continue;
                }
                size_t memberIndex = (x + i + L) % L + L * ((y + j + L) % L);
                scores[k] += field[memberIndex];
            }
        }

        if (field[k] == 0)
        {
            scores[k] = scores[k] * b;
        }
    }
}

void NovakMayTriangularGame::update_field(const std::vector<double> &scores, int time_moment, int perCalFrom, int perCalTill){
    //Field
    std::vector<char> currentField(field);
//    std::copy(field.begin(), field.end(), currentField.begin());

    //Strategy
    for (size_t k = 0; k < L*L; k++) {
        int y = k / L; // Row
        int x = k % L; // Col

        size_t bestStrategyIndex = k;

        for (int i = -1; i <= 1; i++) //Row
        {
            for (int j = -1; j <= 1; j++) //Col
            {
                if(i*j == 1){
                    continue;
                }
                size_t memberIndex = (x + i + L) % L + L * ((y + j + L) % L);

                if (scores[bestStrategyIndex] < scores[memberIndex])
                {
                    bestStrategyIndex = memberIndex;
                }
            }
        }
        if((perCalFrom >= 0)&&(time_moment > perCalFrom)&&(time_moment < perCalTill)&&(k != bestStrategyIndex)){
            unchanged[k] = 0;
        }
        field[k] = currentField[bestStrategyIndex];
    }
    currentField.clear();
}

void MeanTriangularGame::calculate_scores(std::vector<double> &scores){
    scores.assign(L*L, 0);
    double density = densities.back();

    //Payoffs
    for (size_t k = 0; k < L*L; k++) {
        int y = k / L; // Row
        int x = k % L; // Col

        for (int i = -1; i <= 1; i++) //Row
        {
            for (int j = -1; j <= 1; j++) //Col
            {
                if(i*j == 1){
                    continue;
                }
                size_t memberIndex = (x + i + L) % L + L * ((y + j + L) % L);
                if((i == 0)&&(j == 0)){
                    scores[k] += density;
                } else {
                    scores[k] += field[memberIndex];// == 0 ? 1 : 0;
                }
            }
        }

        if (field[k] == 0)
        {
            scores[k] = scores[k] * b;
        }
    }
}



void MeanTriangularGame::update_field(const std::vector<double> &scores, int time_moment, int perCalFrom, int perCalTill){
    //Field
    std::vector<char> currentField(field);
//    std::copy(field.begin(), field.end(), currentField.begin());

    //Strategy
    for (size_t k = 0; k < L*L; k++) {
        int y = k / L; // Row
        int x = k % L; // Col

        size_t bestStrategyIndex = k;

        for (int i = -1; i <= 1; i++) //Row
        {
            for (int j = -1; j <= 1; j++) //Col
            {
                if(i*j == 1){
                    continue;
                }
                size_t memberIndex = (x + i + L) % L + L * ((y + j + L) % L);

                if (scores[bestStrategyIndex] < scores[memberIndex])
                {
                    bestStrategyIndex = memberIndex;
                }
            }
        }
        if((perCalFrom >= 0)&&(time_moment > perCalFrom)&&(time_moment < perCalTill)&&(k != bestStrategyIndex)){
            unchanged[k] = 0;
        }
        field[k] = currentField[bestStrategyIndex];
    }
    currentField.clear();
}
