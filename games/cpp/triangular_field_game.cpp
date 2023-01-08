#include "triangular_field_game.h"
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <fstream>

void TriangularFieldGame::calculate_scores(std::vector<double> &scores)
{
  scores.assign(2 * L * L, 0);
  double density, densitycur;
  // Payoffs
  int offset = 0;
  for (int off = 0; off < 2; ++off) {
    offset = L * L * off;
    density = densities[densities.size() - 1 - off];
    densitycur = densities[densities.size() - 2 + off];

    #pragma omp parallel for
    for (size_t current = 0; current < L * L; current++) {
      int y = current / L; // Row
      int x = current % L; // Col

      for (int i = -1; i <= 1; i++) // Row
      {
        for (int j = -1; j <= 1; j++) // Col
        {
          if(i*j == 1){
            continue;
          }
          size_t memberIndex = (x + i + L) % L + L * ((y + j + L) % L);
          if ((i == 0) && (j == 0)) {
            scores[offset + current] += lam * densitycur + mu * density;
          } else {
            scores[offset + current] += field[offset + memberIndex]; // == 0 ? 1 : 0;
          }
        }
      }

      if (field[offset + current] == 0) {
          scores[offset + current] = scores[offset + current] * (b1 * (1 - off) + b2 * off);
      }
    }
  }
}

void TriangularFieldGame::update_field(
  const std::vector<double> &scores,
  int time_moment,
  int perCalFrom,
  int perCalTill)
{
  // Field
  std::vector<char> currentField(field);
  std::vector<double> rand_w(L * L, 0);

  // Strategy
  int offset = 0;
  for (int off = 0; off < 2; ++off) {
    offset = L * L * off;
    if(K != 0) {
      for (size_t i = 0; i < L * L; i++) {
        rand_w[i] = random();
      }
    }

    #pragma omp parallel for
    for (size_t current = 0; current < L * L; current++) {
      // Highest score of defector and cooperator
      double highest_score[] = {-1., -1.};
      int current_strategy, opposite_strategy;
      double score_diff, W;
      bool changed = false;
      
      int y = current / L; // Row
      int x = current % L; // Col

      for (int i = -1; i <= 1; i++) // Row
      {
        for (int j = -1; j <= 1; j++) // Col
        {
          if(i*j == 1){
            continue;
          }
          size_t memberIndex = (x + i + L) % L + L * ((y + j + L) % L);
          if (highest_score[(int)currentField[offset + memberIndex]] < scores[offset + memberIndex]) {
            highest_score[(int)currentField[offset + memberIndex]] = scores[offset + memberIndex];
          }
        }
      }
      current_strategy = currentField[offset + current];
      opposite_strategy = ((int)currentField[offset + current] == 1) ? 0 : 1;

      score_diff = highest_score[current_strategy] - highest_score[opposite_strategy];
      if (K != 0) {
        W = 1 / (1 + std::exp(score_diff / K));
        changed = rand_w[current] <= W;
      } else {
        changed = score_diff < 0;
      }
      if (changed) {
        field[offset + current] = opposite_strategy;
        if ((perCalFrom >= 0) && (time_moment > perCalFrom) &&
          (time_moment < perCalTill)) {
            unchanged[offset + current] = 0;
          }
      }
    }
    currentField.clear();
  }
}