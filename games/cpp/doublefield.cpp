#include "games.h"
#include <cmath>
#include <iostream>
#include <numeric>

/*
 * initialize game field with vector for densities
 *
 * size - length of field
 * _b - initial payoff parameter
 */
DoubleMeanFieldGame::DoubleMeanFieldGame(size_t size, double _b1, double _b2, double _lam, double _mu)
    : AbstractSpatialGame(size, (_b1 + _b2) / 2) {
  L = size;
  field.assign(2 * L * L, 0);
  unchanged.assign(2 * L * L, 1);
  densities.push_back(0);
  b1 = _b1;
  b2 = _b2;
  lam = _lam;
  mu = _mu;
}

/*
 * return vector of densities for all steps after resetting field
 */
std::vector<double> DoubleMeanFieldGame::get_bs() { return {b1, b2}; }

std::vector<double> DoubleMeanFieldGame::get_koef() { return {lam, mu}; }

/*
 * Set new payoff parameter
 * (doesn't reset the statistic)
 */
void DoubleMeanFieldGame::set_b(double new_b1, double new_b2) {
  b1 = new_b1;
  b2 = new_b2;
}

/*
 * Set new fields
 * (statistic would be reset)
 */
void DoubleMeanFieldGame::set_field(const std::vector<int> &new_field1,
                                    const std::vector<int> &new_field2) {
  if ((new_field1.size() != L * L) || (new_field2.size() != L * L)) {
    throw std::length_error("Wrong size");
  }
  for (size_t i = 0; i < L * L; i++) {
    field[i] = new_field1[i];
  }
  for (size_t i = 0; i < L * L; i++) {
    field[L * L + i] = new_field2[i];
  }
  //    field.assign(new_field.begin(), new_field.end());
  densities.clear();
  unchanged.assign(2 * L * L, 1);
  densities.push_back(
      static_cast<double>(accumulate(new_field1.begin(), new_field1.end(), 0)) /
      new_field1.size());
  densities.push_back(
      static_cast<double>(accumulate(new_field2.begin(), new_field2.end(), 0)) /
      new_field2.size());
}

/*
 * Return persistence
 */
std::vector<double> DoubleMeanFieldGame::get_persistences() {
  return {static_cast<double>(std::accumulate(unchanged.begin(),
                                              unchanged.begin() + L * L, 0)) /
              (L * L),
          static_cast<double>(
              std::accumulate(unchanged.begin() + L * L, unchanged.end(), 0)) /
              (L * L)};
}

void DoubleMeanFieldGame::calculate_scores(std::vector<double> &scores) {
  scores.assign(2 * L * L, 0);
  double density, densitycur;
  // Payoffs
  int offset = 0;
  for (int off = 0; off < 2; ++off) {
    offset = L * L * off;
    density = densities[densities.size() - 1 - off];
    densitycur = densities[densities.size() - 2 + off];
    for (size_t k = 0; k < L * L; k++) {
      int y = k / L; // Row
      int x = k % L; // Col

      for (int i = -1; i <= 1; i++) // Row
      {
        for (int j = -1; j <= 1; j++) // Col
        {
          size_t memberIndex = (x + i + L) % L + L * ((y + j + L) % L);
          if ((i == 0) && (j == 0)) {
            scores[offset + k] += lam * densitycur + mu * density;
          } else {
            scores[offset + k] += field[offset + memberIndex]; // == 0 ? 1 : 0;
          }
        }
      }

      if (field[offset + k] == 0) {
        scores[offset + k] = scores[offset + k] * (b1 * (1 - off) + b2 * off);
      }
    }
  }
}

void DoubleMeanFieldGame::update_field(const std::vector<double> &scores,
                                       int time_moment, int perCalFrom,
                                       int perCalTill) {
  // Field
  std::vector<char> currentField(field);

  // Strategy
  int offset = 0;
  for (int off = 0; off < 2; ++off) {
    offset = L * L * off;
    for (size_t k = 0; k < L * L; k++) {
      int y = k / L; // Row
      int x = k % L; // Col

      size_t bestStrategyIndex = k;

      for (int i = -1; i <= 1; i++) // Row
      {
        for (int j = -1; j <= 1; j++) // Col
        {
          size_t memberIndex = (x + i + L) % L + L * ((y + j + L) % L);

          if (scores[offset + bestStrategyIndex] <
              scores[offset + memberIndex]) {
            bestStrategyIndex = memberIndex;
          }
        }
      }
      if ((perCalFrom >= 0) && (time_moment > perCalFrom) &&
          (time_moment < perCalTill) && (k != bestStrategyIndex)) {
        unchanged[offset + k] = 0;
      }
      field[offset + k] = currentField[offset + bestStrategyIndex];
    }
  }
  currentField.clear();
}

/*
 * Evolve num_steps
 */
void DoubleMeanFieldGame::evolve(int num_steps, int perCalFrom,
                                 int perCalTill) {
  std::vector<double> scores(2 * L * L, 0);

  for (int step = 0, time_moment = densities.size() / 2; step < num_steps;
       step++, time_moment++) {

    // Payoffs
    scores.assign(2 * L * L, 0);
    calculate_scores(scores);

    // Strategy
    update_field(scores, time_moment, perCalFrom, perCalTill);
    densities.push_back(
        (1. * std::accumulate(field.begin(), field.begin() + L * L, 0)) /
        (L * L));
    densities.push_back(
        (1. * std::accumulate(field.begin() + L * L, field.end(), 0)) /
        (L * L));
  }

  scores.clear();
}
