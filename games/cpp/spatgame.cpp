#include "spatgame.h"
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>


/*
 * initialize game field with vector for densities
 *
 * size - length of field
 * _b - initial payoff parameter
 */
AbstractSpatialGame::AbstractSpatialGame(size_t size, double _b1, double _b2, double _lam, double _mu, double _k, int default_seed){
    L = size;
    field.assign(2 * L * L, 0);
    unchanged.assign(2 * L * L, 1);
    densities.push_back(0);
    b1 = _b1;
    b2 = _b2;
    lam = _lam;
    mu = _mu;
    K = _k;
    std::mt19937 rnd_engine (default_seed);
	  std::uniform_real_distribution<double> distribution(0.0, 1.0);
}

/*
 * Generate uniformaly value from [0;1]
 */
double AbstractSpatialGame::random() {
    return distribution(rnd_engine);
}

/*
 * return vector of densities for all steps after resetting field
 */
std::vector<double> AbstractSpatialGame::get_densities(){
    return densities;
}

size_t AbstractSpatialGame::size(){
    return L;
}

/*
 * return vector of densities for all steps after resetting field
 */
std::vector<double> AbstractSpatialGame::get_bs() { return {b1, b2}; }

std::vector<double> AbstractSpatialGame::get_koef() { return {lam, mu}; }

double AbstractSpatialGame::get_K() { return K; }

/*
 * Set new payoff parameter
 * (doesn't reset the statistic)
 */
void AbstractSpatialGame::set_b(double new_b1, double new_b2) {
  b1 = new_b1;
  b2 = new_b2;
}

void AbstractSpatialGame::set_koef(double new_lam, double new_mu){
  lam = new_lam;
  mu = new_mu;
}

/*
 * Set measure of stochastic uncertainties (noise) K
 */
void AbstractSpatialGame::set_K(double new_k) {
  K = new_k;
}

std::vector<int> AbstractSpatialGame::get_field(){
    std::vector<int> ifield;
    for(auto t: field){
        ifield.push_back(static_cast<int>(t));
    }
//    ifield.assign(field.begin(), field.end());
    return ifield;
}

/*
 * Set new field
 * (statistic would be reset)
 */
void AbstractSpatialGame::set_field(const std::vector<int> &new_field1,
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
std::vector<double> AbstractSpatialGame::get_persistences() {
  return {static_cast<double>(std::accumulate(unchanged.begin(),
                                              unchanged.begin() + L * L, 0)) /
              (L * L),
          static_cast<double>(
              std::accumulate(unchanged.begin() + L * L, unchanged.end(), 0)) /
              (L * L)};
}


void AbstractSpatialGame::calculate_scores(std::vector<double> &scores) {
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

void AbstractSpatialGame::update_field(const std::vector<double> &scores,
                                       int time_moment, int perCalFrom,
                                       int perCalTill) {
  // Field
  std::vector<char> currentField(field);

  std::vector<double> rand_w(L * L, 0);

  // Strategy
  int offset = 0;
  for (int off = 0; off < 2; ++off) {
    offset = L * L * off;
    if(K != 0) {
      for (size_t i = 0; i < L * L; i++)
      {
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
        if (rand_w[current] <= W) {
          changed = true;
        } else {
          changed = false;
        }
      } else {
        if (score_diff < 0) {
          changed = true;
        } else {
          changed = false;
        }
      }
      if (changed) {
        field[offset + current] = opposite_strategy;
        if ((perCalFrom >= 0) && (time_moment > perCalFrom) &&
          (time_moment < perCalTill)) {
            unchanged[offset + current] = 0;
          }
      }
    }
  }
  currentField.clear();
}

/*
 * Evolve num_steps
 */
void AbstractSpatialGame::evolve(int num_steps, int perCalFrom,
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

std::vector<int> AbstractSpatialGame::mn_distribution(){
    std::vector<double> score(L*L, 0);
    calculate_scores(score);
    std::vector<int> nmdistr(2*9*9, 0);
    int offset, moffset, m, n, is, x1, x2, x3, y1, y2, y3;
    double n_sc, m_sc;

    #pragma omp parallel for
      for (int off = 0; off < 2; ++off) {
      offset = L * L * off;
      moffset = 81 * off;
      for(size_t x = 0; x < L; ++x){
          for (size_t y = 0; y < L; ++y) {
              n = m = -1;
              n_sc = m_sc = -22;
              for (int i = -1; i < 2; ++i) {
                  for (int j = -1; j < 2; ++j) {
                      is = ((L+y+j)%L)*L+(L + x+i)%L;

                      // Neigbours coordinates
                      x1 = (L+x+i-1)%L;
                      x2 = (L+x+i)%L;
                      x3 = (L+x+i+1)%L;

                      y1 = (L+y+j-1)%L;
                      y2 = (L+y+j)%L;
                      y3 = (L+y+j+1)%L;

                      if((field[offset + is] == 0)&&(score[offset + is] > n_sc)){
                          n = field[offset + y1*L + x1] + field[offset + y1*L + x2] + field[offset + y1*L + x3]
                                  + field[offset + y2*L + x1] +  field[offset + y2*L + x3]
                                  + field[offset + y3*L + x1] + field[offset + y3*L + x2] + field[offset + y3*L + x3];

                          n_sc = score[offset + is];
                      }

                      if((field[offset + is] == 1)&&(score[offset + is] > m_sc)){
                          m = field[offset + y1*L + x1] + field[offset + y1*L + x2] + field[offset + y1*L + x3]
                              + field[offset + y2*L + x1] +  field[offset + y2*L + x3]
                              + field[offset + y3*L + x1] + field[offset + y3*L + x2] + field[offset + y3*L + x3];

                          m_sc = score[offset + is];
                      }
                  }
              }

              if((m >= 0)&&(n >= 0)){
                  nmdistr[moffset + m*9 + n] += 1;
              }
  //            std::cout << n << " " << m << std::endl;
          }
      }
    }
    return nmdistr;
}

/*
 * Methods that simplify NumPy Array creation
 */
char* AbstractSpatialGame::get_field_pointer(){
    return &field[0];
}

int AbstractSpatialGame::get_densities_size(){
    return densities.size();
}

double* AbstractSpatialGame::get_densities_pointer(){
    return &densities[0];
}
