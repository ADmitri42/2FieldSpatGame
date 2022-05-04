//
// Created by Dima on 12.02.2020.
//

#ifndef SPATIAL_EVOLUTIONARY_GAME_SPATGAME_H
#define SPATIAL_EVOLUTIONARY_GAME_SPATGAME_H

#include <sys/types.h>
#include <vector>

class AbstractSpatialGame {
protected:
  std::vector<char> field;
  std::vector<char> unchanged;

  std::vector<double> densities;

  size_t L;
  double b1, b2, lam, mu;
  int perCalFrom;
  int perCalTill;

public:
  AbstractSpatialGame(size_t size, double _b1 = 1.6, double _b2 = 1.6, double _lam=0, double _mu=1);
  virtual ~AbstractSpatialGame(){};

  virtual void calculate_scores(std::vector<double> &scores);
  virtual void update_field(const std::vector<double> &scores, int time_moment,
                            int percfrom = -1, int perctill = -1);
  void evolve(int num_steps = 0, int percfrom = -1, int perctill = -1);

  std::vector<double> get_densities();
  size_t size();
  virtual std::vector<double> get_bs();
  virtual std::vector<double> get_koef();
  void set_b(double new_b1, double new_b2);
  void set_koef(double new_lam, double new_mu);
  std::vector<int> get_field();
  void set_field(const std::vector<int> &new_field1,
                 const std::vector<int> &new_field2);
  std::vector<double> get_persistences();

  std::vector<int> mn_distribution();
  /*
   * Function to create numpy array
   */
  char *get_field_pointer();
  int get_densities_size();
  double *get_densities_pointer();
};

#endif // SPATIAL_EVOLUTIONARY_GAME_SPATGAME_H
