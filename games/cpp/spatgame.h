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
  double b;
  int perCalFrom;
  int perCalTill;

public:
  AbstractSpatialGame(size_t size, double _b = 1.8);
  virtual ~AbstractSpatialGame(){};

  virtual void calculate_scores(std::vector<double> &scores);
  virtual void update_field(const std::vector<double> &scores, int time_moment,
                            int percfrom = -1, int perctill = -1);
  void evolve(int num_steps = 0, int percfrom = -1, int perctill = -1);

  std::vector<double> get_densities();
  size_t size();
  virtual double get_b();
  void set_b(double new_b);
  std::vector<int> get_field();
  void set_field(const std::vector<int> &new_field);
  double get_persistence();
  /*
   * Function to create numpy array
   */
  char *get_field_pointer();
  int get_densities_size();
  double *get_densities_pointer();
};

#endif // SPATIAL_EVOLUTIONARY_GAME_SPATGAME_H
