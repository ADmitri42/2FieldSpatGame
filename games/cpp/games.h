#ifndef __EVOLVE_FIELD__
#define __EVOLVE_FIELD__

#include "spatgame.h"
#include <sys/types.h>
#include <vector>

/*
 * Double field basic game
 *
 * Two field stored one after th other
 * Density for field 1 has indeces 2n
 * For field 2 indeces 2n+1
 */

class DoubleMeanFieldGame : public AbstractSpatialGame {
protected:
  double b1, b2, lam, mu;

public:
  DoubleMeanFieldGame(size_t size, double _b1 = 1.6, double _b2 = 1.6, double _lam=0, double _mu=1);

  virtual void calculate_scores(std::vector<double> &scores);
  virtual void update_field(const std::vector<double> &scores, int time_moment,
                            int percfrom = -1, int perctill = -1);
  void evolve(int num_steps = 0, int percfrom = -1, int perctill = -1);

  std::vector<double> get_bs();
  std::vector<double> get_koef();
  void set_b(double new_b1, double new_b2);
  void set_field(const std::vector<int> &new_field1,
                 const std::vector<int> &new_field2);
  std::vector<double> get_persistences();
};

#endif
