/*
    Class for spatial game on triangular mesh
*/

#ifndef __TRIANGUAL_FIELD__
#define __TRIANGUAL_FIELD__

#include "spatgame.h"
#include <sys/types.h>
#include <vector>

class TriangularFieldGame : public AbstractSpatialGame {
public:
    TriangularFieldGame(size_t size, double _b1 = 1.6, double _b2 = 1.6, double _lam=0, double _mu=1, double _k=0, int default_seed=42) : AbstractSpatialGame(size, _b1, _b2, _lam, _mu, _k, default_seed){};
    void calculate_scores(std::vector<double> &scores);
    void update_field(
        const std::vector<double> &scores,
        int time_moment,
        int percfrom = -1,
        int perctill = -1
    );
};

#endif