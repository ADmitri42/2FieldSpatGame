#include "spatgame.h"
#include <cmath>
#include <iostream>
#include <numeric>


/*
 * initialize game field with vector for densities
 *
 * size - length of field
 * _b - initial payoff parameter
 */
AbstractSpatialGame::AbstractSpatialGame(size_t size, double _b){
    L = size;
    field.assign(L*L, 0);
    unchanged.assign(L*L, 1);
    densities.push_back(0);
    b = _b;
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

double AbstractSpatialGame::get_b(){
    return b;
}

/*
 * Set new payoff parameter
 * (doesn't reset the statistic)
 */
void AbstractSpatialGame::set_b(double new_b){
    b = new_b;
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
void AbstractSpatialGame::set_field(const std::vector<int> &new_field){
    if(new_field.size() != L*L){
        throw std::length_error("Wrong size");
    }
    for(size_t i = 0; i < field.size(); i++){
        field[i] = new_field[i];
    }
//    field.assign(new_field.begin(), new_field.end());
    densities.clear();
    unchanged.assign(L*L, 1);
    densities.push_back(static_cast<double>(accumulate(field.begin(),field.end(),0))/field.size());
}

/*
 * Return persistence
 */
double AbstractSpatialGame::get_persistence(){
    return static_cast<double>(std::accumulate(unchanged.begin(),unchanged.end(),0))/unchanged.size();
}


void AbstractSpatialGame::calculate_scores(std::vector<double> &scores){
    scores.assign(L*L, 0);

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
                    scores[k] += 0;
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


void AbstractSpatialGame::update_field(const std::vector<double> &scores, int time_moment, int perCalFrom, int perCalTill){
    //Field
    std::vector<char> currentField(field);

    //Strategy
    for (size_t k = 0; k < L*L; k++) {
        int y = k / L; // Row
        int x = k % L; // Col

        size_t bestStrategyIndex = k;

        for (int i = -1; i <= 1; i++) //Row
        {
            for (int j = -1; j <= 1; j++) //Col
            {
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

/*
 * Evolve num_steps
 */
void AbstractSpatialGame::evolve(int num_steps, int perCalFrom, int perCalTill)
{
    std::vector<double> scores(L*L, 0);

    for(int step = 0,time_moment = densities.size(); step < num_steps; step++, time_moment++)
    {

        //Payoffs
        calculate_scores(scores);

        //Strategy
        update_field(scores, time_moment, perCalFrom, perCalTill);
        densities.push_back((1.*std::accumulate(field.begin(),field.end(),0))/field.size());
    }

    scores.clear();
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