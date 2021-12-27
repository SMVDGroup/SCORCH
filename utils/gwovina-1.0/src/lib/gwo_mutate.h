
#ifndef VINA_DE_MUTATE_H
#define VINA_DE_MUTATE_H

#include "model.h"
#include "quasi_newton.h"
#include <vector>
#include <cmath>
#include <queue>

// does not set model
				
void gwo_mutate_conf(std::vector<output_type>& populations, int position,output_type alpha, output_type beta, output_type delta, float a, rng& generator);

//void gwo_mutate_conf_team(std::vector<output_type>& populations, int position, std::priority_queue<std::pair<float, int> >& best_team, float a, rng& generator);
void gwo_mutate_conf_team2(std::vector<output_type>& populations, int position, std::vector<output_type>& leaders, float a, rng& generator);

#endif
