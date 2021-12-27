
#ifndef VINA_DE_MUTATE_H
#define VINA_DE_MUTATE_H

#include "model.h"
#include "quasi_newton.h"
#include <vector>
#include <cmath>
// does not set model
/*
void de_mutate_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a,
		  const model& m, rng& generator, const precalculate&,const igrid&, change& ,const vec&, quasi_newton&);

void de_mutate_receptor_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a,
			 const model& m, rng& generator, const precalculate&,const igrid&, change& ,const vec&, quasi_newton&);
			 

void de_mutate_ligand_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a,
		        const model& m, rng& generator, const precalculate&,const igrid&, change& ,const vec&, quasi_newton&);
				*/
				
void gwo_mutate_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a, rng& generator);

void gwo_mutate_sub_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a,
			 rng& generator, std::vector<int> &flex_list);
			 
void gwo_mutate_receptor_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a,
			 rng& generator);
			 
void gwo_mutate_receptor_sub_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a,
			 rng& generator, std::vector<int> &flex_list);

void gwo_mutate_ligand_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a,
		        rng& generator);

#endif
