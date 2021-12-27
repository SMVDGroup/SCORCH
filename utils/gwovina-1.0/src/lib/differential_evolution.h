#ifndef VINA_DIFFERENTIAL_EVOLUTION_H
#define VINA_DIFFERENTIAL_EVOLUTION_H

#include "ssd.h"
#include "incrementable.h"
#include "common.h"
#include "conf.h"
#include <vector>
#include <queue>

struct gwo{
	unsigned num_steps;
	fl temperature;
	vec hunt_cap;
	fl min_rmsd;
	sz num_saved_mins;
	fl mutation_amplitude;
	ssd ssd_par;
	
	//std::vector<output_tpye> populations;
	float F;
	float CR;
	int populations_size;
	gwo() : num_steps(100000), temperature(1.2), hunt_cap(10, 1.5, 10), min_rmsd(0.5), num_saved_mins(50), mutation_amplitude(2), F(0.5), CR(0.8), populations_size(5) {} // T = 600K, R = 2cal/(K*mol) -> temperature = RT = 1.2;  num_steps = 50*lig_atoms = 2500

	output_type operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const;
	// out is sorted
	//void operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator, int num_birds, double w, double c1, double c2) const;
	void operator()(model& m, output_container& out, 
					const precalculate& p, const igrid& ig, 
					const precalculate& p_widened, const igrid& ig_widened, 
					const vec& corner1, const vec& corner2,
					incrementable* increment_me, rng& generator) const;
	//void init_parameters(int size, float F, float CR, int populations_size);
};

#endif

