#include "coords.h"
#include "quasi_newton.h"
#include "gwo_mutate.h"
#include "mutate.h"
#include "gwo.h"
#include <fstream>

output_type grey_wolf_optimizer::operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator, int num_wolves) const {
	output_container tmp;
	this->operator()(m, tmp, p, ig, p_widened, ig_widened, corner1, corner2, increment_me, generator, num_wolves); // call the version that produces the whole container
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}

bool metropolis_accept(fl old_f, fl new_f, fl temperature, rng& generator) {
	if(new_f < old_f) return true;
	const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
	return random_fl(0, 1, generator) < acceptance_probability;
}

// out is sorted
struct Comp{
    bool operator()(const output_type & a, const output_type & b){
        return a.e < b.e;
    }
};

int in_bound(sz size)
{
	int num = size;
	//num /= 2;
	/*
	if (num < 5)
		return 5;
	if (num > 12)
		return 12;
	*/
	return num;
}

int cal_leaders(sz num_wolves)
{
	int num = num_wolves;
	if (num < 12)
		return 3;
	else
	{
		if (num % 4 == 0)
			return num / 4;
		else
			return 1 + num / 4;
	}
}


void generate_population(model& m, const precalculate& p, const igrid& ig, change& g,
						 std::vector<output_type>& populations, int& num_wolves, const vec& corner1, const vec& corner2, 
						 rng& generator, rotamer &rotamers, vec &authentic_v, quasi_newton &quasi_newton_par)
{
	conf_size s = m.get_size();
	output_type tmp(s, 0);
	
	tmp.c.randomize(corner1, corner2, generator);
	for (int i = 0 ; i < tmp.c.flex.size(); i++)
	{
		tmp.c.flex[i].set_to_null();
	}
	quasi_newton_par(m, p, ig, tmp, g, authentic_v);

	populations.push_back(tmp);
	for (int i = 1; i < num_wolves; i++)
	{
		output_type tmp(s, 0);
		bool accept = false;
		tmp.c.randomize(corner1, corner2, generator);
		do
		{
			accept = false;
			for (int j = 0; j < rotamers.res.size(); j++)
			{
				int rot_index = rotamers.res[j].rand_rot();
				rot& r = rotamers.res[j].rots[rot_index];
				for (int k = 0; k < tmp.c.flex[j].torsions.size(); k++)
				{
					float mean = r.chi[k];
					float sig = r.sig[k];
					tmp.c.flex[j].torsions[k] = (random_normal(mean, sig, generator) - rotamers.res[j].ori[k]) / 180 * pi ;
					
					for (int kk = 0; kk < populations.size(); kk++)
					{
						if (abs(populations[kk].c.flex[j].torsions[k] - tmp.c.flex[j].torsions[k]) > 0.87)
						{
							accept = true;
						}
					}
				}
			}
			quasi_newton_par(m, p, ig, tmp, g, authentic_v);

		} while (!accept);
		populations.push_back(tmp);
	}
}

void grey_wolf_optimizer::operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator, int num_wolves) const {
	//std::cout<<"a";   //testing

	vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
	conf_size s = m.get_size();
	change g(s);
	fl best_e = max_fl;
	quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
	
	std::vector<output_type> populations;	//wolf polulations
	std::priority_queue<std::pair<float, int> > best_team;
	std::vector<output_type> p_best;
	std::vector<fl> clash;
	output_container tmp_populations;
	std::vector<atom_range> flex_sc;
	std::vector<atom_range> ligand_range;
	std::vector<std::pair<int, int> > ligand_flex_pairs;
	std::vector<int> clashing_flex;
	rotamer rotamers = m.rotamers;
	
	/*
	int num_wolves = 5 + in_bound(m.ligand_degrees_of_freedom(0));
	int max_steps = 100 + 20 * m.ligand_degrees_of_freedom(0);
	int num_leaders = cal_leaders(num_wolves);
	*/
	
	int max_steps = 100000;
	//int num_wolves = 12;
	int num_leaders = 3;
		
	/*
	int max_steps = num_steps;
	int num_wolves = populations_size;
	int num_leaders = 3;
	*/
	
	if (rotamers.get_size() != 0)
	{
		generate_population(m, p, ig, g, populations, num_wolves, corner1, corner2, generator, rotamers, authentic_v, quasi_newton_par);
		p_best = populations;
	}
	else
	{
		for (int i = 0; i < num_wolves; i++)
		{
			
			output_type tmp(s,0);
			bool pass;	
			int count = 0;
			tmp.c.randomize(corner1, corner2, generator);
			
			quasi_newton_par(m, p, ig, tmp, g, authentic_v);
			populations.push_back(tmp);
			p_best.push_back(tmp);
		}
	}
	
	
	double energy=0;
	int count_step = 0;
	int count=0;
	int generations = 0; // generation counting
	bool better = false;
	bool first = true;
	float improve_energy_each_generation = 0;
	fl best_energy = 10000;
	output_type best(s,0);
	
	for (int i = 0; i < populations.size(); i++)
	{
		best_team.push(std::make_pair(populations[i].e, i));
	}
	while(best_team.size() > num_leaders)
	{

		best_team.pop();
	}
	
	VINA_U_FOR(step, max_steps){
		if(increment_me)
			++(*increment_me);
		generations++;
		
		fl a = 2.0 * (1 - ((double) step / (double) max_steps)); //2.0*
		
		std::vector<output_type> new_populations = populations;
		
		for (int i = 0; i < new_populations.size(); i++)
		{
			fl r = random_fl(0, 1, generator);
			if (r > 0.2)
			{
				std::vector<output_type> leaders;
				std::priority_queue<std::pair<float, int> > best_team_tmp = best_team;
				while(best_team_tmp.size() > 0)
				{
					leaders.push_back(new_populations[best_team_tmp.top().second]);
					best_team_tmp.pop();
				}
				gwo_mutate_conf_team2(new_populations, i, leaders, a, generator);
				
			}
			
			else
			{
				output_type candidate = new_populations[i];
				float F = random_fl(0, 0.2, generator);
				candidate.c.ligands[0].rigid.position += F * random_inside_sphere(generator);
				fl gr = m.gyration_radius(0);
				F = random_fl(0, 0.2, generator);
				if(gr > epsilon_fl) { // FIXME? just doing nothing for 0-radius molecules. do some other mutation?
					vec rotation; 
					rotation = F / gr * random_inside_sphere(generator); 
					quaternion_increment(candidate.c.ligands[0].rigid.orientation, rotation);
				}
				
				for (int j = 0; j < candidate.c.ligands[0].torsions.size(); j++)
				{
					F = random_fl(0, 0.2, generator);
					candidate.c.ligands[0].torsions[j] += F * random_fl(-pi, pi, generator);
				}
				for (int j = 0; j < candidate.c.flex.size(); j++)
				{
					for (int k = 0; k < candidate.c.flex[j].torsions.size(); k++)
					{
						F = random_fl(0, 0.2, generator);
						candidate.c.flex[j].torsions[k] += F*random_fl(-pi, pi, generator);
					}
				}
				new_populations[i] = candidate;
			}
			
			quasi_newton_par(m, p, ig, new_populations[i], g, hunt_cap);
			best_team.push(std::make_pair(new_populations[i].e, i));
			best_team.pop();
			
		}
		
		better = false;
		for (int i = 0; i < new_populations.size(); i++)
		{
			if (new_populations[i].e - p_best[i].e < -0.0001)
			{
				better = true;
				p_best[i] = new_populations[i];
				count = 0;
				// add to the result list every loops				
				output_type tmp = new_populations[i];
				m.set(tmp.c);
				tmp.coords = m.get_heavy_atom_movable_coords();
				add_to_output_container(out, tmp, min_rmsd, num_saved_mins);
				
			}
			
		}
		
		if (!better)
		{
			count++;
			step = step + (int) pow(2.0, count) - 1;
		}
		if (step > max_steps / 2)
		{
			if (first)
			{
				first = false;
				count = 0;
				step = max_steps / 2;
			}

		}
		
		populations = p_best;
	}
	
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order

}
