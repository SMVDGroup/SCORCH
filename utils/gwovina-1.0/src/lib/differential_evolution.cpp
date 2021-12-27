//#include "de_search.h"
#include "coords.h"
#include "quasi_newton.h"
#include "de_mutate.h"
#include "mutate.h"
#include "differential_evolution.h"
#include <fstream>

output_type gwo::operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
	output_container tmp;
	this->operator()(m, tmp, p, ig, p_widened, ig_widened, corner1, corner2, increment_me, generator); // call the version that produces the whole container
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}

bool metropolis_accept(fl old_f, fl new_f, fl temperature, rng& generator) {
	if(new_f < old_f) return true;
	const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
	return random_fl(0, 1, generator) < acceptance_probability;
}
/*
void differential_evolution::init_parameters(int size, float F, float CR, int populations_size)
{

}*/

// out is sorted
struct Comp{
    bool operator()(const output_type & a, const output_type & b){
        return a.e < b.e;
    }
};

int in_bound(sz size)
{
	int num = size;
	num /= 2;
	if (num < 5)
		return 5;
	if (num > 12)
		return 12;
	return num;
}


template<typename T>
atom_range my_get_atom_range(const T& t) {
	atom_range tmp = t.node;
	VINA_FOR_IN(i, t.children) {
		atom_range r = my_get_atom_range(t.children[i]);
		if(tmp.begin > r.begin) tmp.begin = r.begin;
		if(tmp.end   < r.end  ) tmp.end   = r.end;
	}
	return tmp;
}

/*
std::vector<std::pair<int, int> > initialize_ligand_flex_pair(std::vector<atom_range> &ligand_range, std::vector<atom_range> &flex_sc)
{
	std::vector<std::pair<int, int> > ligand_flex_pairs;
	for (int i = 0; i < ligand_range.size(); i++)
	{
		atom_range lig = ligand_range[i];
		for (int j = lig.begin; j < lig.end; j++)
		{
			for (int k = 0; k < flex_sc.size; k++)
			{
				atom_range flex = flex_sc[k];
				
			}
		}
	}
}*/
bool in_range(atom_range ar, int index)
{
	bool result = false;
	if (index >= ar.begin && index < ar.end)
		result = true;
	return result;
}

bool is_clash(const interacting_pair& ip, const atomv &atoms, const vecv &coords) {
	
	const fl r = std::sqrt(vec_distance_sqr(coords[ip.a], coords[ip.b]));
	const fl covalent_r = atoms[ip.a].covalent_radius() + atoms[ip.b].covalent_radius();
	assert(r >= 0);
	assert(covalent_r > epsilon_fl);
	const fl x = r / covalent_r;
	if(x > 2) return false;
	return true;
}

std::vector<bool> find_clash_atom(model &m)
{
	std::vector<bool> clashing_atoms;
	for (int i = 0; i < m.num_movable_atoms(); i++)
	{
		clashing_atoms.push_back(false);
	}
	//std::ofstream outout("outout.dat");
	for (int j = 0; j < m.num_my_flex_pairs(); j++)
	{
		int a = m.my_flex_pairs[j].a;
		int b = m.my_flex_pairs[j].b;
		
		//outout << a << "	" << b << std::endl;
		if (clashing_atoms[a] && clashing_atoms[b])
			continue;
		if (is_clash(m.my_flex_pairs[j], m.atoms, m.coords))
		{
			//std::cout << a << "	" << b << std::endl;
			clashing_atoms[a] = true;
			clashing_atoms[b] = true;
		}
	}
	return clashing_atoms;
}

std::vector<int> find_clash_flex(model &m, std::vector<atom_range> &flex_sc)
{
	std::vector<int> result;
	std::vector<bool> clashing_atoms;
	clashing_atoms = find_clash_atom(m);
	/*
	for (int j = 0; j < clashing_atoms.size(); j++)
	{
		std::cout << clashing_atoms[j] << std::endl;
	}*/
	
	for (int i = 0; i < flex_sc.size(); i++)
	{
		atom_range ar = flex_sc[i];
		for (int j = ar.begin; j < ar.end; j++)
		{
			if (clashing_atoms[j])
			{
				result.push_back(i);
				break;
			}
		}
	}
	return result;
}

output_type gwo_sidechain(output_type &candidate, rng& generator, model& m, const precalculate& p, const igrid& ig, change& g, const vec& hunt_cap, quasi_newton& quasi_newton_par, std::vector<int> &clashing_flex, std::vector<atom_range> &flex_sc)
{
	if (clashing_flex.size() == 0)
		return candidate;
	output_type best = candidate;
	best.e = m.eval(p, ig, hunt_cap, best.c);
	int wolf_size = 5; //5
	int step = 10; //10
	vec authentic_v(1000, 1000, 1000);
	
	/*
	for (int i = 0 ; i < 10000; i++)
	{
		int which_flex = random_int(0, clashing_flex.size()-1, generator);
		{
			int which_torsion = random_int(0, candidate.c.flex[clashing_flex[which_flex]].torsions.size()-1, generator);
			{
				output_type tmp = candidate;
				tmp.c.flex[clashing_flex[which_flex]].torsions[which_torsion] += random_fl(-pi, pi, generator);
				tmp.e = m.eval(p, ig, hunt_cap, tmp.c);
				if (tmp.e < candidate.e)
				{
					candidate = tmp;
				}
			}
		}
	}
	quasi_newton_par(m, p, ig, candidate, g, hunt_cap);
	return candidate;
	*/
	
	std::vector<output_type> populations;
	for (int i = 0; i < wolf_size; i++)
	{
		output_type tmp = candidate;
		/*
		for (int j = 0; j < tmp.c.flex.size(); j++)
		{
			for (int k = 0; k < tmp.c.flex[j].torsions.size(); k++)
			{
						tmp.c.flex[j].torsions[k] = random_fl(-pi, pi, generator);
			}
		}*/
		
		for (int j = 0; j < clashing_flex.size(); j++)
		{
			for (int k = 0; k < tmp.c.flex[clashing_flex[j]].torsions.size(); k++)
			{
				tmp.c.flex[clashing_flex[j]].torsions[k] = random_fl(-pi, pi, generator);
			}
		}
		//quasi_newton_par(m, p, ig, tmp, g, authentic_v);
		tmp.e = 	m.eval(p, ig, authentic_v, tmp.c);
		populations.push_back(tmp);
	}
	
	std::priority_queue<std::pair<float, int> > best_team;
	for (int i = 0; i < step; i++)
	{
		
		for (int j = 0; j < populations.size(); j++)
		{
			best_team.push(std::make_pair(populations[j].e, j));
		}
		while(best_team.size() > 3)
		{
			best_team.pop();
		}
		int delta = best_team.top().second;
		best_team.pop();
		int beta = best_team.top().second;
		best_team.pop();
		int alpha = best_team.top().second;
		best_team.pop();
		std::vector<output_type> new_populations = populations;
		fl a = 2.0 * (1 - ((double) i / (double) step ));
		for (int k = 0; k < new_populations.size(); k++)
		{
			//gwo_mutate_receptor_conf(new_populations, k, new_populations[alpha], new_populations[beta], new_populations[delta],a, generator);
			gwo_mutate_receptor_sub_conf(new_populations, k, new_populations[alpha], new_populations[beta], new_populations[delta],a, generator, clashing_flex);

			new_populations[k].e = 	m.eval(p, ig, hunt_cap, new_populations[k].c);
			//quasi_newton_par(m, p, ig, new_populations[k], g, hunt_cap);
		}
		populations = new_populations;

	}
	for (int i = 0; i < populations.size(); i++)
	{
		if (populations[i].e < best.e)
		{
			best = populations[i];
		}
	}
	
	//quasi_newton_par(m, p, ig, best, g, hunt_cap);
	return best;
}

output_type gwo_ligand(output_type &candidate, rng& generator, model& m, const precalculate& p, const igrid& ig, change& g, const vec& hunt_cap, quasi_newton& quasi_newton_par, const vec& corner1, const vec& corner2)
{
	output_type best = candidate;
	int wolf_size = 5;
	int step = 100000;
	conf_size s = m.get_size();
	vec authentic_v(1000, 1000, 1000);
	bool better = false;
	bool first = true;
	int count = 0;

	std::vector<output_type> populations;
	std::vector<output_type> p_best;
	for (int i = 0; i < wolf_size; i++)
	{
		output_type tmp(s, 0);
		tmp.c.randomize(corner1, corner2, generator);
		tmp.c.flex = candidate.c.flex;
		quasi_newton_par(m, p, ig, tmp, g, authentic_v);
		p_best.push_back(tmp);
		populations.push_back(tmp);
	}
	
	std::priority_queue<std::pair<float, int> > best_team;
	
	int increment = 1;
	for (int i = 0; i < step; i += increment)
	{
		
		for (int j = 0; j < populations.size(); j++)
		{
			best_team.push(std::make_pair(populations[j].e, j));
		}
		while(best_team.size() > 3)
		{
			best_team.pop();
		}
		int delta = best_team.top().second;
		best_team.pop();
		int beta = best_team.top().second;
		best_team.pop();
		int alpha = best_team.top().second;
		best_team.pop();
		std::vector<output_type> new_populations = populations;
		fl a = 2.0 * (1 - ((double) i / (double) step ));
		for (int k = 0; k < new_populations.size(); k++)
		{
			gwo_mutate_ligand_conf(new_populations, k, new_populations[alpha], new_populations[beta], new_populations[delta],a, generator);
			quasi_newton_par(m, p, ig, new_populations[k], g, hunt_cap);
		}
		
		better = false;
		for (int j = 0; j < new_populations.size(); j++)
		{

			if (new_populations[j].e - p_best[j].e < -0.0001)
			{
				better = true;
				p_best[j] = new_populations[j];
				count = 0;
			}
		}
		
		if (!better)
		{
			count++;
		}
		
		increment = (int) pow(2, count);

		if (i > step / 2)
		{
			if (first)
			{
				first = false;
				count = 0;
				i = step / 2;
			}
		}
		
		populations = new_populations;
	}
	
	for (int i = 0; i < p_best.size(); i++)
	{
		if (p_best[i].e < best.e)
		{
			best = p_best[i];
		}
	}

	quasi_newton_par(m, p, ig, best, g, hunt_cap);
	return best;
}

/*std::vector<output_type>*/
output_type gwo_receptors(model& m, const precalculate& p, const vec& corner1, const vec& corner2, rng& generator, const vec& hunt_cap, int num_rec)
{
	std::vector<output_type> receptors;
	vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
	conf_size s = m.get_size();
	change g(s);
	fl best_e = max_fl;
	output_type best(s,0);
	
	//quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
	
	
	std::vector<output_type> populations;
	std::vector<output_type> p_best;
	std::priority_queue<std::pair<float, int> > best_team;
	bool better = false;
	bool first = true;
	int count = 0;

	int num_wolves = num_rec;
	
	for (int i = 0; i < num_wolves; i++)
	{
		output_type tmp(s,0);
		tmp.c.randomize(corner1, corner2, generator);
		tmp.e = m.eval_without_ligand(p, authentic_v, tmp.c);
		//std::cout << tmp.e << std::endl;
		populations.push_back(tmp);
		p_best.push_back(tmp);
	}
	
	int step = 100000;
	int increment = 1;
	for (int i = 0; i < step; i += increment)
	{
		for (int j = 0; j < populations.size(); j++)
		{
			best_team.push(std::make_pair(populations[j].e, j));
		}
		
		while(best_team.size() > 3)
		{
			best_team.pop();
		}
		int delta = best_team.top().second;
		best_team.pop();
		int beta = best_team.top().second;
		best_team.pop();
		int alpha = best_team.top().second;
		best_team.pop();
		fl a = 2.0 * (1 - ((double) i / (double) step));
		
		std::vector<output_type> new_populations = populations;

		for (int j = 0; j < new_populations.size(); j++)
		{
			{
				gwo_mutate_receptor_conf(new_populations, j, new_populations[alpha], new_populations[beta], new_populations[delta],a, generator);
			}

		}
		
		better = false;
		for (int j = 0; j < new_populations.size(); j++)
		{

			new_populations[j].e = m.eval_without_ligand(p, hunt_cap, new_populations[j].c);
			if (new_populations[j].e - p_best[j].e < -0.0001)
			{
				better = true;
				p_best[j] = new_populations[j];
				count = 0;
			}
		}
		
		populations = new_populations;
		
		if (!better)
		{
			count++;
		}
		
		increment = (int) pow(2, count);

		if (i > step / 2)
		{
			if (first)
			{
				first = false;
				count = 0;
				i = step / 2;
			}
		}
	}
	
	for (int i = 0; i < p_best.size(); i++)
	{
		if (p_best[i].e < best_e)
		{
			best_e = p_best[i].e;
			best = p_best[i];
		}
	}
	m.set(best.c);
	m.write_structure("my_structure.pdb");
	//return p_best;
	return best;
}

/*std::vector<output_type>*/
output_type gwo_conf(model& m, output_type &candidate, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, rng& generator, const vec& hunt_cap, quasi_newton& quasi_newton_par2, std::vector<atom_range> &flex_sc, int num_wolves, int flex_rand)
{
	//std::vector<output_type> receptors;
	vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
	conf_size s = m.get_size();
	change g(s);
	fl best_e = max_fl;
	output_type best = candidate;
	ssd ssd_par;
	quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
	
	
	std::vector<output_type> populations;
	std::vector<output_type> p_best;
	std::priority_queue<std::pair<float, int> > best_team;
	bool better = false;
	bool first = true;
	int count = 0;

	//int num_wolves = num_rec;
	//populations.push_back(candidate);
	for (int i = 0; i < num_wolves; i++)
	{
		output_type tmp = candidate;
		tmp.c.ligands[0].randomize(corner1, corner2, generator);
		quasi_newton_par(m, p, ig, tmp, g, hunt_cap);
		populations.push_back(tmp);
		p_best.push_back(tmp);
	}
	
	int step = 100000;
	int increment = 1;
	for (int i = 0; i < step; i += increment)
	{
		for (int j = 0; j < populations.size(); j++)
		{
			best_team.push(std::make_pair(populations[j].e, j));
		}
		
		while(best_team.size() > 3)
		{
			best_team.pop();
		}
		int delta = best_team.top().second;
		best_team.pop();
		int beta = best_team.top().second;
		best_team.pop();
		int alpha = best_team.top().second;
		best_team.pop();
		fl a = 2.0 * (1 - ((double) i / (double) step));
		
		std::vector<output_type> new_populations = populations;

		for (int j = 0; j < new_populations.size(); j++)
		{
			{
				//gwo_mutate_conf(new_populations, j, new_populations[alpha], new_populations[beta], new_populations[delta],a, generator);
				gwo_mutate_ligand_conf(new_populations, j, new_populations[alpha], new_populations[beta], new_populations[delta],a, generator);
			}
			std::vector<int> clashing_flex;
			clashing_flex.push_back(flex_rand);
			new_populations[j] = gwo_sidechain(new_populations[j], generator, m, p, ig,g, hunt_cap, quasi_newton_par, clashing_flex, flex_sc);
		}
		
		better = false;
		for (int j = 0; j < new_populations.size(); j++)
		{
			quasi_newton_par(m, p, ig, new_populations[j], g, hunt_cap);
			if (new_populations[j].e - p_best[j].e < -0.0001)
			{
				better = true;
				p_best[j] = new_populations[j];
				count = 0;
			}
		}
		
		populations = new_populations;
		
		if (!better)
		{
			count++;
		}
		
		increment = (int) pow(2, count);

		if (i > step / 2)
		{
			if (first)
			{
				first = false;
				count = 0;
				i = step / 2;
			}
		}
	}
	
	for (int i = 0; i < p_best.size(); i++)
	{
		if (p_best[i].e < best_e)
		{
			best_e = p_best[i].e;
			best = p_best[i];
		}
	}
	//m.set(best.c);
	//m.write_structure("my_structure.pdb");
	//return p_best;
	return best;
}

std::vector<output_type> gwo_sub_conf(model& m, output_type &candidate, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, rng& generator, const vec& hunt_cap, quasi_newton& quasi_newton_par2, std::vector<atom_range> &flex_sc, int num_wolves, std::vector<int> &flex_list)
{
	//std::vector<output_type> receptors;
	vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
	conf_size s = m.get_size();
	change g(s);
	fl best_e = max_fl;
	output_type best = candidate;
	ssd ssd_par;
	quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
	
	
	std::vector<output_type> populations;
	std::vector<output_type> p_best;
	std::priority_queue<std::pair<float, int> > best_team;
	bool better = false;
	bool first = true;
	int count = 0;

	//int num_wolves = 5;
	//populations.push_back(candidate);
	for (int i = 0; i < num_wolves; i++)
	{
		output_type tmp = candidate;
		tmp.c.ligands[0].randomize(corner1, corner2, generator);
		for (int j = 0; j < flex_list.size(); j++)
		{
			int f = flex_list[j];
			for (int k = 0; k < tmp.c.flex[f].torsions.size(); k++)
			{
				tmp.c.flex[f].torsions[k] = random_fl(-pi, pi, generator);
			}
		}
		quasi_newton_par(m, p, ig, tmp, g, hunt_cap);
		populations.push_back(tmp);
		p_best.push_back(tmp);
	}
	
	int step = 100000;
	int increment = 1;
	for (int i = 0; i < step; i += increment)
	{
		for (int j = 0; j < populations.size(); j++)
		{
			best_team.push(std::make_pair(populations[j].e, j));
		}
		
		while(best_team.size() > 3)
		{
			best_team.pop();
		}
		int delta = best_team.top().second;
		best_team.pop();
		int beta = best_team.top().second;
		best_team.pop();
		int alpha = best_team.top().second;
		best_team.pop();
		fl a = 2.0 * (1 - ((double) i / (double) step));
		
		std::vector<output_type> new_populations = populations;

		for (int j = 0; j < new_populations.size(); j++)
		{
			{
				//gwo_mutate_conf(new_populations, j, new_populations[alpha], new_populations[beta], new_populations[delta],a, generator);
				//gwo_mutate_ligand_conf(new_populations, j, new_populations[alpha], new_populations[beta], new_populations[delta],a, generator);
				gwo_mutate_sub_conf(new_populations, j, new_populations[alpha], new_populations[beta], new_populations[delta],a, generator, flex_list);
			}
			/*
			std::vector<int> clashing_flex;
			clashing_flex.push_back(flex_rand);
			new_populations[j] = gwo_sidechain(new_populations[j], generator, m, p, ig,g, hunt_cap, quasi_newton_par, clashing_flex, flex_sc);
			*/
		}
		
		better = false;
		for (int j = 0; j < new_populations.size(); j++)
		{
			quasi_newton_par(m, p, ig, new_populations[j], g, hunt_cap);
			if (new_populations[j].e - p_best[j].e < -0.0001)
			{
				better = true;
				p_best[j] = new_populations[j];
				count = 0;
			}
		}
		
		populations = new_populations;
		
		if (!better)
		{
			count++;
		}
		
		increment = (int) pow(2, count);

		if (i > step / 2)
		{
			if (first)
			{
				first = false;
				count = 0;
				i = step / 2;
			}
		}
	}
	
	for (int i = 0; i < p_best.size(); i++)
	{
		if (p_best[i].e < best_e)
		{
			best_e = p_best[i].e;
			best = p_best[i];
		}
	}
	//m.set(best.c);
	//m.write_structure("my_structure.pdb");
	return p_best;
	//return best;
}

void gwo::operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
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
	

	
	for(int j=0; j<m.flex.size(); j++)
	{
		const residue& res = m.flex[j];
		atom_range ar = my_get_atom_range(res);
		atom_range tmp(0, 0);
		tmp = ar;
		flex_sc.push_back(tmp);
	}
	
	int num_wolves = in_bound(m.ligand_degrees_of_freedom(0));
	//num_wolves = 12;
	std::cout << num_wolves << std::endl;
	output_type rec_candidate = gwo_receptors(m, p, corner1, corner2, generator, hunt_cap, num_wolves);
	for (int i = 0; i < num_wolves; i++)
	{
		
		output_type tmp(s,0);
		bool pass;
		int count = 0;
		//tmp.c.randomize(corner1, corner2, generator);
		
		do
		{
			pass = true;
			tmp.c.randomize(corner1, corner2, generator);
			m.set(tmp.c);
			tmp.coords = m.get_heavy_atom_movable_coords();
			
			for (int j = 0; j < populations.size(); j++)
			{
				fl dist = rmsd_upper_bound(tmp.coords, populations[j].coords);
				if ( dist < 5.0)
				{
					pass = false;
					count ++;
					break;
				}
			}
			if (count > 100)
				pass = true;

		} while (!pass);
		quasi_newton_par(m, p, ig, tmp, g, authentic_v);
		populations.push_back(tmp);
		p_best.push_back(tmp);
	}
	populations[0] = rec_candidate;
	p_best[0] = rec_candidate;
	
	double energy=0;
	int count_step = 0;
	int count=0;
	int generations = 0; // generation counting
	bool better = false;
	bool first = true;
	float improve_energy_each_generation = 0;
	fl best_energy = 10000;
	output_type best(s,0);
		
	VINA_U_FOR(step, num_steps){
		if(increment_me)
			++(*increment_me);
		generations++;
		for (int i = 0; i < populations.size(); i++)
		{
			best_team.push(std::make_pair(populations[i].e, i));
		}
		
		while(best_team.size() > 3)
		{

			best_team.pop();
		}
		int delta = best_team.top().second;
		best_team.pop();
		int beta = best_team.top().second;
		best_team.pop();
		int alpha = best_team.top().second;
		best_team.pop();
		fl a = 2.0 * (1 - ((double) step / (double) num_steps));
		
		std::vector<output_type> new_populations = populations;

		for (int i = 0; i < new_populations.size(); i++)
		{
			fl r = random_fl(0, 1, generator);
			if (r > 0.2)
			{
				gwo_mutate_conf(new_populations, i, new_populations[alpha], new_populations[beta], new_populations[delta],a, generator);
				//gwo_mutate_ligand_conf(new_populations, i, new_populations[alpha], new_populations[beta], new_populations[delta],a, generator);
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
				/*
				m.set(new_populations[i].c);
				std::vector<int> clashing_flex;
				clashing_flex = find_clash_flex(m, flex_sc);
				new_populations[i] = gwo_sidechain(new_populations[i], generator, m, p, ig,g, hunt_cap, quasi_newton_par, clashing_flex, flex_sc);
				*/
				
		}
		
		better = false;
		for (int i = 0; i < new_populations.size(); i++)
		{
			quasi_newton_par(m, p, ig, new_populations[i], g, hunt_cap);
			if (new_populations[i].e - p_best[i].e < -0.0001)
			{
				better = true;
				p_best[i] = new_populations[i];
				count = 0;
				/*
				output_type tmp = new_populations[i];
				//quasi_newton_par(m, p, ig, tmp, g, authentic_v);
				m.set(tmp.c);
				tmp.coords = m.get_heavy_atom_movable_coords();
				add_to_output_container(out, tmp, min_rmsd, num_saved_mins);
				*/
			}
			
			if (new_populations[i].e < best_energy)
			{
				best_energy = new_populations[i].e;
				best = new_populations[i];
			}

		}
		

		
		if (!better)
		{
			count++;
			step = step + (int) pow(2, count);
		}
		if (step > num_steps / 2)
		{
			if (first)
			{
				first = false;
				count = 0;
				step = num_steps / 2;
			}
		}
		//populations = new_populations;
		populations = p_best;
	}
	
	for (int i = 0; i < p_best.size(); i++)
	{
		//output_type tmp = gwo_ligand(p_best[i], generator, m, p, ig, g, hunt_cap, quasi_newton_par, corner1, corner2);
		
		
		output_type tmp = p_best[i];
		for (int j = 0; j < 10; j++)
		{
			quasi_newton_par(m, p, ig, tmp, g, authentic_v);
		}
		
		if (tmp.e < best_energy)
		{
			best_energy = tmp.e;
			best = tmp;
		}
		m.set(tmp.c);
		tmp.coords = m.get_heavy_atom_movable_coords();
		add_to_output_container(out, tmp, min_rmsd, num_saved_mins);
	}

	std::cout <<	generations << std::endl;
	std::cout << "best energy: " << best_energy << std::endl;
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order

}
