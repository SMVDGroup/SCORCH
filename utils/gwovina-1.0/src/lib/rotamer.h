#ifndef ROTAMER_H
#define ROTAMER_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "random.h"

struct rot{
	float prob;
	std::vector<float> chi;
	std::vector<float> sig;
	
	void print() const
	{
		std::cout << prob << "	";
		for (int i = 0; i < chi.size(); i++)
		{
			std::cout << chi[i] << "	";
		}
		for (int i = 0; i < sig.size(); i++)
		{
			std::cout << sig[i] << "	";
		}
		std::cout << std::endl;
	}
};

struct flex{
	std::vector<float> ori;
	std::vector<rot> rots;
	int rand_rot()
	{
		//double x = (double) rand() / (RAND_MAX + 1.0);
		//double ori_x = x;
		rng generator;
		float x = random_fl(0,1, generator);
		for (int i =0; i < rots.size(); i++)
		{
			if (x < rots[i].prob)
			{
				//std::cout << ori_x << std::endl;
				return i;
			}
			else{
				x -= rots[i].prob;
			}
		}
		//std::cout << ori_x << std::endl;

		return 0;
	}
};
 
class rotamer{
	public:
		rotamer();
		rotamer(std::string filename);
		void read_rotamer(std::string filename);
		int get_size();
		void print();
		void print_res(int id);
		void print_rot(rot r);
		std::vector<flex> res;

	private:
	};

#endif
