#include "rotamer.h"


rotamer::rotamer()
{
	
}

rotamer::rotamer(std::string filename)
{
	read_rotamer(filename);
}

int rotamer::get_size()
{
	return res.size();
}

void rotamer::read_rotamer(std::string filename)
{
	std::fstream in_file(filename.c_str());
	std::string line;
	while (getline(in_file, line))
	{
		int flex_size;
		if (line.substr(0, 8) == "BEGINROT")
		{
			std::string tmp;
			getline(in_file, tmp);	//res_id
			//std::cout << tmp << std::endl;
			getline(in_file, tmp);	//ori angle
			std::vector<float> ori_angle;
			std::vector<rot> rots;
			std::stringstream ss(tmp);
			float tmp_fl;
			while (ss >> tmp_fl)
			{
				ori_angle.push_back(tmp_fl);
				//std::cout << tmp_fl << std::endl;
			}
			flex_size = ori_angle.size();
			
			std::string tmp_line;
			while (getline(in_file, tmp_line))
			{
				
				if (tmp_line.substr(0, 6) == "ENDROT")
				{
					flex tmp_flex;
					tmp_flex.ori = ori_angle;
					tmp_flex.rots = rots;
					res.push_back(tmp_flex);
					break;
				}
				else
				{
					std::stringstream ss(tmp_line);
					std::string rubbish;
					float value;
					for (int i = 0; i < 8; i++)
					{
						ss >> rubbish;
					}
					float p;
					float chi;
					float sig;
					std::vector<float> chi_vec;
					std::vector<float> sig_vec;
					ss >> p;
					for (int i = 0; i < flex_size; i++) // mean
					{
						ss >> chi;
						chi_vec.push_back(chi);
					}
					for (int i = 0; i < 4 - flex_size; i++)
					{
						ss >> rubbish;
					}
					
					for (int i = 0; i < flex_size; i++)  //sig
					{
						ss >> sig;
						sig_vec.push_back(sig);
					}
					for (int i = 0; i < 4 - flex_size; i++)
					{
						ss >> rubbish;
					}
					rot r;
					r.prob = p;
					r.chi = chi_vec;
					r.sig = sig_vec;
					if (p != 0)
						rots.push_back(r);
				}
			}
		}

	}
}

void rotamer::print()
{
	for (int i = 0; i < res.size(); i++)
	{
		print_res(i);
	}
}

void rotamer::print_res(int id)
{
	const flex &tmp = res[id];
	for (int i = 0; i < tmp.ori.size(); i++)
	{
		std::cout << tmp.ori[i] << "	";
	}
	std::cout << std::endl;
	for (int j = 0; j < tmp.rots.size(); j++)
	{
		tmp.rots[j].print();
	}
}
