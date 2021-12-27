/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#include "de_mutate.h"

sz count_mutable_entities(const conf& c) {
	sz counter = 0;
	VINA_FOR_IN(i, c.ligands)
		counter += 2 + c.ligands[i].torsions.size();
	VINA_FOR_IN(i, c.flex)
		counter += c.flex[i].torsions.size();
	return counter;
}

sz count_ligand_mutable_entities(const conf& c) {
	sz counter = 0;
	VINA_FOR_IN(i, c.ligands)
		counter += 2 + c.ligands[i].torsions.size();
	//VINA_FOR_IN(i, c.flex)
	//	counter += c.flex[i].torsions.size();
	return counter;
}

sz count_receptor_mutable_entities(const conf& c) {
	sz counter = 0;
	//VINA_FOR_IN(i, c.ligands)
		//counter += 2 + c.ligands[i].torsions.size();
	VINA_FOR_IN(i, c.flex)
		counter += c.flex[i].torsions.size();
	return counter;
}


void gwo_mutate_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a,
					rng& generator) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion
	
	int changing = 0;
	sz mutable_entities_num = count_mutable_entities(populations[position].c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	
	output_type mutating_candidate = populations[position];
	//output_type difference = mutation(populations, position, F, generator);
	//output_type movement = cal_movement(mutating_candidate, alpha, beta, delta, generator);
	//output_type new_candidate = populations[position];
	//output_type movement;
	/*********************************position****************************/
	fl A1, A2, A3, C1, C2 ,C3, r1,r2;
	//conf_type D1, D2, D3, X1, X2, X3;
	output_type D1 = mutating_candidate;
	output_type D2 = mutating_candidate;
	output_type D3 = mutating_candidate;
	output_type X1 = mutating_candidate;
	output_type X2 = mutating_candidate;
	output_type X3 = mutating_candidate;
	output_type movement = mutating_candidate;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A1 = 2 * a * r1 - a;
	C1 = 2 * r2;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A2 = 2 * a * r1 - a;
	C2 = 2 * r2;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A3 = 2 * a * r1 - a;
	C3 = 2 * r2;
	
	D1.c.ligands[0].rigid.position = C1 * alpha.c.ligands[0].rigid.position - mutating_candidate.c.ligands[0].rigid.position;
	D2.c.ligands[0].rigid.position = C2 * beta.c.ligands[0].rigid.position - mutating_candidate.c.ligands[0].rigid.position;
	D3.c.ligands[0].rigid.position = C3 * delta.c.ligands[0].rigid.position - mutating_candidate.c.ligands[0].rigid.position;
	D1.c.ligands[0].rigid.position[0] = abs(D1.c.ligands[0].rigid.position[0]);
	D1.c.ligands[0].rigid.position[1] = abs(D1.c.ligands[0].rigid.position[1]);
	D1.c.ligands[0].rigid.position[2] = abs(D1.c.ligands[0].rigid.position[2]);
	D2.c.ligands[0].rigid.position[0] = abs(D1.c.ligands[0].rigid.position[0]);
	D2.c.ligands[0].rigid.position[1] = abs(D1.c.ligands[0].rigid.position[1]);
	D2.c.ligands[0].rigid.position[2] = abs(D1.c.ligands[0].rigid.position[2]);
	D3.c.ligands[0].rigid.position[0] = abs(D1.c.ligands[0].rigid.position[0]);
	D3.c.ligands[0].rigid.position[1] = abs(D1.c.ligands[0].rigid.position[1]);
	D3.c.ligands[0].rigid.position[2] = abs(D1.c.ligands[0].rigid.position[2]);
	
	X1.c.ligands[0].rigid.position = alpha.c.ligands[0].rigid.position - A1 * D1.c.ligands[0].rigid.position;
	X2.c.ligands[0].rigid.position = beta.c.ligands[0].rigid.position - A2 * D2.c.ligands[0].rigid.position;
	X3.c.ligands[0].rigid.position = delta.c.ligands[0].rigid.position - A3 * D3.c.ligands[0].rigid.position;
	
	movement.c.ligands[0].rigid.position = (X1.c.ligands[0].rigid.position + X2.c.ligands[0].rigid.position + X3.c.ligands[0].rigid.position) * (1.0/3.0);
	/********************************************orientation****************************************/
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A1 = 2 * a * r1 - a;
	C1 = 2 * r2;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A2 = 2 * a * r1 - a;
	C2 = 2 * r2;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A3 = 2 * a * r1 - a;
	C3 = 2 * r2;
	/*
	D1.c.ligands[0].rigid.orientation = abs(angle_to_quaternion(C1 * quaternion_to_angle(alpha.c.ligands[0].rigid.orientation) - quaternion_to_angle( mutating_candidate.c.ligands[0].rigid.orientation)));
	D2.c.ligands[0].rigid.orientation = abs(angle_to_quaternion(C2 * quaternion_to_angle(beta.c.ligands[0].rigid.orientation) - quaternion_to_angle( mutating_candidate.c.ligands[0].rigid.orientation)));
	D3.c.ligands[0].rigid.orientation = abs(angle_to_quaternion(C3 * quaternion_to_angle(delta.c.ligands[0].rigid.orientation) - quaternion_to_angle( mutating_candidate.c.ligands[0].rigid.orientation)));

	X1.c.ligands[0].rigid.orientation = angle_to_quaternion( quaternion_to_angle(alpha.c.ligands[0].rigid.orientation) - A1 * quaternion_to_angle( D1.c.ligands[0].rigid.orientation));
	X2.c.ligands[0].rigid.orientation = angle_to_quaternion( quaternion_to_angle(beta.c.ligands[0].rigid.orientation) - A2 * quaternion_to_angle( D2.c.ligands[0].rigid.orientation));
	X3.c.ligands[0].rigid.orientation = angle_to_quaternion( quaternion_to_angle(delta.c.ligands[0].rigid.orientation) - A3 * quaternion_to_angle( D3.c.ligands[0].rigid.orientation));
	*/
	vec rotation = quaternion_to_angle(alpha.c.ligands[0].rigid.orientation);
	vec D1_rotation, D2_rotation, D3_rotation, X1_rotation, X2_rotation, X3_rotation;
	D1_rotation = quaternion_difference(my_angle_to_quaternion(rotation, C1), mutating_candidate.c.ligands[0].rigid.orientation);
	rotation = quaternion_to_angle(beta.c.ligands[0].rigid.orientation);
	D2_rotation = quaternion_difference(my_angle_to_quaternion(rotation, C2), mutating_candidate.c.ligands[0].rigid.orientation);
	rotation = quaternion_to_angle(delta.c.ligands[0].rigid.orientation);
	D3_rotation = quaternion_difference(my_angle_to_quaternion(rotation, C3), mutating_candidate.c.ligands[0].rigid.orientation);
	
	X1_rotation = quaternion_difference(alpha.c.ligands[0].rigid.orientation, my_angle_to_quaternion(D1_rotation, A1));
	X2_rotation = quaternion_difference(beta.c.ligands[0].rigid.orientation, my_angle_to_quaternion(D2_rotation, A2));
	X3_rotation = quaternion_difference(delta.c.ligands[0].rigid.orientation, my_angle_to_quaternion(D3_rotation, A3));
	
	qt tmp = my_angle_to_quaternion(X1_rotation, 1/3);
	quaternion_increment(tmp, quaternion_to_angle(my_angle_to_quaternion(X2_rotation, 1/3)));
	quaternion_increment(tmp, quaternion_to_angle(my_angle_to_quaternion(X2_rotation, 1/3)));
	
	//movement.c.ligands[0].rigid.orientation = angle_to_quaternion((quaternion_to_angle(X1.c.ligands[0].rigid.orientation) + quaternion_to_angle(X2.c.ligands[0].rigid.orientation) + quaternion_to_angle(X3.c.ligands[0].rigid.orientation))*(1.0 / 3.0));
	movement.c.ligands[0].rigid.orientation = tmp;
	/***********************************************torsion****************************************/
	for (int i = 0; i < mutating_candidate.c.ligands[0].torsions.size();i++)
	{
		fl A1, A2, A3, C1, C2 ,C3, r1,r2;
		r1 = random_fl(0.0, 1.0, generator);
		r2 = random_fl(0.0, 1.0, generator);
		A1 = 2 * a * r1 - a;
		C1 = 2 * r2;

		r1 = random_fl(0.0, 1.0, generator);
		r2 = random_fl(0.0, 1.0, generator);
		A2 = 2 * a * r1 - a;
		C2 = 2 * r2;

		r1 = random_fl(0.0, 1.0, generator);
		r2 = random_fl(0.0, 1.0, generator);
		A3 = 2 * a * r1 - a;
		C3 = 2 * r2;
		
		D1.c.ligands[0].torsions[i] = abs(C1 * alpha.c.ligands[0].torsions[i] - mutating_candidate.c.ligands[0].torsions[i]);
		D2.c.ligands[0].torsions[i] = abs(C2 * beta.c.ligands[0].torsions[i] - mutating_candidate.c.ligands[0].torsions[i]);
		D3.c.ligands[0].torsions[i] = abs(C3 * delta.c.ligands[0].torsions[i] - mutating_candidate.c.ligands[0].torsions[i]);
		
		X1.c.ligands[0].torsions[i] =  alpha.c.ligands[0].torsions[i] - A1 * D1.c.ligands[0].torsions[i];
		X2.c.ligands[0].torsions[i] =  beta.c.ligands[0].torsions[i] - A2 * D2.c.ligands[0].torsions[i];
		X3.c.ligands[0].torsions[i] =  delta.c.ligands[0].torsions[i] - A3 * D3.c.ligands[0].torsions[i];
		
		movement.c.ligands[0].torsions[i] = (X1.c.ligands[0].torsions[i] + X2.c.ligands[0].torsions[i] + X3.c.ligands[0].torsions[i]) / 3.0;
	}
	/*****************************************************flex**********************************************/
	for (int i = 0; i < mutating_candidate.c.flex.size();i++)
	{
		for (int j = 0; j < mutating_candidate.c.flex[i].torsions.size();j++)
		{
			fl A1, A2, A3, C1, C2 ,C3, r1,r2;
			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A1 = 2 * a * r1 - a;
			C1 = 2 * r2;

			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A2 = 2 * a * r1 - a;
			C2 = 2 * r2;

			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A3 = 2 * a * r1 - a;
			C3 = 2 * r2;
			
			D1.c.flex[i].torsions[j] = abs(C1 * alpha.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j]);
			D2.c.flex[i].torsions[j] = abs(C2 * beta.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j]);
			D3.c.flex[i].torsions[j] = abs(C3 * delta.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j]);
			
			X1.c.flex[i].torsions[j] =  alpha.c.flex[i].torsions[j] - A1 * D1.c.flex[i].torsions[j];
			X2.c.flex[i].torsions[j] =  beta.c.flex[i].torsions[j] - A2 * D2.c.flex[i].torsions[j];
			X3.c.flex[i].torsions[j] =  delta.c.flex[i].torsions[j] - A3 * D3.c.flex[i].torsions[j];
			
			movement.c.flex[i].torsions[j] = (X1.c.flex[i].torsions[j] + X2.c.flex[i].torsions[j] + X3.c.flex[i].torsions[j]) / 3.0;
		}
	}
	populations[position] = movement;
	/*******************************************************************************************************/
}

void gwo_mutate_sub_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a,
			 rng& generator, std::vector<int> &flex_list){
		int changing = 0;
	sz mutable_entities_num = count_mutable_entities(populations[position].c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	
	output_type mutating_candidate = populations[position];
	//output_type difference = mutation(populations, position, F, generator);
	//output_type movement = cal_movement(mutating_candidate, alpha, beta, delta, generator);
	//output_type new_candidate = populations[position];
	//output_type movement;
	/*********************************position****************************/
	fl A1, A2, A3, C1, C2 ,C3, r1,r2;
	//conf_type D1, D2, D3, X1, X2, X3;
	output_type D1 = mutating_candidate;
	output_type D2 = mutating_candidate;
	output_type D3 = mutating_candidate;
	output_type X1 = mutating_candidate;
	output_type X2 = mutating_candidate;
	output_type X3 = mutating_candidate;
	output_type movement = mutating_candidate;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A1 = 2 * a * r1 - a;
	C1 = 2 * r2;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A2 = 2 * a * r1 - a;
	C2 = 2 * r2;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A3 = 2 * a * r1 - a;
	C3 = 2 * r2;
	
	D1.c.ligands[0].rigid.position = C1 * alpha.c.ligands[0].rigid.position - mutating_candidate.c.ligands[0].rigid.position;
	D2.c.ligands[0].rigid.position = C2 * beta.c.ligands[0].rigid.position - mutating_candidate.c.ligands[0].rigid.position;
	D3.c.ligands[0].rigid.position = C3 * delta.c.ligands[0].rigid.position - mutating_candidate.c.ligands[0].rigid.position;
	D1.c.ligands[0].rigid.position[0] = abs(D1.c.ligands[0].rigid.position[0]);
	D1.c.ligands[0].rigid.position[1] = abs(D1.c.ligands[0].rigid.position[1]);
	D1.c.ligands[0].rigid.position[2] = abs(D1.c.ligands[0].rigid.position[2]);
	D2.c.ligands[0].rigid.position[0] = abs(D1.c.ligands[0].rigid.position[0]);
	D2.c.ligands[0].rigid.position[1] = abs(D1.c.ligands[0].rigid.position[1]);
	D2.c.ligands[0].rigid.position[2] = abs(D1.c.ligands[0].rigid.position[2]);
	D3.c.ligands[0].rigid.position[0] = abs(D1.c.ligands[0].rigid.position[0]);
	D3.c.ligands[0].rigid.position[1] = abs(D1.c.ligands[0].rigid.position[1]);
	D3.c.ligands[0].rigid.position[2] = abs(D1.c.ligands[0].rigid.position[2]);
	
	X1.c.ligands[0].rigid.position = alpha.c.ligands[0].rigid.position - A1 * D1.c.ligands[0].rigid.position;
	X2.c.ligands[0].rigid.position = beta.c.ligands[0].rigid.position - A2 * D2.c.ligands[0].rigid.position;
	X3.c.ligands[0].rigid.position = delta.c.ligands[0].rigid.position - A3 * D3.c.ligands[0].rigid.position;
	
	movement.c.ligands[0].rigid.position = (X1.c.ligands[0].rigid.position + X2.c.ligands[0].rigid.position + X3.c.ligands[0].rigid.position) * (1.0/3.0);
	/********************************************orientation****************************************/
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A1 = 2 * a * r1 - a;
	C1 = 2 * r2;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A2 = 2 * a * r1 - a;
	C2 = 2 * r2;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A3 = 2 * a * r1 - a;
	C3 = 2 * r2;
	/*
	D1.c.ligands[0].rigid.orientation = abs(angle_to_quaternion(C1 * quaternion_to_angle(alpha.c.ligands[0].rigid.orientation) - quaternion_to_angle( mutating_candidate.c.ligands[0].rigid.orientation)));
	D2.c.ligands[0].rigid.orientation = abs(angle_to_quaternion(C2 * quaternion_to_angle(beta.c.ligands[0].rigid.orientation) - quaternion_to_angle( mutating_candidate.c.ligands[0].rigid.orientation)));
	D3.c.ligands[0].rigid.orientation = abs(angle_to_quaternion(C3 * quaternion_to_angle(delta.c.ligands[0].rigid.orientation) - quaternion_to_angle( mutating_candidate.c.ligands[0].rigid.orientation)));

	X1.c.ligands[0].rigid.orientation = angle_to_quaternion( quaternion_to_angle(alpha.c.ligands[0].rigid.orientation) - A1 * quaternion_to_angle( D1.c.ligands[0].rigid.orientation));
	X2.c.ligands[0].rigid.orientation = angle_to_quaternion( quaternion_to_angle(beta.c.ligands[0].rigid.orientation) - A2 * quaternion_to_angle( D2.c.ligands[0].rigid.orientation));
	X3.c.ligands[0].rigid.orientation = angle_to_quaternion( quaternion_to_angle(delta.c.ligands[0].rigid.orientation) - A3 * quaternion_to_angle( D3.c.ligands[0].rigid.orientation));
	*/
	vec rotation = quaternion_to_angle(alpha.c.ligands[0].rigid.orientation);
	vec D1_rotation, D2_rotation, D3_rotation, X1_rotation, X2_rotation, X3_rotation;
	D1_rotation = quaternion_difference(my_angle_to_quaternion(rotation, C1), mutating_candidate.c.ligands[0].rigid.orientation);
	rotation = quaternion_to_angle(beta.c.ligands[0].rigid.orientation);
	D2_rotation = quaternion_difference(my_angle_to_quaternion(rotation, C2), mutating_candidate.c.ligands[0].rigid.orientation);
	rotation = quaternion_to_angle(delta.c.ligands[0].rigid.orientation);
	D3_rotation = quaternion_difference(my_angle_to_quaternion(rotation, C3), mutating_candidate.c.ligands[0].rigid.orientation);
	
	X1_rotation = quaternion_difference(alpha.c.ligands[0].rigid.orientation, my_angle_to_quaternion(D1_rotation, A1));
	X2_rotation = quaternion_difference(beta.c.ligands[0].rigid.orientation, my_angle_to_quaternion(D2_rotation, A2));
	X3_rotation = quaternion_difference(delta.c.ligands[0].rigid.orientation, my_angle_to_quaternion(D3_rotation, A3));
	
	qt tmp = my_angle_to_quaternion(X1_rotation, 1/3);
	quaternion_increment(tmp, quaternion_to_angle(my_angle_to_quaternion(X2_rotation, 1/3)));
	quaternion_increment(tmp, quaternion_to_angle(my_angle_to_quaternion(X2_rotation, 1/3)));
	
	//movement.c.ligands[0].rigid.orientation = angle_to_quaternion((quaternion_to_angle(X1.c.ligands[0].rigid.orientation) + quaternion_to_angle(X2.c.ligands[0].rigid.orientation) + quaternion_to_angle(X3.c.ligands[0].rigid.orientation))*(1.0 / 3.0));
	movement.c.ligands[0].rigid.orientation = tmp;
	/***********************************************torsion****************************************/
	for (int i = 0; i < mutating_candidate.c.ligands[0].torsions.size();i++)
	{
		fl A1, A2, A3, C1, C2 ,C3, r1,r2;
		r1 = random_fl(0.0, 1.0, generator);
		r2 = random_fl(0.0, 1.0, generator);
		A1 = 2 * a * r1 - a;
		C1 = 2 * r2;

		r1 = random_fl(0.0, 1.0, generator);
		r2 = random_fl(0.0, 1.0, generator);
		A2 = 2 * a * r1 - a;
		C2 = 2 * r2;

		r1 = random_fl(0.0, 1.0, generator);
		r2 = random_fl(0.0, 1.0, generator);
		A3 = 2 * a * r1 - a;
		C3 = 2 * r2;
		
		D1.c.ligands[0].torsions[i] = abs(C1 * alpha.c.ligands[0].torsions[i] - mutating_candidate.c.ligands[0].torsions[i]);
		D2.c.ligands[0].torsions[i] = abs(C2 * beta.c.ligands[0].torsions[i] - mutating_candidate.c.ligands[0].torsions[i]);
		D3.c.ligands[0].torsions[i] = abs(C3 * delta.c.ligands[0].torsions[i] - mutating_candidate.c.ligands[0].torsions[i]);
		
		X1.c.ligands[0].torsions[i] =  alpha.c.ligands[0].torsions[i] - A1 * D1.c.ligands[0].torsions[i];
		X2.c.ligands[0].torsions[i] =  beta.c.ligands[0].torsions[i] - A2 * D2.c.ligands[0].torsions[i];
		X3.c.ligands[0].torsions[i] =  delta.c.ligands[0].torsions[i] - A3 * D3.c.ligands[0].torsions[i];
		
		movement.c.ligands[0].torsions[i] = (X1.c.ligands[0].torsions[i] + X2.c.ligands[0].torsions[i] + X3.c.ligands[0].torsions[i]) / 3.0;
	}
	/*****************************************************flex**********************************************/
	//for (int i = 0; i < mutating_candidate.c.flex.size();i++)
	for (int k = 0; k < flex_list.size(); k++)
	{
		int i = flex_list[k];
		
		for (int j = 0; j < mutating_candidate.c.flex[i].torsions.size();j++)
		{
			fl A1, A2, A3, C1, C2 ,C3, r1,r2;
			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A1 = 2 * a * r1 - a;
			C1 = 2 * r2;

			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A2 = 2 * a * r1 - a;
			C2 = 2 * r2;

			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A3 = 2 * a * r1 - a;
			C3 = 2 * r2;
			
			D1.c.flex[i].torsions[j] = abs(C1 * alpha.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j]);
			D2.c.flex[i].torsions[j] = abs(C2 * beta.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j]);
			D3.c.flex[i].torsions[j] = abs(C3 * delta.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j]);
			
			X1.c.flex[i].torsions[j] =  alpha.c.flex[i].torsions[j] - A1 * D1.c.flex[i].torsions[j];
			X2.c.flex[i].torsions[j] =  beta.c.flex[i].torsions[j] - A2 * D2.c.flex[i].torsions[j];
			X3.c.flex[i].torsions[j] =  delta.c.flex[i].torsions[j] - A3 * D3.c.flex[i].torsions[j];
			
			movement.c.flex[i].torsions[j] = (X1.c.flex[i].torsions[j] + X2.c.flex[i].torsions[j] + X3.c.flex[i].torsions[j]) / 3.0;
		}
	}
	populations[position] = movement;			 
				 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void gwo_mutate_receptor_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a,
					rng& generator) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion
	
	int changing = 0;
	sz mutable_entities_num = count_mutable_entities(populations[position].c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	
	output_type mutating_candidate = populations[position];
	//output_type difference = mutation(populations, position, F, generator);
	//output_type movement = cal_movement(mutating_candidate, alpha, beta, delta, generator);
	//output_type new_candidate = populations[position];
	//output_type movement;
	/*****************************************************flex**********************************************/
	fl A1, A2, A3, C1, C2 ,C3, r1,r2;
	//conf_type D1, D2, D3, X1, X2, X3;
	output_type D1 = mutating_candidate;
	output_type D2 = mutating_candidate;
	output_type D3 = mutating_candidate;
	output_type X1 = mutating_candidate;
	output_type X2 = mutating_candidate;
	output_type X3 = mutating_candidate;
	output_type movement = mutating_candidate;
	for (int i = 0; i < mutating_candidate.c.flex.size();i++)
	{
		for (int j = 0; j < mutating_candidate.c.flex[i].torsions.size();j++)
		{
			fl A1, A2, A3, C1, C2 ,C3, r1,r2;
			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A1 = 2 * a * r1 - a;
			C1 = 2 * r2;

			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A2 = 2 * a * r1 - a;
			C2 = 2 * r2;

			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A3 = 2 * a * r1 - a;
			C3 = 2 * r2;
			
			D1.c.flex[i].torsions[j] = C1 * alpha.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j];
			D2.c.flex[i].torsions[j] = C2 * beta.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j];
			D3.c.flex[i].torsions[j] = C3 * delta.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j];
			
			X1.c.flex[i].torsions[j] =  alpha.c.flex[i].torsions[j] - A1 * D1.c.flex[i].torsions[j];
			X2.c.flex[i].torsions[j] =  beta.c.flex[i].torsions[j] - A2 * D2.c.flex[i].torsions[j];
			X3.c.flex[i].torsions[j] =  delta.c.flex[i].torsions[j] - A3 * D3.c.flex[i].torsions[j];
			
			movement.c.flex[i].torsions[j] = (X1.c.flex[i].torsions[j] + X2.c.flex[i].torsions[j] + X3.c.flex[i].torsions[j]) / 3;
		}
	}
	populations[position] = movement;
	
}

void gwo_mutate_receptor_sub_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a,rng& generator, std::vector<int> &flex_list){
	
	int changing = 0;
	sz mutable_entities_num = count_mutable_entities(populations[position].c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	
	output_type mutating_candidate = populations[position];
	//output_type difference = mutation(populations, position, F, generator);
	//output_type movement = cal_movement(mutating_candidate, alpha, beta, delta, generator);
	//output_type new_candidate = populations[position];
	//output_type movement;
	/*****************************************************flex**********************************************/
	fl A1, A2, A3, C1, C2 ,C3, r1,r2;
	//conf_type D1, D2, D3, X1, X2, X3;
	output_type D1 = mutating_candidate;
	output_type D2 = mutating_candidate;
	output_type D3 = mutating_candidate;
	output_type X1 = mutating_candidate;
	output_type X2 = mutating_candidate;
	output_type X3 = mutating_candidate;
	output_type movement = mutating_candidate;
	//for (int i = 0; i < mutating_candidate.c.flex.size();i++)
	for (int k = 0; k < flex_list.size(); k++)
	{
		int i = flex_list[k];
		for (int j = 0; j < mutating_candidate.c.flex[i].torsions.size();j++)
		{
			fl A1, A2, A3, C1, C2 ,C3, r1,r2;
			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A1 = 2 * a * r1 - a;
			C1 = 2 * r2;

			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A2 = 2 * a * r1 - a;
			C2 = 2 * r2;

			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A3 = 2 * a * r1 - a;
			C3 = 2 * r2;
			
			D1.c.flex[i].torsions[j] = C1 * alpha.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j];
			D2.c.flex[i].torsions[j] = C2 * beta.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j];
			D3.c.flex[i].torsions[j] = C3 * delta.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j];
			
			X1.c.flex[i].torsions[j] =  alpha.c.flex[i].torsions[j] - A1 * D1.c.flex[i].torsions[j];
			X2.c.flex[i].torsions[j] =  beta.c.flex[i].torsions[j] - A2 * D2.c.flex[i].torsions[j];
			X3.c.flex[i].torsions[j] =  delta.c.flex[i].torsions[j] - A3 * D3.c.flex[i].torsions[j];
			
			movement.c.flex[i].torsions[j] = (X1.c.flex[i].torsions[j] + X2.c.flex[i].torsions[j] + X3.c.flex[i].torsions[j]) / 3;
		}
	}
	populations[position] = movement;

	
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void gwo_mutate_ligand_conf(std::vector<output_type>& populations, int position, output_type alpha, output_type beta, output_type delta, float a,
					rng& generator) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion
	
	int changing = 0;
	sz mutable_entities_num = count_mutable_entities(populations[position].c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	
	output_type mutating_candidate = populations[position];
	//output_type difference = mutation(populations, position, F, generator);
	//output_type movement = cal_movement(mutating_candidate, alpha, beta, delta, generator);
	//output_type new_candidate = populations[position];
	//output_type movement;
	/*********************************position****************************/
	fl A1, A2, A3, C1, C2 ,C3, r1,r2;
	//conf_type D1, D2, D3, X1, X2, X3;
	output_type D1 = mutating_candidate;
	output_type D2 = mutating_candidate;
	output_type D3 = mutating_candidate;
	output_type X1 = mutating_candidate;
	output_type X2 = mutating_candidate;
	output_type X3 = mutating_candidate;
	output_type movement = mutating_candidate;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A1 = 2 * a * r1 - a;
	C1 = 2 * r2;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A2 = 2 * a * r1 - a;
	C2 = 2 * r2;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A3 = 2 * a * r1 - a;
	C3 = 2 * r2;
	
	D1.c.ligands[0].rigid.position = C1 * alpha.c.ligands[0].rigid.position - mutating_candidate.c.ligands[0].rigid.position;
	D2.c.ligands[0].rigid.position = C2 * beta.c.ligands[0].rigid.position - mutating_candidate.c.ligands[0].rigid.position;
	D3.c.ligands[0].rigid.position = C3 * delta.c.ligands[0].rigid.position - mutating_candidate.c.ligands[0].rigid.position;
	D1.c.ligands[0].rigid.position[0] = abs(D1.c.ligands[0].rigid.position[0]);
	D1.c.ligands[0].rigid.position[1] = abs(D1.c.ligands[0].rigid.position[1]);
	D1.c.ligands[0].rigid.position[2] = abs(D1.c.ligands[0].rigid.position[2]);
	D2.c.ligands[0].rigid.position[0] = abs(D1.c.ligands[0].rigid.position[0]);
	D2.c.ligands[0].rigid.position[1] = abs(D1.c.ligands[0].rigid.position[1]);
	D2.c.ligands[0].rigid.position[2] = abs(D1.c.ligands[0].rigid.position[2]);
	D3.c.ligands[0].rigid.position[0] = abs(D1.c.ligands[0].rigid.position[0]);
	D3.c.ligands[0].rigid.position[1] = abs(D1.c.ligands[0].rigid.position[1]);
	D3.c.ligands[0].rigid.position[2] = abs(D1.c.ligands[0].rigid.position[2]);
	
	X1.c.ligands[0].rigid.position = alpha.c.ligands[0].rigid.position - A1 * D1.c.ligands[0].rigid.position;
	X2.c.ligands[0].rigid.position = beta.c.ligands[0].rigid.position - A2 * D2.c.ligands[0].rigid.position;
	X3.c.ligands[0].rigid.position = delta.c.ligands[0].rigid.position - A3 * D3.c.ligands[0].rigid.position;
	
	movement.c.ligands[0].rigid.position = (X1.c.ligands[0].rigid.position + X2.c.ligands[0].rigid.position + X3.c.ligands[0].rigid.position) * (1.0/3.0);
	/********************************************orientation****************************************/
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A1 = 2 * a * r1 - a;
	C1 = 2 * r2;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A2 = 2 * a * r1 - a;
	C2 = 2 * r2;
	
	r1 = random_fl(0.0, 1.0, generator);
	r2 = random_fl(0.0, 1.0, generator);
	A3 = 2 * a * r1 - a;
	C3 = 2 * r2;
	/*
	D1.c.ligands[0].rigid.orientation = abs(angle_to_quaternion(C1 * quaternion_to_angle(alpha.c.ligands[0].rigid.orientation) - quaternion_to_angle( mutating_candidate.c.ligands[0].rigid.orientation)));
	D2.c.ligands[0].rigid.orientation = abs(angle_to_quaternion(C2 * quaternion_to_angle(beta.c.ligands[0].rigid.orientation) - quaternion_to_angle( mutating_candidate.c.ligands[0].rigid.orientation)));
	D3.c.ligands[0].rigid.orientation = abs(angle_to_quaternion(C3 * quaternion_to_angle(delta.c.ligands[0].rigid.orientation) - quaternion_to_angle( mutating_candidate.c.ligands[0].rigid.orientation)));

	X1.c.ligands[0].rigid.orientation = angle_to_quaternion( quaternion_to_angle(alpha.c.ligands[0].rigid.orientation) - A1 * quaternion_to_angle( D1.c.ligands[0].rigid.orientation));
	X2.c.ligands[0].rigid.orientation = angle_to_quaternion( quaternion_to_angle(beta.c.ligands[0].rigid.orientation) - A2 * quaternion_to_angle( D2.c.ligands[0].rigid.orientation));
	X3.c.ligands[0].rigid.orientation = angle_to_quaternion( quaternion_to_angle(delta.c.ligands[0].rigid.orientation) - A3 * quaternion_to_angle( D3.c.ligands[0].rigid.orientation));
	*/
	vec rotation = quaternion_to_angle(alpha.c.ligands[0].rigid.orientation);
	vec D1_rotation, D2_rotation, D3_rotation, X1_rotation, X2_rotation, X3_rotation;
	D1_rotation = quaternion_difference(my_angle_to_quaternion(rotation, C1), mutating_candidate.c.ligands[0].rigid.orientation);
	rotation = quaternion_to_angle(beta.c.ligands[0].rigid.orientation);
	D2_rotation = quaternion_difference(my_angle_to_quaternion(rotation, C2), mutating_candidate.c.ligands[0].rigid.orientation);
	rotation = quaternion_to_angle(delta.c.ligands[0].rigid.orientation);
	D3_rotation = quaternion_difference(my_angle_to_quaternion(rotation, C3), mutating_candidate.c.ligands[0].rigid.orientation);
	
	X1_rotation = quaternion_difference(alpha.c.ligands[0].rigid.orientation, my_angle_to_quaternion(D1_rotation, A1));
	X2_rotation = quaternion_difference(beta.c.ligands[0].rigid.orientation, my_angle_to_quaternion(D2_rotation, A2));
	X3_rotation = quaternion_difference(delta.c.ligands[0].rigid.orientation, my_angle_to_quaternion(D3_rotation, A3));
	
	qt tmp = my_angle_to_quaternion(X1_rotation, 1/3);
	quaternion_increment(tmp, quaternion_to_angle(my_angle_to_quaternion(X2_rotation, 1/3)));
	quaternion_increment(tmp, quaternion_to_angle(my_angle_to_quaternion(X2_rotation, 1/3)));
	
	//movement.c.ligands[0].rigid.orientation = angle_to_quaternion((quaternion_to_angle(X1.c.ligands[0].rigid.orientation) + quaternion_to_angle(X2.c.ligands[0].rigid.orientation) + quaternion_to_angle(X3.c.ligands[0].rigid.orientation))*(1.0 / 3.0));
	movement.c.ligands[0].rigid.orientation = tmp;
	/***********************************************torsion****************************************/
	for (int i = 0; i < mutating_candidate.c.ligands[0].torsions.size();i++)
	{
		fl A1, A2, A3, C1, C2 ,C3, r1,r2;
		r1 = random_fl(0.0, 1.0, generator);
		r2 = random_fl(0.0, 1.0, generator);
		A1 = 2 * a * r1 - a;
		C1 = 2 * r2;

		r1 = random_fl(0.0, 1.0, generator);
		r2 = random_fl(0.0, 1.0, generator);
		A2 = 2 * a * r1 - a;
		C2 = 2 * r2;

		r1 = random_fl(0.0, 1.0, generator);
		r2 = random_fl(0.0, 1.0, generator);
		A3 = 2 * a * r1 - a;
		C3 = 2 * r2;
		
		D1.c.ligands[0].torsions[i] = abs(C1 * alpha.c.ligands[0].torsions[i] - mutating_candidate.c.ligands[0].torsions[i]);
		D2.c.ligands[0].torsions[i] = abs(C2 * beta.c.ligands[0].torsions[i] - mutating_candidate.c.ligands[0].torsions[i]);
		D3.c.ligands[0].torsions[i] = abs(C3 * delta.c.ligands[0].torsions[i] - mutating_candidate.c.ligands[0].torsions[i]);
		
		X1.c.ligands[0].torsions[i] =  alpha.c.ligands[0].torsions[i] - A1 * D1.c.ligands[0].torsions[i];
		X2.c.ligands[0].torsions[i] =  beta.c.ligands[0].torsions[i] - A2 * D2.c.ligands[0].torsions[i];
		X3.c.ligands[0].torsions[i] =  delta.c.ligands[0].torsions[i] - A3 * D3.c.ligands[0].torsions[i];
		
		movement.c.ligands[0].torsions[i] = (X1.c.ligands[0].torsions[i] + X2.c.ligands[0].torsions[i] + X3.c.ligands[0].torsions[i]) / 3.0;
	}
	/*****************************************************flex**********************************************/
	/*
	for (int i = 0; i < mutating_candidate.c.flex.size();i++)
	{
		for (int j = 0; j < mutating_candidate.c.flex[i].torsions.size();j++)
		{
			fl A1, A2, A3, C1, C2 ,C3, r1,r2;
			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A1 = 2 * a * r1 - a;
			C1 = 2 * r2;

			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A2 = 2 * a * r1 - a;
			C2 = 2 * r2;

			r1 = random_fl(0.0, 1.0, generator);
			r2 = random_fl(0.0, 1.0, generator);
			A3 = 2 * a * r1 - a;
			C3 = 2 * r2;
			
			D1.c.flex[i].torsions[j] = abs(C1 * alpha.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j]);
			D2.c.flex[i].torsions[j] = abs(C2 * beta.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j]);
			D3.c.flex[i].torsions[j] = abs(C3 * delta.c.flex[i].torsions[j] - mutating_candidate.c.flex[i].torsions[j]);
			
			X1.c.flex[i].torsions[j] =  alpha.c.flex[i].torsions[j] - A1 * D1.c.flex[i].torsions[j];
			X2.c.flex[i].torsions[j] =  beta.c.flex[i].torsions[j] - A2 * D2.c.flex[i].torsions[j];
			X3.c.flex[i].torsions[j] =  delta.c.flex[i].torsions[j] - A3 * D3.c.flex[i].torsions[j];
			
			movement.c.flex[i].torsions[j] = (X1.c.flex[i].torsions[j] + X2.c.flex[i].torsions[j] + X3.c.flex[i].torsions[j]) / 3.0;
		}
	}*/
	populations[position] = movement;
	
}
