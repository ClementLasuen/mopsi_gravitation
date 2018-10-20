#pragma once
#include<vector>
#include<iostream>
using namespace std;
#include<cmath>

const int nb_planetes =2;
const float h = 0.001;
const int nb_iterations = 100;
const float m [nb_planetes] = {1.0,1.0};
const int G =1;


float norme(float* v);
std::vector<float*>* point_fixe(std::vector<float*>* qn, std::vector<float*>* pn );
std::vector<float *> interaction(std::vector<float *> q);
float ecart(std::vector<float *> q0, std::vector<float *> q1);
