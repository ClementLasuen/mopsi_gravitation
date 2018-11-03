#pragma once
#include<vector>
#include<iostream>
using namespace std;
#include<cmath>
#include <Imagine/Graphics.h>
using namespace Imagine;

const int nb_planetes =2;
const float h = 0.001;
const int nb_iterations = 100;
const float m [nb_planetes] = {1.0,1.0};
const int G =1;



FVector<FVector<float,3>,nb_planetes> f();

// ------------------------------------------- Fonctions pratiques ----------------------------------------

float norme(FVector<float,3> v);
FVector<FVector<float,3>,nb_planetes> interaction(FVector<FVector<float,3>,nb_planetes> q);


FVector<FVector<float,3>,nb_planetes>* point_fixe(FVector<FVector<float,3>,nb_planetes> qn, FVector<FVector<float,3>,nb_planetes> pn );

float ecart(FVector<FVector<float,3>,nb_planetes> q0, FVector<FVector<float,3>,nb_planetes> q1);

//---------------------------------- Methodes d'integration -------------------------------------------

FVector<FVector<float,3>,nb_planetes>* euler_implicite(FVector<FVector<float,3>,nb_planetes> q, FVector<FVector<float,3>,nb_planetes> p0, float masse);
