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
const float m [nb_planetes] = {1.0, 0.00095};
const float G =2.959*0.0001;

// ------------------------------------------- Fonctions pratiques ----------------------------------------

float norme(FVector<float,3> v);

// interaction(q)[i] renvoie la résultante des forces exercées sur la planète i par les autres planètes
// dont la position est donnée par q
FVector<FVector<float,3>,nb_planetes> interaction(FVector<FVector<float,3>,nb_planetes> q);

//---------------------------------- Methodes d'integration -------------------------------------------

// Renvoie une "distance" entre les positions q0 et q1
// Utilisée dans la calcul du point fixe pour euler implicite
float ecart(FVector<FVector<float,3>,nb_planetes> q0, FVector<FVector<float,3>,nb_planetes> q1);

// Euler explicite

FVector<FVector<float,3>,nb_planetes>* euler_explicite(FVector<FVector<float,3>,nb_planetes> q, FVector<FVector<float,3>,nb_planetes> p0, float masse);


// Euler implicite

// Calcule le point fixe à l'itération n
// Ne prend pas en argument la masse car le tableau des masses est donné et CONSTANT
FVector<FVector<float,3>,nb_planetes>* pf_euler_implicite(FVector<FVector<float,3>,nb_planetes> qn, FVector<FVector<float,3>,nb_planetes> pn );

// Renvoie l'ensemble des positions et des quantités de mouvement
FVector<FVector<float,3>,nb_planetes>* euler_implicite(FVector<FVector<float,3>,nb_planetes> q, FVector<FVector<float,3>,nb_planetes> p0);

// Euler symplectique

FVector<FVector<float,3>,nb_planetes>* pf_euler_symplectique(FVector<FVector<float,3>,nb_planetes> qn, FVector<FVector<float,3>,nb_planetes> pn );

FVector<FVector<float,3>,nb_planetes>* euler_symplectique(FVector<FVector<float,3>,nb_planetes> q, FVector<FVector<float,3>,nb_planetes> p0);
