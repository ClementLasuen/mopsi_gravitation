#pragma once
#include<vector>
#include<iostream>
using namespace std;
#include<cmath>
#include <Imagine/Graphics.h>
using namespace Imagine;

/*

  Unites : UA, masse soleil, jour terrestre

*/

const int nb_planetes =2;
const double h = 1.;
const int nb_iterations = 2*pow(10.,4.);
const double m [nb_planetes] = {1.00000597682, 0.00095};
const double G = 2.95912208286*pow(10.,-4.);

// ------------------------------------------- Fonctions pratiques ----------------------------------------

double norme(FVector<double,3> v);

// interaction(q)[i] :
// --> prend en entré la position q (vecteur n_planetes*3) de toute les planètes
// --> renvoie la résultante des forces exercées sur la planète i par les autres planètes
FVector<FVector<double,3>,nb_planetes> interaction(FVector<FVector<double,3>,nb_planetes> q0);

//---------------------------------- Methodes d'integration -------------------------------------------

// Renvoie une "distance" entre les positions q0 et q1
// Utilisée dans la calcul du point fixe pour euler implicite
double ecart(FVector<FVector<double,3>,nb_planetes> q0, FVector<FVector<double,3>,nb_planetes> q1);

// Euler explicite

FVector<FVector<double,3>,nb_planetes>* euler_explicite(FVector<FVector<double,3>,nb_planetes> q0, FVector<FVector<double,3>,nb_planetes> p0);


// Euler implicite

// Calcule le point fixe à l'itération n
// Ne prend pas en argument la masse car le tableau des masses est donné et CONSTANT
FVector<FVector<double,3>,nb_planetes>* pf_euler_implicite(FVector<FVector<double,3>,nb_planetes> qn, FVector<FVector<double,3>,nb_planetes> pn );

// Renvoie l'ensemble des positions et des quantités de mouvement
FVector<FVector<double,3>,nb_planetes>* euler_implicite(FVector<FVector<double,3>,nb_planetes> q0, FVector<FVector<double,3>,nb_planetes> p0);

// Euler symplectique

FVector<FVector<double,3>,nb_planetes>* pf_euler_symplectique(FVector<FVector<double,3>,nb_planetes> qn, FVector<FVector<double,3>,nb_planetes> pn );

FVector<FVector<double,3>,nb_planetes>* euler_symplectique(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p0);
