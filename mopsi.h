#pragma once
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
using namespace std;
#include<cmath>
#include <Imagine/Graphics.h>
using namespace Imagine;

const int nb_planetes =5;
const int nb_iterations = 20000;
const double m [nb_planetes] = {1.00000597682,0.0009548, 0.00029, 0.0000437, 0.0000518};
const double G =2.959122*0.0001;

// -------------------------------------------- Initialisation -------------------------------------------

// Initialisation

// Soleil
const FVector<double,3> q_soleil = {0.0,0.0,0.0};
const FVector<double,3> p_soleil = {0.0,0.0,0.0};

// Jupiter

const FVector<double,3> q_jupiter = {-3.5024 , -3.8170 , -1.5508};
const FVector<double,3> p_jupiter1 = {0.0056543 , -0.0041249, -0.0019059};
const FVector<double,3> p_jupiter = m[1]*p_jupiter1;

// Saturne

const FVector<double,3> q_saturne = {9.1,-3.0,-1.6};
const FVector<double,3> p_saturne1 = {0.0017, 0.0048,0.0019};
const FVector<double,3> p_saturne = m[2]*p_saturne1;

// Uranus
const FVector<double,3> q_uranus = {8.3101, -16.290, -7.2521};
const FVector<double,3> p_uranus1 = {0.0035418, 0.0013710, 0.00055029};
const FVector<double,3> p_uranus = m[3]*p_uranus1;
// Neptune

const FVector<double,3> q_neptune = {-15.539, -25.223, -3.1902};
const FVector<double,3> p_neptune1 = {0.0028893, -0.0017070, -0.0013650};
const FVector<double,3> p_neptune = m[4]*p_neptune1;


// ------------------------------------------- Fonctions pratiques ----------------------------------------

double norme(FVector<double,3> v);

// interaction(q)[i] renvoie la résultante des forces exercées sur la planète i par les autres planètes
// dont la position est donnée par q
FVector<FVector<double,3>,nb_planetes> interaction(FVector<FVector<double,3>,nb_planetes> q);

// ----------------------------------------- Hamiltonien----------------------------------------------

double H(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p);
double H_modifie_ES(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p);
double H_modifie_V(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p);

FVector<FVector<double, 3*nb_planetes>,3*nb_planetes >  Hessienne(FVector<FVector<double,3>,nb_planetes> q);

// Calcul de la Hessienne par differences finies
FVector<FVector<double, 3*nb_planetes>,3*nb_planetes > Hessienne2(FVector<FVector<double, 3>, nb_planetes> q);

//---------------------------------- Methodes d'integration -------------------------------------------

// Renvoie une "distance" entre les positions q0 et q1
// Utilisée dans la calcul du point fixe pour euler implicite
double ecart(FVector<FVector<double,3>,nb_planetes> q0, FVector<FVector<double,3>,nb_planetes> q1);

// Euler explicite

void euler_explicite(double h, bool ecriture = true);
/* ecrit un fichier avec pour premiere ligne le nombre d'iterations puis le nombre de planetes
   En deuxieme ligne il y a la masse de toute les planetes
   Vient ensuite n_iterations lignes avec les coordonnes puis les vitesses de toutes les planetes (donc nb_iterations lignes de 6*nb_planete nombres)
*/

// Euler implicite

// Calcule le point fixe à l'itération n
// Ne prend pas en argument la masse car le tableau des masses est donné et CONSTANT
FVector<FVector<double,3>,nb_planetes>* pf_euler_implicite(double h,FVector<FVector<double,3>,nb_planetes> qn, FVector<FVector<double,3>,nb_planetes> pn);
/* ecrit un fichier avec pour premiere ligne le nombre d'iterations puis le nombre de planetes
   En deuxieme ligne il y a la masse de toute les planetes
   Vient ensuite n_iterations lignes avec les coordonnes puis les vitesses de toutes les planetes (donc nb_iterations lignes de 6*nb_planete nombres)
*/

// Renvoie l'ensemble des positions et des quantités de mouvement
void euler_implicite(double h, bool ecriture = true);

// Euler symplectique

void euler_symplectique(double h, bool ecriture = true);


// Verlet

void changement_variables(FVector<FVector<double,3>,nb_planetes> &p);

void changement_variables_inverse(FVector<FVector<double,3>,nb_planetes> &p);

void verlet(double h, bool ecriture = true);

// Differences finies
// Le potentiel est noté V

double V (FVector<FVector<double, 3>, nb_planetes> q);
FVector<FVector<double,3>,nb_planetes> dV (FVector<FVector<double,3>,nb_planetes> q);
FVector<FVector<double, 3*nb_planetes>,3*nb_planetes > d2V (FVector<FVector<double,3>,nb_planetes> q);
