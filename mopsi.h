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
const double h = 10;
const int nb_iterations = 200000;
const double m [nb_planetes] = {1.00000597682,0.0009548, 0.00029, 0.0000437, 0.0000518};
const double G =2.959122*0.0001;

// ------------------------------------------- Fonctions pratiques ----------------------------------------

double norme(FVector<double,3> v);

// interaction(q)[i] renvoie la résultante des forces exercées sur la planète i par les autres planètes
// dont la position est donnée par q
FVector<FVector<double,3>,nb_planetes> interaction(FVector<FVector<double,3>,nb_planetes> q);

// ----------------------------------------- Hamiltonien----------------------------------------------

double H(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p);
double H_modifie_ES(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p);
double H_modifie_V(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p);

FVector<FVector<double, 3*nb_planetes>,3*nb_planetes >  Hessienne(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p);
//---------------------------------- Methodes d'integration -------------------------------------------

// Renvoie une "distance" entre les positions q0 et q1
// Utilisée dans la calcul du point fixe pour euler implicite
double ecart(FVector<FVector<double,3>,nb_planetes> q0, FVector<FVector<double,3>,nb_planetes> q1);

// Euler explicite

FVector<FVector<double,3>,nb_planetes>* euler_explicite(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p0, bool ecriture = true);
/* ecrit un fichier avec pour premiere ligne le nombre d'iterations puis le nombre de planetes
   En deuxieme ligne il y a la masse de toute les planetes
   Vient ensuite n_iterations lignes avec les coordonnes puis les vitesses de toutes les planetes (donc nb_iterations lignes de 6*nb_planete nombres)
*/

// Euler implicite

// Calcule le point fixe à l'itération n
// Ne prend pas en argument la masse car le tableau des masses est donné et CONSTANT
FVector<FVector<double,3>,nb_planetes>* pf_euler_implicite(FVector<FVector<double,3>,nb_planetes> qn, FVector<FVector<double,3>,nb_planetes> pn);
/* ecrit un fichier avec pour premiere ligne le nombre d'iterations puis le nombre de planetes
   En deuxieme ligne il y a la masse de toute les planetes
   Vient ensuite n_iterations lignes avec les coordonnes puis les vitesses de toutes les planetes (donc nb_iterations lignes de 6*nb_planete nombres)
*/

// Renvoie l'ensemble des positions et des quantités de mouvement
FVector<FVector<double,3>,nb_planetes>* euler_implicite(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p0, bool ecriture = true);

// Euler symplectique

FVector<FVector<double,3>,nb_planetes>* pf_euler_symplectique(FVector<FVector<double,3>,nb_planetes> qn, FVector<FVector<double,3>,nb_planetes> pn );

FVector<FVector<double,3>,nb_planetes>* euler_symplectique(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p0, bool ecriture = true);

// Verlet

void changement_variables(FVector<FVector<double,3>,nb_planetes> &p);
FVector<FVector<double,3>,nb_planetes>* verlet(FVector<FVector<double,3>,nb_planetes> q0, FVector<FVector<double,3>,nb_planetes> p0, bool ecriture = true);
