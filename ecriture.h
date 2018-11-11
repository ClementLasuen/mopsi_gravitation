#pragma once

#include <string>
#include <fstream>

using namespace std;

void ecriture(string file_name, FVector<FVector<double, _Tp2>, nb_planetes> *trajectory);
/* ecrit un fichier avec pour premiere ligne le nombre d'iterations puis le nombre de planetes
   En deuxieme ligne il y a la ligne de toute les planetes
   Vient ensuite n_iterations lignes avec les coordonnes de toutes les planetes les unes a la suite des autres (donc nb_iterations lignes de 3*nb_planete nombres)
   Vient ensuite n_iterations lignes avec les VITESSES de toutes les planetes les unes a la suite des autres (donc nb_iterations lignes de 3*nb_planete nombres)
*/
