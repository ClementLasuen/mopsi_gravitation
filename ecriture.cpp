#include "ecriture.h"
#include "mopsi.h"

void ecriture(string file_name, FVector<FVector<double, 3>, nb_planetes>* trajectory){ // rajouter une variable pour savoirsi on reecrit du debut du fichier ou non
    ofstream fichier(file_name.c_str(), ios::out|ios::app); // On va ecrire a la fin du fichier
    if (fichier){
        fichier << nb_iterations <<" "<<nb_planetes<<endl; // On ecrit les donnees
        for(int j=0;j<nb_planetes;j++) // On ecrit la ligne avec les masses
            fichier << m[j]<<" ";
        fichier <<endl;
        for(int i=0;i<nb_iterations;i++){
            FVector<FVector<double, 3>, nb_planetes> coordonnees = trajectory[i];
            for(int j=0;j<nb_planetes;j++)
                fichier << coordonnees[j][0] << " " << coordonnes[j][1] << " " << coordonnes[j][2] <<endl;; // On ecrit les trajectoires uniquement
        }
        for(int i=0;i<nb_iterations;i++){
            FVector<FVector<double, 3>, nb_planetes> vitesses = trajectory[nb_iterations+i]/m;
            for(int j=0;j<nb_planetes;j++)
                fichier << vitesses[j][0] << " " << vitesses[j][1] << " " << vitesses[j][2] <<endl; // On ecrit les vitesses uniquement
        }
    }
    else{
        cerr<<"Impossible d'ouvrir le fichier"<<endl;
    }
}
