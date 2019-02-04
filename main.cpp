#include "mopsi.h"

#include "stdlib.h"

int main(){

    double h;
    string methode;

    cout << "Quel est le pas de temps ? " << endl;
    cin >> h ;


    cout << "Quelle methode numerique utiliser ? " << endl;
    cout << " - EE : Euler explicite " << endl;
    cout << " - EI : Euler implicite " << endl;
    cout << " - ES : Euler symplectique" << endl;
    cout << " - V : Verlet" << endl;
    cout << "EE, ES, EI ou V ?" << endl;
    cin >> methode ;

    bool choix_methode = true;
    if(methode == "ES") euler_symplectique(h);
    else if(methode == "EE") euler_explicite(h);
    else if(methode == "EI") euler_implicite(h);
    else if(methode == "V")  verlet(h);
    else{
        choix_methode=false;
        cout << "recommencez et choisissez une methode" << endl;
    }
        if(choix_methode)
            system("cd ../mopsi_gravitation;python -c 'from func import *;nb_iterations, nb_planetes,positions,vitesses = get_data();enhanced_plot(nb_iterations, nb_planetes,positions,vitesses)'");

    return 0;

}
