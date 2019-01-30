#include "mopsi.h"

#include "stdlib.h"

int main(){

    double h;
    string methode;

    cout << "Quel est le pas de temps ? " << endl;
    cin >> h ;


    cout << "Quelle methode numerique a utiliser ?" << endl;
    //cout << "explicite, implicite, symplectique1,  symplectique2, Verlet" << endl;
    cout << "EE, ES, EI ou V" << endl;
    cin >> methode ;
    //FVector<FVector<double,3>,nb_planetes>* resu;
    bool choix_methode = true;
    if(methode == "ES") euler_symplectique(h);
    else if(methode == "EE") euler_explicite(h);
    else if(methode == "EI") euler_implicite(h);
    else if(methode == "V")  verlet(h);
    else{
        choix_methode=false;
        cout << "choisis une methode" << endl;
    }

    /*
    Py_Initialize();
    Py_Finalize();
    */

    //system(".\python3.4 test_hess.py");

    return 0;
}
