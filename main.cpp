
#include "mopsi.h"

#include "stdlib.h"

int main(){

    double h;
    string methode;

    cout << "Quel est le pas de temps ? " << endl;
    cin >> h ;


    cout << "Quelle methode numerique a utiliser ?" << endl;
    //cout << "explicite, implicite, symplectique1,  symplectique2, Verlet" << endl;
    cout << "ES ou V" << endl;
    cin >> methode ;
    //FVector<FVector<double,3>,nb_planetes>* resu;
    bool choix_methode = true;
    if(methode == "ES") euler_symplectique(h);
    else if(methode == "V")  verlet(h);
    else{
        choix_methode=false;
        cout << "choisis une methode" << endl;
    }


    //system(".\python3.4 test_hess.py");

    return 0;
}
