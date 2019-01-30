#include "mopsi.h"
#include</home/quentin/anaconda3/include/python3.5m/Python.h>


int main()
{
        double h;
    string methode;

    cout << "Quel est le pas de temps ? " << endl;
    cin >> h ;


    cout << "Quelle methode numerique a utiliser ?" << endl;
    cout << "explicite, implicite, symplectique1,  symplectique2, Verlet" << endl;
    cin >> methode ;
    //FVector<FVector<double,3>,nb_planetes>* resu;
    bool choix_methode = true;
    if(methode == "explicite")  euler_explicite(h);
    /*else if(methode == "implicite")  resu = euler_implicite(h);
    else if(methode == "symplectique1")  resu = euler_symplectique(h);*/
    else if(methode == "symplectique2") euler_symplectique_sans_pf(h);
    else if(methode == "Verlet")  verlet(h);
    else{
        choix_methode=false;
        cout << "choisis une methode" << endl;
    }

    /*
    Py_Initialize();
    Py_Finalize();
    */

    return 0;
}
