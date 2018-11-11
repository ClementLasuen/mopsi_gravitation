#include"mopsi.h"



int main()
{
    // Initialisation

    // Soleil
    FVector<double,3> q_soleil = {0.0,0.0,0.0};
    FVector<double,3> p_soleil = {0.0,0.0,0.0};

    // Jupiter

    FVector<double,3> q_jupiter = {-3.5024 , -3.8170 , -1.5508};
    FVector<double,3> p_jupiter = {0.0056543 , -0.0041249, -0.0019059};
    p_jupiter = m[1]*p_jupiter;

    // Saturne

    FVector<double,3> q_saturne = {9.1,-3.0,-1.6};
    FVector<double,3> p_saturne = {0.0017, 0.0048,0.0019};
    p_saturne = m[2]*p_saturne;

    // Uranus
    FVector<double,3> q_uranus = {8.3101, -16.290, -7.2521};
    FVector<double,3> p_uranus = {0.0035418, 0.0013710, 0.00055029};
    p_uranus = m[3]*p_uranus;
    // Neptune

    FVector<double,3> q_neptune = {-15.539, -25.223, -3.1902};
    FVector<double,3> p_neptune = {0.0028893, -0.0017070, -0.0013650};
    p_neptune = m[4]*p_neptune;

    FVector<FVector<double,3>,nb_planetes> p;// = {p_soleil,p_jupiter, p_saturne, p_uranus,p_neptune};
    FVector<FVector<double,3>,nb_planetes> q;// = {q_soleil,q_jupiter, q_saturne, q_uranus, q_neptune};

    q[0]=q_soleil;
    q[1]=q_jupiter;
    q[2]=q_saturne;
    q[3]=q_uranus;
    q[4]=q_neptune;

    p[0]=p_soleil;
    p[1]=p_jupiter;
    p[2]=p_saturne;
    p[3]=p_uranus;
    p[4]=p_neptune;

    FVector<FVector<double,3>,nb_planetes>* resu = euler_implicite(q,p);
    /*openWindow(500,500);
    // Je projete sur le plan 0xy pour voir si c'est potable

    for(int i =0; i<nb_iterations;i++){
        fillCircle(resu[i][0][0]+200,resu[i][0][1]+200,2,YELLOW);
        fillCircle(resu[i][1][0]+200,resu[i][1][1]+200,2,RED);
        cout << resu[i][1] << endl;
    }

    endGraphics();*/

    ofstream valeur_H("C:/Users/Utilisateur/Downloads/Tp2_Initial/Tennis/test_h.txt");
    if(valeur_H){
        for(int i =0; i<nb_iterations;i++){
            valeur_H << H(resu[i],resu[i+nb_iterations]) << endl;
        }
    }
    else cout << "pb ouverture" << endl;

	return 0;
}
