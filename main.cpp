#include"mopsi.h"
int main()
{
    // Initialisation, on teste juste avec Jupiter
    FVector<float,3> q_soleil = {0.0,0.0,0.0};
    FVector<float,3> q_jupiter = {-3.5,-3.81,-1.55};
    FVector<float,3> p_soleil = {0.0,0.0,0.0};
    FVector<float,3> p_jupiter = {0.0057, -0.004,-0.0019};
    FVector<FVector<float,3>,nb_planetes> p = {p_soleil,p_jupiter};
    FVector<FVector<float,3>,nb_planetes> q = {q_soleil,q_jupiter};

    FVector<FVector<float,3>,nb_planetes>* resu = euler_implicite(q,p);
    openWindow(500,500);
    // Je projete sur le plan 0xy pour voir si c'est potable
    for(int i =0; i<nb_iterations;i++){
        fillCircle(resu[i][0][0]+200,resu[i][0][1]+200,2,YELLOW);
        fillCircle(resu[i][1][0]+200,resu[i][1][1]+200,2,RED);
        //cout << resu[i][0][1] << endl;
    }
    endGraphics();
	return 0;
}
