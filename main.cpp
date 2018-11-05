#include"mopsi.h"
int main()
{
    // Initialisation, on teste juste avec Jupiter
    FVector<double,3> q_soleil = {0.0,0.0,0.0};
    FVector<double,3> q_jupiter = {-3.5023653,-3.8169847,-1.5507963};
    FVector<double,3> p_soleil = {0.0,0.0,0.0};
    FVector<double,3> p_jupiter = {0.00565429, -0.00412490,-0.00190589};
    FVector<FVector<double,3>,nb_planetes> p = {p_soleil,p_jupiter};
    FVector<FVector<double,3>,nb_planetes> q = {q_soleil,q_jupiter};

    FVector<FVector<double,3>,nb_planetes>* resu = euler_explicite(q,p);
    /*
    openWindow(500,500);
    // Je projete sur le plan 0xy pour voir si c'est potable
    for(int i =0; i<nb_iterations;i++){
        fillCircle(resu[i][0][0]+200,resu[i][0][1]+200,2,YELLOW);
        fillCircle(resu[i][1][0]+200,resu[i][1][1]+200,2,RED);
        //cout << resu[i][0][1] << endl;
    }
    endGraphics();
    */
    for(int i=0;i<nb_iterations/1000;i++)
        cout<<resu[i*1000][1]<<" "<<resu[nb_iterations+i*1000][1]<<endl;

	return 0;
}
