#include"mopsi.h"

// --------------------------------Fonctions pratiques---------------------------------

float norme(FVector<float, 3> v){
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

FVector<FVector<float, 3>, nb_planetes> interaction(FVector<FVector<float, 3>, nb_planetes> q){
    FVector<FVector<float, 3>, nb_planetes> resu;
    for(int i = 0;i<nb_planetes;i++){
        for(int j = 0; j<nb_planetes;j++){
            if(j!=i){
                resu[i] += (q[i] - q[j])*G* m[j] /(norme(q[i] - q[j])*norme(q[i] - q[j])*norme(q[i] - q[j]));
            }
        }
    }
    return resu;
}

// ------------------------------------------Methodes d'integration------------------------------

float ecart(FVector<FVector<float, 3>, nb_planetes> q0, FVector<FVector<float, 3>, nb_planetes> q1){
    float resu = 0.0;
    for(int i=0;i<nb_planetes;i++){
        resu += norme(q0[i] - q1[i]);
    }
    return resu;
}


// Euler explicite

FVector<FVector<float,3>,nb_planetes>* euler_explicite(FVector<FVector<float,3>,nb_planetes> q0, FVector<FVector<float,3>,nb_planetes> p0, float masse){

    FVector<FVector<float,3>,nb_planetes>* resu = new FVector<FVector<float,3>,nb_planetes> [2*nb_iterations];
    resu[0]=q0;
    resu[nb_iterations]=p0;

    for(int i =1;i<nb_iterations;i++){
        FVector<FVector<float,3>,nb_planetes>* resu_ = new FVector<FVector<float,3>,nb_planetes> [2];
        for(int j=0;j<nb_planetes;j++){
            resu_[0][j] = resu[i-1][j];
            resu_[1][j] = resu[nb_iterations+i-1][j];
        }
        FVector<FVector<float,3>,nb_planetes> p_point = interaction(resu_[0]);
        for(int j=0;j<nb_planetes;j++){
            resu_[1][j] = resu_[1][j] + h*(interaction(resu_[0])[j]);
            resu_[0][j] = resu_[0][j] + h*resu_[1][j];
        }
        resu[i] = resu_[0];
        resu[nb_iterations+i] = resu_[1];
    }
    return resu;
}

// Euler imlicite


FVector<FVector<float,3>,nb_planetes>* pf_euler_implicite(FVector<FVector<float,3>,nb_planetes> qn,FVector<FVector<float,3>,nb_planetes> pn, float masse ){
    float epsilon = 0.0001; // precision comment choisir epsilon ?
    int compteur =0;
    FVector<FVector<float,3>,nb_planetes>* resu = new FVector<FVector<float,3>,nb_planetes> [2];
    FVector<FVector<float,3>,nb_planetes> q0,q1,p0,p1;
    for(int i =0;i<nb_planetes;i++){  // je pense qu'on peut eviter la boucle for et juste ajouter les vecteurs
            q0[i] = qn[i];
            p0[i] = pn[i];
            q1[i] = qn[i] + h*p0[i]/masse;
            p1[i] = pn[i] + h*interaction(q0)[i];
    }
    while( compteur<1000 && max(ecart(q0,q1),ecart(p0,p1))> epsilon ){
        compteur+=1;
        for(int i =0;i<nb_planetes;i++){ // idem qu'au dessus
                q0[i] = q1[i];
                p0[i] = p1[i];
                q1[i] = qn[i] + h*p0[i]/masse;
                p1[i] = pn[i] + h*interaction(q0)[i];
        }
    }
    resu[0] = q1;
    resu[1] = p1;
    return resu;
}

FVector<FVector<float,3>,nb_planetes>* euler_implicite(FVector<FVector<float,3>,nb_planetes> q0, FVector<FVector<float,3>,nb_planetes> p0, float masse){

    FVector<FVector<float,3>,nb_planetes>* resu = new FVector<FVector<float,3>,nb_planetes> [2*nb_iterations];
    resu[0]=q0;
    resu[nb_iterations]=p0;
    for(int i =1;i<nb_iterations;i++){
        FVector<FVector<float,3>,nb_planetes>* resu_ = new FVector<FVector<float,3>,nb_planetes> [2];
        resu[i] = resu_[0];
        resu[nb_iterations+i] = resu_[1];
    }
    return resu;
}

// Euler symplectique (implicite sur q, explicite sur p)

FVector<FVector<float,3>,nb_planetes>* pf_euler_symplectique(FVector<FVector<float,3>,nb_planetes> qn,FVector<FVector<float,3>,nb_planetes> pn, float masse ){
    float epsilon = 0.0001; // precision Comment choisir epsilon ?
    int compteur =0;
    FVector<FVector<float,3>,nb_planetes>* resu = new FVector<FVector<float,3>,nb_planetes> [2];
    FVector<FVector<float,3>,nb_planetes> q0,q1,p0,p1;
    for(int i =0;i<nb_planetes;i++){
            q0[i] = qn[i];
            p0[i] = pn[i];
            q1[i] = qn[i] + h*pn[i]/masse;
            p1[i] = pn[i] + h*interaction(q0)[i];
    }
    while( compteur<1000 && max(ecart(q0,q1),ecart(p0,p1))> epsilon ){
        compteur+=1;
        for(int i =0;i<nb_planetes;i++){
                q0[i] = q1[i];
                //p0[i] = p1[i];
                q1[i] = qn[i] + h*pn[i]/masse;
                p1[i] = pn[i] + h*interaction(q0)[i];
        }
    }
    resu[0] = q1;
    resu[1] = p1;
    return resu;
}

FVector<FVector<float,3>,nb_planetes>* euler_symplectique(FVector<FVector<float,3>,nb_planetes> q0, FVector<FVector<float,3>,nb_planetes> p0, float masse){

    FVector<FVector<float,3>,nb_planetes>* resu = new FVector<FVector<float,3>,nb_planetes> [2*nb_iterations];
    resu[0]=q0;
    resu[nb_iterations]=p0;
    for(int i =1;i<nb_iterations;i++){
        FVector<FVector<float,3>,nb_planetes>* resu_ = new FVector<FVector<float,3>,nb_planetes> [2];
        resu_ = pf_euler_implicite(resu[i-1],resu[nb_iterations + i-1], masse);
        resu[i] = resu_[0];
        resu[nb_iterations+i] = resu_[1];
    }
    return resu;
}
