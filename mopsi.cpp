#include"mopsi.h"

// --------------------------------Fonctions pratiques---------------------------------

double norme(FVector<double, 3> v){
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

FVector<FVector<double, 3>, nb_planetes> interaction(FVector<FVector<double, 3>, nb_planetes> q0){
    FVector<FVector<double, 3>, nb_planetes> force;
    for(int i = 0;i<nb_planetes;i++){
        force[i] = {0.0,0.0,0.0};
        for(int j = 0; j<nb_planetes;j++){
            if(j!=i){
                force[i] += -(q0[i] - q0[j])*G*m[j]*m[i]/pow(norme(q0[i] - q0[j]),3.);
            }
        }
    }
    return force;
}

// ------------------------------------------Methodes d'integration------------------------------

double ecart(FVector<FVector<double, 3>, nb_planetes> q0, FVector<FVector<double, 3>, nb_planetes> q1){
    double resu = 0.0;
    for(int i=0;i<nb_planetes;i++){
        resu += norme(q0[i] - q1[i]);
    }
    return resu;
}


// Euler explicite

FVector<FVector<double,3>,nb_planetes>* euler_explicite(FVector<FVector<double,3>,nb_planetes> q0, FVector<FVector<double,3>,nb_planetes> p0){

    FVector<FVector<double,3>,nb_planetes>* resu = new FVector<FVector<double,3>,nb_planetes> [2*nb_iterations];
    resu[0]=q0;
    resu[nb_iterations]=p0;

    for(int i=1;i<nb_iterations;i++){
        FVector<FVector<double,3>,nb_planetes>* resu_ = new FVector<FVector<double,3>,nb_planetes> [2];
        for(int j=0;j<nb_planetes;j++){
            resu_[0][j] = resu[i-1][j];
            resu_[1][j] = resu[nb_iterations+i-1][j];
        }
        FVector<FVector<double,3>,nb_planetes> p_point = interaction(resu_[0]);
        for(int j=0;j<nb_planetes;j++){
            resu_[0][j] = resu_[0][j] + h*resu_[1][j]/m[j];
            resu_[1][j] = resu_[1][j] + h*p_point[j];
        }
        resu[i] = resu_[0];
        resu[nb_iterations+i] = resu_[1];
    }
    return resu;
}

// Euler imlicite


FVector<FVector<double,3>,nb_planetes>* pf_euler_implicite(FVector<FVector<double,3>,nb_planetes> qn,FVector<FVector<double,3>,nb_planetes> pn ){
    double epsilon = 0.00001; // precision comment choisir epsilon ?
    int compteur =0;
    FVector<FVector<double,3>,nb_planetes>* resu = new FVector<FVector<double,3>,nb_planetes> [2];
    FVector<FVector<double,3>,nb_planetes> q0,q1,p0,p1;
    q0 = qn;
    p0 = pn;
    FVector<FVector<double, 3>, nb_planetes> force = interaction(q0);
    for(int i =0;i<nb_planetes;i++){  // je pense qu'on peut eviter la boucle for et juste ajouter les vecteurs
            q1[i] = q0[i] + h*p0[i]/m[i];
            p1[i] = p0[i] + h*force[i]; // A CETTE ETAPE Q0 N'EST PAS TOTALEMENT CALCULEEEEE
    }
    while( compteur<10000 && max(ecart(q0,q1),ecart(p0,p1))> epsilon ){

        compteur+=1;
        q0=q1;
        p0=p1;
        force = interaction(q0);
        for(int i =0;i<nb_planetes;i++){
                q1[i] = qn[i] + h*p0[i]/m[i];
                p1[i] = pn[i] + h*force[i];
        }
    }
    resu[0] = q1;
    resu[1] = p1;
    return resu;
}

FVector<FVector<double,3>,nb_planetes>* euler_implicite(FVector<FVector<double,3>,nb_planetes> q0, FVector<FVector<double,3>,nb_planetes> p0){

    FVector<FVector<double,3>,nb_planetes>* resu = new FVector<FVector<double,3>,nb_planetes> [2*nb_iterations];
    // resu[i] donne les positions de toutes les planètes à l'iteration i;
    // resu[i + nb_iterations] donne les quantites de mouvement
    resu[0]=q0;
    resu[nb_iterations]=p0;
    FVector<FVector<double,3>,nb_planetes>* point_fixe = new FVector<FVector<double,3>,nb_planetes> [2];
    for(int i =1;i<nb_iterations;i++){
        point_fixe = pf_euler_implicite(resu[i-1],resu[i-1 + nb_iterations]);
        resu[i] = point_fixe[0];
        resu[nb_iterations+i] = point_fixe[1];
        //cout << resu[i+nb_iterations][1] << endl;
    }
    return resu;
}

// Euler symplectique (implicite sur q, explicite sur p)

FVector<FVector<double,3>,nb_planetes>* pf_euler_symplectique(FVector<FVector<double,3>,nb_planetes> qn,FVector<FVector<double,3>,nb_planetes> pn ){
    double epsilon = 0.0001; // precision Comment choisir epsilon ?
    int compteur =0;
    FVector<FVector<double,3>,nb_planetes>* resu = new FVector<FVector<double,3>,nb_planetes> [2];
    FVector<FVector<double,3>,nb_planetes> q0,q1,p0,p1;
    for(int i =0;i<nb_planetes;i++){
            q0[i] = qn[i];
            p0[i] = pn[i];
            q1[i] = qn[i] + h*pn[i]/m[i];
            p1[i] = pn[i] + h*interaction(q0)[i];
    }
    while( compteur<1000 && max(ecart(q0,q1),ecart(p0,p1))> epsilon ){
        compteur+=1;
        for(int i =0;i<nb_planetes;i++){
                q0[i] = q1[i];
                p0[i] = p1[i];
                q1[i] = qn[i] + h*pn[i]/m[i];
                p1[i] = pn[i] + h*interaction(q0)[i];
        }
    }
    resu[0] = q1;
    resu[1] = p1;
    return resu;
}

FVector<FVector<double,3>,nb_planetes>* euler_symplectique(FVector<FVector<double,3>,nb_planetes> q0, FVector<FVector<double,3>,nb_planetes> p0){

    FVector<FVector<double,3>,nb_planetes>* resu = new FVector<FVector<double,3>,nb_planetes> [2*nb_iterations];
    resu[0]=q0;
    resu[nb_iterations]=p0;
    for(int i =1;i<nb_iterations;i++){
        FVector<FVector<double,3>,nb_planetes>* resu_ = new FVector<FVector<double,3>,nb_planetes> [2];
        resu_ = pf_euler_symplectique(resu[i-1],resu[nb_iterations + i-1]);
        resu[i] = resu_[0];
        resu[nb_iterations+i] = resu_[1];
    }
    return resu;
}
