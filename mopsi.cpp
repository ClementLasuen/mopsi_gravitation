#include"mopsi.h"

// --------------------------------Fonctions pratiques---------------------------------

double norme(FVector<double, 3> v){
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

FVector<FVector<double, 3>, nb_planetes> interaction(FVector<FVector<double, 3>, nb_planetes> q){
    FVector<FVector<double,3>,nb_planetes> resu;
    double r;
    for(int i = 0;i<nb_planetes;i++){
        resu[i] = {0.0,0.0,0.0};
        for(int j = 0; j<nb_planetes;j++){
            if(j!=i){
                r = norme(q[i] - q[j])*norme(q[i] - q[j])*norme(q[i] - q[j]);
                resu[i] += - (q[i] - q[j])*G*m[i]* m[j] *1./r;
            }
        }
    }
    return resu;
}

double H(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p){
    double resu = norme(p[0])*norme(p[0])/(2*m[0]) ;
    for(int i =1;i<nb_planetes;i++){
        resu += norme(p[i])*norme(p[i])/(2*m[i]);
        for(int j=0;j<i-1;j++){
            resu -= G*m[i]*m[j]/norme(q[i]-q[j]);
        }
    }
    return resu;
}

// ------------------------------------------Methodes d'integration------------------------------

double ecart(FVector<FVector<double, 3>, nb_planetes> q0, FVector<FVector<double,3>, nb_planetes> q1){
    double resu = 0.0;
    for(int i=0;i<nb_planetes;i++){
        resu += norme(q0[i] - q1[i]);
    }
    return resu;
}


// Euler explicite

FVector<FVector<double,3>,nb_planetes>* euler_explicite(FVector<FVector<double,3>,nb_planetes> q0, FVector<FVector<double,3>,nb_planetes> p0, bool ecriture){

    string file_name = string("../mopsi_gravitation/Datas/euler_explicite.txt"); //+string<int>(nb_iterations)+string("_")+string<int>(h);
    ofstream fichier(file_name.c_str(), ios::out|ios::trunc); // On va ecrire a la fin du fichier
    if (fichier){
        if(ecriture){
            cout <<"OK"<<endl;
            fichier << nb_iterations <<" "<<nb_planetes<<endl; // On ecrit les donnees
            for(int j=0;j<nb_planetes;j++) // On ecrit la ligne avec les masses
                fichier << m[j]<<" ";
            fichier <<endl;
            for(int j=0;j<nb_planetes;j++)
                fichier << q0[j][0] << " " << q0[j][1] << " " << q0[j][2] << p0[j][0] << " " << p0[j][1] << " " << p0[j][2] << endl;; // On ecrit les conditions initiales
        }
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
            if(ecriture){
                for(int j=0;j<nb_planetes;j++)
                    fichier << resu_[0][j][0] << " " << resu_[0][j][1] << " " << resu_[0][j][2] << resu_[1][j][0]/m[j] << " " << resu_[1][j][1]/m[j] << " " << resu_[1][j][2]/m[j] << endl; // On ecrit les conditions initiales
            }
            delete[] resu_;
        }
        fichier.close();
        return resu;
    }
    else{
        cerr<<"Impossible d'ouvrir le fichier"<<endl;
    }
}

// Euler imlicite


FVector<FVector<double, 3>, nb_planetes> *pf_euler_implicite(FVector<FVector<double, 3>, nb_planetes> qn, FVector<FVector<double, 3>, nb_planetes> pn ){
    double epsilon = 0.000001; // precision comment choisir epsilon ?
    int compteur =0;
    FVector<FVector<double,3>,nb_planetes>* resu = new FVector<FVector<double,3>,nb_planetes> [2];
    FVector<FVector<double,3>,nb_planetes> q0,q1,p0,p1;
    q0 = qn;
    p0 = pn;
    FVector<FVector<double, 3>, nb_planetes> force = interaction(q0);
    for(int i =0;i<nb_planetes;i++){  // je pense qu'on peut eviter la boucle for et juste ajouter les vecteurs
            q1[i] = qn[i] + h*p0[i]/m[i];
            p1[i] = pn[i] + h*force[i];
    }
    while( compteur<100 && max(ecart(q0,q1),ecart(p0,p1))> epsilon ){

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

    }
    return resu;
}

// Euler symplectique (implicite sur q, explicite sur p)

FVector<FVector<double, 3>, nb_planetes> *pf_euler_symplectique(FVector<FVector<double, 3>, nb_planetes> qn, FVector<FVector<double, 3>, nb_planetes> pn ){
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

FVector<FVector<double, 3>, nb_planetes> *euler_symplectique(FVector<FVector<double, 3>, nb_planetes> q0, FVector<FVector<double, 3>, nb_planetes> p0){

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

// Verlet

void changement_variables(FVector<FVector<double, 3>, nb_planetes> &p){
    FVector<FVector<double,3>,nb_planetes> p_ = p;
    for(int i=0;i<nb_planetes;i++){
        p[i]= p_[i]/m[i];
    }
}

FVector<FVector<double, 3>, nb_planetes> *verlet(FVector<FVector<double, 3>, nb_planetes> q0, FVector<FVector<double, 3>, nb_planetes> p0){
    FVector<FVector<double,3>,nb_planetes>* resu = new FVector<FVector<double,3>,nb_planetes> [3*nb_iterations];

    // resu[i] donne les positions de toutes les planètes à l'iteration i;
    // resu[i + nb_iterations] donne les quantites de mouvement
    // resu[i + nb_iterations] donne p_n+1/2

    changement_variables(p0);
    resu[0]=q0;
    resu[nb_iterations]=p0;

    FVector<FVector<double, 3>, nb_planetes> force,force_1;
    force_1 = interaction(q0);

    for(int i = 0;i<nb_iterations-1;i++){
        force = force_1;

        //resu[2*nb_iterations + i+1] = resu[nb_iterations+i] + force*h/2; NE MARCHE PAS
        //resu[i+1] = resu[i] + h*resu[2*nb_iterations + i]; NE MARCHE PAS

        for(int j=0;j<nb_planetes;j++){

            resu[2*nb_iterations + i+1][j] = resu[nb_iterations+i][j] + force[j]*h/(2*m[j]);
            resu[i+1][j] = resu[i][j] + h*resu[2*nb_iterations + i][j];

        }
        force_1 = interaction(resu[i+1]);
        for(int j =0;j<nb_planetes;j++){
            resu[nb_iterations + i+1][j] = resu[2*nb_iterations + i+1][j] + h*force_1[j]/(2*m[j]);
        }
    }
    return resu;
}

// Ne pas oublier de le faire lorsque on appelle interaction -> force/masse




