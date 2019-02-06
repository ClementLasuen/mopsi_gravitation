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
                resu[i] += - (q[i] - q[j])*G*m[i]* m[j]/r;
            }
        }
    }
    return resu;
}

// --------------------------------------- Hamiltonien ------------------------------------------

double H(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p){
    double resu = 0 ;
    for(int i =0;i<nb_planetes;i++){
        resu += norme(p[i])*norme(p[i])/(2*m[i]);
        for(int j=0;j<i;j++){
            resu -= G*m[i]*m[j]/norme(q[i]-q[j]);
        }
    }
    return resu;
}
// Hamiltonien modifié pour euler symplectique (implicite sur p) / il s'agit de H3 (voir article) -> *h pour H

double H_modifie_ES(FVector<FVector<double,3>,nb_planetes> q, FVector<FVector<double,3>,nb_planetes> p){
    double resu = 0;
    FVector<FVector<double, 3>, nb_planetes> delta_V = interaction(q);
    for(int i=0;i<nb_planetes;i++){
        for(int j=0;j<3;j++) {
            resu +=  -p[i][j]* delta_V[i][j]/(2*m[i]);
        }
    }
    return resu;
}

// Hamiltonien modifié pour Verlet / il s'agit de H3 (voir article) -> *h² pour H

double H_modifie_V(FVector<FVector<double, 3>, nb_planetes> q, FVector<FVector<double, 3>, nb_planetes> p){
    double resu = 0.0;
    FVector<FVector<double, 3>, nb_planetes> delta_V = dV(q);  //interaction(q);
    FVector<FVector<double, 3>, nb_planetes> delta_V_M = delta_V;
    changement_variables(delta_V_M);
    for(int i=0;i<nb_planetes;i++){
        for(int coord_i=0;coord_i<3;coord_i++)   resu += -delta_V[i][coord_i]*delta_V_M[i][coord_i]/24.0;
    }

    // Calcul de la Hessienne

    FVector<FVector<double, 3*nb_planetes>,3*nb_planetes > hess = Hessienne(q);// Hessienne(q);
    for(int i=0;i<nb_planetes;i++){
        for(int j=0;j<nb_planetes;j++){
            for(int coord_i=0;coord_i <3;coord_i++){
                for(int coord_j=0;coord_j <3;coord_j++) resu += p[j][coord_j]*hess[j*3 +coord_j][i*3+coord_i]*p[i][coord_i] /12.0;
            }
        }
    }
    return resu;
}

double egalite_hessiennes(FVector<FVector<double, 3>, nb_planetes> q){
    double resu=0;
    FVector<FVector<double, 3*nb_planetes>,3*nb_planetes > H1 = Hessienne(q);
    //FVector<FVector<double, 3*nb_planetes>,3*nb_planetes > H2 = d2V(q);
    for(int i=0;i<nb_planetes*3;i++){
        for(int j=0;j<3*nb_planetes;j++){
            //resu += (H1[i][j] - H2[i][j])*(H1[i][j] - H2[i][j]);
            resu += H1[i][j]*H1[i][j];
        }
    }
    return sqrt(resu);
}


// Calcul la hessienne du potentiel selon q

FVector<FVector<double, 3*nb_planetes>,3*nb_planetes > Hessienne(FVector<FVector<double, 3>, nb_planetes> q){

    FVector<FVector<double, 3*nb_planetes>,3*nb_planetes > result;
    // coeff yi xi , i pour une meme planete

    double dxi_dyi=0;
    double dxi_dyj=0;


    for(int i=0;i<nb_planetes;i++){
        for(int j=0;j<nb_planetes;j++){
            if(j==i){

                for(int coord =0; coord<3; coord++){
                    for(int coord2=0; coord2<3; coord2++){
                        dxi_dyi=0;
                        for(int k=0;k<nb_planetes;k++){
                            if(k!=i){
                                dxi_dyi -= 3*G*m[i]*m[k]*(q[i][coord] -q[k][coord])*(q[i][coord2]-q[k][coord2])/pow(norme(q[i] - q[k]),5.);
                                if (coord==coord2)
                                    dxi_dyi += G*m[i]*m[k]/pow(norme(q[i] - q[k]),3.);
                            }
                        }
                        result[i*3+coord][j*3 +coord2] = dxi_dyi;
                    }
                }
            }
            else{
                for(int coord =0;coord<3;coord++){
                    for(int coord2=0; coord2<3;coord2++){
                        dxi_dyj =  3*G*m[i]*m[j]*(q[i][coord]-q[j][coord])*(q[i][coord2]-q[j][coord2])/pow(norme(q[i] - q[j]),5.);
                        if (coord==coord2)
                            dxi_dyj -= G*m[i]*m[j]/pow(norme(q[i] - q[j]),3.);
                        result[i*3+coord][j*3 +coord2]= dxi_dyj;
                    }
                }
            }
        }
    }
    return result;
}


// ------------------------------------------Methodes d'integration------------------------------

double ecart(FVector<FVector<double, 3>, nb_planetes> q0, FVector<FVector<double,3>, nb_planetes> q1){
    double resu = 0.0;
    for(int i=0;i<nb_planetes;i++) resu += norme(q0[i] - q1[i]);
    return resu;
}


// Euler explicite

void euler_explicite(double h, bool ecriture){
    // --------------------------------- Initialistation ---------------------------------------------------------

    FVector<FVector<double,3>,nb_planetes> p0,q0,p,q;

    q0[0]=q_soleil;
    q0[1]=q_jupiter;
    q0[2]=q_saturne;
    q0[3]=q_uranus;
    q0[4]=q_neptune;

    p0[0]=p_soleil;
    p0[1]=p_jupiter;
    p0[2]=p_saturne;
    p0[3]=p_uranus;
    p0[4]=p_neptune;

    p=p0;
    q=q0;


    //-------------------------------------------------------------------------------------------------------------
    string file_name = string("../mopsi_gravitation/Datas/coord.txt");
    //+string<int>(nb_iterations)+string("_")+string<int>(h);
    string file_name_H = string("../mopsi_gravitation/Datas/H.txt");
    ofstream fichier(file_name.c_str(), ios::out|ios::trunc);
    ofstream fichier_H(file_name_H.c_str(), ios::out|ios::trunc);// On va ecrire a la fin du fichier

    if (fichier){
        if(ecriture){
            cout <<"calcul des trajectoires"<<endl;
            fichier << nb_iterations <<" "<<nb_planetes<<endl; // On ecrit les donnees
            for(int j=0;j<nb_planetes;j++) // On ecrit la ligne avec les masses
                fichier << m[j]<<" ";
            fichier <<endl;
            for(int j=0;j<nb_planetes;j++){
                fichier << q0[j][0] << " " << q0[j][1] << " " << q0[j][2] << " " << p0[j][0]/m[j] << " " << p0[j][1]/m[j] << " " << p0[j][2]/m[j] << " "; // On ecrit les conditions initiales

            }
            fichier << endl;
        }


        for(int i=1;i<nb_iterations;i++){
            if (i%(nb_iterations/100)==0)               // On affiche l'avancée de l'ecriture
                cout << int(i/int(nb_iterations/100)) << endl;

            FVector<FVector<double,3>,nb_planetes> p_point = interaction(q0);
            for(int j=0;j<nb_planetes;j++){

                q[j] = q[j] + h*p0[j]/m[j];
                p[j] = p[j]+ h*p_point[j];
            }

            if(ecriture){
                for(int j=0;j<nb_planetes;j++){
                    fichier << q0[j][0] << " " << q0[j][1] << " " << q0[j][2] << " " << p0[j][0]/m[j] << " " << p0[j][1]/m[j] << " " << p0[j][2]/m[j]<< " "; // On ecrit les positions puis vitesses d'une planete
            }
            fichier << endl;
            fichier_H << H(q,p);
            fichier_H << endl;

            }
            p0=p;
            q0=q;


        }
        fichier.close();
        fichier_H.close();

    }
    else{
        cerr<<"Impossible d'ouvrir le fichier"<<endl;
    }
}

// Euler implicite


FVector<FVector<double, 3>, nb_planetes> *pf_euler_implicite(double h,FVector<FVector<double, 3>, nb_planetes> qn, FVector<FVector<double, 3>, nb_planetes> pn ){
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

void euler_implicite(double h, bool ecriture){


    // --------------------------------- Initialistation ---------------------------------------------------------

    FVector<FVector<double,3>,nb_planetes> p0,q0,p,q;

    q0[0]=q_soleil;
    q0[1]=q_jupiter;
    q0[2]=q_saturne;
    q0[3]=q_uranus;
    q0[4]=q_neptune;

    p0[0]=p_soleil;
    p0[1]=p_jupiter;
    p0[2]=p_saturne;
    p0[3]=p_uranus;
    p0[4]=p_neptune;


    p=p0;
    q=q0;
    //-------------------------------------------------------------------------------------------------------------

    string file_name = string("../mopsi_gravitation/Datas/coord.txt"); //+string<int>(nb_iterations)+string("_")+string<int>(h);
    ofstream fichier(file_name.c_str(), ios::out|ios::trunc); // On va ecrire a la fin du fichier

    string file_name_H = string("../mopsi_gravitation/Datas/H.txt");
    ofstream fichier_H(file_name_H.c_str(), ios::out|ios::trunc);// On va ecrire a la fin du fichier

    if (fichier){
        if(ecriture){
            cout <<"calcul des trajectoires"<<endl;
            fichier << nb_iterations <<" "<<nb_planetes<<endl; // On ecrit les donnees
            for(int j=0;j<nb_planetes;j++) // On ecrit la ligne avec les masses
                fichier << m[j]<<" ";
            fichier <<endl;
            for(int j=0;j<nb_planetes;j++)
                fichier << q0[j][0] << " " << q0[j][1] << " " << q0[j][2] << " " << p0[j][0]/m[j] << " " << p0[j][1]/m[j] << " " << p0[j][2]/m[j] << " "; // On ecrit les conditions initiales
            fichier << endl;
        }

        FVector<FVector<double,3>,nb_planetes>* point_fixe = new FVector<FVector<double,3>,nb_planetes> [2];
        for(int i =1;i<nb_iterations;i++){
            if (i%(nb_iterations/100)==0)               // On affiche l'avancée de l'ecriture
                cout << int(i/int(nb_iterations/100)) << endl;
            point_fixe = pf_euler_implicite(h,q,p);
            q = point_fixe[0];
            p = point_fixe[1];


            if(ecriture){
                for(int j=0;j<nb_planetes;j++)
                    fichier << q[j][0] << " " << q[j][1] << " " << q[j][2] << " " << p[j][0]/m[j] << " " << p[j][1]/m[j] << " " << p[j][2]/m[j]<< " "; // On ecrit les positions puis vitesses d'une planete
            fichier << endl;
            }
        fichier_H << H(q,p);
        fichier_H << endl;
        }
    }
    else{
        cerr<<"Impossible d'ouvrir le fichier"<<endl;
    }
}

// Euler symplectique (implicite sur q, explicite sur p)

void euler_symplectique(double h, bool ecriture){

    // --------------------------------- Initialistation ---------------------------------------------------------

    FVector<FVector<double,3>,nb_planetes> p0,q0,p,q;

    q0[0]=q_soleil;
    q0[1]=q_jupiter;
    q0[2]=q_saturne;
    q0[3]=q_uranus;
    q0[4]=q_neptune;

    p0[0]=p_soleil;
    p0[1]=p_jupiter;
    p0[2]=p_saturne;
    p0[3]=p_uranus;
    p0[4]=p_neptune;

    p=p0;
    q=q0;

    //-------------------------------------------------------------------------------------------------------------

    string file_name = string("../mopsi_gravitation/Datas/coord.txt"); //+string<int>(nb_iterations)+string("_")+string<int>(h);
    string file_name_H = string("../mopsi_gravitation/Datas/H.txt");
    ofstream fichier(file_name.c_str(), ios::out|ios::trunc); // On va ecrire a la fin du fichier
    ofstream fichier_H(file_name_H.c_str(), ios::out|ios::trunc);

    string file_name_H_modifie = string("../mopsi_gravitation/Datas/H_modifie.txt");
    ofstream fichier_H_modifie(file_name_H_modifie.c_str(), ios::out|ios::trunc);


    if (fichier){
        if(ecriture){
            cout <<"calcul des trajectoires"<<endl;
            fichier << nb_iterations <<" "<<nb_planetes<<endl; // On ecrit les donnees
            for(int j=0;j<nb_planetes;j++) // On ecrit la ligne avec les masses
                fichier << m[j]<<" ";
            fichier <<endl;
            for(int j=0;j<nb_planetes;j++)
                fichier << q0[j][0] << " " << q0[j][1] << " " << q0[j][2] << " " << p0[j][0]/m[j] << " " << p0[j][1]/m[j] << " " << p0[j][2]/m[j] << " "; // On ecrit les conditions initiales
            fichier << endl;
        }

        for(int i=1;i<nb_iterations;i++){
            if (i%(nb_iterations/100)==0)               // On affiche l'avancée de l'ecriture
                cout << int(i/int(nb_iterations/100)) << endl;
            for(int j=0; j<nb_planetes;j++) q[j] += h*p0[j]/m[j];

            FVector<FVector<double,3>,nb_planetes> delta_V = interaction(q);
            for(int j=0; j<nb_planetes;j++) p[j] += h*delta_V[j];

            if(ecriture){
                for(int j=0;j<nb_planetes;j++)
                    fichier << q[j][0] << " " << q[j][1] << " " << q[j][2] << " " << p[j][0]/m[j] << " " << p[j][1]/m[j] << " " << p[j][2]/m[j]<< " "; // On ecrit les positions puis vitesses d'une planete
            fichier << endl;
            fichier_H << H(q,p);
            fichier_H << endl;
            fichier_H_modifie << H(q,p) + h*H_modifie_ES(q,p);
            fichier_H_modifie << endl;
            p0=p;
            q0=q;
            }
        }
        fichier.close();
        fichier_H.close();
        fichier_H_modifie.close();
    }
    else{
        cerr<<"Impossible d'ouvrir le fichier"<<endl;
    }
}




// Verlet

void changement_variables(FVector<FVector<double, 3>, nb_planetes> &p){
    for(int i=0;i<nb_planetes;i++) p[i]= p[i]/m[i];
}

void changement_variables_inverse(FVector<FVector<double, 3>, nb_planetes> &p){
    for(int i=0;i<nb_planetes;i++) p[i]= p[i]*m[i];
}

void verlet(double h,bool ecriture){

    // --------------------------------- Initialistation ---------------------------------------------------------

    FVector<FVector<double,3>,nb_planetes> p0,q0,p,q;

    q0[0]=q_soleil;
    q0[1]=q_jupiter;
    q0[2]=q_saturne;
    q0[3]=q_uranus;
    q0[4]=q_neptune;

    p0[0]=p_soleil;
    p0[1]=p_jupiter;
    p0[2]=p_saturne;
    p0[3]=p_uranus;
    p0[4]=p_neptune;

    p = p0;
    q = q0;
    double hamiltonien, hamiltonien_modifie;
    double difference_grad = 0.0;
    FVector<FVector<double, 3>, nb_planetes> grad_diff_finies;

    FVector<FVector<double, 3>, nb_planetes> grad_initial = dV(q);

    //-------------------------------------------------------------------------------------------------------------

    string file_name = string("../mopsi_gravitation/Datas/coord.txt");
    //+string<int>(nb_iterations)+string("_")+string<int>(h);
    string file_name_H = string("../mopsi_gravitation/Datas/H.txt");
    ofstream fichier(file_name.c_str(), ios::out|ios::trunc); // On va ecrire a la fin du fichier
    ofstream fichier_H(file_name_H.c_str(), ios::out|ios::trunc);

    string file_name_H_modifie = string("../mopsi_gravitation/Datas/H_modifie.txt");
    ofstream fichier_H_modifie(file_name_H_modifie.c_str(), ios::out|ios::trunc);
    if (fichier){
        if(ecriture){
            cout <<"calcul des trajectoires"<<endl;
            fichier << nb_iterations <<" "<<nb_planetes<<endl; // On ecrit les donnees
            for(int j=0;j<nb_planetes;j++) // On ecrit la ligne avec les masses
                fichier << m[j]<<" ";
            fichier <<endl;
            for(int j=0;j<nb_planetes;j++)
                fichier << q0[j][0] << " " << q0[j][1] << " " << q0[j][2] << " " << p0[j][0]/m[j] << " " << p0[j][1]/m[j] << " " << p0[j][2]/m[j] << " "; // On ecrit les conditions initiales
            fichier << endl;
        }
        changement_variables(p0);
        changement_variables(p);

        FVector<FVector<double, 3>, nb_planetes> force = interaction(q);

        for(int i = 0;i<nb_iterations-1;i++){

            if (i%int(nb_iterations/100)==0)
                cout << int(i/int(nb_iterations/100)) << endl;


            for(int j=0;j<nb_planetes;j++){
                p[j] += force[j]*h/(2.0*m[j]);
                q[j] += h*p[j];
            }

            force = interaction(q);

            for(int j =0;j<nb_planetes;j++)  p[j] += h*force[j]/(2.0*m[j]);

            if(ecriture){
                for(int j=0;j<nb_planetes;j++)
                    fichier << q[j][0] << " " << q[j][1] << " " << q[j][2] << " " <<p[j][0] << " " << p[j][1] << " " << p[j][2]<< " "; // On ecrit les positions puis vitesses d'une planete
            fichier << endl;

            changement_variables_inverse(p);
            hamiltonien = H(q,p);
            changement_variables(p);

            fichier_H << hamiltonien;
            fichier_H << endl;

            fichier_H_modifie << hamiltonien + h*h*H_modifie_V(q,p);
            //fichier_H_modifie << h*h*H_modifie_V(q,p);
            //fichier_H_modifie << egalite_hessiennes(q);

            // Test on retrouve le gradient ?

            /*
            difference_grad = 0;
            grad_diff_finies = dV(q); // Attention, - dV = force
            for(int i=0;i<nb_planetes;i++){   //difference_grad += (force[i] + grad_diff_finies[i])*(force[i] + grad_diff_finies[i]);
                for(int coord =0;coord<3;coord++){
                    difference_grad += (grad_diff_finies[i][coord] - grad_initial[i][coord])*(grad_diff_finies[i][coord] - grad_initial[i][coord]);
                }
            }

            fichier_H_modifie << sqrt(difference_grad);

            */

            //changement_variables(p);
            fichier_H_modifie << endl;
            }
            p0=p;
            q0=q;
        }
        fichier.close();
        fichier_H.close();
        fichier_H_modifie.close();
    }
    else{
        cerr<<"Impossible d'ouvrir le fichier"<<endl;
    }
}

// Calcul du gradient et de la Hessienne en differences finies

double V(FVector<FVector<double, 3>, nb_planetes> q){
    double resu = 0.0;
    for(int i =0;i<nb_planetes;i++){
        for(int j=0;j<i;j++) resu -= G*m[i]*m[j]/norme(q[i]-q[j]);
    }
    return resu;
}

FVector<FVector<double, 3>, nb_planetes> dV(FVector<FVector<double, 3>, nb_planetes> q){
    double epsilon = 0.0001;
    FVector<FVector<double, 3>, nb_planetes> grad , q_plus_epsilon , q_moins_epsilon ;

    for(int i=0;i<nb_planetes;i++){
        for(int coord = 0;coord<3;coord++){
            q_plus_epsilon = q;
            q_moins_epsilon = q;
            q_plus_epsilon[i][coord] += epsilon;
            q_moins_epsilon[i][coord] -= epsilon;

            grad[i][coord] = (V(q_plus_epsilon)/epsilon - V(q_moins_epsilon)/epsilon)/2.0;
        }
    }
    return grad;
}

FVector<FVector<double, 3*nb_planetes>,3*nb_planetes > d2V(FVector<FVector<double, 3>, nb_planetes> q){
    double epsilon = 0.0001;
    FVector<FVector<double, 3*nb_planetes>,3*nb_planetes > resu;

    FVector<FVector<double, 3>, nb_planetes> q_plus_epsilon, q_moins_epsilon;

    for(int i =0;i<nb_planetes;i++){
        for(int j=0;j<nb_planetes;j++){
            for(int coord_i =0;coord_i<3;coord_i++){
                for(int coord_j=0;coord_j<3;coord_j++){

                    q_plus_epsilon = q;
                    q_moins_epsilon = q;
                    q_plus_epsilon[j][coord_j] += epsilon;
                    q_moins_epsilon[j][coord_j] += - epsilon;

                    resu[3*i + coord_i][3*j + coord_j] = (dV(q_plus_epsilon)[i][coord_i]/epsilon - dV(q_moins_epsilon)[i][coord_i]/epsilon)/2.0;
                }
            }
        }
    }
    return resu;
}
