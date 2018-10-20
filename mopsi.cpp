
#include"mopsi.h"



std::vector<float*>* q = new std::vector< float*> [nb_planetes];
std::vector<float*>* p = new std::vector< float*> [nb_planetes];

float norme(float* v){
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

std::vector<float*> interaction(std::vector<float*> q){
    std::vector<float*> resu;
    for(int i = 0;i<nb_planetes;i++){
        for(int j = 0; j<nb_planetes;j++){
            if(j!=i){
                resu[i] += (q[i] - q[j])*G* m[j] /(norme(q[i] - q[j])*norme(q[i] - q[j])*norme(q[i] - q[j]))
            }
        }
    }
    return resu;
}

float ecart(std::vector<float *> q0, std::vector<float *> q1){
    float resu = 0.0;
    for(int i=0;i<nb_planetes;i++){
        resu += norme(q0[i] - q1[i]);
    }
    return resu;
}

std::vector<float*>* point_fixe(std::vector<float*>* qn, std::vector<float*>* pn, float masse ){
    epsilon = 0.00001; // precision

    std::vector<float*>* q1 = new std::vector<float*> [nb_planetes];
    std::vector<float*>* p1 = new std::vector<float*> [nb_planetes];
    std::vector<float*>* q0 = new std::vector<float*> [nb_planetes];
    std::vector<float*>* p0 = new std::vector<float*> [nb_planetes];
    for(int i =0;i<nb_planetes;i++){
        q0[i] = new float [3];
        p0[i] = new float [3];
        q1[i] = new float [3];
        p1[i] = new float [3];
        /*for(int j=0;j<3;j++){
            q0[i][j] = qn[i][j];
            p0[i][j] = pn[i][j];
            q1[i][j] = qn[i][j] + h*p0[i][j]/masse;
            p1[i][j] = pn[i][j] + h*interaction(p[0]);
        }*/
        q0=qn;
        p0=pn;
        q1 = qn + p0*h/masse;
        p1 = pn + h*interaction(q0);
    }
    while( compteur<1000){ // rajouter condition sur abs(q1 - q0), abs(p1 - p0)
        compteur+=1;
        p0 = p1;
        q0 = q1;
        q1 = qn + p0*h/masse;
        p1 = pn + h*interaction(q0);
    }
}
