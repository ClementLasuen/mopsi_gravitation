# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""
from func import *

#%%

######################################### TEEEEEEEEST !!!!!!!!

def norme(q):
    return sqrt(q[0]**2 + q[1]**2 + q[2]**2)

def H(q,v):
    
    #passer de la vitesse à la quantité de mouvement
    for i in range(nb_planetes):
        for j in range(3):
            v[i,j] = liste_m[i]*v[i,j]
    resu=0
    for i in range(nb_planetes):
        resu += norme(v[i,:])*2 /liste_m[i]
        for j in range(i):
            resu -= G*liste_m[i]*liste_m[j]/norme(positions[i,:] - positions[j,:])
    
    return resu;
    
    
    
import unittest

class cas_tests(unittest.TestCase):
   
     def test_BOUNDED(self,verlet):
         C = 1000
         # JE CALCULE COMME CA CAR JE N AI PAS ACCES AU PAS DE TEMPS
         h = int ((positions[1,0,0] - positions[0,0,0])/vitesses[1,0,0] )
         if(verlet):
             h = h*h
         for i in range(nb_iterations) :
             self.assertTrue(abs(H(position[i,:,:],vitesses[i,:,:])) <= h*C)

    
#if __name__ == "__main__":
#    unittest.main()

######################################################################################"


#%%
             
file_path = "Datas/euler_symplectique_sans_pf.txt"

nb_iterations, nb_planetes,positions,vitesses = get_data(file_path)

enhanced_plot(nb_iterations, nb_planetes,positions,vitesses)
