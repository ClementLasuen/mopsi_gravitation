#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from func import *

#%%
    
file_path = "Datas/coord.txt"

nb_iterations,nb_planetes,positions,vitesses = get_data(file_path)
    

#%%

#first_plot(nb_iterations, nb_planetes, positions)


#%%

H = "H.txt"
H_modifie = "H_modifie.txt"

print(nb_iterations)

plot_H(nb_iterations, H, H_modifie)




