#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from func import *

#%%
    
file_path = "Datas/euler_symplectique_sans_pf.txt"

nb_iterations,nb_planetes,positions,vitesses = get_data(file_path)
    

#%%

first_plot(nb_iterations, nb_planetes, positions)


#%%

file_name = "test_h_tot_ES_sspf.txt"

plot_H(nb_iterations, file_name)





