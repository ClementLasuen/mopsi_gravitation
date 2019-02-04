# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""
from func import *

             
file_path = "Datas/coord.txt"

nb_iterations, nb_planetes,positions,vitesses = get_data(file_path)

enhanced_plot(nb_iterations, nb_planetes,positions,vitesses)
