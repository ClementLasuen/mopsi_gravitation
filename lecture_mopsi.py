#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si
import scipy.optimize as so
from mpl_toolkits.mplot3d import Axes3D

def get_data(file_path):
    path = "/media/OS/Users/Quentin/Documents/ENPC/2A/MOPSI/mopsi_gravitation/"+file_path
    # prend le fichier renvoyé par le code C++ pour une trajectoire
    # renvoie le nombre d'iteration, le nombre de planetes, la liste des position de toutes les planetes a chaque temps, la liste des vitesses de toutes les planetes a chaque temps
    with open(path) as datas:
        lines = list(map(str.rstrip, datas.readlines()))
        metadatas = lines[0].split(' ')
        nb_iterations = int(metadatas[0])
        nb_planetes = int(metadatas[1])
        liste_m = np.array([0. for i in range(nb_planetes)])
        masses = lines[1].split(' ')
        for i in range(nb_planetes):
            liste_m[i] = float(masses[i])
        positions = np.zeros((nb_iterations,nb_planetes,3), dtype=float)
        vitesses = np.zeros((nb_iterations,nb_planetes,3), dtype=float)
        for i in range(nb_iterations):
            coordonnees = lines[2+i].split(' ')
            for j in range(nb_planetes):
                positions[i,j,0] = float(coordonnees[6*j+0])
                positions[i,j,1] = float(coordonnees[6*j+1])
                positions[i,j,2] = float(coordonnees[6*j+2])
                vitesses[i,j,0] = float(coordonnees[6*j+3])
                vitesses[i,j,1] = float(coordonnees[6*j+4])
                vitesses[i,j,2] = float(coordonnees[6*j+5])
        return(nb_iterations,nb_planetes,positions,vitesses)

def plot_H(nb_iterations):
    resu=[]
    with open("Datas/test_h.txt") as f:
        for i in range(nb_iterations):
            texte = f.readline()
            resu.append(float(texte))
    print(len(resu))
    X=[i for i in range(nb_iterations)]
    plt.plot(X,resu)
    plt.show()
    

#%% Affichage trajectoires

file_path = "Datas/euler_explicite.txt"

nb_iterations,nb_planetes,positions,vitesses = get_data(file_path)

ax = plt.gca(projection='3d')
ax.plot([0.],[0.],[0.],'o',label='Soleil')
ax.plot(positions[:,1,0],positions[:,1,1],positions[:,1,2],label='Jupiter')
ax.plot(positions[:,2,0],positions[:,2,1],positions[:,2,2],label='Jupiter')
ax.plot(positions[:,3,0],positions[:,3,1],positions[:,3,2],label='Jupiter')
ax.plot(positions[:,4,0],positions[:,4,1],positions[:,4,2],label='Jupiter')
ax.plot(positions[:,5,0],positions[:,5,1],positions[:,5,2],label='Jupiter')
ax.plot(positions[:,6,0],positions[:,6,1],positions[:,6,2],label='Jupiter')
plt.legend()
plt.axis('equal')
plt.show()