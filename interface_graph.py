# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""
from vpython import *

import numpy as np
#terre = sphere(pos=vector(0,0,0), radius=10, color=color.blue)
#lune = sphere()
#vaisseau = sphere(pos=vector(0,30,0), radius=3, color=color.red)
#
#angle = 0.0
#distance = 40.0
#i=0
#while i<1000:
#   pos1=vector(distance * cos(angle),distance * sin(angle),0) 
#   lune.pos=pos1
#   angle = angle + pi/6
#   rate(10)
  
  #lune.x = distance * cos(angle)
  #lune.y = distance * sin(angle)
  #vaisseau.x = vaisseau.x + random() - 0.5
  #vaisseau.y = vaisseau.y + random() - 0.5


def get_data():
    # prend le fichier renvoyé par le code C++ pour une trajectoire
    # renvoie le nombre d'iteration, le nombre de planetes, la liste des position de toutes les planetes a chaque temps, la liste des vitesses de toutes les planetes a chaque temps
    with open("/media/OS/Users/Quentin/Documents/ENPC/2A/MOPSI/mopsi_gravitation/Datas/euler_explicite.txt") as datas:
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
    

nb_iterations, nb_planetes,positions,vitesse = get_data()
agrandissement =5
pos1=vector(positions[0][0][0],positions[0][0][1],positions[0][0][2])
soleil = sphere(pos=pos1,radius =4,color=color.yellow, make_trail = True, trail_type="curve",
              interval=10)#, materials = materials.wood)

pos1=vector(positions[0][1][0],positions[0][1][1],positions[0][1][2])

jupiter = sphere(pos=pos1,radius =2,color=color.red,make_trail = True, trail_type="curve",
              interval=10)

pos1=vector(positions[0][2][0],positions[0][2][1],positions[0][2][2])
saturne = sphere(pos=pos1,radius =2,color=color.white,make_trail = True, trail_type="curve",
              interval=10)

pos1=vector(positions[0][3][0],positions[0][3][1],positions[0][3][2])
uranus = sphere(pos=pos1,radius =2,color=color.green,make_trail = True, trail_type="curve",
              interval=10)

pos1=vector(positions[0][4][0],positions[0][4][1],positions[0][4][2])
neptune = sphere(pos=pos1,radius =2,color=color.blue,make_trail = True, trail_type="curve",
              interval=10)

for i in range(nb_iterations):
    
    pos1=vector(positions[i][0][0],positions[i][0][1],positions[i][0][2])
    soleil.pos=pos1
    
    pos1=vector(positions[i][1][0],positions[i][1][1],positions[i][1][2])
    jupiter.pos=pos1
    
    pos1=vector(positions[i][2][0],positions[i][2][1],positions[i][2][2])
    saturne.pos=pos1
    
    pos1=vector(positions[i][3][0],positions[i][3][1],positions[i][3][2])
    uranus.pos=pos1
    
    pos1=vector(positions[i][4][0],positions[i][4][1],positions[i][4][2])
    neptune.pos=pos1

    rate(100)
