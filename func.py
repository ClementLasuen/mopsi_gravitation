#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 15:25:10 2019

@author: quentin
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si
import scipy.optimize as so
from mpl_toolkits.mplot3d import Axes3D
from vpython import *

# Cette fonction stocke les données le fichier coord.txt dans les vecteurs positions et vitesses
# position[i,j,x] est la coordonnée x de la planète j à l'itération i (idem pour les vitesses).
def get_data(file_path = "Datas/coord.txt"):
    # Lit le fichier C++ avec les données
    path = file_path
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

# Cette fonction trace l'énergie (hamiltonien) du système en fonction des itérations.
# Les valeurs de H sont calculées lors du calcul des positions et de vitesses puis stockées
# dans le fichier H.txt ou H_modifie.txt
def plot_H(nb_iterations, file_name, file_name_modifie):
    # Lit le fichier C++ avec les valeurs de H, et les affiches
    resu=[]
    resu1=[]
    with open("Datas/"+file_name) as f:
        for i in range(nb_iterations-1):
            texte = f.readline()
            resu.append(float(texte))
    with open("Datas/"+file_name_modifie) as f1:
        for i in range(nb_iterations-1):
            texte = f1.readline()
            resu1.append(float(texte))
    X=[i for i in range(nb_iterations-1)]
    plt.figure()
    plt.plot(X,resu, label = "normal")
    plt.plot(X,resu1, label = "modifie")
    plt.legend()
    plt.show()
    return()
    
def first_plot(nb_iterations, nb_planetes, positions):
    # Fait un affichage sommaire des trajectoires
    ax = plt.gca(projection='3d')
    ax.plot([0.],[0.],[0.],'o',label='Soleil')
    for i in range(1,nb_planetes):
        ax.plot(positions[:,i,0],positions[:,i,1],positions[:,i,2])
    plt.axis('equal')
    plt.show()
    
    return()
    
#affichage avancé
    
pause_var = False
    
def enhanced_plot(nb_iterations, nb_planetes, positions, vitesses):
    
    scene.title = "Enhanced 3D of surfaces using bump maps"
    scene.caption = "Drag the single light with the left button, rotate with the right button. \n\n"
    scene.fullscreen = 1
    
    scene.camera.pos = vector(0,0,0)
    
    
    # Initialisation des planètes
    
    pos1=vector(positions[0][0][0],positions[0][0][1],positions[0][0][2])
    soleil = sphere(pos=pos1,radius =4,texture = 'sun.jpg', emissive =True,
                    interval=10)
    
    pos1=vector(positions[0][1][0],positions[0][1][1],positions[0][1][2])
    
    jupiter = sphere(pos=pos1,radius =2,texture = 'jupiter.jpg',
                     interval=10)
    
    pos1=vector(positions[0][2][0],positions[0][2][1],positions[0][2][2])
    saturne = sphere(pos=pos1,radius =2,texture = 'saturn.png',
                     interval=10)
    
    pos1=vector(positions[0][3][0],positions[0][3][1],positions[0][3][2])
    uranus = sphere(pos=pos1,radius =2,texture = 'uranus.jpg',
                    interval=10)
    
    pos1=vector(positions[0][4][0],positions[0][4][1],positions[0][4][2])
    neptune = sphere(pos=pos1,radius =2,texture = 'neptune.jpg',
                     interval=10)
    
    # Affichage de la trajectoire
    
    jupiter.trail = curve(color = vector(1,0.8,0.65), radius=0.05)
    saturne.trail = curve(color = vector(1,0.8,0.5), radius=0.05)
    uranus.trail = curve(color = color.cyan, radius=0.05)
    neptune.trail = curve(color = color.blue, radius=0.05)
    
    # Affichage de l'anneau de Saturne 
    precious = ring(pos=saturne.pos, axis = saturne.pos,radius = 3, thickness=0.1)
    
    # Permet que l'echelle soit mise de facon a voir toutes les planetes
    scene.autoscale=False
    
    #Le soleil est source de lumiere
    lamp = local_light(pos=vector(0,0,0), color=color.white)
    
    # On ajoute un décor
    decor = sphere(pos=vector(0,0,0),radius = 80,texture='ciel.jpg')
    
    # vecteur de rotation des planètes sur elles-mêmes
    rotation = vector(-0.5,1,0)
    
    
    
    #boutton de pause / play
    
    c = controls 
    
    def pause():
        global pause_var
        pause_var = not pause_var
        return()
        
    i = 0
    pas_affichage = 10
    

    def rate_modif():
        value+=1
        return()
    
    # Permet de gerer la vitesse d'affichage
    # Le temps entre deux affichages est value millisecondes
    s = slider(bind = rate_modif,min=1, max=150, value=20)
    scene.append_to_caption('')
    b = button(bind = pause, action = pause, text='play/pause')

    while i<nb_iterations:
        
        rate(s.value)
        
        if pause_var:
            continue
            
        # mise à jour de la position des planètes, de leur trace, et de leur rotation
        
        pos1=vector(positions[i][0][1],positions[i][0][2],positions[i][0][0])
        soleil.pos=pos1
        
        pos1=vector(positions[i][1][1],positions[i][1][2],positions[i][1][0])
        jupiter.pos=pos1
        jupiter.trail.append(pos=jupiter.pos)
        jupiter.rotate(angle=0.2, axis=rotation)
        
        pos1=vector(positions[i][2][1],positions[i][2][2],positions[i][2][0])
        saturne.pos=pos1
        saturne.trail.append(pos=saturne.pos)
        saturne.rotate(angle=0.2, axis=rotation)
        
        
        precious.pos = saturne.pos
        
        pos1=vector(positions[i][3][1],positions[i][3][2],positions[i][3][0])
        uranus.pos=pos1
        uranus.trail.append(pos=uranus.pos)
        saturne.rotate(angle=0.2, axis=rotation)
        
        pos1=vector(positions[i][4][1],positions[i][4][2],positions[i][4][0])
        neptune.pos=pos1
        neptune.trail.append(pos=neptune.pos)
        neptune.rotate(angle=0.2, axis=rotation)
        
        i+=pas_affichage

    return()
    
