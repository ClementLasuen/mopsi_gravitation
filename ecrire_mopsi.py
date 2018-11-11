# with open("test_h.txt") as f:
#     texte = f.readline()
#     texte2 = f.readline()
# 
# print(texte)
# print(texte2)

import matplotlib.pyplot as plt

def print_H(nb_iterations):
    resu=[]
    with open("C:/Users/Utilisateur/Downloads/Tp2_Initial/Tennis/test_h.txt") as f:
        for i in range(nb_iterations):
            texte = f.readline()
            resu.append(float(texte))
    print(len(resu))
    X=[i for i in range(nb_iterations)]
    plt.plot(X,resu)
    plt.show()