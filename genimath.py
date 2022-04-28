# -*- coding: utf-8 -*-
"""
TP2 méthode LU sur des matrice tridiag
Antonin Henriques 
"""
import numpy as np


def TRIreduite(A):
    """
    

    Parameters
    ----------
    A : Matrix 
        DESCRIPTION : Matrice tridiag de forme complete.

    Returns
    -------
    B : Matrix
        DESCRIPTION : Matrice tridiag de forme reduite

    """
    n,m=np.shape(A) # recupère la forme de matrice
    b= np.diag(A,-1) # recupère la sous diagonale en forme de vecteur
    c= np.diag(A,1) # recupère la  diagonale en forme de vecteur
    a= np.diag(A,0) # recupère la sur diagonale en forme de vecteur
   
    B = np.zeros((n,3))
    B[1:,0] = b # place la sous diagonale dans la première colonne 
    B[:,1] = a # place la  diagonale dans la deuxième colonne 
    B[:n-1,2] = c # place la sur diagonale dans la troisième colonne 
    return B 



def TRIcomplete (B) :
    """
    

    Parameters
    ----------
    B : Matrix
        DESCRIPTION : Matrice tridiag de forme reduite 

    Returns
    -------
    A : Matrix 
        DESCRIPTION : Matrice tridiag de forme complete

    """
    n,m=np.shape(B)
    b= np.diag(B[1:,0],-1) # créé une matrice ayant pour sous diagonale la première colonne de la matrice réduite
    c= np.diag(B[:n-1,2],1) # créé une matrice ayant pour sur diagonale la troisième colonne de la matrice réduite
    a= np.diag(B[:,1]) # créé une matrice ayant pour  diagonale la deuxième colonne colonne de la matrice réduite
    A=a+b+c # créé la matrice complète issue de la concaténation des trois matrices précédentes
    return A


def produitTRIvect (A,x) :
    """
    

    Parameters
    ----------
    A : Matrix
        DESCRIPTION : Matrice tridiag de forme reduite 
    x : Vector
        DESCRIPTION : vecteur quelconque

    Returns
    -------
    Ax : Vector
        DESCRIPTION : resultat du produit 

    """
    n1,m1=np.shape(A)
    n2=np.size(x)
    
    if n2 != n1 : # vérification de la possibilité du produit
        print(' produit impossible')
        return
    
    Ax=np.zeros((n2,1))
    
    for i in range (n1) :
        for j in range(m1):
            if i+j-1>=0 and i+j-1<=n1-1 :
                    Ax[i]+=A[i,j]*x[i+j-1] # calcul de chaque terme du vecteur (selon la définitioon formelle d'un produit matriciel)
                    
    return Ax


def DecompositionLU(A):  
    """
    

    Parameters
    ----------
    A : Matrix
        DESCRIPTION : Matrice tridiag de forme complete ( on peut mettre une matrice quelconque il 
                      faudra juste modifier return U et L au lieu de B)

    Returns
    -------
    B : Matrix
        DESCRIPTION :  Matrice contenent la sous diag de L ,la diag et la sur diag de U dans cet ordre et en colonne

    """
    n,n = A.shape
    U = np.array(A,float)
    L = np.eye(n)
    for i in range(n):
        if U[i,i]==0: # vérification de la compatibilité avec la méthode de Gauss
            print('Un pivot est nul')
            return U
        for j in range(i+1,n):
            g = U[j,i]/U[i,i] # calcul du coéfficient
            U[j,:]=U[j,:]-g*U[i,:] # changement de la ligne suivante 
            L[j,i]=g # insertion du coefficient dans L
    a=np.diag(L,-1) # extraction de la sous diagonale de L
    b=np.diag(U)  # extraction de la  diagonale de U
    c=np.diag(U,1) # extraction de la sur diagonale de U
    
    B = np.zeros((n,3))
    B[1:,0] = a # première colonne de M
    B[:,1] = b # deuxième colonne de M
    B[:n-1,2] = c # troisième colonne de M
    
    return B
 
 
def triLU (A) :
    """
    

    Parameters
    ----------
    A : Matrix
        DESCRIPTION : Matrice tridiag de forme reduite

    Returns
    -------
    M : Matrix
        DESCRIPTION :  Matrice contenent la sous diag de L ,la diag et la sur diag de U dans cet ordre et en colonne

    """
    n,m=np.shape(A)
    M=np.array(A,float) 

    for i in range (1,n) :
            g=A[i,0]/M[i-1,1] # calcul du coefficient 
            M[i,0]= g # sous diag de L
            M[i,1]=M[i,1]-g*M[i-1,2] # on fais l'opération sur la ligne suivante pour obtenir la diagonale de U
            # la sur diag est dejà mise quand on a copié A dans M au debut 
           
    
            
    return M

    
def TriLUResol(M,b):
    """
    

    Parameters
    ----------
    M : Matrix
        DESCRIPTION : Matrice contenent la sous diag de L ,la diag et la sur diag de U dans cet ordre et en colonne
    b : vector
        DESCRIPTION : parti de droite du système de cramel

    Returns
    -------
    x : Vector
        DESCRIPTION :  Solution du système

    """
    n = len(M)
    y = np.zeros((n,1)) # initialistion y
    x = np.zeros((n,1)) # initialistion x
    y[0] = b[0] # premeir terme de y
    for i in range(1,n):
        y[i] = b[i] - M[i,0]*y[i-1] # calcul des termes de y par redescente
    x[n-1] = y[n-1]/M[n-1,1] # Dernier terme de x
    for i in range(2,n+1):
        x[n-i] = (y[n-i]-x[n-i+1]*M[n-i,2]) / M[n-i,1] # calcul de x par remonter
    return x


def TriResol(A,b):
    """
    

    Parameters
    ----------
    A : Matrix
        DESCRIPTION : matrice tridiagonale sous forme reduite
    b : Vector
        DESCRIPTION : parti de droite du système de cramel

    Returns
    -------
    x : Vector
        DESCRIPTION : Solution du système

    """
    M = triLU(A) # decomposition LU
    x = TriLUResol(M,b) # resolution par redescente puis remonter
    return x







    
#---------------------test---------------------

# création d'un système tri diagonale
nb=4 #nombre de variable
Ar=np.random.random(size=(nb,3)) # matrice tridiagonale reduite
b = np.random.random(size=(nb,))

x=TriResol(Ar, b)

print('la solution trouver avec la méthode de frome reduite :' )
print(x)

print (b)
print(produitTRIvect(Ar, x)) # on utilise la fonction produit trivect pour multiplier 
# notre matrice reduite au solution trouver pour retrouver b









