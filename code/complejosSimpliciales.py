# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 15:24:22 2020

@author: Alejandro
"""

from itertools import combinations, chain

def powerset(iterable):
    """
    Optiene un chain con todos los subconjuntos del iterable
    """
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

class Complejo():
    def __init__(self, carasMaximales):
        """
        Constructor del complejo simplicial abstracto a partir de sus caras maximales
        carasMaximales: set of tuples
        """
        #Concatenamos el los conjuntos obtenidos de cada cara maximal
        self.caras = set()
        for cara in carasMaximales:
            self.caras |= set(powerset(cara)) 
        #Quitamos el conjunto vacio
        self.caras -= {()}
   
    def __str__(self):
        return "Caras: " + str(self.caras)
    
    def getCaras(self):
        """
        Devuelve el conjunto de todas las caras del complejo simplicial
        """
        return self.caras
    
    def dim(self):
        """
        Devuelve la dimensión del complejo simplicial
        """
        return max([len(caras) for caras in self.caras]) - 1
    
    def getCarasN(self, dimension):
        """
        Devuelve el conjunto de todas las caras de dimension dada
        dimension: int
        """
        return set(c for c in self.caras if len(c) == dimension + 1)
    
    def st(self, v):
        """
        Calcular la estrella del simplice v
        v: set
        """
        return set(c for c in self.caras if v.issubset(c))
    
    def lk(self, v):
        """
        Calcular el de un simplice v
        v: set
        """
        #Calculamos la estrella de v
        st = self.st(v)
        
        #Calculamos la estrella cerrada de v
        st_ = set()
        for cara in st:
            if cara not in st_:
                st_ |= set(powerset(cara)) 
        #Quitamos el conjunto vacio
        st_ -= {()}
        
        #Devolvemos el link de v
        return st_ - st
        

if __name__ == "__main__":
    
    comp1 = Complejo([{0,1,2},{2,3},{3,4}])
    print(comp1)
    
    #Dimension del complejo
    dimComp1 = comp1.dim()
    print("Dimension:", dimComp1)
    
    #Todas las caras del complejo
    print("Todas las caras: ", comp1.getCaras())
    
    #Caras de cierta dimension
    for n in range(dimComp1 + 1):
        print("Todas las caras de dim", n , ": ", comp1.getCarasN(n))
    #Estrella
    print("Estrella de {2}: ", comp1.st({2}))
    
    #Link
    print("Link de {2}: ", comp1.lk({2}))
    
    
