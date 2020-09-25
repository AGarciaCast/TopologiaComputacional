# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 15:24:22 2020

@author: Alejandro
"""

from itertools import combinations, chain
import networkx as nx
import matplotlib.pyplot as plt

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
        carasMaximales: list of tuples
        """
        self.carasMaximales = set(carasMaximales)
        
        #Concatenamos el los conjuntos obtenidos de cada cara maximal
        self.caras = set()
        for cara in carasMaximales:
            self.caras |= set(powerset(cara)) 
        #Quitamos el conjunto vacio
        self.caras -= {()}
        
        #Añadimos peso
        self.caras = set([(cara, 0.0) for cara in self.caras])
   
    def setCaras(self, carasNuevas, peso):
        """
        Insertar nuevas caras y sus correspondientes subconjuntos con un peso dado
        """
        #Concatenamos el los conjuntos obtenidos de cada cara maximal
        for cara in carasNuevas:
            #Quitamos las caras ya añadidas previamente para actualizar el peso
            self.caras -= set([caraAnt for caraAnt in self.caras if caraAnt[0] in powerset(cara)])
            #Añadimos las nuevas caras con su correspondiente peso
            self.caras |= set([(c, peso) for c in powerset(cara)]) 
        
        self.caras -= {((), peso)}
    
    def getCaras(self):
        """
        Devuelve el conjunto de todas las caras del complejo simplicial
        """
        return set([cara[0] for cara in self.caras])
    
    def dim(self):
        """
        Devuelve la dimensión del complejo simplicial
        """
        return max([len(caras) for caras in self.carasMaximales]) - 1
    
    def getCarasN(self, dimension):
        """
        Devuelve el conjunto de todas las caras de dimension dada
        dimension: int
        """
        return set(c for c in self.getCaras() if len(c) == dimension + 1)
    
    def st(self, v):
        """
        Calcular la estrella del simplice v
        v: set
        """
        return set(c for c in self.getCaras() if v.issubset(c))
    
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
        
    def isconnected(self):
        """
        Comprobar la conexion de un complejo simplicial
        """
        #Para ello comprobamos que su 1-esqueleto sea conexo
        k1Graph = nx.Graph()
        k1Graph.add_edges_from(self.getCarasN(1))
        return nx.is_connected(k1Graph)
    
    def drawK1(self):
        k1Graph = nx.Graph()
        k1Graph.add_edges_from(self.getCarasN(1))
        plt.figure().add_subplot(111)
        nx.draw_networkx(k1Graph, with_labels=True, pos= nx.kamada_kawai_layout(k1Graph))

    def caractEuler(self):
        return sum([(-1)**k * len(self.getCarasN(k)) for k in range(self.dim() + 1)])
            
    
    def __str__(self):
        return "Caras: " + str(self.caras)

if __name__ == "__main__":
    
    comp1 = Complejo([(0,1,2),(2,3),(3,4)])
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
    
    
    #Representacion 1-esqueleto
    print("1-esqueleto:")
    comp1.drawK1()
    
    #Insertar nuevas caras
    comp1.setCaras([(1,5,6,7), (8,9)], 0.1)
    print(comp1)
    comp1.drawK1()
    
    #Conexion
    print("Es conexo?: ", comp1.isconnected())
    
    #Característica de Euler
    print("Característica de Euler: ", comp1.caractEuler())
    
    