# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 15:24:22 2020.

@author: Alejandro
"""

from itertools import combinations, chain
import networkx as nx
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay, Voronoi, voronoi_plot_2d
import matplotlib.colors
import numpy as np
from numpy.linalg import matrix_rank
import imageio
import sympy as sy

t = sy.symbols('t', real=True)


def puntosCurvaRuido(curva, t, t0, t1, numPuntos=10, mu=0, sigma=0.1):
    """
    Obtener conjunto discretos de puntos de una curva con ruido.

    curva: list.
    t: t sympy symbol
    t0: float
        Inicio intervalo.
    t1: float
        Final intervalo.
    numPuntos: int. Por defecto 10.
    mu: float. Por defecto 0
        Media para la distribución normal.
    sigma: float. Por defecto 0.1
        Desviación típica para la distribución normal.
    """
    valores = np.linspace(t0, t1, num=numPuntos)
    puntosCurva = np.array([[x.subs(t, v) for x in curva] for v in valores], dtype=np.float64)
    ruido = np.random.normal(mu, sigma, [numPuntos, len(curva)])

    return puntosCurva + ruido


def pivotar(M, k, m):
    """
    Pivota en el elemento (k,m) intercambiando filas y columnas.

    M: np.array.
    k: int.
    m: int.
    """
    if M[k, m] != 1:
        encontrado = False
        i = k
        j = m + 1
        while not encontrado and i < M.shape[0] and j < M.shape[1]:
            if M[i, j] == 1:
                # Intercambio columnas
                M[:, [m, j]] = M[:, [j, m]]
                # Intercambio filas
                M[[k, i], :] = M[[i, k], :]
                encontrado = True

            j += 1
            if j == M.shape[1]:
                j = m
                i += 1
    else:
        encontrado = True

    return encontrado


def sumaFilZ2(M, i, j):
    """
    Suma: filai + filaj (sobre la j).

    M: np.array.
    i: int.
    j: int.
    """
    M[j, :] = (M[i, :] + M[j, :]) % 2


def sumaColZ2(M, i, j):
    """
    Suma: coli + colj (sobre la j).

    M: np.array.
    i: int.
    j: int.
    """
    M[:, j] = (M[:, i] + M[:, j]) % 2


def normSmithZ2(M):
    """
    Obtener la forma normal de Smith de una matriz con coefs en Z2.

    M: np.array.
    """
    n = 0
    cols = M.shape[1]
    fils = M.shape[0]
    while n < cols and n < fils and pivotar(M, n, n):
        # Recorrer fila
        for j in range(n + 1, cols):
            if M[n, j] == 1:
                sumaColZ2(M, n, j)
        # Recorrer columna
        for i in range(n + 1, fils):
            if M[i, n] == 1:
                sumaFilZ2(M, n, i)
        n += 1

    return M


def powerset(iterable):
    """
    Optiene un chain con todos los subconjuntos del iterable.

    iterable: iterable.
    """
    # "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def ordCaras(cara):
    """
    Relacion de orden de las caras de una filtracion.

    cara: tuple().
    """
    return (cara[1], len(cara[0]) - 1, cara[0])


def distancia(p1, p2):
    """
    Distancia euclídea entre los puntos p1 y p2.

    p1: list.
    p2: list.
    """
    p1Np = np.array(p1)
    p2Np = np.array(p2)
    return np.sqrt(np.dot(p1Np - p2Np, p1Np - p2Np))


def radioCircunscrita(p1, p2, p3):
    """
    Dados los vertices de un triangulo obtenemos el radio de la cirunferencia circuncentra.

    p1: tuple.
    p2: tuple.
    p3: tuple.
    """
    a = distancia(p1, p2)
    b = distancia(p1, p3)
    c = distancia(p2, p3)

    s = (a + b + c) / 2

    return (a * b * c) / (4 * np.sqrt(s * (s - a) * (s - b) * (s - c)))


def analisisComplejo(comp, simplice):
    """
    Alisis de las propiedades del complejo simplicial.

    comp: Complejo.
    simplice: set(tuple)
        simplice como ref para ejemplos.
    """
    # Todas las caras del complejo
    print(f"Todas las caras: {comp.getCaras()}")

    # Dimension del complejo
    dimComp = comp.dim()
    print(f"Dimension: {dimComp}")

    # Caras de cierta dimension
    for n in range(dimComp + 1):
        print(f"Todas las caras de dim {n}: {comp.getCarasN(n)}")

    # Característica de Euler
    print(f"Característica de Euler: {comp.caractEuler()}")

    # Estrella
    print(f"Estrella de {simplice}: {comp.st(simplice)}")

    # Link
    print(f"Link de {simplice}: {comp.lk(simplice)}")

    # Componentes conexas
    print(f"Numero de componentes conexas: {comp.compConexas()}")

    # 1-esqueleto
    print(f"El 1-esqueleto es: {comp.k_esqueleto(1)}")


def drawVor(puntos):
    """
    Representacion de las celdas de Voronoi de una nube de puntos.

    puntos: np.array.
    """
    vor = Voronoi(puntos)
    voronoi_plot_2d(vor, show_vertices=False, line_width=2,
                    line_colors='blue', line_alpha=0.6)
    plt.plot(puntos[:, 0], puntos[:, 1], 'ko')
    return vor


def delaunay(puntos):
    """
    Generar el triangulacion de Delaunay y su representacion junto a las celdas de Voronoi.

    puntos: np.array.
    """
    drawVor(puntos)

    Del = Delaunay(puntos)
    c = np.ones(len(puntos))
    cmap = matplotlib.colors.ListedColormap("limegreen")
    plt.tripcolor(puntos[:, 0], puntos[:, 1], Del.simplices, c, edgecolor="k", lw=2,
                  cmap=cmap)
    plt.plot(puntos[:, 0], puntos[:, 1], 'ko')
    plt.show()

    return Complejo([tuple(sorted(triangulo)) for triangulo in Del.simplices])


def alfaComplejo(puntos):
    """
    Genera la filtracion de alfa complejos de la triangulacion de Delaunay.

    puntos: np.array.
    """
    Del = delaunay(puntos)
    # Introducimos los 0-simplices
    alfa = Complejo(Del.getCarasN(0))

    # Introducimos los 2-simplices
    traingNuev = ((t, radioCircunscrita(puntos[t[0]], puntos[t[1]], puntos[t[2]]))
                  for t in Del.getCarasN(2))

    for t in traingNuev:
        alfa.setCaras([t[0]], t[1])

    # Introducimos los 1-simplices
    for arista in Del.getCarasN(1):
        # print(arista)
        p1 = puntos[arista[0]]
        p2 = puntos[arista[1]]
        pMedio = ((p2[0] + p1[0]) / 2, (p2[1] + p1[1]) / 2)
        d = distancia(p1, p2) / 2

        pesoTriangMin = -1
        for triang in Del.getCarasN(2):
            # print(arista, triang)
            difTriangArista = set(triang) - set(arista)

            if len(difTriangArista) == 1 and distancia(puntos[difTriangArista.pop()], pMedio) < d:
                pesoTriang = alfa.umbral(triang)
                if pesoTriangMin < 0 or pesoTriang < pesoTriangMin:
                    pesoTriangMin = pesoTriang

        alfa.setCaras([arista], d if pesoTriangMin < 0 else pesoTriangMin)

    return alfa


def plotalpha(puntos, K):
    """
    Representar el alpha complejo del complejo K.

    puntos: np.array.
    K: Complejo.
    """
    dim = K.dim()

    if dim > 1:
        c = np.ones(len(puntos))
        cmap = matplotlib.colors.ListedColormap("limegreen")
        plt.tripcolor(puntos[:, 0], puntos[:, 1], list(K.getCarasN(2)), c, edgecolor="k", lw=2,
                      cmap=cmap)

    plt.plot(puntos[:, 0], puntos[:, 1], 'ko')

    if dim > 0:
        for arista in K.getCarasN(1):
            p1 = puntos[arista[0]]
            p2 = puntos[arista[1]]
            plt.plot([p1[0], p2[0]], [p1[1], p2[1]], 'k')

    # plt.show()


def vietorisRips(puntos):
    """
    Calculo del complejo Vietoris Rips de una nube de puntos.

    puntos: np.array.
    """
    nsimplex = Complejo([tuple(range(len(puntos)))])
    VR = Complejo(list(nsimplex.getCarasN(0)))

    for arista in nsimplex.getCarasN(1):
        VR.setCaras([arista], 0.5 * distancia(puntos[arista[0]], puntos[arista[1]]))

    for i in range(2, len(puntos)):
        for simplex in nsimplex.getCarasN(i):
            lista = []

            for arista in combinations(simplex, 2):
                lista.append(0.5 * distancia(puntos[arista[0]], puntos[arista[1]]))
            VR.insert([simplex], max(lista))

    return VR


class Complejo():
    """Clase del complejo simplicial."""

    def __init__(self, carasMaximales=[]):
        """
        Complejo simplicial abstracto a partir de sus caras maximales.

        carasMaximales: list(tuple). Por defecto [].
        """
        # Concatenamos el los conjuntos obtenidos de cada cara maximal
        self.caras = set()
        for cara in carasMaximales:
            if cara not in self.caras:
                self.caras |= set(powerset(cara))
        # Quitamos el conjunto vacio
        self.caras -= {()}

        # Añadimos peso
        self.caras = set([(cara, 0.0) for cara in self.caras])

        self.carasOrd = sorted(list(self.caras), key=ordCaras)

    def setCaras(self, carasNuevas, peso=0.0):
        """
        Insertar nuevas caras y sus correspondientes subconjuntos con un peso dado.

        carasNuevas: list(tuple).
        peso: float. Por defecto 0.0.
        """
        for cara in carasNuevas:
            for caraGen in powerset(cara):
                if caraGen == tuple():
                    continue

                encontrado = False
                for caraAnt in self.caras:
                    if caraGen == caraAnt[0]:
                        encontrado = True
                        if caraAnt[1] > peso:
                            self.caras -= {caraAnt}
                            self.caras |= {(caraGen, peso)}

                        break

                if not encontrado:
                    self.caras |= {(caraGen, peso)}

        self.carasOrd = sorted(self.caras, key=ordCaras)

    def getCaras(self):
        """Devuelve el conjunto de todas las caras del complejo simplicial."""
        return set([cara[0] for cara in self.caras])

    def getCarasOrd(self):
        """Devuelve el conjunto de las caras ordenadas segun su filtracion."""
        return [cara[0] for cara in self.carasOrd]

    def umbrales(self):
        """Devuelve el conjunto de las umbrales ordenados segun la filtracion."""
        return list(dict.fromkeys([cara[1] for cara in self.carasOrd]))

    def umbral(self, cara):
        """
        Obtiene el umbral de una cara dada.

        cara: tuple.
        """
        index = 0
        encontrado = False
        while index < len(self.carasOrd) and not encontrado:
            encontrado = self.carasOrd[index][0] == cara
            index += int(not encontrado)

        return self.carasOrd[index][1] if encontrado else None

    def dim(self):
        """Devuelve la dimensión del complejo simplicial."""
        return max([len(caras[0]) for caras in self.caras]) - 1

    def getCarasN(self, dimension):
        """
        Devuelve el conjunto de todas las caras de dimension dada.

        dimension: int.
        """
        return set(c for c in self.getCaras() if len(c) == dimension + 1)

    def st(self, v):
        """
        Calcular la estrella del simplice v.

        v: set.
        """
        return set(c for c in self.getCaras() if v.issubset(c))

    def lk(self, v):
        """
        Calcular el de un simplice v.

        v: set.
        """
        # Calculamos la estrella de v
        st = self.st(v)

        # Calculamos la estrella cerrada de v
        st_ = set()
        for cara in st:
            if cara not in st_:
                st_ |= set(powerset(cara))
        # Quitamos el conjunto vacio
        st_ -= {()}

        # Devolvemos el link de v
        return st_ - st

    def compConexas(self):
        """Comprobar la conexion de un complejo simplicial."""
        # Para ello comprobamos que su 1-esqueleto sea conexo
        k1Graph = nx.Graph()
        k1Graph.add_nodes_from([vertice[0] for vertice in self.getCarasN(0)])
        k1Graph.add_edges_from(self.getCarasN(1))
        return nx.number_connected_components(k1Graph)

    def k_esqueleto(self, k):
        """
        Calcular el k-esqueleto de un complejo simplicial.

        k: int.
        """
        return set(c for c in self.getCaras() if len(c) <= k + 1)

    def drawK1(self):
        """Representación gráfica del 1-esqueleto."""
        k1Graph = nx.Graph()
        k1Graph.add_nodes_from([vertice[0] for vertice in self.getCarasN(0)])
        k1Graph.add_edges_from(self.getCarasN(1))
        plt.figure().add_subplot(111)
        nx.draw_networkx(k1Graph, with_labels=True)

    def caractEuler(self):
        """Obtención de la característica de Euler."""
        return sum([(-1)**k * len(self.getCarasN(k)) for k in range(self.dim() + 1)])

    def filtracion(self, a):
        """
        Obtener las caras con peso menor o igual que un valor.

        a: float.
        """
        i = 0
        caras = list()
        while i < len(self.carasOrd) and self.carasOrd[i][1] <= a:
            caras.append(self.carasOrd[i])
            i += 1

        result = Complejo()
        for cara, peso in caras:
            result.setCaras([cara], peso)

        return result

    def matrizBorde(self, p):
        """
        Calculo de la matriz borde de dimensión dada.

        p: int.
        """
        if p < 0:
            return None

        carasP = sorted(list(self.getCarasN(p)))

        if p == 0:
            m = np.zeros((1, len(carasP)), dtype=int)
        else:
            carasP_1 = sorted(list(self.getCarasN(p - 1)))
            d = self.dim()
            if p == d + 1:
                m = np.zeros((len(carasP_1), 1), dtype=int)
            elif p > d:
                m = None
            else:
                m = np.zeros((len(carasP_1), len(carasP)), dtype=int)
                for j in range(len(carasP)):
                    caraP = set(carasP[j])
                    for i in range(len(carasP_1)):
                        m[i, j] = int(set(carasP_1[i]).issubset(caraP))

        return m

    def betti(self, p):
        """
        Calculo del número de p de Betti.

        p: int.
        """
        if p == 0:
            Zp = len(self.getCarasN(0))
        else:
            Mp = normSmithZ2(self.matrizBorde(p))
            Zp = Mp.shape[1] - matrix_rank(Mp)

        Bp = matrix_rank(normSmithZ2(self.matrizBorde(p + 1)))

        return Zp - Bp

    def allBettis(self):
        """Calculo de todos los números de Betti."""
        Zps = np.array([len(self.getCarasN(0))], dtype=int)
        Bps = np.array([], dtype=int)
        d = self.dim()
        for p in range(1, d + 2):
            Mp = normSmithZ2(self.matrizBorde(p))

            if p <= d:
                Zps = np.append(Zps, Mp.shape[1] - matrix_rank(Mp))

            Bps = np.append(Bps, matrix_rank(Mp))

        return list(Zps - Bps)

    def __str__(self):
        """El toString del complejo."""
        return "Caras: " + str(self.caras)


if __name__ == "__main__":
    """
    comp1 = Complejo([(0, 1, 2, 3)])
    print("-------------COMP1-------------")
    analisisComplejo(comp1, set((0, 1)))

    comp2 = Complejo(list(comp1.k_esqueleto(2)))
    print("\n-------------COMP2-------------")
    analisisComplejo(comp2, set((0,)))

    comp3 = Complejo([(0, 1), (1, 2, 3, 4), (4, 5), (5, 6), (4, 6), (6, 7, 8), (8, 9)])
    print("\n-------------COMP3-------------")
    analisisComplejo(comp3, set((4,)))

    comp4 = Complejo(list(comp3.k_esqueleto(1)))
    print("\n-------------COMP4-------------")
    analisisComplejo(comp4, set((4,)))

    comp5 = Complejo([(0, 1, 2), (2, 3), (3, 4)])
    print("\n-------------COMP5-------------")
    analisisComplejo(comp5, set((2,)))

    comp6 = Complejo([(1, 2, 4), (1, 3, 6), (1, 4, 6), (2, 3, 5), (2, 4, 5), (3, 5, 6)])
    print("\n-------------COMP6-------------")
    analisisComplejo(comp6, set((1, 4)))

    comp7 = Complejo(list(comp6.k_esqueleto(1)))
    print("\n-------------COMP7-------------")
    analisisComplejo(comp7, set((1, 4)))

    comp8 = Complejo([(1, 2, 4), (2, 4, 5), (2, 3, 5), (3, 5, 6), (1, 3, 6), (1, 4, 6),
                      (4, 5, 7), (5, 7, 8), (5, 6, 8), (6, 8, 9), (4, 6, 9), (4, 7, 9),
                      (1, 7, 8), (1, 2, 8), (2, 8, 9), (2, 3, 9), (3, 7, 9), (1, 3, 7)])
    print("\n-------------COMP8-------------")
    analisisComplejo(comp8, set((1,)))

    comp9 = Complejo(list(comp8.k_esqueleto(1)))
    print("\n-------------COMP9-------------")
    analisisComplejo(comp9, set((1,)))

    comp10 = Complejo([(1, 2, 6), (2, 3, 4), (1, 3, 4), (1, 2, 5), (2, 3, 5), (1, 3, 6),
                       (2, 4, 6), (1, 4, 5), (3, 5, 6), (4, 5, 6)])
    print("\n-------------COMP10-------------")
    analisisComplejo(comp10, set((1,)))

    comp11 = Complejo(list(comp10.k_esqueleto(1)))
    print("\n-------------COMP11-------------")
    analisisComplejo(comp11, set((1,)))

    comp12 = Complejo([(0,), (1,), (2, 3), (4, 5), (5, 6), (4, 6), (6, 7, 8, 9)])
    print("\n-------------COMP12-------------")
    analisisComplejo(comp12, set((6,)))

    # Ejemplo filtracion
    print("\n-------------COMP13-------------")
    comp13 = Complejo()
    comp13.setCaras([(0, 1)], 1.0)

    comp13.setCaras([(1, 2), (2, 3), (2, 4)], 2.0)
    comp13.setCaras([(3, 4)], 3.0)
    comp13.setCaras([(2, 3, 4)], 4.0)

    # Todas las caras del complejo
    print(f"Todas las caras: {comp13.getCaras()}")

    # Umbral
    print(f"Umbral de {{3}}: {comp13.umbral((3,))}")

    # Filtraciones
    K1 = comp13.filtracion(1.0)
    K2 = comp13.filtracion(2.0)
    K3 = comp13.filtracion(3.0)
    K4 = comp13.filtracion(4.0)

    # Todas las caras de las filtraciones
    print(f"Todas las caras de K1: {K1.getCaras()}")
    print(f"Todas las caras de K2: {K2.getCaras()}")
    print(f"Todas las caras de K3: {K3.getCaras()}")
    print(f"Todas las caras de K4: {K4.getCaras()}")

    # Caras ordenedas por filtracion
    print(f"Caras ordenadas segun las filtraciones: {comp13.getCarasOrd()}")

    """

    """
    points = np.array([(0.38021546727456423, 0.46419202339598786),
                       (0.7951628297672293, 0.49263630135869474),
                       (0.566623772375203, 0.038325621649018426),
                       (0.3369306814864865, 0.7103735061134965),
                       (0.08272837815822842, 0.2263273314352896),
                       (0.5180166301873989, 0.6271769943824689),
                       (0.33691411899985035, 0.8402045183219995),
                       (0.33244488399729255, 0.4524636520475205),
                       (0.11778991601260325, 0.6657734204021165),
                       (0.9384303415747769, 0.2313873874340855)])

    points = np.array([[0.8957641450573793, 0.2950833519989374],
                       [0.028621391963087994, 0.9440875759025237],
                       [0.517621505875702, 0.1236620161847416],
                       [0.7871047164191424, 0.7777474116014623],
                       [0.21869796914805273, 0.7233589914276723],
                       [0.9891035292480995, 0.6032186214942837],
                       [0.30113764052453484, 0.613321425324272],
                       [0.18407448222466916, 0.7868606964403773],
                       [0.4496777667376678, 0.874366215574117],
                       [0.08225571534539433, 0.616710205071694]])

    curva1 = [4 * sy.sin(t), 9 * sy.cos(t)]
    curva2 = [1 + 3 * t**2, t**3 - 2 * t]
    points = puntosCurvaRuido(curva2, t, -2, 2, numPuntos=30)
    print(points)

    plt.plot(points[:, 0], points[:, 1], 'ko')
    plt.show()

    vor = drawVor(points)

    delaunay(points)

    alpha = alfaComplejo(points)
    print(alpha)
    i = 0
    images = []
    for valor in alpha.umbrales():
        # print(valor)
        K = alpha.filtracion(valor)
        fig = voronoi_plot_2d(vor, show_vertices=False, line_width=2, line_colors='blue', lines_alpha=0.6)
        plotalpha(points, K)
        plt.title("Radio: " + str(valor))
        fig.savefig(f"imgTemp/im{i}.png")
        images.append(imageio.imread(f"imgTemp/im{i}.png"))
        i += 1
        plt.show()

    imageio.mimsave('alphaGif/alpha.gif', images)
    """
