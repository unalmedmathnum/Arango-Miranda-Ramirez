from sympy.matrices import Matrix, eye
from sympy import roots, shape
from sympy.abc import (
    lamda
)  # Se usa 'lamda' porque lambda es una expresión reservada en Python


def metodo_polinomio_caracteristico(matriz, find_eigenvectors=False):
    """Lee una matriz cuadrada y retorna sus autovalores y el polinomio caracteristico asociado
    Params:
        matriz: matriz cuadrada que se debe expresar como un array de arrays
        e.g [[2, 1], [3, 4]].
        find_eigenvectors (opcional): si es True, intenta calcular el espacio
        nulo para hallar los autovectores asociados. No funciona con todas
        las matrices.
    Returns:
        diccionario que contiene:
        "matriz" matriz de Sympy ingresada.
        "autovalores" valores propios hallados.
        "autovectores" vectores propios hallados. Su posición está relacionada con la
            de los autovalores.
        "polinomio_caracteristico" polinomio asociado a la matriz en términos de lambda.
    """
    M = Matrix(matriz)  # Generamos una matriz con los datos ingresados
    n = shape(M)[0]  # n será la dim de M
    mi = M - (lamda * eye(n))  # mi es (A-lambdaI)
    mi_det = mi.det()  # hallar el determinante de mi para obtener el char. poly.
    mi_roots = roots(mi_det)  # hallar raíces de det(A-lambdaI) = 0
    mi_autovalor = []  # lista de autovalores
    mi_autovector = []  # lista de autovectores
    for eigenvalue in mi_roots.keys():
        mi_autovalor.append(eigenvalue)  # guardamos el valor del valor propio
        mi_1 = M - eigenvalue * eye(n)  # Reemplazamos iterativamente lambda en (A-lambdaI)x = 0 por cada uno
                                        # de nuestros valores propios
        if (find_eigenvectors == True): # Hallar el espacio nulo para encontrar el vector propio
            mi_autovector.append(mi_1.nullspace())  

    return {  # retorna la matriz, sus autovalores y vectores, y el polinomio característico
        "matriz": M,
        "autovalores": mi_autovalor,
        "autovectores": mi_autovector,
        "polinomio_caracteristico": mi_det,
    }
