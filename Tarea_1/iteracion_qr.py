# -*- coding: utf-8 -*-
"""iteracion_QR.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1MOuw4DeRtZRuMmqxjdBtp2hHX_xHTpJI
"""

import numpy as np
from numpy.linalg import qr

def QR_Iteracion(A, tol=1e-6, max_iter=1000):
    #Definimos la función con una matriz inicial cuadrada, una tolerancia estándar y un número máximo de iteraciones estándar

    # Esta función va a retornar a la matriz A después de las iteraciones y cuántas iteraciones fueron necesarias

    # Copia inicial de la matriz
    A_actual = np.array(A, dtype=float)

    for iteracion in range(max_iter):
        # Hacemos el primer paso de descomposición QR, incluido en el paquete de Python
        Q, R = np.linalg.qr(A_actual)

        # Actualizamos la matriz usando Q y R multiplicándolas en el orden inverso para generar la siguiente iteración
        A_nueva = R @ Q

        # Verificamos la norma de Froebenius de la matriz obtenida, norma que nos permite comparar todas las entradas de las matrices, y aumentamos la iteración
        if np.linalg.norm(A_nueva - A_actual, ord='fro') < tol:
            return A_nueva, iteracion + 1

        # Actualizar A para la siguiente iteración
        A_actual = A_nueva

    # Si no converge, devolver el resultado actual; si converge, devolver los elementos de la diagonal de la matriz de convergencia, el número de iteraciones y los valores y vectores propios de dicha matriz
    return np.diag(A_actual), max_iter, np.linalg.eig(A_actual)