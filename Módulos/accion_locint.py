
from scipy.integrate import nquad

def Accion_LocInt(T, f_test_, n, limites):
    """
    Calcula la acción de una distribución localmente integrable T sobre una función test f_test_,
    usando cuadratura gaussiana en n dimensiones de la biblioteca scipy.

    Args:
        T (function): La distribución localmente integrable, representada como una función matemática.
        f_test_ (function): La función test con soporte compacto.
        n (int): Número de dimensiones del espacio.
        limites (list): Lista de límites de integración para cada dimensión.
                        Ejemplo: [[xmin, xmax], [ymin, ymax], ...].
                        
    Returns:
        float: El resultado de la acción integral.
    """
    #definir la función producto de T y f_test_
    def integrando(*args):
        return T(*args) * f_test_(*args)

    #calcular la integral múltiple
    resultado, error = nquad(integrando, limites)
    return resultado, error