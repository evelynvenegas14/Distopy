
from Módulos.deltadirac_1d import DeltaDirac1D
from sympy import symbols, Abs, sqrt, Expr, lambdify


class FuncionTestCartesiana:
    """
    Función test en coordenadas cartesianas para el caso de una distribución de carga puntual.
    En el contexto de la electrostática, la función test principal es la que usamos 
    para medir la distancia desde una distribución de carga hasta el observador. Es decir,

    :math:`f(x)=1/|x - x_0|`

    donde :math:`x_0` es la posición de la distribución de carga y :math:`x` la posición del
    observador.
    Args:
                x0 (float, array o symbol): posición de la distribución de carga, puede recibir
                un valor numérico o simbólico.

    Returns:
                function: Una función que representa la distancia desde el punto x0.
    """
    def __init__(self, x0):#inicializa clase
        self.x0 = x0
        self.expr, self.vars = self._build_expression(x0)

    def build(self, x0):
        """
        Construye la expresión simbólica de la función test en función de la posición x0.
        Args:
            x0 (float, array o symbol): Posición de la distribución de carga

        Returns:
            tuple: Una tupla que contiene la expresión simbólica y las variables utilizadas.
        """
        if isinstance(x0, (int, float)): #este if verifica la dimensión de x0
            #caso unidimensional numérico
            x = symbols('x')
            expr = 1 / Abs(x - x0)
            return expr, (x,)

        elif isinstance(x0, Expr):
            #caso unidimensional simbólico
            x = symbols('x')
            expr = 1 / Abs(x - x0)
            return expr, (x,)

        elif isinstance(x0, (list, tuple)):
            #caso multidimensional
            n = len(x0) #dimensión 2 o 3, no más
            if n == 2:
                x, y = symbols('x y')
                expr = 1 / sqrt((x - x0[0])**2 + (y - x0[1])**2)
                return expr, (x, y)
            elif n == 3:
                x, y, z = symbols('x y z')
                expr = 1 / sqrt((x - x0[0])**2 + (y - x0[1])**2 + (z - x0[2])**2)
                return expr, (x, y, z)
            else:
                raise ValueError("Solo se permiten 1, 2 o 3 dimensiones.")
        else:
            raise TypeError("Tipo no válido para x0.")

    def simbolico(self): 
        """
        Devuelve la expresión simbólica de la función test.
        Returns:
            sympy.Expr: La expresión simbólica de la función test.
        """
        return self.expr

    def evaluar(self, **kwargs):
        """
        Evalúa la función test con los valores proporcionados.
        Args:
            **kwargs: Argumentos clave-valor donde las claves son las variables de la función test
                      y los valores son los valores a evaluar.
        Returns:
            float: El resultado de la evaluación de la función test.
        """
        return self.expr.subs(kwargs).evalf()

    def como_funcion(self):
        """
        Convierte la expresión simbólica de la función test en una función de Python.
        Returns:
            function: Una función anónima que toma los valores de las variables 
                      y devuelve el resultado de la evaluación de la función test.
        """
        return lambdify(self.vars, self.expr)

#ejemplos de uso  
f = FuncionTestCartesiana([3,8,5]) #crea el objetote como tal

#y dsp llamas a cada metodo, q ya esta hecho
print(f.symbolic())              # ➜ 1/|x - 3|
print(f.evaluate(x=5))           # ➜ 0.5
func = f.as_function()
print(func(3,1,4))                  # ➜ 1/7 = 0.142857...

