
from scipy.misc import derivative
import sympy as sp
from sympy import Float, sympify
from scipy.misc import derivative
from sympy import symbols, Function, Derivative

class DeltaDirac:
    """
    Clase para representar la delta de Dirac en el sentido de distribuciones.
    """

    def __init__(self, f, x, x_0):
        """
        Inicializa la clase DeltaDirac.

        :param f: La función test.
        :param x: Un vector (lista de variables o variable unidimensional).
        :param x_0: Un vector (lista de valores simbólicos) o un valor único.
        """
        self.f = f
        self.x = x
        self.x_0 = x_0

    def accion(self):
        """
        Evalúa en el punto x_0.

        Returns:
            sympy.Expr: Resultado de la evaluación de manera simbólica.
        """
        f = self.f  # Copia de la función test para trabajar
        return f.subs(self.x, self.x_0)


    def derivada(self, n):
        """
        Calcula la n-ésima derivada de la delta de Dirac en sentido de distribuciones
        La función test tiene que ser clase C1 y tener derivada analítica para esta función.

        :param n: Orden de la derivada.
        :return: La n-ésima derivada de la delta de Dirac en acción con la función test.
        """
        if n == 0:
            return self.accion()
        else:
            try:
                derivada_simb = self.f.diff(self.x, n).subs(self.x, self.x_0)
                return (-1)**n * derivada_simb
            except:
                def func(x_val):
                    return self.f.subs(self.x, x_val)

                derivada_num = derivative(func, float(self.x_0), dx=1e-5, n=n)
                resultado_simbolico = Float((-1)**n * derivada_num) #conversion a float, pero el de sympy!
                return resultado_simbolico
            
    class DerivadaDiscontinua:
        """
        Clase para calcular la n-ésima derivada de la delta de Dirac en sentido de distribuciones
        para el caso donde la función test es discontinua en x_0.
        """

        def __init__(self, f, x, discontinuidades):
            """
            Inicializa la clase DerivadaDiscontinua.

            :param f: La función test discontinua.
            :param x: Un vector (lista de variables o variable unidimensional).
            :param x_0: Un vector (lista de valores simbólicos) o un valor único.
            """
            self.f = f
            self.x = x
            self.discontinuidades = discontinuidades

        def salto_funcion(self, x_0):
            """
            Calcula el salto de la función [[f]] = f(x_0^+) - f(x_0^-), para una sola
            discontinuidad puntual en x_0
            Returns:
                sympy.Expr: El salto de la función en x_0.
            """
            # Límite por la derecha (x_0+)
            limite_derecha = sp.limit(self.f, self.x, x_0, dir='+')
            # Límite por la izquierda (x_0-)
            limite_izquierda = sp.limit(self.f, self.x, x_0, dir='-')
            # Salto: [[f]] = f(x_0^+) - f(x_0^-)
            return limite_derecha - limite_izquierda

        def derivada_sin_precaucion_n(self, n):
            """
            Calcula la derivada n-ésima de la función simbólica f con respecto a x 
            de forma recursiva.
            Intenta hacerlo simbólicamente; si falla, lo hace numéricamente.

            Args:
                n (int): orden de la derivada.

            Returns:
                sympy.Expr o float: expresión simbólica o aproximación numérica.
            """
            if 0 <= n <= 2:
                return self.f #regresa la función original
            try:
                derivada_simbolica = sp.diff(self.f, self.x, n)
                return derivada_simbolica #deriva las veces necesarias usando la función diff
            except Exception:
                def func(x_val):
                    return float(self.f.subs(self.x, x_val)) #convierte a float para evitar problemas con sympy
                return derivative(func, self.x_0, dx=1e-5, n=n) #tolerancia de 1e-5
            else:
                raise NotImplementedError("El orden de la derivada debe ser 0, 1 o 2 para este método.")


        def derivada_discontinua_1(self, n, x_0):
            """
            Calcula la n-ésima derivada de la delta de Dirac en sentido de distribuciones
            para el caso donde la función test es discontinua en x_0, con una sola discontinuidad.

            Por ahora solo funciona para n=0, 1, 2

            :param n: Orden de la derivada.
            :param x_0: Punto de discontinuidad.
            :return: Expresión simbólica de la derivada distribuida.
            """
            delta = Function('δ')  # símbolo de la delta de Dirac
            if n == 0:
                return self.accion()
            
            elif n == 1:
                # f' = {f'} + [[f]]·δ(x₀)
                derivada_regular = self.derivada_sin_precaucion(1)
                salto = self.salto_funcion(x_0)
                termino_delta = salto * delta(x_0)
                return derivada_regular + termino_delta

            elif n == 2:
                # f'' = {f''} + [[f']]·δ(x₀) + [[f]]·δ'(x₀)
                derivada_regular = self.derivada_sin_precaucion(2)
                salto_derivada = self.salto_funcion(x_0)
                salto_funcion = self.salto_funcion(x_0)  # salto de la función en x_0

                delta_prima = DeltaDirac.derivada(delta(x_0), self.x)  # derivada simbólica de la delta
                
                termino_delta = salto_derivada * delta(x_0)
                termino_delta_prima = salto_funcion * delta_prima

                return derivada_regular + termino_delta + termino_delta_prima

            else:
                raise NotImplementedError("Solo se implementa para n=0, 1 o 2 en este momento.")
            
        def derivada_discontinua_n(self, n):
            """
            Calcula la primera derivada de la delta de Dirac en sentido de distribuciones para
            el caso donde la función test es discontinua en x_0_i, con múltiples discontinuidades.

            :param n: Orden de la derivada.
            :return: La primera derivada de la función en sentido de distribuciones.
            """
            #f' = {f'} + Σ [[f]](a_i)·δ(x₀ᵢ)

            delta = Function('δ')  # símbolo de la delta de Dirac
            derivada_regular = self.derivada_sin_precaucion_n(1)
            suma_terminos = 0
            for x_0 in self.discontinuidades:
                salto = self.salto_funcion(self.f, x_0)
                termino_delta = salto * delta(x_0)
                suma_terminos += termino_delta
            
            return derivada_regular + suma_terminos


            
