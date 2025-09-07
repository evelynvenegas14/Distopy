
from funcion_test_main import FuncionTest
import sympy as sp
from deltadirac_1d import DeltaDirac

class OperadoresDiferenciales:
    """
    Clase para representar un operador diferencial.
    """
    
    def __init__(self, operador):
        self.operador = operador

    def vector_normal(superficie_discontinuidad):
        """
        Calcula el vector normal a una superficie de discontinuidad. Recibe la
        parametrización de la superficie en cartesianas como función simbólica y devuelve 
        el vector ortogonal normalizado.
        """
        x,y,z=sp.symbols('x y z')

        N = sp.CoordSys3D('N') # Definimos un sistema de coordenadas 3D

        #el vector normal se calcula como
        #n= grad(phi)/||grad(phi)||

        # Calculamos el gradiente (vector normal)
        grad_sup = sp.gradient(superficie_discontinuidad, N)
        # Normalizamos el vector normal
        vector_normal = grad_sup / sp.sqrt(grad_sup.dot(grad_sup))
        return vector_normal #lo devuelve en sus coodenadas x,y,z
    
    def salto_funcion(self):
        """
        Calcula el salto [[f]] = f(x + eps*n) - f(x - eps*n)
        como función simbólica de (x, y, z), para funciones escalares o vectoriales.
        """
        eps = sp.Symbol('eps', real=True, positive=True)
        x, y, z = sp.symbols('x y z')
        pos = sp.Matrix([x, y, z])
        
        # Vector normal puede depender de x, y, z simbólicamente
        n = self.vector_normal  # vector normal unitario
        
        # Desplazamiento sobre la normal
        pos_plus = pos + eps * n
        pos_minus = pos - eps * n
        
        # Si la función es escalar
        if isinstance(self.funcion, (sp.Basic, sp.Expr)):
            f_plus = self.funcion.subs({x: pos_plus[0], y: pos_plus[1], z: pos_plus[2]})
            f_minus = self.funcion.subs({x: pos_minus[0], y: pos_minus[1], z: pos_minus[2]})
            salto = sp.simplify(sp.limit(f_plus, eps, 0, dir='+') - sp.limit(f_minus, eps, 0, dir='+'))
        
        # Si la función es vectorial (tipo Matrix)
        elif isinstance(self.funcion, sp.Matrix):
            salto = sp.Matrix([
                sp.simplify(
                    sp.limit(comp.subs({x: pos_plus[0], y: pos_plus[1], z: pos_plus[2]}), eps, 0, dir='+') -
                    sp.limit(comp.subs({x: pos_minus[0], y: pos_minus[1], z: pos_minus[2]}), eps, 0, dir='+')
                )
                for comp in self.funcion
            ])
        
        else:
            raise TypeError("La función debe ser escalar (Expr) o vectorial (Matrix)")
        
        return salto



    class Gradiente:
        """
        Calcula el gradiente de una función discontinua en sentido de distribuciones, usando
        la fórmula:
        nabla f = {nabla f} + n_out [[f]] delta_s.
        """
        #nabla f = {nabla f}+n_out [[f]]delta_s
        
        def __init__(self, funcion, vector_normal, salto_funcion):
            self.funcion = funcion
            self.vector_normal = vector_normal
            self.salto_funcion= salto_funcion
        
        def gradiente_sin_precaucion(self):
            """
            Calcula el gradiente de la función sin considerar la discontinuidad.
            """
            # Ordena los símbolos para consistencia
            variables = sorted(self.funcion.free_symbols, key=lambda s: s.name)
            return sp.Matrix([sp.diff(self.funcion, var) for var in variables])
        
        def gradiente_isod(self):
            """
            Junta el cálculo completo para tener el gradiente en sentido de distribuciones.
            """
            gradiente = self.gradiente_sin_precaucion() #gradiente sin considerar la discontinuidad
            salto = self.salto_funcion() if callable(self.saltofuncion) else self.salto_funcion
            #devuelve el gradiente con el término de discontinuidad
            return gradiente + self.vector_normal * salto * FuncionTest.delta_s()


    class Divergencia:
        """
        Calcula el gradiente de una función discontinua en sentido de distribuciones, usando
        la fórmula:
        nabla dot f = {nabla dot f} + n_out [[f]] delta_s.
        """
        
        def __init__(self, funcion, vector_normal, salto_funcion):
            self.funcion = funcion
            self.vector_normal = vector_normal
            self.salto_funcion = salto_funcion

        def divergencia_sin_precaucion(self):
            """
            Calcula la divergencia de un campo vectorial sin considerar la discontinuidad.
            """
            variables = sorted(self.funcion.free_symbols, key=lambda s: s.name)
            return sum(sp.diff(self.funcion[i], variables[i]) for i in range(len(variables)))
            
        def divergencia_isod(self):
            """
            Junta el cálculo completo para tener la divergencia en sentido de distribuciones.
            """
            divergencia = self.divergencia_sin_precaucion()
            salto = self.salto_funcion() if callable(self.salto_funcion) else self.salto_funcion
            return divergencia + self.vector_normal * salto * FuncionTest.delta_s()

    class Rotacional:
        """
        Clase para representar el operador rotacional en sentido de distribuciones, usando la fórmula:
        curl f = {curl f} + n_out x [[f]] delta_s.
        """
        
        def __init__(self, vector, vector_normal, salto_funcion):
            self.vector = vector
            self.vector_normal=vector_normal
            self.salto_funcion=salto_funcion

        def rotacional_sin_precaucion(self):
            """
            Calcula el rotacional de un campo vectorial sin considerar la discontinuidad.
            """
            x, y, z = sp.symbols('x y z')
            variables = [x, y, z]
            return sp.Matrix([
                sp.diff(self.vector[2], variables[1]) - sp.diff(self.vector[1], variables[2]),
                sp.diff(self.vector[0], variables[2]) - sp.diff(self.vector[2], variables[0]),
                sp.diff(self.vector[1], variables[0]) - sp.diff(self.vector[0], variables[1])
            ])
        def rotacional_isod(self):
            """
            Junta el cálculo completo para tener el rotacional en sentido de distribuciones.
            """
            rotacional = self.rotacional_sin_precaucion()
            salto = self.salto_funcion() if callable(self.salto_funcion) else self.salto_funcion
            return rotacional + self.vector_normal.cross(salto) * FuncionTest.delta_s()
        

    class Laplaciano:
        """
        Clase para representar el operador laplaciano en sentido de distribuciones, usando la fórmula
        nabla^2 f = {nabla^2 f} + [[\partial f / \partial n]] delta_s + nabla \dot (n [[f]] delta_s).
        """
        
        def __init__(self, funcion, vector_normal, salto_funcion, salto_derivada_normal):
            #va a recibir el salto de la derivada normal, no supe hacer que lo calcule
            self.funcion = funcion
            self.vector_normal = vector_normal
            self.salto_funcion = salto_funcion
            self.salto_derivada_normal = salto_derivada_normal

        def laplaciano_sin_precaucion(self):
            """
            Calcula el laplaciano de una función sin considerar la discontinuidad.
            """
            variables = sorted(self.funcion.free_symbols, key=lambda s: s.name)
            return sum(sp.diff(self.funcion, var, 2) for var in variables) #segunda derivada
        
        def laplaciano_isod(self):
            """
            Junta el cálculo completo para tener el laplaciano en sentido de distribuciones.
            """
            laplaciano = self.laplaciano_sin_precaucion()
            salto = self.salto_funcion() if callable(self.salto_funcion) else self.salto_funcion
            
            return (laplaciano + self.salto_derivada_normal * FuncionTest.delta_s() +
                    sp.div(self.vector_normal * salto * FuncionTest.delta_s()))

