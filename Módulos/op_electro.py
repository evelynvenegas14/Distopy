
class OperadoresElectrostaticos:
    """
    Clase para representar operadores electrostáticos en un sistema de coordenadas cartesianas.
    """

    def __init__(self, operador):
        self.operador = operador

    class CampoElectrico:
        """
        Clase para representar el campo eléctrico.
        """
        
        def __init__(self, vector):
            self.vector = vector

    class PotencialElectrico:
        """
        Clase para representar el potencial eléctrico.
        """
        
        def __init__(self, funcion):
            self.funcion = funcion

    class FuerzaElectrica:
        """
        Clase para representar la fuerza eléctrica.
        """
        
        def __init__(self, vector):
            self.vector = vector

    class TorcaElectrica:
        """
        Clase para representar la torca eléctrica.
        """
        
        def __init__(self, vector):
            self.vector = vector