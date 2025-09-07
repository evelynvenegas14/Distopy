from setuptools import setup

setup(
    name='DistoPy',
    version='0.1',
    packages=['DistoPy'],
    install_requires=[
        'sympy',
        'numpy'],
    description='Esta librería contiene las funciones necesarias para usar la teoría de '
    'distribuciones para la teoría Electrostática en Python de forma simbólica',
    author='Evelyn Venegas Agustín // Sherezade García Otero',
    author_email= 'evelyn_venegas@ciencias.unam.mx // shere.o@ciencias.unam.mx'
)