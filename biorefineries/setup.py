# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 16:17:00 2017

@author: Yoel Cortes-Pena
"""
from setuptools import setup
#from Cython.Build import cythonize
#import numpy

setup(
    name='biorefineries',
    packages=['biorefineries'],
    license='MIT',
    version='0.1.2',
    description="Biorefinery models in BioSTEAM",
    long_description=open('README.rst').read(),
    #ext_modules=cythonize('biosteam/equilibrium/unifac.pyx'),
    #include_dirs=[numpy.get_include()],
    author='Yoel Cortes-Pena',
    install_requires=['biosteam==1.0.7'],
    python_requires=">=3.6",
    package_data=
        {'biorefineries': ['biorefineries/*',
                           'lipidcane/*', 
                           'lipidcane/utils/*', 
                           'lipidcane/species/*',
                           'lipidcane/species/tripalmitin_liquid.xlsx', 
                           'lipidcane/system/*', 
                           'sugarcane/*',
                           'cornstover/*', 
                           'cornstover/_humbird2011.xlsx',
                      ]},
    platforms=['Windows', 'Mac', 'Linux'],
    author_email='yoelcortes@gmail.com',
    url='https://github.com/yoelcortes/Bioindustrial-Park/tree/master/biorefineries',
    download_url='https://github.com/yoelcortes/Bioindustrial-Park/tree/master/biorefineries',
    classifiers=['Development Status :: 3 - Alpha',
                 'Environment :: Console',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3.6',
                 'Programming Language :: Python :: 3.7',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Mathematics'],
    keywords='chemical process simmulation bioprocess engineering mass CABBI biorefinery biofuel bioproducts lipid-cane corn stover',
)