# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 16:17:00 2017

@author: Yoel Cortes-Pena
"""
from setuptools import setup

setup(
    name='biorefineries',
    packages=['biorefineries'],
    license='MIT',
    version='2.8.0',
    description="Biorefinery models in BioSTEAM",
    long_description=open('README.rst').read(),
    author='Yoel Cortes-Pena',
    install_requires=['biosteam>=2.9.0',
                      'lazypkg==1.4'],
    python_requires=">=3.6",
    package_data=
        {'biorefineries': ['biorefineries/*',
                           'lipidcane/*', 
                           'lipidcane/utils/*', 
                           'cornstover/*', 
                           'cornstover/_humbird2011.xlsx',
                      ]},
    platforms=['Windows', 'Mac', 'Linux'],
    author_email='yoelcortes@gmail.com',
    url='https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/BioSTEAM 2.x.x/biorefineries',
    download_url='https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/BioSTEAM 2.x.x/biorefineries',
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