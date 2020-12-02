# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
from setuptools import setup

setup(
    name='biorefineries',
    packages=['biorefineries'],
    license='MIT',
    version='2.16.5',
    description="Biorefinery models in BioSTEAM",
    long_description=open('README.rst').read(),
    author='Yoel Cortes-Pena',
    install_requires=['biosteam>=2.21.8'],
    python_requires=">=3.6",
    package_data=
        {'biorefineries': ['biorefineries/*',
                           'lipidcane/*', 
                           'lipidcane/utils/*', 
                           'cornstover/*', 
                           'corn/*',
                           'cornstover/_humbird2011.xlsx',
                           'sugarcane/*',
                           'fattyalcohols/*',
                           'fattyalcohols/units/*',
                           'LAOs/*',
                           'LAOs/units/*',
                           'tests/*',
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