#!/usr/bin/env python
from setuptools import setup, find_packages

build_exe_options = {
    "includes": [
        "sys",
        "os",
    ],
    "excludes": [

    ],
}

setup(
    name="deliqc",
    version="1.0",
    description="DNA encoded library information quality control",
    url="https://github.com/Gillingham-Lab/deliqc",

    author="Basilius Sauter",
    author_email="basilius.sauter@unibas.ch",

    python_requires=">3.9.0",
    packages=find_packages(),

    install_requires=[
        "argh",
        "colorama",
        "numpy",
        "scipy",
        "BioPython",
        "matplotlib"
    ],

    entry_points={
        "console_scripts": [
            'deliqc=deliqc.cli.main:main',
        ]
    }
)
