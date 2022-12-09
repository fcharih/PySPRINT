from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    name="pysprint",
    version="1.0",
    rust_extensions=[RustExtension("rsprint.rsprint", binding=Binding.PyO3)],
    packages=["rsprint"],
    # rust extensions are not zip safe, just like C-extensions.
    zip_safe=False,
    entry_points={
    'console_scripts': [
        'pysprint = rsprint.entrypoints:main'
     ],
    },
    install_requires=[
        "biopython",
        "numpy",
        "mpi4py",
        "matplotlib",
        "loguru",
        "pandas"
    ]
)
