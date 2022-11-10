from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    name="rsprint",
    version="1.0",
    rust_extensions=[RustExtension("rsprint.rsprint", binding=Binding.PyO3)],
    packages=["rsprint"],
    # rust extensions are not zip safe, just like C-extensions.
    zip_safe=False,
    entry_points={
    'console_scripts': [
        'rsprint = rsprint.entrypoints:main'
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