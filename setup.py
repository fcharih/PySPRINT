from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    name="sprint",
    version="1.0",
    rust_extensions=[RustExtension("sprint.sprint", binding=Binding.PyO3)],
    packages=["sprint"],
    # rust extensions are not zip safe, just like C-extensions.
    zip_safe=False,
    entry_points={
    'console_scripts': [
        'sprint = sprint.entrypoints:main'
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
