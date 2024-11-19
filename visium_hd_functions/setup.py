from setuptools import setup, find_packages

setup(
    name="visium_hd_functions",  # Name of the package
    version="0.1.0",  # Version of the package
    packages=find_packages(),  # Automatically find all submodules in this directory
    install_requires=[  # Any dependencies specific to this package
        "numpy",
        "pandas",
        "matplotlib",
        "shapely",
        "scanpy",
        "geopandas",
    ],
    description="A collection of functions for Visium HD data processing",
    author="Selim Romero",  # Replace with your name or organization
    author_email="ssromerogon@tamu.edu",  # Replace with your email
    url="https://github.com/sromero1992/visium_hd_processing",  # URL for your GitHub repo or package
)

