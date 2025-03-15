import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()
DESCRIPTION = "A collection of functions for Visium HD data processing",

# Find packages inside the 'qfeatures' folder
PACKAGES = find_packages(where="visium_hd_functions")

KEYWORDS = [
    "Visium HD",
    "Nuclei segmentation",
    "Stardist",
    "Image recognition",
    "single-cell"                                                                                                                                                                               
]

INSTALL_REQUIRES=[  # Any dependencies specific to this package
        "numpy",
        "pandas",
        "matplotlib",
        "shapely",
        "scanpy",
        "geopandas",
        ] 

# Verbose output: print the list of packages found
print(f"Packages found: {PACKAGES}")

# Print the install_requires for visibility
print(f"Install requires: {INSTALL_REQUIRES}")


setup(
    name="visium_hd_functions",  # Name of the package
    version="0.1.0",  # Version of the package
    long_description=README,
    description=DESCRIPTION,
    author='Selim Romero',
    author_email="ssromerogon@tamu.edu",  # Replace with your email
    license="MIT"
    url="https://github.com/sromero1992/visium_hd_processing",  # URL for your GitHub repo or package
    packages=PACKAGES,  # Look for packages inside the 'qfeatures' folder                                                                                                                       
    keywords=KEYWORDS,                                                                                                                                                                          
    package_dir={'visium_hd_functions': 'visium_hd_functions'},  # Adjust to match your structure                                                                                                                   
    install_requires=INSTALL_REQUIRES,  # Add dependencies                                                                                                                                      
    python_requires='>=3.9,<3.10',  # Specify Python version compatibility for any 3.9.x version                                                                                                
    description="A collection of functions for Visium HD data processing",
)

