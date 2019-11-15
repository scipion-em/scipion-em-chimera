"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-em-chimera',  # Required
    version='2.0.1',  # Required
    description='Plugin to use chimera programs within the Scipion framework',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/scipion-em/scipion-em-chimera',  # Optional
    author='I2PC, Marta Martinez, Roberto Marabini',  # Optional
    author_email='scipion@cnb.csic.es',  # Optional
    keywords='scipion cryoem imageprocessing scipion-2.0',  # Optional
    packages=find_packages(),
    include_package_data=True,
    package_data={  # Optional
       'chimera': [''],
    }
)
