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
from chimera import __version__
here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()


# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-em-chimera',  # Required
    version=__version__,  # Required
    description='Plugin to use chimera programs within the Scipion framework',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/scipion-em/scipion-em-chimera',  # Optional
    author='Marta Martinez, Roberto Marabini',  # Optional
    author_email='scipion@cnb.csic.es',  # Optional
    keywords='scipion cryoem imageprocessing scipion-3.0',  # Optional
    packages=find_packages(),
    include_package_data=True,
    package_data={  # Optional
        'chimera': [''],
    },
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3'
    ],
    install_requires=[requirements], # Optional
    entry_points={
            'pyworkflow.plugin': 'chimera = chimera'
        },

)
