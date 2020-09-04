from setuptools import setup, find_packages
import fastroot
from os import walk, listdir
from os.path import join, normpath, isfile

param = {
    'name': fastroot.PROGRAM_NAME,
    'version': fastroot.PROGRAM_VERSION,
    'author': fastroot.PROGRAM_AUTHOR,
    'license': fastroot.PROGRAM_LICENSE,
    'packages': find_packages(),
    'include_package_data': True,
    'scripts': ['FastRoot.py','compute_RTT.py','compute_variance.py'],
    'zip_safe': True,
    'install_requires': ['treeswift>=1.1.14', 'cvxopt>=1.2.5', 'numpy>=1.19.0'],
    'keywords': 'Phylogenetics Evolution Biology',
    'classifiers': ["Environment :: Console",
                    "Intended Audience :: Developers",
                    "Intended Audience :: Science/Research",
                    "Natural Language :: English",
                    "Operating System :: OS Independent",
                    "Programming Language :: Python",
                    "Topic :: Scientific/Engineering :: Bio-Informatics",
                    ],
}

setup(**param)
