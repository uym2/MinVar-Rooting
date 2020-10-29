from setuptools import setup, find_packages
import fastroot
from os import walk, listdir
from os.path import join, normpath, isfile

def recursive_list_dir(path):
    listing=[]
    for x in walk(path):
        if isfile(x[0]):
            listing.append(x[0].split(path+'/')[1])
        for y in listdir(x[0]):
            z = normpath(join(x[0],y))
            if isfile(z):
                listing.append(z.split(path+'/')[1])
    return listing

param = {
    'name': fastroot.PROGRAM_NAME,
    'version': fastroot.PROGRAM_VERSION,
    'author': fastroot.PROGRAM_AUTHOR,
    'license': fastroot.PROGRAM_LICENSE,
    'package_data':{'':recursive_list_dir('fastroot_tests')},
    'packages': find_packages(),
    'include_package_data': True,
    'scripts': ['FastRoot.py','compute_RTT.py','compute_variance.py','FastRoot_tests.py'],
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
    #'test_suite': 'test_suite',
    #'tests_require': ['nose'],                
}

setup(**param)
