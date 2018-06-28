try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'PyMeeus',
    'version': '0.0.1',
    'description': 'Python implementation of Jean Meeus astronomical routines',
    'keywords': 'Meeus astronomy module library',
    'license': 'LGPLv3',
    'author': 'Dagoberto Salazar',
    'author_email': 'dagoberto.salazar@gmail.com',
    'url': 'https://github.com/architest/pymeeus',
    'download_url': 'https://github.com/architest/pymeeus',
    'install_requires': ['nose'],
    'packages': ['pymeeus'],
    'scripts': ['example.py'],
    'py_modules': ['base'],
    'classifiers': [
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3',
        'Operating System :: OS Independent',
        'Programming Language :: Pythoni :: 2.7',
        'Programming Language :: Pythoni :: 3.6',
        'Topic :: Scientific/Engineering :: Astronomy'
    ]
}

setup(**config)
