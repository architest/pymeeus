try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from os import path
import pypandoc


# Converts Markdown to reStructured
z = pypandoc.convert('README.md', 'rst', format='markdown')

# Writes converted file
with open('README.rst', 'w') as outfile:
    outfile.write(z)


here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.md')) as f:
    long_description = f.read()


config = {
    'name': 'PyMeeus',
    'version': '0.0.1',
    'description': 'Python implementation of Jean Meeus astronomical routines',
    'long_description': long_description,
    'keywords': 'Meeus astronomy module library',
    'license': 'LGPLv3',
    'author': 'Dagoberto Salazar',
    'author_email': 'dagoberto.salazar@gmail.com',
    'url': 'https://github.com/architest/pymeeus',
    'download_url': 'https://github.com/architest/pymeeus',
    'install_requires': ['nose'],
    'packages': ['pymeeus'],            # 'scripts': ['example.py'],
    'py_modules': ['base'],
    'classifiers': [
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Astronomy'
    ]
}

setup(**config)
