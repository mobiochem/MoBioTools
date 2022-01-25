import os
import io
from setuptools import setup, find_packages

#def read(fname):
#    with io.open(os.path.join(os.path.dirname(__file__), fname), encoding='utf-8') as f:
#        return f.read()


setup(
        name = 'mobiotools',
        version = '0.0.0',
        author = 'Gustavo Cardenas, Juan J. Nogueira',
        author_email = 'gustavo.cardenas@uam.es',
        package_dir = {'': 'src'},
        packages = find_packages('src'),
        scripts = ['src/qminputs/main_qminputs.py', 
                   'src/qminputs/nwchem_parser.py', 
                   'src/qminputs/molcas_parser.py', 
                   'src/qminputs/gau_parser.py', 
                   'src/qminputs/common_parser.py', 
                   'src/qminputs/constants.py',
                   'src/pyoverlaps/cart2sph.py',
                   'src/pyoverlaps/ovlp_wrapper.py',
                   'src/pyoverlaps/permutations.py',
                   'src/pyoverlaps/align.py',
                   'src/pyoverlaps/matrixio.py',
                   'src/pyoverlaps/parse_molden.py',
                   'src/pyoverlaps/pyoverlaps.py',
                   'scripts/get_num_geometries.py',
                   'src/pyoverlaps/intwrap.so'],
        python_requires = '>=3.6',
        install_requires = ["numpy",
                            "pytraj"])
