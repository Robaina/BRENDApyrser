from setuptools import find_packages, setup
import importlib


def check_package(package):
    try:
        importlib.import_module(package)
    except Exception:
        raise ImportError("Requires {package} to "
                          "be installed before running setup.py (pip install {package})"
                          .format(package=package))
    finally:
        pass


_ = [check_package(p) for p in ['numpy', 'pandas']]

setup(
    name='brendapyrser',
    version='0.1.0',
    description='Tools to parse de BRENDA database',
    # long_description=open('DESCRIPTION.rst').read(),
    url='https://github.com/Robaina/BRENDA_database',
    author='Semidan Robaina Estevez',
    author_email='semidan.robaina@gmail.com',
    license='Creative Commons Attribution 4.0 International',
    packages=find_packages(),
    zip_safe=False,
    install_requires=['numpy', 'pandas'],
    include_package_data=False,  # to include data files in installation
    # package_data={'': ['data/*.csv']},
)
