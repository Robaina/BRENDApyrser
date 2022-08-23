from setuptools import setup
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), 'r', encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='brendapyrser',
    version='0.0.2',
    description='Tools to parse the BRENDA database',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/robaina/BRENDA_database',
    donwload_url='https://github.com/robaina/BRENDA_database',
    author='Semidán Robaina Estévez, 2020-2022',
    author_email='srobaina@ull.edu.es',
    maintainer='Semidán Robaina Estévez',
    maintainer_email='srobaina@ull.edu.es',
    license='BSD-3-Clause license',
    install_requires=['numpy', 'pandas', 'importlib-metadata >= 1.0 ; python_version < "3.8"'],
    packages=['brendapyrser']
)