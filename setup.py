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
    author='Semidán Robaina Estévez',
    author_email='srobaina@ull.edu.es',
    maintainer='Semidán Robaina Estévez',
    maintainer_email='srobaina@ull.edu.es',
    license='Creative Commons Attribution 4.0 International',
    install_requires=['numpy', 'pandas'],
    packages=['brendapyrser']
)
