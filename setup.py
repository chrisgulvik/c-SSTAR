import glob
import os
from setuptools import setup, find_packages


def get_version(version='Major.Minor.Patch'):
    version_file = os.path.join('csstar', 'VERSION')
    if os.path.exists(version_file):
        with open(version_file) as f:
            version = f.readline().strip()
    return version


def readme():
    with open('README.md', 'r') as f:
        return f.read()


setup(
    name='c-SSTAR',
    version=get_version(),
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, <4',
    author='Christopher A. Gulvik',
    packages=find_packages(),
    entry_points={'console_scripts': ['csstar = csstar.__main__:main']},
    scripts=glob.glob('scripts/*'),
    package_data={'csstar': ['VERSION', 'data/*', 'data/*/*']},
    include_package_data=True,
    zip_safe=False,
    url='https://github.com/chrisgulvik/c-SSTAR',
    license='Apache 2.0',
    description=(
        'Antibiotic resistance gene detection from contigs in FastA'
        ' or GenBank format'),
    long_description=readme(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=['biopython >= 1.68'])
