from setuptools import find_packages, setup

with open('README.md', 'r') as readme:
    long_description = readme.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='ICESat-2-sea-ice-thickness',
    version='1.0',
    author='Alek Petty',
    author_email='alek.a.petty@nasa.gov',
    description='Code repository for producing sea ice thickness estimates from ICESat-2 freeboard data',
    license='MIT',
    long_description=long_description,
    url='https://github.com/akpetty/ICESat-2-sea-ice-thickness',
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 5 - Production/Stable',
        'Topic :: Scientific/Engineering'
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3 :: Only',
    ],
    python_requires='>=3,<3.7',
    install_requires=['basemap>=1.1.0', 'numpy>=1.13', 'proj4>=5.2.0', 'pandas>=0.22', 'geos>=3.6.2', 'h5p>=2.8.0', 'matplotlib>=2.2.2', 'scipy>=0.19.1', 'xarray>=0.12.3', 'dask>=0.17.5', 'netcdf4>=1.3.1'],
    keywords='seaice'
)