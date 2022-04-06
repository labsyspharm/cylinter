from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    LONG_DESCRIPTION = fh.read()

requires = [
    'numba==0.53.1',
    'PyYAML==5.4',
    'pyarrow==5.0.0',
    'pandas>=1.1.5',
    'hdbscan==0.8.26',
    'joblib==0.17.0',  # joblib==1.0.0 parallel processing fails on larger datasets
    'matplotlib==3.3.3',
    'seaborn==0.11.1',
    'scikit-learn==0.24.2',
    'scikit-image==0.18.1',
    'napari==0.4.10',
    'PyQt5==5.15.2',  # must include if not installing napari[all]
    'zarr==2.6.1',
    'natsort==7.1.0',
    'umap-learn==0.5.1',
    'hurry.filesize==0.9',
    'synapseclient==2.4.0',
    'pyparsing==2.0.3',  # version compatible with hdbscan and matplotlib
    'PyOpenGL-accelerate',  # suppresses "No OpenGL_accelerate..." message
    'cellcutter>=0.2.4',
]

VERSION = '0.0.40'
DESCRIPTION = 'CyLinter'
AUTHOR = 'Gregory J. Baker'
AUTHOR_EMAIL = 'gregory_baker2@hms.harvard.edu'
LICENSE = 'MIT License'
HOMEPAGE = 'https://github.com/labsyspharm/cylinter'

setup(
    name='cylinter',
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    setup_requires=[],
    packages=find_packages(),
    install_requires=requires,
    data_files=[('', ['cylinter/config.yml']),
                ('', ['cylinter/prep_subprocess.sh']),
                ],
    entry_points={
        'console_scripts': [
            'cylinter=cylinter.cylinter:main',
            'prep=cylinter.prep:main',
        ]
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: %s' % LICENSE,
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Visualization'
    ],
    python_requires='>=3.6.0',
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    url=HOMEPAGE,
    download_url='%s/archive/v%s.tar.gz' % (HOMEPAGE, VERSION),
    keywords='scripts single cell data science',
    zip_safe=False,
)
