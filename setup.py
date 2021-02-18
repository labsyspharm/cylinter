from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    LONG_DESCRIPTION = fh.read()

requires = [
    'PyYAML==5.4',
    'pyarrow==2.0.0',
    'pandas==1.2.0',
    'hdbscan==0.8.26',
    'joblib==0.17.0',  # newly-released joblib==1.0.0 results in issues with parallel processing with larger datasets
    'matplotlib==3.3.3',
    'seaborn==0.11.1',
    'scikit-image==0.18.1',
    'napari==0.3.6',  # points layer doesn't work with napari[all] (v0.4.3)
    'PyQt5==5.15.2',  # must include if not installing napari[all]
    'zarr==2.6.1',
    'natsort==7.1.0',
    'rpy2==3.4.2',
    'bridson==0.1.0',
    'hurry.filesize==0.9',
    'synapseclient==2.0.0',  # v2.2.2 (current version as of 01/25/21 fails to transfer)
]

VERSION = '0.0.17'
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
    setup_requires=['cython>=0.29.21', 'numpy>=1.19.5'],
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
    python_requires='==3.8.5',
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    url=HOMEPAGE,
    download_url='%s/archive/v%s.tar.gz' % (HOMEPAGE, VERSION),
    keywords='scripts single cell data science',
    zip_safe=False,
)
