from setuptools import setup, find_packages

setup(
    name='multiblast',
    version='2.1.0',
    packages=find_packages(where='bin'),
    package_dir={'': 'bin'},
    include_package_data=True,
    install_requires=[
        'pandas',
        'biopython'
    ],
    entry_point={
        'console_scripts': [
            'multiblast = multiblast:main'
        ]
    },
    author='Elijah R. Bring Horvath',
    description='A BLAST-based genomics toolkit',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/ERBringHorvath/multiBLAST',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent'
    ],
    python_requires='>=3.8'
)