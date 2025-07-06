from setuptools import setup, find_packages

setup(
    name='seqforge',
    version='0.0.1',
    packages=find_packages(where='bin'),
    package_dir={'': 'bin'},
    include_package_data=True,
    install_requires=[
        'pandas',
        'biopython',
        'tqdm',
        'logomaker',
        'numpy',
        'seaborn',
        'matplotlib'
    ],
    entry_point={
        'console_scripts': [
            'seqforge = seqforge:main'
        ]
    },
    author='Elijah R. Bring Horvath',
    description='A multi-use genomics toolkit',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/ERBringHorvath/SeqForge',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT',
        'Operating System :: OS Independent'
    ],
    python_requires='>=3.10'
)
