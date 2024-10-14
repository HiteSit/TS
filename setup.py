from setuptools import setup, find_packages

setup(
    name='TS',
    version='1',
    description='Thomson Sampling',
    author='Pat Walters',
    author_email='your.email@example.com',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.25',
        'pandas>=2.0',
        'rdkit>=2022.03',
        'tqdm>=4.65',
        'useful-rdkit-utils>=0.2.7',
        'pillow>=9.5',
        'notebook',
        'scikit-learn',
        'seaborn',
        "datamol"
    ],
    include_package_data=True,
)