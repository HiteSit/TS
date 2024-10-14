from setuptools import setup, find_packages

setup(
    name='TS',
    version='1.0',
    description='Thomson Sampling',
    author='Pat Walters',
    author_email='your.email@example.com',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'rdkit>=2022.03',
        'tqdm',
        'useful-rdkit-utils>=0.2.7',
        'pillow>=9.5',
        'notebook',
        'scikit-learn',
        'seaborn',
        "datamol"
    ],
    include_package_data=True,
)