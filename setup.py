from setuptools import setup, find_packages

# setup.py

setup(
    name='algepy', 
    version='0.0.1',      
    packages=find_packages(where="src"),
    install_requires=[
    ],     
    description='A Python library for efficient number theory and abstract algebra.',  
    python_requires='>=3.7',  
    url='https://github.com/sumaddury/algepy',
    author='Sucheer Maddury',
    author_email='sm2939@cornell.edu',
    package_dir={"": "src"},
)