from setuptools import setup, find_packages

# setup.py

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name='algepy-tools', 
    version='0.1.2',      
    packages=find_packages(),
    install_requires=requirements,     
    description='A Python library for efficient number theory and abstract algebra.',  
    python_requires='>=3.7',  
    url='https://github.com/sumaddury/algepy',
    author='Sucheer Maddury',
    author_email='sm2939@cornell.edu',
)