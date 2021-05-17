# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py

from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.readlines()[1].strip()

setup(
    name='pyphase',
    version='1.0.1',
    description='An open source phase retrieval code',
    long_description=readme,
    long_description_content_type='text/x-rst',
    author='Max Langer',
    author_email='max.langer@creatis.insa-lyon.fr',
    url='https://gitlab.in2p3.fr/mlanger/pyPhase/',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    install_requires=required
)


# -----------------------------------------------------------------------------
# Uploading the package on pypi

# Steps
# 1 - change version in setup.py file
# 2 - commit, tag. git push --tags
# 3 - setup: python3 setup.py sdist bdist_wheel
# 4 - twine: see below

# On TEST pypi:
# twine upload --repository-url https://test.pypi.org/legacy/ dist/*
# On REAL pyip
# twine upload  dist/*

# test with
# pip uninstall pyphase
# pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple pyphase
# https://test.pypi.org/project/pyphase/
