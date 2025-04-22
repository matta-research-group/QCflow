from setuptools import find_packages, setup

# this can be more complicated and clever, but for now it's minimal
setup(
    name='qcflow',
    version='0.4.0',
    packages=find_packages(include=['QCflow', 'QCflow.*']),
)
