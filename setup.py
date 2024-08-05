from setuptools import setup, find_packages

packages = [
        'pandas',
        'matplotlib',
        'icecream',
        'numpy',
        'scipy',
]

setup(
    name='junction_rivers_stage2',
    version='1.0.0',
    author='Windlab',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=packages,
)
