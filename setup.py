#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('CHANGELOG.md') as changelog_file:
    changelog = changelog_file.read()

requirements = ['click>=7.0', 'numpy', 'pandas', 'pyfaidx', 'h5py']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

setup(
    author="Sebastian Röner",
    author_email='sebastian.roener@charite.de',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Script that encodes genomic regions from bed files to one hot encoded arrays",
    entry_points={
        'console_scripts': [
            'bed_to_1hot=bed_to_1hot.__main__:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + changelog,
    include_package_data=True,
    keywords='bed_to_1hot',
    name='bed_to_1hot',
    packages=find_packages(include=['bed_to_1hot', 'bed_to_1hot.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/sroener/bed_to_1hot',
    version='0.1.0',
    zip_safe=False,
)