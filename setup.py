from setuptools import setup, find_packages

with open('README.rst') as file:
    long_description = file.read()

setup(
    name = "pHcalc",
    version = "0.1.2",

    description = "Systemtic pH calculation package for Python.",
    long_description = long_description,
    url = "https://github.com/rnelsonchem/pHcalc",

    author = "Ryan Nelson",
    author_email = "rnelsonchem@gmail.com",

    license = "BSD",
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],

    keywords = "pH systematic distribution titration strong weak acid base " +\
            "speciation",

    packages = find_packages(),
    install_requires = [
        'numpy>=1.10.0',
        'scipy>=0.17.0',
    ],

)

