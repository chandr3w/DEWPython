import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='DEWPython',
    version='1.0.0',
    author='chandr3w',
    author_email='amchan@caltech.edu',
    description='Python-Implemented Deep Earth Water Model',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/chandr3w/DEWPython',
    packages=setuptools.find_packages(),
    #package_dir={'': 'DEWPython'},
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ],
    package_data={'': ['LICENSE.txt', 'mineralDictionary.txt','aqueousLst.txt','gasLst.txt','input.csv','dielectric.csv','Wat_den.csv','water_gibbs.csv']},
    include_package_data=True
)
