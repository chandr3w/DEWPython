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
    download_url = https://github.com/chandr3w/DEWPython/archive/1.0.0.tar.gz,
    install_requires=[
        'numpy',
        'pandas',
        'sys',
        'threading',
        'subprocess',
        'os',
        'json',
        'matplotlib',
        'mpl_toolkits',
        'collections',
      ],
    #package_dir={'': 'DEWPython'},
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ],
    include_package_data=True
)
