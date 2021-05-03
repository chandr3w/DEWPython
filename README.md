# DEWPython

[![License: MIT License](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT)
[![Twitter Follow](https://img.shields.io/twitter/follow/chandr3w.svg?style=flat-square&logo=twitter&label=Follow)](https://twitter.com/chandr3w)

V 1.3.7

The DEW package allows the user to compute the thermodynamic and elastic properties of various aqueous inputs for a general range of 100-1200C and a pressure range of 1.0-60 Kb. It is based on the [DEW spreadsheet](http://www.dewcommunity.org/) and behaves similarly.The DEW package additionally provides integrated support for [SUPCRTBL](https://models.earth.indiana.edu/supcrtbl.php) and can be used to directly import and compare species between the two models.

## Getting Started

This section provide a basic example on the running of DEW. Because the input fields are interactive, it relies on user input in order to function (which cannot be demonstrated here). [Full documentation](https://chandr3w.github.io/DEW_amchan/) for the class is available externally.

### Download
Download and install using using 
```
pip install DEWPython
```

### Running DEW
Import DEWPython and DEWPython.DEWModel
```
import DEWPython
from DEWPython import DEWModel as dm
```
You are now ready to use the DEW module! This command will execute the imports of packages and the initialization of the model. If you wish to change base parameters please see the documentation. The only non-standard package that DEW is dependent on is the "pandas" package, however a full list of dependencies is included below.

### Using the Model
DEW is an object-oriented class dependent on the DEWEquations class. To run it, first initialize a DEW object:
```
reaction = dm.DEW()
```
If you run
```
reaction.options()
```
It will give you a brief overview of the documentation that follows.

From here, set the inputs, outputs, and preferences interactively. *Throughout every stage in the model, parameters can be queried for debugging. See the documentation for more details.* 

```
reaction.set_inputs()
# Initialize the LHS of the reaction
reaction.set_outputs()
# Initialize the RHS of the reaction
reaction.set_preferences()
# By default, the reaction runs at Psat conditions (see documentation). The set_preferences() method interactively changes this.
```
If you are unsure what the name of your inputs/outputs are in slop16, you can run a search on a string of components/part of the name that your input or output contains
```
dm.search(string)
```
If custom options are selected in the preferences, the relevant CSV files must be updated. Otherwise, you must run the following command:
```
reaction.set_TPRho()
```
This command will interactively initialize temperature and pressure arrays. Finally, run
```
reaction.calculate()
```
You can now query
```
reaction.delG
# Gibbs of Formation
reaction.delV
# Volume change of Formation
reaction.logK
# The log K value of the reaction
reaction.make_plots()
# Automatically constructs the plots for the model.
```
You can now also export your plots to CSV format using the command
```
reaction.export_to_csv()
```

### Running supcrt
##### NOTE: TO RUN SUPCRT YOU *MUST* HAVE THE .CON/.DAT FILES IN YOUR WORKING DIRECTORY

Included in the DEW_Folder is supcrt96 and SUPCRTBL, a similar program to calculate the properties of species at different temeperature and pressure conditions. The defauly is set to supcrt96.
```
reaction.run_supcrt()
```
This command will interactively run SUPCRTBL inline and create output files. It automatically updates "reaction.supcrtFile", a variable that stores the most recently produced Supcrt file.

Now run:
```
reaction.calculate_supcrt()
```
This function takes one optional argument of a SUPCRTBL output file. For the input file or the file previously calculated it will calculate "reaction.supcrtOut", a dictionary that can be queried for 'delG', 'delV', 'LogK', 'delH', 'delS', 'delCp', 'DH2O', 'Temperature', and 'Pressure'. 

Finally, run:
```
reaction.make_supcrt_plots()
```
To produce the same plots as DEW.

### Automatic Input
Added in version 1.3.1, you can now run supcrt/DEW without running the input loops. For DEW you can use  
```
reaction.run(pt_arr, min_inp =[], aq_inp = [], g_inp = [], h2o_inp = 0, min_out = [],aq_out =[], g_out = [],h2o_out = 0, 
        ptInp = 'Psat', rhoWat = 'Z&D 2005', forceBool = False, dieEQ = 'Supcrt', forceSC = True, 
        WFEQ ='D&H 1978', dsV = True, pdsV = True, DV = True, EQ = 1, dEQ = 1, pst = True, mWn = 1, makeP = False))
```
Where pt_arr is either of shape ((tmin, tmax, tstep)) for Psat inputs or (((tmin, tmax, tstep),(pmin, pmax, pstep))) for non-psat.

Similarly, for supcrt you can run
```
def supcrt_inp(rxn_lst, reaction_type = 'psat')
```
Where the first value is a list of reaction arrays where each array is split into tuples of reactants formatted akin to supcrt (e.g., (-1, Ca+2)) and the reaction type determines which of the pre-established reaction files are run.

### Range of validity
Certain equations within DEW are valid to certain values (as are the properties of specific mineral species). For more information please see the [DEW website](http://www.dewcommunity.org/)

### Dependencies

* Pandas
* Numpy
* Sys
* Threading
* Subprocess
* Os
* Json
* Matplotlib
* Matplotlib Toolkits
* Collections

## References
* Huang, F., & Sverjensky, D. A. (2019). Extended Deep Earth Water Model for predicting major element mantle metasomatism. Geochimica et Cosmochimica Acta, 254, 192-230.
* Facq, S., Daniel, I., Montagnac, G., Cardon, H., & Sverjensky, D. A. (2016). Carbon speciation in saline solutions in equilibrium with aragonite at high pressure. Chemical Geology, 431, 44-53.
* Facq, S., Daniel, I., Montagnac, G., Cardon, H., & Sverjensky, D. A. (2014). In situ Raman study and thermodynamic model of aqueous carbonate speciation in equilibrium with aragonite under subduction zone conditions. Geochimica et Cosmochimica Acta, 132, 375-390.
* Johnson, J.W., Oelkers, E.H. and Helgeson, H.C. (1992) SUPCRT92 - A software package for calculating the standard molal thermodynamic properties of minerals, gases, aqueous species, and reactions from 1-bar to 5000-bar and 0C to 1000C. Computer and Geosciences 18:899-947.
* Sverjensky, D. A., Harrison, B., & Azzolini, D. (2014). Water in the deep Earth: The dielectric constant and the solubilities of quartz and corundum to 60kb and 1200 C. Geochimica et Cosmochimica Acta, 129, 125-145.
* Pan, D., Spanu, L., Harrison, B., Sverjensky, D. A., & Galli, G. (2013). Dielectric properties of water under extreme conditions and transport of carbonates in the deep Earth. Proceedings of the National Academy of Sciences, 110(17), 6646-6650.
* Zimmer, K., Zhang, Y.L., Lu, P., Chen, Y.Y., Zhang, G.R., Dalkilic, M. and Zhu, C. (2016) SUPCRTBL: A revised and extended thermodynamic dataset and software package of SUPCRT92. Computer and Geosciences 90:97-111. 

## Authors

* **Andrew Chan** - *Div. of Geological and Planetary Sciences, California Institute of Technology, Pasadena, CA, USA 91125* 
* **Mohit Melwani Daswani** - *Group 3226 (Planetary Interiors and Geophysics). NASA Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA 91109*
* **Steven Vance** - *Group 3226 (Planetary Interiors and Geophysics). NASA Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA 91109*

## Change log

### Changes since 1.2.0
- Updated stored species to include slop16
- Integrated SUPCRTBL/supcrt16
- Updated functionality (plots, automatic feeding, search function)

### Planned updates

## License

This project is licensed under the MIT License :

Copyright (c) 2021, A. Chan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files, to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Acknowledgments

This work was produced with the financial support provided by the NASA Jet Propulsion Laboratory, the California Institute of Technology Summer Undergraduate Research Fellowship program, and the Caltech Assocaites.
