# SOPA - SHARAD Open Processing Attempt

SOPA is thought to be a suite of open source algorithms for SHAllow RADar (SHARAD) low-level data processing written in python v3.

At the time of writing the suite includes:
- the possibility to load CO-SHARPS raw data (both echoes and full ancillary dataset);
- DRSSim, a dummy radar sounder (RS) simulator.

The project SOFA (SHARAD Open Focusing Attempt, http://af-projects.it/sofa) will be included here in the future in order to share the same interfaces for raw data reading.

Unfortunately, almost no documentation is provided. An example program showing how to use the main functions can be found in the root of the project.

Feel free to contact me for any help (contacts here: http://af-projects.it/contacts).

***NOTE:** This project has not been developed during my work at the University of Trento, nor it has been supported by the University of Trento in any way.*

### Requirements
SOPA has been tested only in a Linux environment using python v3.5 and the following modules:
- sys
- math
- os
- glob
- datetime
- gzip (needed only if txt files related to CO-SHARPS ancillary data have been downloaded gzipped)
- numpy 1.17.0
- scipy 1.3.1
- spectral 0.19
- ray 0.7.3 (needed only for parallel computation, not mandatory)

### Current version
0.001

### Disclaimer

The author does not guarantee that this software will always provide correct results nor that it will not crash your hardware. In any case, any use of SOPA is ONLY user responsibility. The use of SOPA or part of it for the creation of any sub-product (e.g., scientific papers, posters, images, other softwares) must be acknowledged.
As you may imagine, the code is not optimized for fast running.