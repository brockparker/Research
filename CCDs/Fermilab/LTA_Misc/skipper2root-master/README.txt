This program process the raw Skipper CCD data.
It computes overscan mean for each sample and subtracts it line by line.
The output file will be a ROOT image containing a TTree with the pixels 
values averaged over all the samples (after subtraction of the corresponding 
overscan value). It's also possible to save an additional TTree that will 
contain the individual values of all the samples.

The script "readHeaderVarsFromRootFile.C" shows how to read the header data
from the root file.
