### Please cite Milo as:
Milo, Revision 1.0.3, M. S. Teynor, N. Wohlgemuth, L. Carlson, J. Huang, S. L. Pugh, B. O. Grant, R. S. Hamilton, R. Carlsen, and D. H. Ess, Brigham Young University, Provo UT, 2021.

### Requirements
#### Python:
Milo has only been tested against Python 3.8. It is expected to work with Python 3.6+, but will not work with Python 3.5 or older.

#### Gaussian:
Milo can interface with both Gaussian 16 and Gaussian 09 to perform force calculation. Milo expects Gaussian to be loaded as a module, as either "g16" for Gaussian 16 or "g09" for Gaussian 09.

#### Operating System:
Milo is only intended to work on a Linux based operating system through the command line.

### Installation Guide
1. Unzip milo-1.0.3.zip in your home directory  
2. Add the following to your .bashrc:  
  `export PYTHONPATH=$PYTHONPATH:$HOME/milo-1.0.3`  
3. (Optional) Also add the following to your .bashrc, to make calling scripts easier:  
	`export PATH=$PATH:$HOME/milo-1.0.3/milo_1_0_3/tools`  
	`module load python/3.8`  
4. Source your .bashrc or restart your session.  

You are now ready to run your first Milo job.  

### Using Milo
To run a Milo job with Gaussian 16:  
	`module load python/3.8`  
	`module load g16`  
	`python -m milo_1_0_3 < job.in > job.out`  

To run a Milo job with Gaussian 09:  
	`module load python/3.8`  
	`module load g09`  
	`python -m milo_1_0_3 < job.in > job.out`  
