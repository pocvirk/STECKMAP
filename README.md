IMPORTANT:
A new environment variable must be defined in your .bashrc or .cshrc or whatever...
STECKMAPROOTDIR 
is the directory where the STECKAMP.tar package has been untarred, and hence the 
directory from which the STECKMAP directory is visible.
Note that STECKMAPROOTDIR should end with a /

pop_paths.i is now the first file to be included, as it redefines all the paths, using
include, "FULLPATH/pop_paths.i"
it recovers the environment variable STECKMAPROOTDIR and uses it to redefine the paths 
relevant for STECKMAP.
then, sfit can be included as usual:
include, "sfit.i"

and the rest is just as usual..

