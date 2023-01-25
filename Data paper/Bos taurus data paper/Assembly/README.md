# SeqOccIn - A Bos taurus data paper assembly scripts

The scripts located in this folder enable to reproduce the different assemblies presented in the data paper. 

ftp.sh and ascp.sh enable to download all the data files from the ENA. If you intend to use ascp.sh you have to modify the script in order to add the correct location of your asperaweb_id_dsa.openssh file.

Once the data has been downloaded you can check the header of the run_assembly.sh file in which you have all the software to be install as well as the versions used for the paper. 

Once the software packages are available in you path you can launch the run_assembly.sh script which will process the data an produce the assemblies. The computer on which you run the script should have at least 500Go RAM and 24 CPUs.
