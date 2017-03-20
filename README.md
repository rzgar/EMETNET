# EMETNET
## Installation of the software
Note: It is strongly recommended that the following steps be carried out on a virtual
machine with freshly installed Ubuntu 14.04 LTS. Unfortunately, the software is not easily installable on newer versions of Ubuntu.
1. Download emetnet_install_files.tar.gz from the following URL:
https://github.com/rzgar/EMETNET/blob/master/emetnet_install_files.tar.gz
2. Change directory to the one where the archive has been downloaded
3. Open the archive using the following command:
tar zxf emetnet_install_files.tar.gz
4. Change directory to the opened archive as follows:
cd emetnet_install_files/
5. Run the installation script as follows:
./install_emetnet.sh

## The organization of the codes and the data
In the "emetnet_install_files/emetnet" folder you can find part of the C++ codes that I wrote in the first stages of my PhD. The results of these analyses are published in the following two papers:
1. Hosseini SR, Barve A, Wagner A. Exhaustive Analysis of a Genotype Space Comprising 10 15 Central Carbon Metabolisms Reveals an Organization Conducive to Metabolic Innovation; 2015; PLOS Computational biology
2. Hosseini SR, Martin OC, Wagner A. Phenotypic innovation through recombination in genome-scale metabolic networks; 2016; Proceedings of the Royal Society B.

Please note that most of the functions in erandomwalk class were previously written by a former PhD student in our lab (Joao Rodriguez). Furthermore, he also made his own library (emetnet_install_files/eutils) for this project. 

In the "DATA" folder, you can find the following 3 files:
1. universe-gcs.net

2. donor.dat
3. recipient.dat

In the "ENVS" folder, there are 50 .flx files each of which defines a specific environment with a given carbon source.

Finally, in the folder "EXAMPLE", I have provided two fully commented C++ codes that you can run as follows:

## Example 1: MCMC Sampling of parental genotypes with genotypic distance=300
## example 2: Recombination between parental genotypes to create 100 recombinant metabolic networks
