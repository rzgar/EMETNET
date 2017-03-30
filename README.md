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
### emetnet_install_files/emetnet
In the "emetnet_install_files/emetnet" folder you can find part of the C++ programs that are described in more detail in the following two publications.

1. Hosseini SR, Barve A, Wagner A. Exhaustive Analysis of a Genotype Space Comprising 10 15 Central Carbon Metabolisms Reveals an Organization Conducive to Metabolic Innovation; 2015; PLOS Computational biology
2. Hosseini SR, Martin OC, Wagner A. Phenotypic innovation through recombination in genome-scale metabolic networks; 2016; Proceedings of the Royal Society B.

Most of the functions in the ‘erandomwalk’ class that is used in these programs were previously written by Joao Rodriguez, a former PhD student in the lab of Andreas Wagner at the University of Zurich. The library ‘emetnet_install_files/eutils’ was also written by Joao Rodriguez. The research that resulted from this work was published in the following paper:

3. Matias Rodrigues JF, Wagner A. Evolutionary plasticity and innovations in complex metabolic reaction networks. PLoS Comput Biol. 2009;5:e1000613.

### ENVS
The "ENVS" folder contains 50 .flx files, each of which defines a specific chemical environment with a unique carbon source.
### DATA
In the "DATA" folder, you can find the following 3 files:
1. universe-gcs.net
A text file with 6588 lines, each representing a biochemical reaction known to exist in prokaryotic organisms.
2. donor.dat
This is an example of a donor genotype, which is exclusively viable on glucose. This file includes a binary vector of length 6588. Each entry of this vector corresponds to a reaction in the prokaryotic reaction universe. An entry equals one if the reaction exists in the donor genotype, and zero otherwise.
3. recipient.dat
Analogous to donor.dat, this file encodes a recipient genotype that is only viable on glucose. The files are used in the following two examples.
### EXAMPLES
The folder "EXAMPLES” contains two fully commented C++ programs that can be executed as follows:

## Example 1: MCMC Sampling of parental genotypes with genotypic distance=300
### sampling_phenotypedist_metropolis.cpp
This program samples from a vast genotype space a pair of genotypes with the same phenotype and a given genotypic distance.
The input files to the program are 
1) The prokaryotic universe of reactions "./DATA/universe-gcs.net"
2) The initial donor genotype "./DATA/donor.dat"
3) The initial recipient genotype "./DATA/recipient.dat"
4) All environment files "./ENVS/*.flx"
5) delta: is the desired genotypic distance between the sampled genotype pair.

The output of the program is a file with extension dat (e.g. "./donor_recipient.dat"), which includes the genotype vector of the sampled pair of donor-recipient genotypes with genotypic distance delta.

You can execute the code as follows:

sampling_phenotypedist_metropolis  ./DATA/universe-gcs.net ./DATA/donor.dat ./DATA/recipient.dat ./donor_recipient.dat --delta 300 ./ENVS/*.flx

## example 2: Recombination between parental genotypes to create 100 recombinant metabolic networks
### recombination_genomescale_phen_dist.cpp
This program simulates recombination in metabolic networks.
Starting from a pair of genotypes with a given genotypic distance (as created by the program described in example 1), this program creates a given number of recombinant genotypes by adding a given number of randomly selected reactions from the donor genotype to the recipient one, and deleting a given number of randomly selected genotypes from the recipient genotype.

The input files to the program are:
1) The prokaryotic universe of reactions "./DATA/universe-gcs.net"
2) The donor-recipient genotype pair "./DATA/donor_recipient.dat" as created by the program described in example1
3) All environment files "./ENVS/.flx"
4) num1: number of reactions present in the recipient (or donor) genotype
5) num2: number of reactions that is present in the donor genotype, but absent in the recipient genotype (corresponding to half of the genotypic distance between donor and recipient).
6) num3: number of reactions that are transferred from the donor genotype to the recipient one.

The program produces the following output files:

1) a file with extension dat (e.g. "./recombinants.dat", which includes the information about how many recombinants lost viablity on one or more carbon sources, and how many gained viability on new carbons sources 
2) a file with extension dat_phen (e.g. ./recombinants.dat_phen), which includes the phenotype of each recombinant genotype.

You can execute the code as follows:

recombination_genomescale_phen_dist ./DATA/universe-gcs.net ./DATA/donor_recipient.dat --outnet ./recombinants.dat  --num1 2079 --num2 150 --num3 5  --iter 100 ./ENVS/*.flx
