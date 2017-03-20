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
In the "emetnet_install_files/emetnet" folder you can find part of the C++ codes that I wrote in the first stages of my PhD. The results of these analyses are published in the following two papers:
1. Hosseini SR, Barve A, Wagner A. Exhaustive Analysis of a Genotype Space Comprising 10 15 Central Carbon Metabolisms Reveals an Organization Conducive to Metabolic Innovation; 2015; PLOS Computational biology
2. Hosseini SR, Martin OC, Wagner A. Phenotypic innovation through recombination in genome-scale metabolic networks; 2016; Proceedings of the Royal Society B.

Please note that most of the functions in erandomwalk class were previously written by a former PhD student in our lab (Joao Rodriguez). Furthermore, he also made his own library (emetnet_install_files/eutils) for this project. 

### ENVS
In the "ENVS" folder, there are 50 .flx files each of which defines a specific environment with a given carbon source.
### DATA
In the "DATA" folder, you can find the following 3 files:
1. universe-gcs.net
This is a text file including 6588 lines each of which represents a given reaction known to exist among prokaryotic organisms.
2. donor.dat
This is an example of a donor genotype, which is exclusively viable on glucose. This file includes a binary vector of length 6588 each elemts of which corresponds to a reaction in the universe, and it is 1 if the reaction exists in the donor genotype and 0 otherwise.
3. recipient.dat
similar to donor.dat this is a genotype that is only viable on glucose. We use these files in our following two examples.
### EXAMPLES
Finally, in the folder "EXAMPLES", I have provided two fully commented C++ codes that you can run as follows:

## Example 1: MCMC Sampling of parental genotypes with genotypic distance=300
### sampling_phenotypedist_metropolis.cpp
This program samples from the vast genotype space a pair of genotypes with the same phenotype and a given genotypic distance.
The inputs of the program are 
1) The universe of the reactions "./DATA/universe-gcs.net"
2) The initial donor genotype "./DATA/donor.dat"
3) The initial recipient genotype "./DATA/recipient.dat"
4) All the environment files "./ENVS/*.flx"
5) delta: is the desired genotypic distance between the sampled genotype pair.

and the output of the program is "./donor_recipient.dat", which includes the genotype of vector of the finally sampled pair of donor-recipient genotypes with the specified genotypic distance.

You can run the code as follows:
sampling_phenotypedist_metropolis  ./DATA/universe-gcs.net ./DATA/donor.dat ./DATA/recipient.dat ./donor_recipient.dat --delta 300 ./ENVS/*.flx

## example 2: Recombination between parental genotypes to create 100 recombinant metabolic networks
### recombination_genomescale_phen_dist.cpp
This program simulates recombination in metabolic networks.
Starting from a pair of genotypes with a given genotypic distance (as created by the program described in example 1), this program creates a given number of recombinant genotypes by adding a given number of randomly selected reactions from the donor genotype to the recipient one, and deleting a given number of randomly selected genotypes from the recipient genotype.

The inputs of the program are:
1) The universe of the reactions "./DATA/universe-gcs.net"
2) The donor-recipient genotype pair "./DATA/donor_recipient.dat" as created by the program described in example1
3) All the environment files "./ENVS/.flx"
4) num1: number of reactions present in the recipient (or donor) genotype
5) num2: number of reactions that is present in the donor genotype, but absent in the recipient genotype (i.e. half of the genotypic distance).
6) num3: number of reactions that are transferred from the donor genotype to the recipient one.

You can run the code as follows:
recombination_genomescale_phen_dist ./DATA/universe-gcs.net ./DATA/donor_recipient.dat --outnet ./recombinants.dat  --num1 2079 --num2 150 --num3 5  --iter 100 ./ENVS/*.flx
