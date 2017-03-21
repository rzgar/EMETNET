#!/usr/bin/env bash

export DEBIAN_FRONTEND=noninteractive

###############################################################################
#						Install Compilation Tools							  #
###############################################################################


# G++ compiler
echo ''
echo '@@@@@@@@@@ INSTALLING G++ @@@@@@@@@@'
echo ''
sudo apt-get --assume-yes install g++

# Automake
echo ''
echo '@@@@@@@@@@ INSTALLING AUTOMAKE @@@@@@@@@@'
echo ''
sudo apt-get --assume-yes install automake

# M4
echo ''
echo '@@@@@@@@@@ INSTALLING M4 @@@@@@@@@@'
echo ''
sudo apt-get --assume-yes install m4

###############################################################################
#						Install Dependecies									  #
###############################################################################

# Zilb
echo ''
echo '@@@@@@@@@@ INSTALLING ZLIB @@@@@@@@@@'
echo ''
wget http://zlib.net/zlib-1.2.8.tar.gz
tar -xvf zlib-1.2.8.tar.gz
cd zlib-1.2.8
./configure
make
sudo make install
cd ..
rm zlib-1.2.8.tar.gz

# Gawk
echo ''
echo '@@@@@@@@@@ INSTALLING GAWK @@@@@@@@@@'
echo ''
sudo apt-get --assume-yes install gawk

# GSL (GNU Scientific Library)
echo ''
echo '@@@@@@@@@@ INSTALLING GSL @@@@@@@@@@'
echo ''
sudo apt-get --assume-yes install libgsl0-dev

# Coin-CLP
echo ''
echo '@@@@@@@@@@ INSTALLING COIN-CLP @@@@@@@@@@'
echo ''
wget https://launchpad.net/ubuntu/+archive/primary/+files/coinor-libclp-dev_1.15.5-1ubuntu3_amd64.deb
sudo dpkg --force-confnew -i coinor-libclp-dev_1.15.5-1ubuntu3_amd64.deb
sudo apt-get --assume-yes install -f
rm coinor-libclp-dev_1.15.5-1ubuntu3_amd64.deb

# Igraph
echo ''
echo '@@@@@@@@@@ INSTALLING IGRAPH @@@@@@@@@@'
echo ''
wget https://launchpad.net/ubuntu/+archive/primary/+files/libigraph0-dev_0.6.5-5_amd64.deb
sudo dpkg --force-confnew -i libigraph0-dev_0.6.5-5_amd64.deb
sudo apt-get --assume-yes install -f
rm libigraph0-dev_0.6.5-5_amd64.deb

# EUtils
echo ''
echo '@@@@@@@@@@ INSTALLING EUTILS @@@@@@@@@@'
echo ''
cd eutils
./configure
make
sudo make install
cd ..

###############################################################################
#						Install Emetnet										  #
###############################################################################

echo ''
echo '@@@@@@@@@@ INSTALLING EMETNET @@@@@@@@@@'
echo ''

cd emetnet
./configure
make clean
make clean
sed -i -e "s/LIBS = -lClp -lCoinUtils -lpthread -lz -leutils/LIBS = -lClp -lCoinUtils -leutils -lgsl -lgslcblas -ligraph -lz -lpthread -lreadline/" Makefile # added -lreadline
make
sudo make install
cd ..

###############################################################################
#						Test Emetnet										  #
###############################################################################

echo ''
echo '@@@@@@@@@@ RUNNING EMETNET WITH TEST DATA @@@@@@@@@@'
echo ''

genotyping test_data/universe-gcs.net test_data/ecoli-gcs.net test_data/ecoli-gcs.dat
