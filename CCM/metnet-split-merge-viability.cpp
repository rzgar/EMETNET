#include "enet.h"
#include "erandomwalk.h"
#include <eutils/logger.h>
#include <eutils/emain.h>
#include <signal.h>
#include <eutils/estrarrayof.h>
#include <eutils/eregexp.h>
#include <eutils/eheap.h>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
/////////////////////////////////////////////////////////////////////////////////////////////////
estr solver="esolver_clp";
int netsize=-1;
int strict=0;
int force_away=0;
int periphery_only=0;
int mutate_transport=0;
int internal_secretion=0;
int only_viable=0;
int s1=0;//The number of viable genotypes in inputfile1
int s2=0;//The number of viable genotypes in inputfile2
enet net;
/////////////////////////////////////////////////////////////////////////////////////////////////
erandomWalk *prw=0x00;
/////////////////////////////////////////////////////////////////////////////////////////////////
int emain() {

  ldieif(argvc<7,"syntax: metnet-split-merge-viability universe.net <inputfilename1.dat> <inputfilename2.dat> <outfilename> <env.flx> --s1 [>1] --s2 [>1] ");	
  epregister(solver);
  epregister(s1);
  epregister(s2);
  eparseArgs(argvc,argv);

  net.load(argv[1]);//Loading the reaction universe (universe.net)
  efile file1;//The file containing the viable genotypes of a merged block from previous steps (inputfilename1.dat)
  efile file2;//The file containing the viable genotypes of a merged block from previous steps (inputfilename2.dat)
  efile fileout;//The outfile storing the viable genotypes after merging the block in inputfile1 and inputfile2 

  file1.open(argv[2],"r");
  file2.open(argv[3],"r");
  fileout.open(argv[4],"w");

 	
  net.correct_malformed();
  erandomWalk rw(net,solver,strict);
  prw=&rw;
  rw.periphery_only=periphery_only;
  rw.mutate_transport=mutate_transport;
  rw.internal_secretion=internal_secretion;
  rw.only_viable=only_viable;
  rw.setRSize(netsize);

  rw.getEnv(argvc,argv);//getting the environment files and read them
  rw.load(net); // load network in object "net" into the erandomwalk 
  rw.calcPhenotype(); //calculate phenotype
  rw.viablePhenotype=rw.phenotype;//checking viability

  int viables1[s1][45];
  int viables2[s2][45];
  int i;
  int linecounter=0;
  estr str1;
  estr str2;
  estrarray parts;

  ////////////// Reading files /////////////////////////////////////////
  linecounter=0;
  while (file1.readln(str1)) {//read inputfile1 line by line
	 parts=str1.explode(" ");
	 for (i=0; i<45; i++) {viables1[linecounter][i]=0;}
         for (i=0; i<parts.size(); ++i) {viables1[linecounter][i]=parts[i].i();}//storing the deleted reactions of the viable genotypes of the inputfile2
	 linecounter++;
  }
  file1.close();

  linecounter=0;
  while (file2.readln(str2)) {//read inputfile2 line by line
	 parts=str2.explode(" ");
	 for (i=0; i<45; i++) {viables2[linecounter][i]=0;}
         for (i=0; i<parts.size(); ++i) {viables2[linecounter][i]=parts[i].i();}//storing the deleted reactions of the viable genotypes of the inputfile2
	 linecounter++;
  }
  file2.close();


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////// Merging step //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<s1;i++){
       for (int j=0;j<s2;j++){
            eintarray finalarray;
            for (int l=0;l<45;l++){if (viables1[i][l]>0){int tmp=(viables1[i][l]+26);finalarray.add(tmp);}}      
            for (int l=0;l<45;l++){if (viables2[j][l]>0){int tmp=(viables2[j][l]+26);finalarray.add(tmp);}}    
            for (int k=0; k<finalarray.size(); k++) {rw.disable(finalarray[k]);}
            rw.calcPhenotype();
            if (rw.isViable()) {
                eintarray sortedindex, tmparr;
	        sortedindex = sort(finalarray);// sort returns the sorted index list
	        for (int p=0;p<sortedindex.size();p++) {tmparr.add(finalarray[sortedindex[p]]);}
	        finalarray = tmparr;
                for (int q=0;q<tmparr.size();q++) {tmparr[q] -= 26;}
                estr zeroloc = intarr2str2(tmparr);
                fileout.write(zeroloc+"\n");
            }
            for (int n=0; n<finalarray.size(); ++n) {rw.activate(finalarray[n]);}
       }
  }
  fileout.close();
  return(0);
}
