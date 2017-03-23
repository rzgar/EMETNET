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
int s1=0;//The index of the viable genotype in inputfile1
int s2=0;//The number of viable genotypes in inputfile2
enet net;
/////////////////////////////////////////////////////////////////////////////////////////////////
erandomWalk *prw=0x00;
/////////////////////////////////////////////////////////////////////////////////////////////////
int emain() {

  ldieif(argvc<7,"syntax: metnet-split-merge-memoryfriendly-viability universe.net <inputfilename1.dat> <inputfilename2.dat> <outfilename> <env.flx> --s1 [>1] --s2 [>1] ");	
  epregister(solver);
  epregister(s1);
  epregister(s2);
  eparseArgs(argvc,argv);

  net.load(argv[1]);//Loading the reaction universe (universe.net)
  efile file1;//The file containing the viable genotypes of a single block (inputfilename1.dat)
  efile file2;//The file containing the viable genotypes of a 4 merged block (inputfilename2.dat)
  efile fileout;//The outfile storing the viable genotypes after merging the block in inputfile1 and the merged blocks in inputfile2 

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


  int viables2[s2][45];
  int viable1[1][45];
  int i;
  int linecounter=0;
  estr str1;
  estr str2;
  estrarray parts;
  cout<<"alak1"<<endl;
  ////////////// Reading files /////////////////////////////////////////
  linecounter=0;
  while (file1.readln(str1)) {//read inputfile1 
	 linecounter++;
         if (linecounter==s1){
	     parts=str1.explode(" ");
	     for (i=0; i<45; i++) {viable1[0][i]=0;}
             for (i=0; i<parts.size(); ++i) {viable1[0][i]=parts[i].i();}//storing the deleted reactions of the viable genotype of the inputfile1
         }
  }
  file1.close();
  cout<<"alak2"<<endl;
  linecounter=0;
  while (file2.readln(str2)) {//read inputfile2 line by line
	 parts=str2.explode(" ");
	 for (i=0; i<45; i++) {viables2[linecounter][i]=0;}
         for (i=0; i<parts.size(); ++i) {viables2[linecounter][i]=parts[i].i();}//storing the deleted reactions of the viable genotypes of the inputfile2
	 linecounter++;
  }
  file2.close();
  cout<<"alak3"<<endl;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////// Merging step //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for (i=0;i<s2;i++){
       int temparray[90];
       cout<<"alak4.1"<<endl;
       for (int j=0;j<90;j++){temparray[j]=0;}
       int count=0;
       for (int j=0;j<45;j++){
            if (viables2[i][j]>0){int tmp=(viables2[i][j]+26);temparray[count]=tmp;count++;}
            if (viable1[0][j]>0){int tmp=(viable1[0][j]+26);temparray[count]=tmp;count++;}
       }
       cout<<"alak4.2"<<endl;
       eintarray finalarray;
       eintarray sortedindex;
       eintarray tmparr;
       for (int k=0;k<count;k++){finalarray.add(0);tmparr.add(0);}
       for (int k=0;k<count;k++){finalarray[k]=temparray[k];}
       for (int k=0;k<count;k++){rw.disable(finalarray[k]);}
       cout<<intarr2str2(finalarray)<<endl;
       rw.calcPhenotype();
       cout<<"alak4.3"<<endl;
       if (rw.isViable()) {
           cout<<"alak4.4"<<endl;
	   sortedindex = sort(finalarray);// sort returns the sorted index list
	   for (int p=0;p<sortedindex.size();p++) {tmparr[p]=finalarray[sortedindex[p]];}
	   finalarray = tmparr;
           for (int q=0;q<tmparr.size();q++) {tmparr[q] -= 26;}
           estr zeroloc = intarr2str2(tmparr);
           fileout.write(zeroloc+"\n");
           cout<<"alak4.5"<<endl;
       }
       for (int n=0; n<finalarray.size(); ++n) {rw.activate(finalarray[n]);}
       cout<<"alak4.6"<<endl;
  }

  fileout.close();
  return(0);
}
