#include "enet.h"
#include "erandomwalk.h"
#include <eutils/logger.h>
#include <eutils/emain.h>
#include <signal.h>
#include <eutils/eheap.h>
#include <algorithm>
#include <eutils/ernd.h>

estr solver="esolver_clp";

int strict=0;
int periphery_only=0;
int mutate_transport=0;
int internal_secretion=0;
int only_viable=0;
enet net;
erandomWalk *prw=0x00;

int num1=0; //number of reactions in the parental networks
double num2=0.0; //genotypic distance between parental genotypes
int num3=0; //number of exchanged reactions
estr outnet="out.net";
int iter=0;

int emain()
{ ldieif(argvc<4,"syntax: ./recombination_genomescale_phen_dist <universe.net> <file.dat> --outnet  --num1 --num2 --num3  --iter <fluxbounds.flx>");  
  epregister(num1);
  epregister(num2);
  epregister(num3);
  epregister(iter);
  epregister(outnet);
  eparseArgs(argvc,argv);
  epregister(solver);
  epregister(strict);
  epregister(periphery_only);
  epregister(mutate_transport);
  epregister(internal_secretion);
  epregister(only_viable);
  int netsize=num1-1;  

  //////////////////////////////////////////////////////////Genotyping////////////////////////////////////////////////////
  enet net;
  net.load(argv[1]);//loading universe
  net.correct_malformed();
  erandomWalk rw(net,solver,strict);
  prw=&rw; 
  rw.periphery_only=periphery_only;
  rw.mutate_transport=mutate_transport;
  rw.internal_secretion=internal_secretion;
  rw.only_viable=only_viable;
  rw.getEnv(argvc,argv);
  rw.load(net);  
  net.correct_malformed();
  efile fu;
  fu.open(argv[2],"r");//opens the file containing the genotypes of the donor and the recipient parents
  efile fw;
  fw.open(outnet,"a");//opens the file that will store the information about innovation and robustness of the recombinant offspring
  estr sttr;
  estr outnet2=outnet+"_phen";//opens the file that will store the information about the phenotype of the recombinant offspring
  efile fw2;
  fw2.open(outnet2,"a");
  //////////////////////
  int id=1;
  eintarray gen1;//genotype of the donor
  eintarray gen2;//genotype of the recipient
  eintarray GEN1;//copy of the genotype of the donor
  eintarray GEN2;//copy of the genotype of the recipient
  eintarray phen1;//phenotype of the donor
  eintarray phen2;//phenotype of the recipient
  int tmp=0;
  while (fu.readln(sttr)) {//read the file containg the genotype of the parents
	    estrarray parts;
        parts=sttr.explode(" ");
        if (id==1){//The fist line of the genotype file
            for (int i=0;i<6588;i++){gen1.add(tmp);}
            for (int i=0;i<6588;i++){gen1[i]=parts[i].i();}//convert string to integer [store the genotype as a binary vector]
            for (int i=0;i<6588;i++){if (gen1[i]==0){rw.disable(i);}}//disable reactions in the universe that are not present in the donor genotype
            rw.periphery_only=periphery_only;
            rw.mutate_transport=mutate_transport;
            rw.internal_secretion=internal_secretion;
            rw.only_viable=only_viable;
            rw.setRSize(netsize);
            rw.calcPhenotype();//calculate the phenotype of the donor
            phen1=rw.phenotype;//store the phenotype of the donor

            fw2.write(intarr2str2(phen1)+"\n");//write the phenotype of the donor to the phenotype file
            GEN1=gen1;//make a copy of the donor genotype
            for (int i=0;i<6588;i++){rw.activate(i);}//reactivate all the reactions in the universe
        }
        else if (id==2){//The second line of the genotype file
                 for (int i=0;i<6588;i++){gen2.add(tmp);}
                 for (int i=0;i<6588;i++){gen2[i]=parts[i].i();}//convert string to integer [store the genotype as a binary vector]
                 for (int i=0;i<6588;i++){if (gen2[i]==0){rw.disable(i);}}//disable reactions in the universe that are not present in the recipient genotype
                 cout<<"alaki6"<<endl;
                 rw.periphery_only=periphery_only;
                 rw.mutate_transport=mutate_transport;
                 rw.internal_secretion=internal_secretion;
                 rw.only_viable=only_viable;
                 rw.setRSize(netsize);
                 rw.calcPhenotype();//calculate the phenotype of the recipient
                 phen2=rw.phenotype;//store the phenotype of the recipient
                 fw2.write(intarr2str2(phen2)+"\n");//write the phenotype of the recipient to the phenotype file
                 GEN2=gen2;//make a copy of the recipient genotype
                 for (int i=0;i<6588;i++){rw.activate(i);}//reactivate all the reactions in the universe
        }
        else {}
        id++;
  }
  for (int countering=0;countering<iter;countering++){//countering is the number of recombinant offspring [sorry for the wierd naming. countering???]
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       ////////////////////////////////////////////////////////      Recombination     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       eintarray gen2=GEN2;//the recipient genotype
      
       ////////////////////////////////////////////////////////      Donor     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       cout<<"#############################################"<<endl;
       eintarray pop;//an array storing randomly shuffled integers from 0 to (num2-1);
       for (int i=0; i<num2; ++i) {pop.add(i);}
       int tmp1; 
       for (int i=(num2-1); i>=0; --i) {// to shuffle the initially ordered array
            int r = (int)(ernd.uniform()*i);
            tmp1 = pop[r];
            pop[r] = pop[i];
            pop[i]=tmp1;
       }
       eintarray pop1;
       for (int j=0;j<num3;j++){//Select the first num3 out of num2 elements from pop array
            int x=pop[j];
            pop1.add(x);
       }
       eintarray sorted1=sort(pop1);//getting the index of the sorted pop1
       eintarray selected1;//finally selected elements
       for (int j=0;j<num3;j++){
            int y=sorted1[j];
            selected1.add(pop1[y]);
       }
      
      
      
      ////////////////////////////////////////////////////////     Recipient     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       eintarray popp;//The vector storing the indices of the reactions that are present in the recipient genotype
       for (int i=682; i<6588; ++i) {
            if (gen2[i]==1){popp.add(i);}
       }
 
       for (int i=(num1-682); i>=0; --i) {// Randomly shuffle the popp array. note 682 are the number of transport reactions thgat we kept constant in this analysis.
            int r = (int)(ernd.uniform()*i);
            tmp1 = popp[r];
            popp[r] = popp[i];
            popp[i]=tmp1;
       }
       eintarray pop2;
       for (int j=0;j<num3;j++){//Select the first num3 out of (num1-682) elements from pop array
           int x=popp[j];
           pop2.add(x);
       }
       eintarray sorted2=sort(pop2);//getting the index of the sorted pop2
       eintarray selected2;//finally selected elements
       for (int j=0;j<num3;j++){
           int y=sorted2[j];
           selected2.add(pop2[y]);
       }
       selected1.add(0);
       selected2.add(0);    
       ////////////////////////////final indexing//////////////////////////
       eintarray final1;// to determine the index of the genes that are present in the donor genotype and are absent in the recipient genotype.
       for (int i=0;i<6588;i++){
           if (gen1[i]==1){
               if (gen2[i]==0){
                  final1.add(i);
               }
           }
       }
      /////////////////Recombinants and outputing////////////////////////
      int count1=0;
      int count2=0;
      for (int i=0;i<6588;i++){
        if (gen2[i]==1){  
           if (selected2[count2]==i){gen2[i]=0;count2++;}//to deletle the num3 selected reactions from the recipient (changing 1 to 0).
           else {}
        }

        else if (gen2[i]==0) {
             if (final1[selected1[count1]]==i){
                 gen2[i]=1;//to add the num3 selected reactions from the donor to the recipient (changing 0 to 1).
                 count1++;
             }
        }
        else {}
      }	
      ///////////////////////////////Phenotyping////////////////////////
      for (int i=0;i<6588;i++){if (gen2[i]==0){rw.disable(i);}}//disable the absent reactions in the offspring genotype
      rw.calcPhenotype();//calculate phenotype
      eintarray phen3 = rw.phenotype;//store the phenotype of the offspring.
      for (int i=0;i<6588;i++){if (gen2[i]==0){rw.activate(i);}}//reactivate the absent reactions in the universe
      double gain1=0;//first type of gain of phenotype: gaining viability on a carbon source on which neither of the parents are viable on.
      double gain2=0;//second type of gain of phenotype: gaining viability on a carbon source on which only the donor genotype is viable on.
      double loss1=0;//first type of the loss of phenotype:
      double loss2=0;//second type of the loss of phenotype:
      for (int k=0;k<50;k++){//considering all 50 carbon sources
          if (phen1[k]==1){
              if (phen2[k]==1){if (phen3[k]==0){loss1=loss1+1;}}
              else {if (phen3[k]==1){gain2=gain2+1;}}
          }
          else {
              if (phen2[k]==1){if (phen3[k]==0){loss2=loss2+1;}}
              else {if (phen3[k]==1){gain1=gain1+1;}}
          }
      }
      estr g1=gain1;
      estr g2=gain2;
      estr l1=loss1;
      estr l2=loss2;
      estr intstr=g1+" "+g2+" "+l1+" "+l2;
      fw.write(intstr+"\n");//writing to the first file
      fw2.write(intarr2str2(phen3)+"\n"); //writing to the phenotype file
  }
  fw.close();
  fw2.close();
  return(0);
}


