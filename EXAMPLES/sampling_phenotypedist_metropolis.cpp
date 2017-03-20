#include "enet.h"
#include "erandomwalk.h"
#include <eutils/logger.h>
#include <eutils/emain.h>
#include <signal.h>
#include <eutils/eheap.h>
#include <algorithm>
#include <eutils/ernd.h>
#include <math.h>
///////////////////////////////////////////////////////////////
estr solver="esolver_clp";
int netsize=-1;
int strict=0;
int force_away=0;
int periphery_only=0;
int mutate_transport=0;
int internal_secretion=0;
int only_viable=0;
enet net;
double delta=0;//is the desired genotypic distance between the resulting parental genotypes
////////////////////////////////////////////////////////////////
erandomWalk *prw=0x00;
////////////////////////////////////////////////////////////////
int emain()
{ ldieif(argvc<5,"syntax: ./sampling_phenotypedist_metropolis  <kegg.net> <genotype1.dat> <genotype2.dat> <outnet.dat> --delta <fluxbounds.flx>");

  epregister(solver);
  epregister(netsize);
  epregister(strict);
  epregister(force_away);
  epregister(periphery_only);
  epregister(mutate_transport);
  epregister(internal_secretion);
  epregister(only_viable);
  epregister(delta);
   eparseArgs(argvc,argv);

  net.load(argv[1]);//Loading the reaction universe (kegg.net)
  efile filein1;//The file containing the initial donor genotype
  efile filein2;//The file containing the initial recipient genotype
  efile fileout;//The output file that will store the final parental genotypes (donor and recipient) with the specified genotypic distance (i.e. delta)

  filein1.open(argv[2],"r");
  filein2.open(argv[3],"r");

  estr sttr1;
  estr sttr2;
  estr sttr3;
  eintarray gen1;//initial donor genotype
  while (filein1.readln(sttr1)) {//reading the file containing the initial donor genotype
        estrarray parts1;
        parts1=sttr1.explode(" ");
        int tmp=0;
        for (int i=0;i<6588;i++){gen1.add(tmp);}
        for (int i=0;i<6588;i++){gen1[i]=parts1[i].i();}
  }

  eintarray gen2;//initial recipient genotype
  while (filein2.readln(sttr2)) {//reading the file containing the initial recipient genotype
        estrarray parts2;
        parts2=sttr2.explode(" ");
        int tmp=0;
        for (int i=0;i<6588;i++){gen2.add(tmp);}
        for (int i=0;i<6588;i++){gen2[i]=parts2[i].i();}
  }
  

  filein1.close();
  filein2.close();

  erandomWalk rw(net,solver,strict);
  prw=&rw; 

  rw.periphery_only=periphery_only;
  rw.mutate_transport=mutate_transport;
  rw.internal_secretion=internal_secretion;
  rw.only_viable=only_viable;
  rw.setRSize(netsize);
  rw.getEnv(argvc,argv);//finding the environment files and read them

  for (int i=0;i<6588;i++){if (gen1[i]==0){rw.disable(i);}}//disable absent reaction in the donor genotype
  rw.calcPhenotype();//calculate phenotype
  eintarray phen1=rw.phenotype;//store the phenotype
  for (int i=0;i<6588;i++){rw.activate(i);}//reactivate all the reactions in the universe

  for (int i=0;i<6588;i++){if (gen2[i]==0){rw.disable(i);}}//disable absent reaction in the recipient genotype
  rw.calcPhenotype();//calculate phenotype
  eintarray phen2=rw.phenotype;//store the phenotype
  for (int i=0;i<6588;i++){rw.activate(i);}//reactivate all the reactions in the universe
  
  eintarray GEN1=gen1;
  eintarray GEN2=gen2;
  eintarray PHEN1=phen1;
  eintarray PHEN2=phen2;


  double gendist=gendistance2(gen1,gen2,&rw);//calculate the genotypic distance between the initial donor and recipient genotypes

  double dist_new=10000.0;
  double dist_old=0.0;

  int present1_size=0;
  int present2_size=0;
  int absent1_size=0;
  int absent2_size=0;

  int gen1_present[6588];for (int i=0;i<6588;i++){gen1_present[i]=0;}
  int gen2_present[6588];for (int i=0;i<6588;i++){gen2_present[i]=0;}
  int gen1_absent[6588]; for (int i=0;i<6588;i++){gen1_absent[i]=0;}
  int gen2_absent[6588]; for (int i=0;i<6588;i++){gen2_absent[i]=0;}

  int d1=0;
  int d2=0;

  int to_delete1=0;
  int to_add1=0;
  int to_delete2=0;
  int to_add2=0;
  int terminator=0;
    

  while (terminator==0) {

        present1_size=0;
        present2_size=0;
        absent1_size=0;
        absent2_size=0;

        for (int i=0;i<6588;i++){
             if (gen1[i]==1&&gen2[i]==0){gen1_present[present1_size]=i;present1_size++;}//stores the set of reactions that are present in the donor genotype and absent in the recipient genotype.
             if (gen1[i]==0&&gen2[i]==1){gen2_present[present2_size]=i;present2_size++;}//stores the set of reactions that are present in the recipient genotype and absent in the donor genotype.

             if (gen1[i]==0){gen1_absent[absent1_size]=i;absent1_size++;}//stores the set of reactions that are absent in the donor genotype
             if (gen2[i]==0){gen2_absent[absent2_size]=i;absent2_size++;} //stores the set of reaction sthat are absent in the recipient genotype.
        }


         dist_old=gendistance2(gen1,gen2,&rw);//old genotypic distance between the donor and the recipient genotypes
         cout<<"###############################################################"<<endl;
         cout<<"Old: "<<dist_old<<endl;
         cout<<"###############################################################"<<endl;
         d1=0;
         to_delete1=0;
         to_add1=0;
         to_delete2=0;
         to_add2=0;

         GEN1=gen1;//copy of the donor genotype
         while (d1==0){
               ///////////// Proposed genetic changes in the donor genotype ///////////////
               to_delete1=(int)(ernd.uniform()*present1_size);//randomly select a present reaction in the donor genotype to delete
               to_add1=(int)(ernd.uniform()*absent1_size); //randomly select an absent reaction in the donor genotype to add
               gen1[(gen1_present[to_delete1])]=0;//change the donor genotype accordingly
               gen1[(gen1_absent[to_add1])]=1;//change the donor genotype accordingly
             
               for (int i=0;i<6588;i++){if (gen1[i]==0){rw.disable(i);}}//disable the absent reactions in the donor genotype
               rw.calcPhenotype();//calculate phenotype
               eintarray phen1=rw.phenotype;//store the phenotype
               if (phen1==PHEN1){d1=1;}//we check whether the phenotype of modified donor is the same as the original donor; if yes d=1, the change is accepted and and the loop terminates.
               else {gen1=GEN1;}//otherwise we do not accept the change and the donor genotype is changed back to the previous donor genotype (GEN1 is the copy of the previous donor genotype).
               for (int i=0;i<6588;i++){rw.activate(i);}//reactivate the deleted genes in the universe.
         }
         d2=0;
         GEN2=gen2;
         while (d2==0){
               ///////////// Proposed genetic changes in the recipient genotype ///////////////
               to_delete2=(int)(ernd.uniform()*present2_size);//randomly select a present reaction in the recipient genotype to delete
               to_add2=(int)(ernd.uniform()*absent2_size);//randomly select an absent reaction in the recipient genotype to add
               gen2[(gen2_present[to_delete2])]=0;//change the recipient genotype accordingly
               gen2[(gen2_absent[to_add2])]=1;//change the recipient genotype accordingly

               for (int i=0;i<6588;i++){if (gen2[i]==0){rw.disable(i);}}//disable the absent reactions in the recipient genotype
               rw.calcPhenotype();//calculate phenotype
               eintarray phen2=rw.phenotype;//store the phenotype
               if (phen2==PHEN2){d2=1;}//we check whether the phenotype of modified recipient is the same as the original recipient; if yes d=1, the change is accepted and and the loop terminates.
               else {gen2=GEN2;}//otherwise we do not accept the change and the recipient genotype is changed back to the previous recipient genotype (GEN2 is the copy of the previous recipient genotype).
               for (int i=0;i<6588;i++){rw.activate(i);}//reactivate the deleted genes in the universe.
        }
        dist_new=gendistance2(gen1,gen2,&rw);//new genotypic distance between the new donor and the new recipient genotypes
        cout<<"###############################################################"<<endl;
        cout<<"New: "<<dist_new<<endl;
        cout<<"###############################################################"<<endl;

        if (dist_new>delta) {if (dist_new>dist_old){gen1=GEN1;gen2=GEN2;}}// If the new genotypic distance is larger than the desired final genotypic distance (delta) and is also laregr than the old distance the change is not accepted and both the donor and recipient genotypes are changed backj to the previous genotypes.
        else if (dist_new<delta){if (dist_new<dist_old){gen1=GEN1;gen2=GEN2;}}//If the new genotypic distance is smaller than the desired final genotypic distance (delta) and is also smaller than the old distance the change is not accepted and both the donor and recipient genotypes are changed backj to the previous genotypes.
        else if (dist_new==delta){terminator=1;break;}//If the new genotypic distance is equal to the desired genotypic distance (delta), the loop terminates and the current genotypes are considered as our final donor and recipient genotypes that have the same phenotype as the original genotypes but with a specified *i.e. delta) genotypic distance).
        else {}
  }


  fileout.open(argv[4],"w");//open the output file

  fileout.write(intarr2str2(gen1)+"\n");//write the donor genotype to the output file.
  fileout.write(intarr2str2(gen2)+"\n");//write the recipient genotype to the output file.
  
  fileout.close();//close the output file

  return(0);
}
