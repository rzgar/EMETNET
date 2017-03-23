#include "enet.h"
#include "erandomwalk.h"
#include <eutils/logger.h>
#include <eutils/emain.h>
#include <signal.h>
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
estr solver="esolver_clp";
int netsize=-1;
int strict=0;
int force_away=0;
int periphery_only=0;
int mutate_transport=0;
int internal_secretion=0;
int only_viable=0;
enet net;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
erandomWalk *prw=0x00;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int emain()
{
  ldieif(argvc<5,"syntax: ./metnet-exgen <universe.net> <inputfilename.dat> <fluxbounds.flx> <outputfilename>");

  epregister(solver);
  eparseArgs(argvc,argv);

  net.load(argv[1]);
  efile inputfile;
  efile outputfile;

  inputfile.open(argv[2],"r");
  outputfile.open(argv[4],"a");

  net.correct_malformed();
  erandomWalk rw(net,solver,strict);
  prw=&rw; 
  rw.periphery_only=periphery_only;
  rw.mutate_transport=mutate_transport;
  rw.internal_secretion=internal_secretion;
  rw.only_viable=only_viable;
  rw.setRSize(netsize);
  
  rw.getEnv(argvc,argv);
  rw.load(net);
  rw.calcPhenotype();
  rw.viablePhenotype=rw.phenotype;

  estr str;
  estrarray parts;

  while (inputfile.readln(str)) {			
  	parts=str.explode(" ");
	eintarray numarray, writearr;
	int tmp = 0;
        int tmpo= 0;
	for (int i=0; i<parts.size(); ++i){
      	     tmp = parts[i].i()+26;
             tmpo= parts[i].i();
             writearr.add(tmpo);						
	     rw.disable(tmp);
             numarray.add(tmp);
	}
	rw.calcPhenotype();
	if (rw.isViable()) {estr zeroloc = intarr2str2(writearr);outputfile.write(zeroloc+"\n");}
	for (int i=0; i<numarray.size(); ++i){rw.activate(numarray[i]);}
   }
   inputfile.close();
   return(0);
}

