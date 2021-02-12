//#include "/home/li/bin/coreMolecularMC.h"
//#include "coreMolecularMC.h"
//#include "/home/li/bin/coreMolecularMC.cpp"//only one instance
#include "coreMolecularMC.h"//only one instance

int main(int argc,char *argv[])
{    /**************  random seed  *******************************/
SeedByTime();

char grname[50];
char lowname[50];char midname[50];char highname[50];
char* itarea;

if(argc>=2) itarea=(argv[1]);
double outer_define1 = atof(argv[2]);

sprintf(grname, "chr2-%s.LJgr", itarea);
sprintf(lowname, "chr2-%s.lowconfig", itarea);
sprintf(midname, "chr2-%s.midconfig", itarea);
sprintf(highname, "chr2-%s.highconfig", itarea);

char buffer[50]; //only storage for 256 characters.
sprintf(buffer, "chr2-0LJHCM2/potSAM2%f.LJpot",outer_define1); 

HCMethod2(grname,lowname,buffer);
    
return 0; 
} 



















