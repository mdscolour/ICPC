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
double outer_define = atof(argv[2]);

sprintf(grname, "chr2-%s.midgr", itarea);
sprintf(lowname, "chr2-%s.lowconfig", itarea);
sprintf(midname, "chr2-%s.midconfig", itarea);
sprintf(highname, "chr2-%s.highconfig", itarea);

char buffer[50]; //only storage for 256 characters.
sprintf(buffer, "potES1%.5f.midpot",outer_define); 

ClassErSearch er1 = ClassErSearch();
er1.doErSearch(outer_define,grname,lowname,buffer);

return 0; 
} 



















