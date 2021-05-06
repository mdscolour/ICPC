//#include "/home/li/bin/coreMolecularMC.h"
//#include "coreMolecularMC.h"
//#include "/home/li/bin/coreMolecularMC.cpp"//only one instance
#include "coreMolecularMC.h"//only one instance

int main(int argc,char *argv[])
{    /**************  random seed  *******************************/
SeedByTime();

double p1=atof(argv[1]);
double p2=atof(argv[2]);
double p3=atof(argv[3]);
double p4=atof(argv[4]);
int ngrpoint=atoi(argv[5]);
double dxt=atof(argv[6]);

coreDirectRun::assignAndWriteV9(p1,p2,p3,p4,ngrpoint,dxt,argv[7],argv[8]);

return 0; 
} 



















