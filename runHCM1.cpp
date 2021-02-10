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
int abindex = atoi(argv[2]);
int abcount=0;
double outer_define1,outer_define2;

for(double a=2;a<19;a++)
for(double b=1;b<a;b++)
{
    //cout<<a<<"  "<<b<<endl;
    outer_define1 = a;
    outer_define2 = b;
    if(abcount==abindex)break; //total 152, if up to 18
    abcount++;
}

sprintf(grname, "chr2-%s.LJgr", itarea);
sprintf(lowname, "chr2-%s.lowconfig", itarea);
sprintf(midname, "chr2-%s.midconfig", itarea);
sprintf(highname, "chr2-%s.highconfig", itarea);

char buffer[50]; //only storage for 256 characters.
sprintf(buffer, "potHCM1%.5f_%.5f.LJpot",outer_define1,outer_define2); 

HCMethod1(outer_define1,outer_define2,grname,midname,buffer);
    
return 0; 
} 



















