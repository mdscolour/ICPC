//#include "/home/li/bin/coreMolecularMC.h"
//#include "coreMolecularMC.h"
//#include "/home/li/bin/coreMolecularMC.cpp"//only one instance
#include "coreMolecularMC.h"//only one instance


double LJPotential(double r,double epsilon=1,double sigma=1) {
    return 4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));}
double QuarticBond(double r,double k0,double k1,double R){
    return 0.5*k0*pow((r-R),2)+0.25*k1*pow(r-R,4);}
double FeneBond(double r,double k,double drm,double r0){
    return -0.5*k*drm*drm*log(1-pow((r-r0)/drm,2));}
//double potential(double r){return LJPotential(r,1,1);}

double potential(double r){
    //return FeneBond(r,1,30,3)+QuarticBond(r,1,1,10)+FeneBond(r,-1,30,8);
    return LJPotential(r,3,3);
}

int main(int argc,char *argv[])
{    /**************  random seed  *******************************/
SeedByTime();

char grname[50];
char lowname[50];
char highname[50];
char* itarea = new char[1];
itarea[0]='0';

sprintf(grname, "chr2-%s.LJgr", itarea);
sprintf(lowname, "chr2-%s.lowconfig", itarea);
//sprintf(highname, "chr2-%s.highconfig", itarea);

coreDirectRun::loadTarget(grname);
coreDirectRun::readConfig(lowname);
coreDirectRun::readyToRun();

double energy = coreDirectRun::assignAndRun(162,3,12,6);

vector<double> gr = coreDirectRun::savegr;

int ngr = int(coreDirectRun::targetX.size());
double dx = coreDirectRun::dx;
int lgr = int(double(ngr) * dx);
        
write2DVector("chr2-0.testLJgr", arange(0, lgr, dx), gr); 
    
    return 0; 
} 



















