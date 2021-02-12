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
sprintf(grname, "../chr2-%s.LJgr", itarea);
sprintf(lowname, "../chr2-%s.lowconfig", itarea);
sprintf(midname, "../chr2-%s.midconfig", itarea);
sprintf(highname, "../chr2-%s.highconfig", itarea);

char canname[50];
sprintf(canname, "pot.can");
FILE *fptr = fopen(canname, "r");
double tmp,p1,p2,p3,p4,eng;
double oldenergy=99999,energy=99999;
vector<double> best;

coreDirectRun::loadTarget(grname);
coreDirectRun::readConfig(midname);
coreDirectRun::readyToRun();

FILE *fptrwrite;
fptrwrite = fopen("pot.res", "w");

int runcount = 0;
while (fscanf(fptr, "%lf %lf %lf %lf %lf %lf\n", &tmp,&p1,&p2,&p3,&p4,&eng) == 6)
{
    cout<<(runcount)<<"    ";
    if(runcount>=100000)break;
    energy = coreDirectRun::assignAndRun(p1,p2,p3,p4);
    printf("%lf %lf %lf %lf %lf %lf\n",tmp,p1,p2,p3,p4,energy);
    fprintf(fptrwrite, "%f %.5f %.5f %.5f %.5f %.5f\n",0,p1,p2,p3,p4,energy);
    if(energy<oldenergy)
    {
        oldenergy = energy;
        best.clear();
        best.push_back(p1);
        best.push_back(p2);
        best.push_back(p3);
        best.push_back(p4);
        cout<<oldenergy<<"    ";
        printfVector(best);
    }
    runcount++;
}
fprintf(fptrwrite, "best is:\n");
fprintf(fptrwrite, "%f %.5f %.5f %.5f %.5f %.5f\n",0,best[0],best[1],best[2],best[3],energy);

return 0; 
} 



















