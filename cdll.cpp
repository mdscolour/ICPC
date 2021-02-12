//#include "/home/li/bin/coreMolecularMC.h"
//#include "coreMolecularMC.h"
//#include "/home/li/bin/coreMolecularMC.cpp"//only one instance
#include "coreMolecularMC.h"//only one instance

extern "C" { 
//     void c_xAdd(double x){ ClassErSearch::ersearchx.push_back(x); }  
//     void c_yAdd(double y){ ClassErSearch::ersearchy.push_back(y); } 
//     void c_xClear(){ ClassErSearch::ersearchx.clear(); } 
//     void c_yClear(){ ClassErSearch::ersearchy.clear(); } 
//     bool c_xyEqual(){return (ClassErSearch::ersearchx.size()==ClassErSearch::ersearchy.size());}
//     
//     void c_runPotInput(double t,double scal,double minvol, double tbox_l,double tbox_r,int npart,int num_configs,int num_steps,int init_step,int ngr,double dx){
//         ClassErSearch::runPotInput(t,scal,minvol,tbox_l,tbox_r,npart,
//                     num_configs,num_steps,init_step,ngr,dx);}
//     double c_getBest(int pos){return ClassErSearch::best[pos];}
//     
//     void c_SeedByTime(){SeedByTime();}
//     void c_loadTarget(){coreDoubleRMC::loadTarget("targetGr.gr");}
//     
//     void c_xPrint(){printfVector(ClassErSearch::ersearchx);}
//     void c_yPrint(){printfVector(ClassErSearch::ersearchy);}
//     void c_bestPrint(){printfVector(ClassErSearch::best);}
    
    void c_loadTarget(char* tag){coreDirectRun::loadTarget(tag);}
    void c_readConfig(char* cog){coreDirectRun::readConfig(cog);}
    void c_readyToRun(){coreDirectRun::readyToRun();}
    double c_assignAndRun(double p1,double p2,double p3,double p4)
    {return coreDirectRun::assignAndRun(p1,p2,p3,p4);}
    double c_getBestGr(int pos){return coreDirectRun::savegr[pos];}
}




