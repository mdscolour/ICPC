//#include "/home/li/bin/coreMolecularMC.h"
//#include "coreMolecularMC.h"
//#include "/home/li/bin/coreMolecularMC.cpp"//only one instance
#include "coreMolecularMC.cpp"//only one instance


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

ClassErSearch::ersearchx = arange(100,300,1);
ClassErSearch::ersearchy.clear();
for(int i=0;i<ClassErSearch::ersearchx.size();i++)
    ClassErSearch::ersearchy.push_back(
        LJlike(ClassErSearch::ersearchx[i],167,6,13,11));

coreDirectRun::loadTarget("standardLJTest.gr");
int ngr = int(coreDirectRun::targetX.size());
double dx=coreDirectRun::dx;
int lgr = int(double(ngr) * dx);
ClassErSearch::runPotInput(1,40,100,0,15000,74,10000,10,500,ngr,dx);
write2DVector("167LJ100kTest.gr",arange(0,lgr,dx),ClassErSearch::best);
cout<<"167:"<<coreDirectRun::getDiff(ClassErSearch::best)<<endl;

ClassErSearch::ersearchx = arange(100,300,1);
ClassErSearch::ersearchy.clear();
for(int i=0;i<ClassErSearch::ersearchx.size();i++)
    ClassErSearch::ersearchy.push_back(
        LJlike(ClassErSearch::ersearchx[i],162,3,12,6));

ClassErSearch::runPotInput(1,40,100,0,15000,74,10000,10,500,ngr,dx);
write2DVector("162LJ100kTest.gr",arange(0,lgr,dx),ClassErSearch::best);
cout<<"162:"<<coreDirectRun::getDiff(ClassErSearch::best)<<endl;

ClassErSearch::ersearchx = arange(100,300,1);
ClassErSearch::ersearchy.clear();
for(int i=0;i<ClassErSearch::ersearchx.size();i++)
    ClassErSearch::ersearchy.push_back(
        LJlike(ClassErSearch::ersearchx[i],165,7,11,9));
ClassErSearch::runPotInput(1,40,100,0,15000,74,10000,10,500,ngr,dx);
write2DVector("165LJ100kTest.gr",arange(0,lgr,dx),ClassErSearch::best);
cout<<"165:"<<coreDirectRun::getDiff(ClassErSearch::best)<<endl;

ClassErSearch::ersearchx = arange(100,300,1);
ClassErSearch::ersearchy.clear();
for(int i=0;i<ClassErSearch::ersearchx.size();i++)
    ClassErSearch::ersearchy.push_back(
        LJlike(ClassErSearch::ersearchx[i],161,3,12,6));

ClassErSearch::runPotInput(1,40,100,0,15000,74,10000,10,500,ngr,dx);
write2DVector("161LJ100kTest.gr",arange(0,lgr,dx),ClassErSearch::best);
cout<<"161:"<<coreDirectRun::getDiff(ClassErSearch::best)<<endl;

//double scaling=0.45;
//if(argc>=2) scaling=atof(argv[1]);

/*    // function of potential, temperature, step size
	coreMolecularMC MMC=coreMolecularMC(&potential,1,0.45);
    
    const char* saveGrName="targetGr.gr";
    int num = 200;
    //double spacing = MMC.getPotMin(30,0.01);

    //MMC.setBox(0,(num+1)*spacing);
    MMC.setBox(0,800);
    //for(int i=0;i<num;i++)MMC.addPart(spacing*i);
    
    MMC.addAllPartEvenly(num);
    //MMC.addAllPartByMinPot(num,30,0.01);
    MMC.readyEverything();
    
    int num_configs=100;
    int num_steps=1;
    int init_step = 500;
    //int each_step = 1;
    MMC.run(init_step);
    
    int npart = MMC.getNpart();
    int nsteps = 0;
    
    //FILE *fptrinit;
    //fptrinit = fopen("inittest.dat", "w");
    
    int lgr = 35;
    double dx=0.01;
    int ngr = int(double(lgr)/dx);
    vector<double> gr(ngr,0.0);//=MMC.getGr(ngr,dx);
    //printfVector(gr);
    
    FILE *fptr;
    fptr = fopen("test.out", "w");
    for(int i=0;i<num_configs;i++)
    {
        //MMC.run(num_steps);
        nsteps+=num_steps;
        VectorAddedByVector(gr,MMC.getGr(ngr,dx));
                
        //spacing = MMC.getSpacing();
        //if(spacing!=spaold)printf("%d   %.8f    %.8f \n",i,MMC.part[num-1],MMC.part[0]);
        //spaold=spacing;
        //fprintf(fptrinit, "%d %.5f\n",i,spacing/num);
        //printf("initializaiton %d steps....done, %.6f spacing\n",nsteps,spacing/num);
        //printf("configuration %d in %d....done.\n",i,num_configs);
        
        fprintf(fptr, "ITEM: TIMESTEP\n%d\n",nsteps);
        fprintf(fptr, "ITEM: NUMBER OF ATOMS\n%d\n",npart);
        fprintf(fptr, "ITEM: BOX BOUNDS pp pp pp\n");
        fprintf(fptr, "0 1\n");
        fprintf(fptr, "0 1\n");
        fprintf(fptr, "0 1\n");
        fprintf(fptr, "ITEM: ATOMS id xs ys zs\n");
        for(int i=0;i<npart;i++){
            fprintf(fptr, "%d %.5f %.5f %.5f\n",i, MMC.getPart(i), 0.0, 0.0);}

        MMC.run(num_steps);
    }
    fclose(fptr);
    //fclose(fptrinit);
    VectorDividedByScalar(gr,num_configs);
    write2DVector(saveGrName,arange(0,lgr,dx),gr);
*/
/*    typedef coreDoubleRMC rmc; 
    
    double rndstep = 10;
    if(argc>=2) rndstep=atof(argv[1]);
    
    vector<double> arrX;// = {0,10,20,30};arrX.push_back(100);
    arrX = linspace(0,30,31);arrX.push_back(1000);
    //coreDoubleRMC* rmc = new coreDoubleRMC(arrX);
    rmc::initPotentialX(arrX);
    //static void asignConfig(double talpha,int tnum_configs,int tnum_steps,
    //                    int tinit_step,int tnum_particles,
    //                    int tmove_scale,int tmove_cut)
    
    rmc::asignConfig(0.5,100,10,5000,200,2*rndstep,rndstep);
    rmc::readyToRun();
    rmc::run(1000);*/
    
//     for(int i=0;i<10;i++)
//     {
//         rmc::run(10000);
//         
//         char buffer[50]; // <- danger, only storage for 256 characters.
//         sprintf(buffer, "savedpotential%d_10000_5k.dat", i);
//         rmc::writePot(buffer);
//     }
    
//     arrX = linspace(0,30,10);arrX.push_back(100);
//     rmc::updatePotX(arrX);
//     rmc::run(4);
//     
//     arrX = linspace(0,30,31);arrX.push_back(100);
//     rmc::updatePotX(arrX);
//     rmc::run(2);
    
    
    
    return 0; 
} 



















