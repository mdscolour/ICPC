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

double global_sigma=0;
double global_epsilon=0;
double global_a=0;
double global_b=0;
double potential(double r){
    if(global_sigma==0)throw range_error("no initialization");
    //return FeneBond(r,1,30,3)+QuarticBond(r,1,1,10)+FeneBond(r,-1,30,8);
//     return LJlike(r,150,8,10,9);
    return LJlike(r,global_sigma,global_epsilon,global_a,global_b);
}

int main(int argc,char *argv[])
{    /**************  random seed  *******************************/
SeedByTime();

char grname[50];
char savegrname[50];
char lowname[50];
char highname[50];
char conforname[50];
char finresname[50];

sprintf(grname, "chr%s-%s.midgr", argv[2],argv[1]);
sprintf(lowname, "chr%s-%s.lowconfig",argv[2],argv[1]);
sprintf(conforname, "chr%s-%s.confor",argv[2],argv[1]);
sprintf(finresname, "chr%s-%s.finrestmp",argv[2],argv[1]);

FILE *fptrpara = fopen(finresname, "r");
while(fscanf(fptrpara, "%lf %lf %lf %lf\n", 
    &global_sigma,&global_epsilon, &global_a, &global_b) != 4){}
fclose(fptrpara);
// global_sigma=164;
// global_epsilon=1;
// global_a=14;
// global_b=13;
printf("%lf %lf %lf %lf\n", global_sigma,global_epsilon, global_a, global_b);

coreDirectRun::loadTarget(grname);//not actually needed
coreDirectRun::readConfig(lowname);
coreDirectRun::readyToRun();
// write2DVector("chr2-0.testLJgr", arange(0, lgr, dx), gr); 
    
	coreMolecularMC MMC=coreMolecularMC(&potential,coreDirectRun::global_temperature, 
            coreDirectRun::global_scaling, coreDirectRun::global_minvol);

    MMC.setBox(coreDirectRun::global_boxl,coreDirectRun::global_boxr);
    MMC.addAllPartEvenly(coreDirectRun::global_npart);
    //cout<<MMC.part[0]<<endl;
    //MMC.part[0]=76;
    MMC.readyEverything();
    
    //const char* saveGrName="res/chr2-0.gr";
    //double spacing;
    
    int num_configs=15000;//coreDirectRun::global_num_configs;
    int num_steps=10;//coreDirectRun::global_num_steps;
    int init_step=5000;//coreDirectRun::global_init_step;
    MMC.run(init_step);
    
    int npart = MMC.getNpart();
    int nsteps = init_step;

    
    FILE *fptr;
    fptr = fopen(conforname, "w");
    for(int i=0;i<num_configs;i++)
    {
        MMC.run(num_steps);
        nsteps+=num_steps;
                
        if(i%1000==0)
        printf("configuration %d in %d....done.\n",i,num_configs);

        //fprintf(fptr, "%d\n\n",npart);
        for(int i=0;i<npart;i++){
            fprintf(fptr, "%.5f ", MMC.getPart(i));}
        fprintf(fptr, "\n");
    }
    fclose(fptr);
    
    return 0; 
} 



















