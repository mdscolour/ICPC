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
char* itarea = new char[1];

//itarea[0]='0';
itarea=(argv[1]);

sprintf(grname, "chr2-%s.midgr", itarea);
sprintf(lowname, "chr2-%s.lowconfig", itarea);
//printf(highname, "chr2-%s.highconfig", itarea);
sprintf(conforname, "res/chr2-%s.pml", itarea);
sprintf(finresname, "res/chr2-%s.finres", itarea);
sprintf(savegrname, "res/chr2-%s.gr", itarea);

// double t1,t2;
// FILE *fptrpara = fopen(finresname, "r");
// while (fscanf(fptrpara, "%lf %lf %lf %lf %lf %lf\n", &t1, &global_sigma,
//     &global_epsilon, &global_a, &global_b, &t2) == 6){}
// fclose(fptrpara);
global_sigma=164;
global_epsilon=1;
global_a=14;
global_b=13;
//printf("%lf %lf %lf %lf %lf %lf\n", t1, global_sigma,
//    global_epsilon, global_a, global_b, t2);

coreDirectRun::loadTarget(grname);
coreDirectRun::readConfig(lowname);
coreDirectRun::readyToRun();
// 
// double energy = coreDirectRun::assignAndRun(162,3,12,6);
// 
// vector<double> gr = coreDirectRun::savegr;
// 
int ngr = int(coreDirectRun::targetX.size());
double dx = coreDirectRun::dx;
int lgr = int(double(ngr) * dx);
//         
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
    
    int num_configs=5000;//coreDirectRun::global_num_configs;
    int num_steps=10;//coreDirectRun::global_num_steps;
    int init_step=5000;//coreDirectRun::global_init_step;
    MMC.run(init_step);
    
    int npart = MMC.getNpart();
    int nsteps = init_step;
    
//    FILE *fptrinit;
//    fptrinit = fopen("inittest.dat", "w");
//    spacing=0.0;double spaold=0;
/*    for(int i=0;i<10000;i++)
    {
        MMC.run(each_step);
        spacing = MMC.part[num-1]-MMC.part[0];
        nsteps+=each_step;
        //printf("initializaiton %d steps....done, %.6f spacing\n",nsteps,spacing/num);
        fprintf(fptrinit, "%d %.8f\n",i,spacing/num);
    }*/
    //nsteps += MMC.runTilThermalBySpacing(each_step,false);//true for small seperation in initial position
    
    //std::cout << "particle\n";
    //MMC.printfPart();
    //std::cout << "distance\n";
    //MMC.printfDisAll();
//     double spacing=0.0,sum_spacing=0.0,sum_spacing_old=0;
//     for(int j=0;j<1000;j++)
//     while(true)
//     {
//         for(int i=0;i<10;i++)
//         {
//             MMC.run(init_step);
//             spacing = MMC.part[num-1]-MMC.part[0];
//             sum_spacing+=spacing;
//             nsteps+=init_step;
//             printf("initializaiton %d mil. steps....done, %.6f spacing\n",nsteps/1000000,spacing/num);
//         }
//         if(sum_spacing<=sum_spacing_old && sum_spacing_old!=0)break;
//         else {sum_spacing_old=sum_spacing;sum_spacing=0;}
//     }
    
    
//    int lgr = 160;
//    int ngr = 320;
//    double dx=0.5;
    vector<double> gr=MMC.getGr(ngr,dx);
    //printfVector(gr);
    
    //FILE *fptr;int st,ed;double pc;
    //fptr = fopen(conforname, "w");
    for(int i=0;i<num_configs;i++)
    {
        MMC.run(num_steps);
        nsteps+=num_steps;
        VectorAddedByVector(gr,MMC.getGr(ngr,dx));
                
        if(i%1000==0)
        printf("configuration %d in %d....done.\n",i,num_configs);
        
        //fprintf(fptr, "ITEM: TIMESTEP\n%d\n",nsteps);
        //fprintf(fptr, "ITEM: NUMBER OF ATOMS\n%d\n",npart);
        //fprintf(fptr, "ITEM: BOX BOUNDS pp pp pp\n");
        //fprintf(fptr, "0 1\n");
        //fprintf(fptr, "0 1\n");
        //fprintf(fptr, "0 1\n");
        //fprintf(fptr, "ITEM: ATOMS id xs ys zs\n");
        
//         st = 0;ed = 9999999;
//         for(int i=0;i<npart;i++)
//         {
//             pc = MMC.getPart(i);
//             if(pc>=st && pc<=ed)fprintf(fptr, "%.5f  ", pc);
//         }
//         fprintf(fptr, "\n");
        
//         fprintf(fptr, "%d\n\n",npart);
//         for(int i=0;i<npart;i++){
//             fprintf(fptr, "%s %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n","H", MMC.getPart(i), 0.0, 0.0, 200.,0.,100.,0.,0.707,0.,0.707);}
//         fclose(fptr);
//         exit(0);
        
        //cmd.load_cgo( [ 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], "cylinder1" )
        //char spnamebuffer[50];
        //sprintf(spnamebuffer, "%d", 9999);
//         vector<double> v = MMC.part;
//         sort(v.begin(),v.end());
//         fprintf(fptr, "cmd.load_cgo( [ 9.0, 0, 0, 0, %f, 0, 0, 5, 55, 0, 0, 55, 0, 0], \"line%d\" )\n", v[0]-50.,9999);
//         
//         for(int i=0;i<10;i++){
//         fprintf(fptr, "cmd.load_cgo( [ 9.0, %f, 25, 0, %f, -25, 0, 200, 14,73,92, 56,86,96], \"cylinder%d\" )\n", v[i]-50.,v[i]+50.,i);
// 
//         fprintf(fptr, "cmd.load_cgo( [ 9.0, %f, 0, 0, %f, 0, 0, 10, 55, 0, 0, 55, 0, 0], \"line%d\" )\n", v[i]+50.,v[i+1]-50.,i);}
//         
//         fprintf(fptr,"cmd.bg_color(\"white\")\n");
//         fclose(fptr);
//         exit(0);
    }
    //fclose(fptr);
    
    //fclose(fptrinit);
    VectorDividedByScalar(gr,num_configs);
    write2DVector(savegrname,arange(0,lgr,dx),gr);

    return 0; 
} 



















