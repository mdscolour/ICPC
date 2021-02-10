#include "MyLib.h"
#define PARA_FUNC getGr
#define ROLLING_SIZE 5000

typedef double (*type_pot)(double);

double LJlike(double r, double sigma, double epsilon, double a, double b)
{
    return 4 * epsilon * (pow(sigma / r, a) - pow(sigma / r, b));
}
////// class declaration
/// class of coreMolecularMC, which will map from potential function to radial distribution function by simulating the configurations
class coreMolecularMC
{
private:
    type_pot funcNam;
    vector<double> disFor;
    vector<double> disBck;
    double minvol;

public:
    double alpha;
    vector<double> part;
    vector<int> ptype;
    int npart;
    double scaling;
    double box_l;
    double box_r;
    double box_s;

    coreMolecularMC(type_pot tfuncNam, double talpha, double tscaling = 0.45,
                    double tminvol = PRECISION)
    {
        funcNam = tfuncNam;
        alpha = talpha;
        minvol = tminvol;
        scaling = tscaling;
    }
    double getSpacing()
    {
        double s = part[npart - 1] - part[0];
        if (s < 0)
            return s + box_s;
        return s;
    }
    void setBox(double tboxl, double tboxr)
    {
        box_l = tboxl;
        box_r = tboxr;
        box_s = tboxr - tboxl;
    }
    bool passExcludedVolume(int ipart, double step);
    void addPart(double pos) { part.push_back(pos); }
    int getNpart() { return npart; }
    double getPart(int i) { return part[i]; }
    vector<double> getPartAll() { return part; }
    void printfDisAll();
    void printfPart()
    {
        for (int i = 0; i < part.size(); i++)
        {
            std::cout << part[i] << " ";
        }
        std::cout << '\n';
    }
    vector<double> getGr(int ngr = 160, double dx = 1);
    double rollingCount(int pos, int rolsiz);
    vector<double> GetAutocorrelation(int nu, double t_intersteps = 1);
    bool thermalOrNotByLastThree(vector<double> &v, double spacing, int sizV = 10);
    int runTilThermalByLastThree(int each_step);
    int runTilThermalBySpacing(int each_step, bool upgoing, int sizV = 10);
    double getPotMin(double maxX = 30, double dx = 0.01);
    void addAllPartByMinPot(int num, double potRange = 50, double potPrecision = 0.01);
    void addAllPartEvenly(int num);

    bool illedPot(double maxX=5,double sep=0.1,double tol=0.1)
{
    vector<double> x = arange(minvol, maxX+minvol,sep);
    int numRange = int(x.size());
    double potminy = 1e99;
    int potminx = -1;
    for (int i = 0; i < numRange; i++)
    {
        double tmpy = funcNam(x[i]);
        if (tmpy <= potminy)
        {
            potminy = tmpy;
            potminx = i;
        }
    }
    if(x[potminx]-minvol<=tol)return true;
    else return false;
}


    bool readyEverything()
    {
        updateAllDis();
        return true;
    }
    void updateAllDis();
    double proposal(int ipart);
    inline bool acceptOrNot(double oldp, double newp);
    bool runOne();
    int run(int num);
    double getOriPot(int ipart);
    double getNewPot(int ipart, double prop);
    void updateOne(int i, double prop);
};

class coreDoubleRMC
{
public:
    static vector<double> targetX;
    static vector<double> targetY;
    static vector<double> potX;
    static vector<double> curinpY;
    static vector<double> propinpY;
    static coreMolecularMC *MMC;
    static double dx;
    static bool potFlag;
    static int numX;
    static double oldp;

    static double alpha;
    static int num_configs;
    static int num_steps;
    static int init_step;
    static int num_particles;
    static double move_scale;
    static double move_cut;

    coreDoubleRMC(vector<double> tpotX)
    {
        loadTarget("targetGr.gr");
        potX = tpotX;
    }

    static void newMMC();

    static void loadTarget(const char *name);
    static void initGusTarget();
    static void initGusZero();
    static double MMCPotential(double r);
    static void updatePotX(vector<double> arrX);
    static inline int proposal();
    static inline int proposalGlobal();
    static inline void proposalGaussian();
    static double getDiff(vector<double> gr);
    static inline bool acceptOrNot(double oldp, double newp);

    static bool runOne();
    static bool runOneGlobal();
    static vector<double> runMMC();

    //below are all inline function
    static void initPotentialX(vector<double> tpotX) //used if no instance
    {
        loadTarget("targetGr.gr");
        potX = tpotX;
    }
    static void readyToRun()
    {
        potFlag = true;
        initGusTarget();
        //initGusZero();
        newMMC();
        oldp = getDiff(runMMC());
    } //initilize old energy
    static void asignConfig(double talpha, int tnum_configs, int tnum_steps,
                            int tinit_step, int tnum_particles,
                            double tmove_scale, double tmove_cut)
    {
        alpha = talpha;
        num_configs = tnum_configs;
        num_steps = tnum_steps;
        init_step = tinit_step;
        num_particles = tnum_particles;
        move_scale = tmove_scale;
        move_cut = tmove_cut;
    }
    static void deleteMMC() { delete MMC; };
    static double getCurrentEnergy()
    {
        potFlag = true;
        oldp = getDiff(runMMC());
        return oldp;
    }
    static void updateOne(int i, double tnewp)
    {
        curinpY[i] = propinpY[i];
        oldp = tnewp;
    }
    static bool run(int N)
    {
        //printf("%d,%d,%d\n",potX.size(),curinpY.size(),propinpY.size());
        int resCount = 0;
        bool resFlag;
        for (int i = 0; i < N; i++)
        {
            //resFlag=runOne();
            //cout<<i<<"    ";
            resFlag = runOneGlobal();
            if (resFlag)
                resCount++;
        }
        printf("%d  #N:%d, step size: %f, success rate: %f\n", init_step, N, move_cut, double(resCount) / double(N));
        return resFlag;
    }
    static void writePot(const char *name = "savedpotential.dat") { write2DVector(name, potX, curinpY); }
};

//all inline and static except run()
class ClassErSearch
{
public:
    static vector<double> ersearchy;
    static vector<double> ersearchx;
    static vector<double> ersearchymid;
    static vector<double> ersearchymodel;
    static vector<double> vecval;
    static vector<double> best;
    static vector<double> bestpara;
    static double bestp;
    static int vnum;

    static double global_scaling;
    static double global_boxl;
    static double global_boxr;
    static int global_npart;
    static double global_minvol;
    static double global_temperature;
    static int global_num_configs;
    static int global_num_steps;
    static int global_init_step;
    //static coreMolecularMC* MMC;

    ClassErSearch() {}
    static void readConfig(const char *name = "config.config")
    {
        FILE *fptr = fopen(name, "r");
        fscanf(fptr, "global_boxl:%lf\n", &global_boxl);
        fscanf(fptr, "global_boxr:%lf\n", &global_boxr);
        fscanf(fptr, "global_npart:%d\n", &global_npart);
        fscanf(fptr, "global_minvol:%lf\n", &global_minvol);
        fscanf(fptr, "global_temperature:%lf\n", &global_temperature);
        fscanf(fptr, "global_scaling:%lf\n", &global_scaling);
        fscanf(fptr, "global_num_configs:%d\n", &global_num_configs);
        fscanf(fptr, "global_num_steps:%d\n", &global_num_steps);
        fscanf(fptr, "global_init_step:%d\n", &global_init_step);
        fclose(fptr);
    }
    static void run(char *argv[], vector<double> tersearchx)
    {
        ersearchx = tersearchx;
        int numx = ersearchx.size();

        double outer_define = atof(argv[1]);
        ersearchymodel = {outer_define, -99, -99, -99}; //-99 refer to ersearchymid
        ersearchymid = {-99, -99, -99};
        vecval = arange(1, 20, 1);

        //coreDoubleRMC::loadTarget("targetGr.gr");
        best = {};
        bestp = 1e99;
        vnum = (int)vecval.size();

        getCom(int(ersearchymid.size()) - 1);

        best.insert(best.begin(), outer_define);
        printfVector(best);
        char buffer[50]; //only storage for 256 characters.
        sprintf(buffer, "%lf.erpotout", bestp);
        write2DVector(buffer, arange(0, best.size()), best);
    }
    static void updateErSearchy()
    {
        ersearchy.clear();
        //         int midcount = 0;
        //         for(double n : ersearchymodel)
        //         {
        //             if(n!=-99) ersearchy.push_back(n);
        //             else
        //             {
        //                 ersearchy.push_back(ersearchymid[midcount]);
        //                 midcount++;
        //             }
        //         }
        int numx = ersearchx.size();
        for (int j = 0; j < numx; j++)
            ersearchy.push_back(LJlike(ersearchx[j],
                                       ersearchymodel[0], ersearchymid[0], ersearchymid[1], ersearchymid[2]));
        if (ersearchy.size() != numx)
            printf("input of ersearchy error\n");
        //ersearchy.push_back(100);
        //for(double n : ersearchymid)ersearchy.push_back(n);
        //ersearchy.push_back(0);
    }

    static double erPotential(double r)
    {
        return interpolate(ersearchx, ersearchy, r, true);
    }
    static double runErSearch()
    {
        updateErSearchy();

        coreMolecularMC MMC = coreMolecularMC(&erPotential,
                                              global_temperature, global_scaling, global_minvol);

        MMC.setBox(global_boxl, global_boxr);
        MMC.addAllPartEvenly(global_npart);
        MMC.readyEverything();

        int num_configs = global_num_configs;
        int num_steps = global_num_steps;
        int init_step = global_init_step;
        MMC.run(init_step);

        int ngr = int(coreDoubleRMC::targetX.size());
        double dx = coreDoubleRMC::dx;
        int lgr = int(double(ngr) * dx);
        vector<double> gr = MMC.PARA_FUNC(ngr, dx);

        for (int i = 0; i < num_configs; i++)
        {
            MMC.run(num_steps);
            VectorAddedByVector(gr, MMC.PARA_FUNC(ngr, dx));
        }
        VectorDividedByScalar(gr, num_configs);

        double newp = coreDoubleRMC::getDiff(gr);
        return newp;
    }
    static double runSpecificPot(bool writeGr = false)
    {
        //set ersearchx and ersearchy in advance

        coreMolecularMC MMC = coreMolecularMC(&erPotential, 1, 0.45);

        MMC.setBox(0, 800);
        int npart = 200;
        MMC.addAllPartEvenly(npart);
        MMC.readyEverything();

        int num_configs = 500;
        int num_steps = 1;
        int init_step = 500;
        MMC.run(init_step);

        //int lgr = 35;
        //double dx=0.01;
        //int ngr = int(double(lgr)/dx);
        int ngr = int(coreDoubleRMC::targetX.size());
        double dx = coreDoubleRMC::dx;
        int lgr = int(double(ngr) * dx);
        vector<double> gr = MMC.PARA_FUNC(ngr, dx);

        FILE *fptr;
        if (writeGr)
            fptr = fopen("test.out", "w");
        for (int i = 0; i < num_configs; i++)
        {
            MMC.run(num_steps);
            VectorAddedByVector(gr, MMC.PARA_FUNC(ngr, dx));

            if (writeGr)
            {
                fprintf(fptr, "ITEM: TIMESTEP\n%d\n", init_step + i * num_steps);
                fprintf(fptr, "ITEM: NUMBER OF ATOMS\n%d\n", npart);
                fprintf(fptr, "ITEM: BOX BOUNDS pp pp pp\n");
                fprintf(fptr, "0 1\n");
                fprintf(fptr, "0 1\n");
                fprintf(fptr, "0 1\n");
                fprintf(fptr, "ITEM: ATOMS id xs ys zs\n");
                for (int i = 0; i < npart; i++)
                {
                    fprintf(fptr, "%d %.5f %.5f %.5f\n", i, MMC.getPart(i), 0.0, 0.0);
                }
            }
        }
        VectorDividedByScalar(gr, num_configs);
        double newp = coreDoubleRMC::getDiff(gr);
        //cout<<newp<<endl;

        if (writeGr)
            fclose(fptr);
        if (writeGr)
            write2DVector("MC result.gr", arange(0, lgr, dx), gr);
        return newp;
    }
    static void runPotInput(double t, double scal, double minvol, double tbox_l, double tbox_r, int npart, int num_configs, int num_steps, int init_step, int ngr, double dx)
    {
        //set ersearchx and ersearchy in advance

        coreMolecularMC MMC = coreMolecularMC(&erPotential, t, scal, minvol);

        MMC.setBox(tbox_l, tbox_r);
        MMC.addAllPartEvenly(npart);
        MMC.readyEverything();

        MMC.run(init_step);

        //int lgr = int(double(ngr)*dx);
        vector<double> gr = MMC.getGr(ngr, dx);

        for (int i = 0; i < num_configs; i++)
        {
            MMC.run(num_steps);
            VectorAddedByVector(gr, MMC.getGr(ngr, dx));
        }
        VectorDividedByScalar(gr, num_configs);
        best.clear();
        best = gr;
        //double newp = coreDoubleRMC::getDiff(gr);
        //cout<<newp<<endl;
    }
    static void tolRunErSearch()
    {
        double newp = runErSearch();
        if (newp <= bestp)
        {
            best = ersearchymid;
            bestp = newp;
            printfVector(ersearchymid);
        }
    }

    //void getCom(int lv);
    static void getCom(int lv)
    {
        if (lv == 0)
            for (int i = 0; i < vnum; i++)
            {
                ersearchymid[lv] = vecval[i];
                tolRunErSearch();
            }
        else
        {
            for (int i = 0; i < vnum; i++)
            {
                ersearchymid[lv] = vecval[i];
                getCom(lv - 1);
            }
        }
    }

    static void randomSearch()
    {
        //set ersearchx in advance
        int numx = int(ersearchx.size());
        vector<double> para(4, 0.0);
        //for(int i=0;i<100000;i++)
        for (int i = 0; i < 10000; i++)
        {
            ersearchy.clear();
            //ersearchy.push_back(1e20);
            //for(int j=0;j<numx-2;j++)ersearchy.push_back(-10*RNG_NAME()+5);
            for (int k = 0; k < 4; k++)
                para[k] = 19 * RNG_NAME() + 1;
            for (int j = 0; j < numx; j++)
                ersearchy.push_back(LJlike(ersearchx[j],
                                           para[0], para[1], para[2], para[3]));
            //ersearchy.push_back(0);
            if (numx != ersearchy.size())
                cout << "error 111" << endl;

            double newp = runSpecificPot();
            if (newp <= bestp)
            {
                bestpara = para;
                best = ersearchy;
                bestp = newp;
                cout << bestp << "   ";
                printfVector(bestpara);
            }
        }
        //ersearchy.clear();ersearchy=best;
        //runSpecificPot(true);
        //write2DVector("erPotential.rndpot",ersearchx,best);
    }
};

class coreSampleMCMC
{
public:
    static vector<double> targetX;
    static vector<double> targetY;
    static vector<double> oldpara;
    static vector<double> curpara;
    static coreMolecularMC *MMC;
    static double dx;
    static double oldp;
    static double rndSize;
    static double alpha;

    static double global_LJpara;
    static double global_scaling;
    static double global_boxl;
    static double global_boxr;
    static int global_npart;
    static double global_minvol;
    static double global_temperature;
    static int global_num_configs;
    static int global_num_steps;
    static int global_init_step;

    static double modelPot(double r)
    {
        //return LJlike(r, 162,curpara[1],9,7);
        return LJlike(r, curpara[0], curpara[1], curpara[2], curpara[3]);
    }

    coreSampleMCMC()
    {
    }

    static void newMMC()
    {
        MMC = new coreMolecularMC(&modelPot, global_temperature, 
            global_scaling, global_minvol);
        MMC->setBox(global_boxl,global_boxr);
        MMC->addAllPartEvenly(global_npart);
        MMC->readyEverything();
    }

    static void loadTarget(const char *name)
    {
        FILE *fptr = fopen(name, "r");
        double x;
        double y;

        while (fscanf(fptr, "%lf %lf\n", &x, &y) == 2)
        {
            targetX.push_back(x);
            targetY.push_back(y);
        }
        dx = targetX[1] - targetX[0];
        //cout<<"load targetGr.gr ~success~"<<endl;
    }

    static void proposal()
    {
        //for(int i=0;i<curpara.size();i++)
        //curpara[0] =fPeriodRestrict(curpara[0] + int(RNG_NAME() * (2*rndSize+1)-rndSize-0.5),0,20);
        curpara[1] =fPeriodRestrict(curpara[1] + int(RNG_NAME() * (2*rndSize+1)-rndSize-0.5),0,20);
        //curpara[2] =fPeriodRestrict(curpara[2] + int(RNG_NAME() * (2*rndSize+1)-rndSize-0.5),0,20);
        //curpara[3] =fPeriodRestrict(curpara[3] + int(RNG_NAME() * (2*rndSize+1)-rndSize-0.5),0,20);
     }
    static double getDiff(vector<double> gr)
    {
        int numGr = int(gr.size());
        if (numGr != targetY.size())
            printf("Gr size not equal:%d %d\n", numGr, int(targetY.size()));

        double sum = 0;
        for (int i = 0; i < numGr; i++)
            sum += pow(gr[i] - targetY[i], 2);
        //for(int i=0;i<numGr;i++)sum+=(gr[i]-targetY[i])*(gr[i]-targetY[i]);
        return sum / (double(numGr));
        //return sum;
    }
    static inline bool acceptOrNot(double oldp, double newp)
    {
        if (newp <= oldp)
            return true;
        if (RNG_NAME() < exp(alpha * (oldp - newp)))
            return true;
        return false;
    }

    static bool runOne()
    {
        proposal();

        double spacing = MMC->getPotMin(30, 0.01);
        if (spacing <= PRECISION)
        { /*printf("trival step\n");*/
            return false;
        }

        double newp = getDiff(runMMC());
        //cout<<" old_SSD:"<<oldp<<" new_SSD:"<<newp;
        bool resFlag = acceptOrNot(oldp, newp);
        if(resFlag){oldpara.clear();oldpara.assign(curpara.begin(), curpara.end());oldp = newp;}
        else {curpara.clear();curpara.assign(oldpara.begin(), oldpara.end());} 
        //cout<<"  "<<resFlag<<endl;
        return resFlag;
    }
    static vector<double> runMMC()
    {
        deleteMMC();
        newMMC();

        MMC->run(global_init_step);
        int ngrpoint = int(targetX.size());
        vector<double> gr = MMC->getGr(ngrpoint, dx);
        vector<double> tgr;
        double countgr = 1;
        for (int i = 0; i < global_num_configs; i++)
        {
            MMC->run(global_num_steps);
            tgr = MMC->getGr(ngrpoint, dx);
            VectorAddedByVector(gr, tgr);
            countgr++;
        }
        VectorDividedByScalar(gr, countgr);
        return gr;
    }
    static void initGus()
    {
        curpara.clear();
        curpara.push_back(1);
        curpara.push_back(1);
        curpara.push_back(1);
        curpara.push_back(1);
        oldpara = curpara;
    }
    static void readyToRun()
    {
        initGus();
        newMMC();
        oldp = getDiff(runMMC());//initilize old energy
    } 
    static void asignRndSize(double trndSize,double talpha=1)
    {
        rndSize = trndSize;
        alpha=talpha;
    }
    static void deleteMMC() { delete MMC; };
    static double getCurrentEnergy()
    {
        oldp = getDiff(runMMC());
        return oldp;
    }
    static bool run(int N,bool printFlag=false)
    {
        //printf("%d,%d,%d\n",potX.size(),curinpY.size(),propinpY.size());
        int resCount = 0;
        bool resFlag;
        for (int i = 0; i < N; i++)
        {
            resFlag=runOne();
            if (resFlag)
                resCount++;
        }
        if(printFlag)printf("#init:%d  #steps:%d, step size: %f, success rate: %f\n", global_init_step, N, rndSize, double(resCount) / double(N));
        return resFlag;
    }
    static void readConfig(const char *name = "config.config")
    {
        FILE *fptr = fopen(name, "r");
        fscanf(fptr, "global_boxl:%lf\n", &global_boxl);
        fscanf(fptr, "global_boxr:%lf\n", &global_boxr);
        fscanf(fptr, "global_npart:%d\n", &global_npart);
        fscanf(fptr, "global_minvol:%lf\n", &global_minvol);
        fscanf(fptr, "global_temperature:%lf\n", &global_temperature);
        fscanf(fptr, "global_scaling:%lf\n", &global_scaling);
        fscanf(fptr, "global_num_configs:%d\n", &global_num_configs);
        fscanf(fptr, "global_num_steps:%d\n", &global_num_steps);
        fscanf(fptr, "global_init_step:%d\n", &global_init_step);
        fclose(fptr);
    }

    static void writePot(const char *name = "potpara.MCpot") { write2DVector(name, arange(0,int(curpara.size()),1), curpara); }
};

class coreDirectRun
{
public:
    static vector<double> targetX;
    static vector<double> targetY;
    static vector<double> curpara;
    static coreMolecularMC *MMC;
    static double dx;

    static double global_scaling;
    static double global_boxl;
    static double global_boxr;
    static int global_npart;
    static double global_minvol;
    static double global_temperature;
    static int global_num_configs;
    static int global_num_steps;
    static int global_init_step;

    static double modelPot(double r)
    {
        //return LJlike(r, 162,curpara[1],9,7);
        return LJlike(r, curpara[0], curpara[1], curpara[2], curpara[3]);
    }
    static double assignAndRun(double p1,double p2,double p3,double p4)
    {
        assignPara(p1,p2,p3,p4);
        return runOne();
    }

    static void assignPara(double p1,double p2,double p3,double p4)
    {
        //if(int(curpara.size())!=4)
        //{
            curpara.clear();
            curpara.push_back(p1);
            curpara.push_back(p2);
            curpara.push_back(p3);
            curpara.push_back(p4);
        //}
        //else
        //{
       // curpara[0]=p1;
       // curpara[1]=p2;
       // curpara[2]=p3;
       // curpara[3]=p4;
       // }
    }
    coreDirectRun()
    {
    }

    static void newMMC()
    {
        MMC = new coreMolecularMC(&modelPot, global_temperature, 
            global_scaling, global_minvol);
        MMC->setBox(global_boxl,global_boxr);
        MMC->addAllPartEvenly(global_npart);
        MMC->readyEverything();
    }

    static void loadTarget(const char *name)
    {
        FILE *fptr = fopen(name, "r");
        double x;
        double y;
        targetX.clear();
        targetY.clear();

        while (fscanf(fptr, "%lf %lf\n", &x, &y) == 2)
        {
            targetX.push_back(x);
            targetY.push_back(y);
        }
        dx = targetX[1] - targetX[0];
        //cout<<"load targetGr.gr ~success~"<<endl;
    }
    
    static double getDiff(vector<double> gr)
    {
        int numGr = int(gr.size());
        if (numGr != targetY.size())
            printf("Gr size not equal:%d %d\n", numGr, int(targetY.size()));

        double sum = 0;
        for (int i = 0; i < numGr; i++)
            sum += pow(gr[i] - targetY[i], 2);
        //for(int i=0;i<numGr;i++)sum+=(gr[i]-targetY[i])*(gr[i]-targetY[i]);
        return sum / (double(numGr));
        //return sum;
    }

    static double runOne()
    {
        double illedFlag = MMC->illedPot();
        if (illedFlag)
        {   
            printf("illed potential\n");
            return 99999;
        }

        return getDiff(runMMC());
    }
    static vector<double> runMMC()
    {
        deleteMMC();
        newMMC();

        MMC->run(global_init_step);
        int ngrpoint = int(targetX.size());
        vector<double> gr = MMC->getGr(ngrpoint, dx);
        vector<double> tgr;
        double countgr = 1;
        for (int i = 0; i < global_num_configs; i++)
        {
            MMC->run(global_num_steps);
            tgr = MMC->getGr(ngrpoint, dx);
            VectorAddedByVector(gr, tgr);
            countgr++;
        }
        VectorDividedByScalar(gr, countgr);
        return gr;
    }
   
    static void readyToRun()
    {
        newMMC();
    } 
    
    static void deleteMMC() { delete MMC; };
    static void readConfig(const char *name = "config.config")
    {
        FILE *fptr = fopen(name, "r");
        fscanf(fptr, "global_boxl:%lf\n", &global_boxl);
        fscanf(fptr, "global_boxr:%lf\n", &global_boxr);
        fscanf(fptr, "global_npart:%d\n", &global_npart);
        fscanf(fptr, "global_minvol:%lf\n", &global_minvol);
        fscanf(fptr, "global_temperature:%lf\n", &global_temperature);
        fscanf(fptr, "global_scaling:%lf\n", &global_scaling);
        fscanf(fptr, "global_num_configs:%d\n", &global_num_configs);
        fscanf(fptr, "global_num_steps:%d\n", &global_num_steps);
        fscanf(fptr, "global_init_step:%d\n", &global_init_step);
        fclose(fptr);
    }
};

class classAnneal
{
    public:
        vector<double> oldpara;
        vector<double> newpara;
        vector<double> vRndSize;
        double energy;
        void readyToRun(const char* tag,const char* cog,vector<double> toldpara,vector<double> tvRndSize)
        {
            coreDirectRun::loadTarget(tag);
            coreDirectRun::readConfig(cog);
            coreDirectRun::readyToRun();
            oldpara.clear();oldpara.assign(toldpara.begin(), toldpara.end());
            newpara = oldpara;
            energy = 99999;
            vRndSize.clear();vRndSize.assign(tvRndSize.begin(), tvRndSize.end());
        }

        vector<double> anneal(double T,double T_min,double lambda,int SAcount_max,const char* name="potparaAll.SApot")
        {
            energy = coreDirectRun::assignAndRun(oldpara[0],oldpara[1],oldpara[2],oldpara[3]);

            //double T = 100;
            //double T_min = 0.0006;
            //double lambda = 0.992;
            
            int SAcount = 0;
            int i;
            double new_E;
            bool apFlag;
    FILE *fptr;
    fptr = fopen(name, "w");
            while(T > T_min && SAcount<SAcount_max)
            {
                i = 1;
                while(i <= 1)
                {
                    newpara = updateNew(oldpara);
                    new_E = coreDirectRun::assignAndRun(newpara[0],newpara[1],newpara[2],newpara[3]);
                    apFlag = acceptOrNot(energy, new_E, T);
                    if(apFlag)
                    {
                        oldpara = newpara;
                        energy = new_E;
                    }
                    i++;
                    SAcount++;
                }
                T *= lambda;

    fprintf(fptr, "%e %.5f %.5f %.5f %.5f %.5f\n",T,oldpara[0],oldpara[1],oldpara[2],oldpara[3],energy);
    //printf("%e %.5f %.5f %.5f %.5f,%.5f\n",T,oldpara[0],oldpara[1],oldpara[2],oldpara[3],energy);
            }
            fclose(fptr);
            return oldpara;
        }
        bool acceptOrNot(double oldp, double newp, double T)
        {
            if (newp <= oldp)return true;
            if (RNG_NAME() < exp((oldp - newp)/T))return true;
            return false;
        }
        vector<double> updateNew(vector<double> v1)
        {
            //// parameter domains set manually
            vector<double> v2={};
            if(vRndSize[0]==0)v2.push_back(v1[0]);
            else v2.push_back( fPeriodRestrict(v1[0] + round(RNG_NAME() * (2*vRndSize[0])-vRndSize[0]),100,350));
            
            if(vRndSize[1]==0)v2.push_back(v1[1]);
            else v2.push_back( fPeriodRestrict(v1[1] + round(RNG_NAME() * (2*vRndSize[1])-vRndSize[1]),1,30));
            
            if(vRndSize[2]==0)v2.push_back(v1[2]);
            else v2.push_back( fPeriodRestrict(v1[2] + round(RNG_NAME() * (2*vRndSize[2])-vRndSize[2]),1,20));
            
            if(vRndSize[3]==0)v2.push_back(v1[3]);
            else v2.push_back( fPeriodRestrict(v1[3] + round(RNG_NAME() * (2*vRndSize[3])-vRndSize[3]),1,20));
            // v2.push_back(12);
            // v2.push_back(6);
            return v2;
        }

};


//void readyToRun(const char* tag,const char* cog,vector<double> toldpara,vector<double> tvRndSize)
//vector<double> anneal(double T,double T_min,double lambda,int SAcount_max,const char* name="potparaAll.SApot")
void SAMethod1(const char* tagnam="targetGr.gr",const char* lowname="config.lowconfig",const char* highname="config.highconfig")
{
    classAnneal* ann;
    ann = new classAnneal();
    char buffer[50];
    vector<double> told,tsiz;
    vector<vector<double>> recPara = {};
    vector<double> tmpPara = {};

    for(double a=2;a<19;a++)
    for(double b=1;b<a;b++)
    {
        sprintf(buffer, "potpara%d_%d.SAM1exppot",int(a),int(b));
        delete ann;
        ann = new classAnneal();
        told = {100.,1.0,a,b};
        tsiz = {30.,3.,0.,0.};
        
        ann->readyToRun(tagnam,lowname,told,tsiz);
        tmpPara = ann->anneal(100,0.001,0.94,150,buffer);
        
        recPara.push_back(tmpPara);
        printf("a=%d b=%d finished.\n",a,b);
    }

coreDirectRun::loadTarget(tagnam);
coreDirectRun::readConfig(highname);
coreDirectRun::readyToRun();

int nrec = int(recPara.size());
double oldenergy = 99999;
double energy;
vector<double> best={};
double p1=0,p2=0,p3=0,p4=0;
FILE *fptr;
//char buffer[50];
sprintf(buffer, "potparaRes.SAM1exppot"); 

fptr = fopen(buffer, "w");
for(double ii=0;ii<nrec;ii+=1)
{
    p1=recPara[ii][0];
    p2=recPara[ii][1];
    p3=recPara[ii][2];
    p4=recPara[ii][3];
    energy = coreDirectRun::assignAndRun(p1,p2,p3,p4);
    fprintf(fptr, "%.5f %.5f %.5f %.5f %.5f\n",p1,p2,p3,p4,energy);
    printf("%.5f %.5f %.5f %.5f,%.5f\n",p1,p2,p3,p4,energy);
    if(energy<oldenergy){best=recPara[ii];oldenergy=energy;}
}
cout<<"best fit:";
printfVector(best);
}

void SAMethod2()
{
    classAnneal* ann;
    ann = new classAnneal();
    char buffer[50];
    vector<double> told,tsiz;
    vector<vector<double>> recPara = {};
    vector<double> tmpPara = {};

    for(double nrep=0;nrep<20;nrep++)
    {
        sprintf(buffer, "potpara%d.SAM2pot",int(nrep));
        delete ann;
        ann = new classAnneal();
        tsiz = {30.,3.,3.,3.};
        told.clear();
        told.push_back(int((250-150)*RNG_NAME()+150));
        told.push_back(int((20-1)*RNG_NAME()+1));
        told.push_back(int((20-1)*RNG_NAME()+1));
        told.push_back(int((20-1)*RNG_NAME()+1));
        
        ann->readyToRun("standardLJTest.gr","low.config",told,tsiz);
        tmpPara = ann->anneal(100,0.001,0.993,1500,buffer);
        
        recPara.push_back(tmpPara);
        printf("#%d finished.\n",nrep);
    }

coreDirectRun::loadTarget("standardLJTest.gr");
coreDirectRun::readConfig("high.config");
coreDirectRun::readyToRun();

int nrec = int(recPara.size());
double oldenergy = 99999;
double energy;
vector<double> best={};
double p1=0,p2=0,p3=0,p4=0;
FILE *fptr;
//char buffer[50];
sprintf(buffer, "potparaRes.SAM2pot"); 

fptr = fopen(buffer, "w");
for(double ii=0;ii<nrec;ii+=1)
{
    p1=recPara[ii][0];
    p2=recPara[ii][1];
    p3=recPara[ii][2];
    p4=recPara[ii][3];
    energy = coreDirectRun::assignAndRun(p1,p2,p3,p4);
    fprintf(fptr, "%.5f %.5f %.5f %.5f %.5f\n",p1,p2,p3,p4,energy);
    printf("%.5f %.5f %.5f %.5f,%.5f\n",p1,p2,p3,p4,energy);
    if(energy<oldenergy){best=recPara[ii];oldenergy=energy;}
}
cout<<"best fit:";
printfVector(best);
}

/*coreDoubleRMC::loadTarget("targetGrshort.gr");
vector<double> best;
double bestp=1e99;
    
for(double i1=-3;i1<=0.5;i1+=0.5)
    for(double i2=-3;i2<=0.5;i2+=0.5)
        for(double i3=-3;i3<=0.5;i3+=0.5)
            for(double i4=-3;i4<=0.5;i4+=0.5)
                for(double i5=-3;i5<=0.5;i5+=0.5)
                    for(double i6=-3;i6<=0.5;i6+=0.5)
                {
                    ersearchy.clear();
                    ersearchy.push_back(100000);
                    ersearchy.push_back(i1);
                    ersearchy.push_back(i2);
                    ersearchy.push_back(i3);
                    ersearchy.push_back(i4);
                    ersearchy.push_back(i5);
                    ersearchy.push_back(i6);
                    ersearchy.push_back(0);
                    double newp=runErSearch();
                    if(newp<=bestp){best=ersearchy;bestp=newp;cout<<bestp<<endl;}
                }
write2DVector("erPotentialshort20.pot",ersearchx,best);*/

/*vector<double> ersearchy;
vector<double> ersearchx = {0,2,3,4,5,6,7,100};
vector<double> ersearchymid = {2,3,4,5,6,7};
vector<double> vecval = arange(-3,0.55,1);
double erPotential(double r)
{
    return interpolate(ersearchx,ersearchy,r,true);
}
double runErSearch()
{
    ersearchy.clear();ersearchy.push_back(10000);
    for(double n : ersearchymid)ersearchy.push_back(n);
    ersearchy.push_back(0);
    
    coreMolecularMC MMC=coreMolecularMC(&erPotential,1,0.45);

    MMC.setBox(0,800);
    MMC.addAllPartEvenly(20);
    MMC.readyEverything();
    
    int num_configs=100;
    int num_steps=10;
    int init_step = 1000;
    MMC.run(init_step);
    
    int lgr = 35;
    int ngr = 70;
    double dx=0.5;
    vector<double> gr=MMC.getGr(ngr,dx);

    for(int i=0;i<num_configs;i++)
    {
        MMC.run(num_steps);
        VectorAddedByVector(gr,MMC.getGr(ngr,dx));
    }
    VectorDividedByScalar(gr,num_configs);
    
    double newp = coreDoubleRMC::getDiff(gr);
    return newp;
}

vector<double> best;
double bestp=1e99;
int vnum = (int)vecval.size();
void tolRun()
{
    double newp=runErSearch();
    if(newp<=bestp){best=ersearchy;bestp=newp;printfVector(ersearchy);}
}

//void getCom(int lv);
void getCom(int lv)
{
    if(lv==0)for(int i=0;i<vnum;i++){ersearchymid[lv]=vecval[i];tolRun();}
    else
    {
        for(int i=0;i<vnum;i++){ersearchymid[lv]=vecval[i];getCom(lv-1);}
    }
}*/
