////////////////////////////////
/////File to generate the schematic figure according to the real density
/////the file address or the number of chromosome need to be changed manually
/////use as:
/////$:g++ thiscpp.cpp
/////$:a.out 0    (0 is section index)
/////$:povray chr2-0.pov -W10240 -H7680
////////////////////////////////

//#include "/home/li/bin/coreMolecularMC.h"
//#include "coreMolecularMC.h"
//#include "/home/li/bin/coreMolecularMC.cpp"//only one instance
#include "./coreMolecularMC.h"//only one instance


int main(int argc,char *argv[])
{    /**************  random seed  *******************************/
SeedByTime();

char conforname[50];
char finresname[50];
char* itarea = new char[1];

//itarea[0]='0';
itarea=(argv[1]);
////////
//change number of chromosome i.e. target file here
sprintf(conforname, "chr6/chr6-%s.pov", itarea);
sprintf(finresname, "chr6/chr6-%s.confor", itarea);
///////
    std::string s; //string
    std::fstream f; //file stream
    f.open(finresname); //open your word list
    std::getline(f, s); 
    std::getline(f, s);//second configuration
    //std::cout << s << "\n"; //output string
    std::istringstream iss(s);
    std::vector<double> v = std::vector<double>{std::istream_iterator<double>(iss),
                         std::istream_iterator<double>()};  
    int npart = int(v.size());
    cout<<"Read successful! Number of nucleosome: "<<v.size()<<endl;

        FILE *fptr;
        fptr = fopen(conforname, "w");
    
        //vector<double> v;// = MMC.part;
        sort(v.begin(),v.end());
        int stp=0;int edp=npart;
        double camerarng = v[edp-1]-v[stp]-100*(edp-stp-1)+8000;
        //double camerarng = v[edp-1]-v[stp]-100*(edp-stp-1);
        double cameramid = v[stp]-4000+camerarng/2.0;
        cout<<cameramid<<endl;
        cout<<camerarng<<endl;
        cout<<camerarng/0.0875/1.414<<endl;
        
        
fprintf(fptr,"//Files with predefined colors and textures\n");
fprintf(fptr,"#include \"colors.inc\"\n");
fprintf(fptr,"#include \"textures.inc\"\n");
fprintf(fptr,"#include \"glass.inc\"\n");
fprintf(fptr,"#include \"shapes.inc\"\n");
fprintf(fptr,"#include \"shapes2.inc\"\n");
fprintf(fptr,"#include \"shapes3.inc\"\n");
fprintf(fptr,"#include \"functions.inc\"\n");
fprintf(fptr,"#include \"math.inc\"\n");
fprintf(fptr,"#include \"transforms.inc\"\n");
fprintf(fptr,"#include \"golds.inc\"\n");
fprintf(fptr,"#include \"metals.inc\"\n");
fprintf(fptr,"#include \"stones.inc\"\n");
fprintf(fptr,"#include \"woods.inc\"\n");
fprintf(fptr,"\n");
fprintf(fptr,"#local nuscal = 0.1;\n");
fprintf(fptr,"//Place the camera\n");
fprintf(fptr,"camera {\n");
fprintf(fptr,"  location 1.05*<nuscal*%.3f,nuscal*%.3f,-nuscal*%.3f> //Camera location\n",
        cameramid,camerarng/0.0875/1.414,camerarng/0.0875/1.414);
fprintf(fptr,"  look_at <nuscal*%.3f,0,-3>     //Where camera is pointing\n",cameramid);
fprintf(fptr,"  angle 5      //Angle of the view--increase to see more, decrease to see less\n");
fprintf(fptr,"}\n");
fprintf(fptr,"\n");
fprintf(fptr,"//Ambient light to \"brighten up\" darker pictures\n");
fprintf(fptr,"global_settings { ambient_light White }\n");
fprintf(fptr,"\n");
fprintf(fptr,"//Place a light--you can have more than one!\n");
fprintf(fptr,"light_source {\n");
fprintf(fptr,"  <nuscal*%.3f,nuscal*%.3f,-nuscal*%.3f> //Camera location\n",cameramid,camerarng/0.0875/1.414,camerarng/0.0875/1.414);
//fprintf(fptr,"  White         //Multiplying by 2 doubles the brightness\n");
fprintf(fptr,"  rgb <255/255,250/255,240/255>*0.7\n");
fprintf(fptr,"  parallel\n");
fprintf(fptr,"}\n");
fprintf(fptr,"\n");
fprintf(fptr,"//Set a background color\n");
fprintf(fptr,"background { color White }\n");
fprintf(fptr,"\n");
fprintf(fptr,"//Create a \"floor\"\n");
fprintf(fptr,"//plane {\n");
fprintf(fptr,"//  <0,0,1>, 0            //This represents the plane 0x+0y+z=0\n");
fprintf(fptr,"//  texture { T_Silver_3A }       //The texture comes from the file \"metals.inc\"\n");
fprintf(fptr,"//}\n");
fprintf(fptr,"\n");
fprintf(fptr,"//Sphere with specified center point and radius\n");
fprintf(fptr,"//The texture comes from the file \"stones.inc\"\n");
fprintf(fptr,"//cylinder { <-50, 0, -20>, <50, 0, 20>, 200 \n");
fprintf(fptr,"//pigment {Blue}\n");
fprintf(fptr,"//finish { ambient .4}\n");
fprintf(fptr,"//}\n");
fprintf(fptr,"\n");
fprintf(fptr,"\n");
fprintf(fptr,"/*object {\n");
fprintf(fptr,"  Round_Cylinder\n");
fprintf(fptr,"   (<0,0,0>,<0,0.3001,0>,0.5,0.15,0)\n");
fprintf(fptr,"  texture{\n");
fprintf(fptr,"    pigment{ color rgb<0.75,0.6,1>}\n");
fprintf(fptr,"    finish { phong 1}\n");
fprintf(fptr,"  } // end of texture\n");
fprintf(fptr,"  scale<1,1,1>\n");
fprintf(fptr,"  rotate<0,0,0>\n");
fprintf(fptr,"  translate<3,0, 0>\n");
fprintf(fptr,"} //------------------------------------*/\n");
fprintf(fptr,"\n");
fprintf(fptr,"#declare Wire_Texture = \n");
fprintf(fptr,"         texture {  pigment{ color rgb<202/255.,202/255.,202/255.> }\n");
fprintf(fptr,"                    finish {\n");
fprintf(fptr,"                    ambient .1\n");
fprintf(fptr,"                    phong .10\n");
fprintf(fptr,"                    phong_size 5}\n");
fprintf(fptr,"                 } // end of texture\n");
fprintf(fptr,"\n");
fprintf(fptr,"#declare Cylinder_Texture = \n");
fprintf(fptr,"texture{\n");
fprintf(fptr,"    //pigment{ color rgb<0.75,0.6,1>}\n");
fprintf(fptr,"   pigment{ color rgb<245/255.,26/255.,15/255.>}\n");
fprintf(fptr,"   finish { phong 1}\n");
fprintf(fptr,"}\n");
fprintf(fptr,"//#macro nucleo(xnu,ynu,znu)\n");
fprintf(fptr,"//--------------------------------------------------------------------------\n");
fprintf(fptr,"//--------------------------------------------------------------------------\n");
fprintf(fptr,"#local D = 0.00001;\n");
fprintf(fptr,"//-------------------------------------------\n");
fprintf(fptr,"#local WR   = 0.9;// Wire radius\n");
fprintf(fptr,"#local WD   = 1.99; // vertical wire distance: min 2*WR  (pitch)\n");
fprintf(fptr,"#local BR   = 7.96; // column base radius (to wire center)\n");
fprintf(fptr,"#local Rmaj = 7.96; // corner radius major of wire \n");
fprintf(fptr,"#local Revolutions = 4 ;// 0.25, 0.5, 0.75, 1, .... \n");
fprintf(fptr,"\n");
fprintf(fptr,"#local rotdegree   = 45;// rotation degrees\n");
fprintf(fptr,"#local nustx = WR/3;\n");
fprintf(fptr,"#local nuedx = 2*WD*sin(rotdegree*pi/180)-WR/3;\n");
fprintf(fptr,"#local nuedy = -2*WD*cos(rotdegree*pi/180);\n");
fprintf(fptr,"//--------------------------------------------------------------------------\n");
fprintf(fptr,"//--------------------------------------------------------------------------\n");
fprintf(fptr,"#local A  = < 0.00, 0.00, -BR >;\n");
fprintf(fptr,"#local B  = <   BR, WD/4, 0.00>;\n");
fprintf(fptr,"#local Co = <   BR, WD/8, -BR >;//corner\n");
fprintf(fptr,"\n");
fprintf(fptr,"//----------------------------------- \n");
fprintf(fptr,"#if(Rmaj <= WR) #local Rmaj = wR+D; #end\n");
fprintf(fptr,"#if(Rmaj >= BR) #local Rmaj = BR-D; #end\n");
fprintf(fptr,"//--------------------------------------------------------------------------\n");
fprintf(fptr,"\n");
fprintf(fptr,"#local Corner_Angle = VAngleD( Co-A,B-Co ); // -> angle in degrees!\n");
fprintf(fptr,"#local Co_Len = Rmaj*tan(radians(Corner_Angle/2)); // corner linear length\n");
fprintf(fptr,"#local Len_X = vlength(Co-A)-Co_Len;  // linear wire part \n");
fprintf(fptr,"#local Wire_Angle = VAngleD( Co-A, <1,0,0>);  // angle of vector ACo against xz plane\n");
fprintf(fptr,"\n");
fprintf(fptr,"#local Th = WD/8*cos(radians(Wire_Angle)); \n");
fprintf(fptr,"#local Inner_Angle = degrees(atan2(Th,BR)); // rotation angle of 2nd part CoB\n");
fprintf(fptr,"\n");
fprintf(fptr,"//--------------------------------------------------------------------------\n");
fprintf(fptr,"#local W_Corner =   \n");
fprintf(fptr,"object{ Segment_of_Torus( Rmaj, // radius major,\n");
fprintf(fptr,"                          WR,   // radius minor,\n");
fprintf(fptr,"                          -Corner_Angle // segment angle\n");
fprintf(fptr,"                        ) //----------------------------\n");
fprintf(fptr,"        rotate<0,90,0> \n");
fprintf(fptr,"        translate<-0,0,+Rmaj>\n");
fprintf(fptr,"      } // -----------------------------------------------------------------\n");
fprintf(fptr,"//--------------------------------------------------------------------------\n");
fprintf(fptr,"#local Quarter =   \n");
fprintf(fptr,"union{ \n");
fprintf(fptr,"     cylinder{ <0,0,0>,< Len_X,0,0>, WR  } \n");
fprintf(fptr,"     object  { W_Corner  rotate<-Inner_Angle,0,0>  translate< Len_X,0,0> }  \n");
fprintf(fptr,"     cylinder{ <0,0,0>,< Len_X,0,0>, WR \n");
fprintf(fptr,"               translate<0,0,-Rmaj> rotate<0,-Corner_Angle,0> translate<0,0,Rmaj>\n");
fprintf(fptr,"               rotate<-Inner_Angle,0,0> \n");
fprintf(fptr,"               translate<Len_X,0,0> \n");
fprintf(fptr,"             } //-------------------- \n");
fprintf(fptr,"  translate<0,0,-BR>\n");
fprintf(fptr,"  rotate<0,0,Wire_Angle>\n");
fprintf(fptr,"} // end of union  ----------------------------------------------------------\n");
fprintf(fptr,"//---------------------------------------------------------------------------\n");
fprintf(fptr,"\n");
fprintf(fptr,"\n");
fprintf(fptr,"// Winding Wire\n");
fprintf(fptr,"#local Nucleo = \n");
fprintf(fptr,"union{ //-----------------------------------\n");
fprintf(fptr,"\n");
fprintf(fptr,"\n");
fprintf(fptr," #local Nr = 0;                // start\n");
fprintf(fptr," #local EndNr = 2*Revolutions; // end\n");
fprintf(fptr," #while (Nr< EndNr)            // loop\n");
fprintf(fptr,"   object{ Quarter\n");
fprintf(fptr,"           translate<0, Nr*WD/4,0>\n");
fprintf(fptr,"           rotate<0,-Nr*90,0>\n");
fprintf(fptr,"           texture{ Wire_Texture } \n");
fprintf(fptr,"         } //----------------\n");
fprintf(fptr,"\n");
fprintf(fptr," #local Nr = Nr + 1; // next Nr\n");
fprintf(fptr," #end // ---------------// end of loop\n");
fprintf(fptr,"// translate<0,WR ,0>\n");
fprintf(fptr,"object {\n");
fprintf(fptr,"  Round_Cylinder\n");
fprintf(fptr,"   (<0,-1.99,0>,<0,5.97,0>,7.96,2.0,0)\n");
fprintf(fptr,"  texture{ Cylinder_Texture }\n");
fprintf(fptr,"} //------------------------------------\n");
fprintf(fptr," scale<1,-1,1>\n");
fprintf(fptr," rotate<0,0,rotdegree> \n");
fprintf(fptr," //translate<xnu,ynu,znu>\n");
fprintf(fptr,"} // end of union -----------------------------------------------------------\n");
fprintf(fptr,"//---------------------------------------------------------------------------\n");
fprintf(fptr,"//--------------------------------------------------------------------------- \n");
fprintf(fptr,"//}\n");
fprintf(fptr,"//#end\n");

        double curpos = v[stp];
        double oldpos = v[stp];
        
        fprintf(fptr,"cylinder{ <nuscal*%.3f+nustx-200,0,-BR>,<nuscal*%.3f+nustx,0,-BR>, WR  \n",curpos,curpos);
        fprintf(fptr,"texture{ Wire_Texture } }\n");
            
        fprintf(fptr,"object{ Nucleo\n");
        fprintf(fptr,"translate<nuscal*%.3f,0,0>}\n",curpos);
        //cout<<v[stp]<<endl;
        for(int i=stp+1;i<edp;i++){
            //cout<<v[i]<<endl;
            
            curpos += v[i]-v[i-1]-100;
            fprintf(fptr,"cylinder{ <nuscal*%.3f+nuedx,nuedy,-BR>,<nuscal*%.3f+nustx,0,-BR>, WR  \n",oldpos,curpos);
            fprintf(fptr,"texture{ Wire_Texture } }\n");
            
            fprintf(fptr,"object{ Nucleo\n");
            fprintf(fptr,"translate<nuscal*%.3f,0,0>}\n",curpos);
            oldpos = curpos;
        }
            fprintf(fptr,"cylinder{ <nuscal*%.3f+nuedx,nuedy,-BR>,<nuscal*%.3f+nuedx+200,0,-BR>, WR  \n",curpos,curpos);
            fprintf(fptr,"texture{ Wire_Texture } }\n");
        

    fclose(fptr);    
    //fclose(fptrinit);
    //VectorDividedByScalar(gr,num_configs);
    //write2DVector(savegrname,arange(0,lgr,dx),gr);

    return 0; 
} 




















