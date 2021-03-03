#PBS -l nodes=1:ppn=1:medium_buster
#PBS -q medium_buster
#PBS -e ncxxx_yyy.errout
#PBS -o ncxxx_yyy.resout
#PBS -l walltime=12:00:00
#PBS -l mem=500mb,vmem=500mb

echo "START_TIME           = `date +'%y-%m-%d %H:%M:%S %s'`"
##ichr=
tdr=SAM2chr2-xxx_program1_yyy  # at least one xxx remain
underWS=/remote/pi310b/li/MolecularMC/v8ES1/
# need to time walltime

cptarget1=${underWS}/*.cpp
cptarget2=${underWS}/*.h
cptarget3=${underWS}/*config
cptarget4=${underWS}/*gr
cptarget5=${underWS}/*.py
#cptarget6=${underWS}/chr4-${ichr}/pot.can

cpback1=${underWS}/chr2-xxx
#tdr=xxx_yyy_prolong

# in case intermediate interruption
cd /tmp
if [ -d li_${tdr} ]
then
    rm -rf li_${tdr}
fi
mkdir li_${tdr}

if [ -f li_${tdr}.tar.gz ]
then
    rm li_${tdr}.tar.gz
fi
#cd ../
#cd /tmp/chu_${tdr}
cd li_${tdr}

##### start command
cp ${cptarget1} .
cp ${cptarget2} .
cp ${cptarget3} .
cp ${cptarget4} .
cp ${cptarget5} .
#cp ${cptarget6} .

g++ -fPIC -shared -o cdll.so cdll.cpp
./runMolecularMC.py SAM2 xxx yyy
#g++ runES1.cpp
#./a.out ${ichr} xxx

##### data transfer back
#cp *.like_bed $cpback1
#cp *_ctcfbin $cpback1
#cp *_ctcfbin47* $cpback1
cp *pot $cpback1
#cp *201res $cpback1
#cp *201gr $cpback1
#cp _chrxxx.like_wig $cpback1
#cp chrxxx $underWS
#cp ../*.out $cpback1

cd ../
#tar -czvf li_${tdr}.tar.gz li_${tdr}/*
rm -rf ./li_${tdr}
#cp ./li_${tdr}.tar.gz $underWS
#rm ./li_${tdr}.tar.gz
### cd out, important in some cases
cd ../


echo "END_TIME             = `date +'%y-%m-%d %H:%M:%S %s'`"
