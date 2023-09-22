##!/bin/bash
#-------------------------
#AUTHOR: John Chavis 
#DATE: 2020-MAY-24
#DESCRIPTION: This sciprt is my attempt to utilzie MC-STAN with our Parallel approach 
#				Essentially MC-STAN will conduct within model samples (gold standard chains)
#				We then use our C program for RJMCMC Crawler 
#
#TO DO
# Find a way 
#
#MIT License
#
#Copyright (c) 2020 John T. Chavis III
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#------------------------


############################################
#
#
#			READ CONFIG FILE 
#
###########################################
#Code snpit from https://stackoverflow.com/questions/16571739/parsing-variables-from-config-file-in-bash

#Useful for delting trailign spaces
shopt -s extglob
configfile="setup.config" # set the actual path name of your (DOS or Unix) config file
#deletes DOS carriage return.
tr -d '\r' < $configfile > $configfile.unix
while IFS='= ' read -r lhs rhs
do
		#Ignore empty lines
    if [[ ! $lhs =~ ^\ *# && -n $lhs ]]; then
        rhs="${rhs%%\#*}"    # Del in line right comments
        rhs="${rhs%%*( )}"   # Del trailing spaces
        rhs="${rhs%\"*}"     # Del opening string quotes 
        rhs="${rhs#\"*}"     # Del closing string quotes 
        rhs="$(echo "${rhs}" | tr -d '[:space:]')" #Remove white space
        declare $lhs="$rhs"
    fi
done < $configfile.unix

rm $configfile.unix
#Echo Results so user can See

printf '=%.0s' {1..80}
echo "cmdstan path is $CMDSTANPATH" 
echo "Number of Cores for cmdstan is $NCORES"
echo "Number of Burn-in (warm-up) samples is $NBURN"
echo "Number of Samples to Collect is  $NLEN"
echo "Number of NCHAIN $NCHAIN"

if $PLOT; then 
    echo "Plotting set to True"
else
    echo "Plotting set to False"
fi


printf '=%.0s' {1..80}


if [ -d  goldStandardChains/ ]; then 
  rm -r goldStandardChains/ 
fi	 	
 	mkdir goldStandardChains/ 


#File name for timing 
TIMEFILE=timings.txt

#if [ -f $TIMEFILE ]; then 
#    rm $TIMEFILE 
#fi

#Grab Current Direcotry 
parWD=`pwd`

#Find number of models 
NMODEL=`find models/ -name "*.stan" -type f | wc -l | xargs`

#Create a directory to store results gold standard chains   
for ((i=0; i<$NMODEL; i++))
do
	MODELDIR="goldStandardChains/model${i}"
	if [ -d  "$MODELDIR" ]; then 
	  rm -r "$MODELDIR" 
	fi	 	
   	mkdir "$MODELDIR"
done


#Make a Dicectory to store proccessed models. Proccessed models take 
#The user defined models (whatever name) and createa Model index for these models
#I.e. 
#Model Index | Model Name 
# 0 				 | Model A 
# 1 				 | Model B 


if [ -d  proccessedModel/ ]; then 
	rm -r proccessedModel/
fi

mkdir proccessedModel/


#Move into the user defined models
cd models/

#Create a dictionary for model index
#Rename files
#I.e. 
#Model A  	0 
#Model B	  1
#....
#Model N 		N 

#STart time stamping 
totalStart=`perl -MTime::HiRes -e 'printf("%.0f\n",Time::HiRes::time()*1000)'`

echo "Proccessing Models..."
MODCOUNT=0
echo "Name,Index" > ../modelIndex.txt


#Find all files
#Keep track of index for renaming
#This help accoutn for file names with spaces
SAVEIFS=$IFS
IFS=$(echo -en "\n\b")
 
MODELARRAY=()
for i in `ls -v *.stan`; do
	#Get name of stan file 
    echo "$i"
    [ -f "$i" ] || break
    #Add model to Index
    echo "${i%.*},${MODCOUNT}" >> ../modelIndex.txt

    #Copy model file over
   	cp "$i" "../proccessedModel/model${MODCOUNT}.stan"

   	#Write to this file 
    #Add to array 
    MODELARRAY+=("${i%.*}")
	  MODCOUNT=$((MODCOUNT+1))
done

#Revet back 
IFS=$SAVEIFS
#echo "${MODELARRAY[@]}"
#Return back to main directoy 
cd ../

###################################
#
#						MC-STAN
#
##################################
#First Buiid Mc-Stan - See cmdSTAN documentation 
cd $CMDSTANPATH
make clean-all
make build -j"$NCORES"



#Build Each of the Users Stand Models 
# This is done in a parallel on cmd cores
echo "Building Models...."
for ((i=0; i<$NMODEL; i++))
do
    ((core=core%$NCORES)); ((core++==0)) && wait
	make "$parWD/proccessedModel/model$i" > /dev/null & 
	echo "Model $i Built."
done


#Make sure all models are built
wait 


#Time stamp of how long it took to build model 
totalBuilt=`perl -MTime::HiRes -e 'printf("%.0f\n",Time::HiRes::time()*1000)'`


cd "$parWD"

#Now run each Model. Done in Parallel 
#RUn for N Bathces based on cmd cores
echo 'Running Models'
cd proccessedModel/
for ((k=0; k< $NMODEL; k++))
do 
	for ((chain=1; chain<=$NCHAIN; chain++))
	do 
	    ((core=core%$NCORES)); ((core++==0)) && wait
	    echo "Model $k Chain $chain Runnning...."
		"./model$k" sample num_samples="$NLEN" num_warmup="$NBURN" \
		random seed=12345 id="$chain" \
		data file=../data/data.json output file="../goldStandardChains/model${k}/model${k}StanResults${chain}.csv" > /dev/null &
		#pids[${chain}]=$!
		#wait		
	done

done

# wait for all pids
#for pid in ${pids[*]}; do
#    wait $pid
#done
wait 


cd ../goldStandardChains/

#Now grab the information we need. For now this is hard coded. We just need paramters after Energy....
#Assuming Current cmd STAN result structure. 
for ((k=0; k< $NMODEL; k++))
do 
	(cd "model${k}/" 
	#Grab header for files
	#Always 1 chain so we make this work 
	grep lp__ "model${k}StanResults1.csv" > columnNames.txt

	#Find where energy is. The rest contian the info we need. 
	COLCUT=`head -n1 columnNames.txt | tr "," "\n" | grep -nx energy__ | cut -d":" -f1`


	#Now write results 
	for ((chain=1; chain<=$NCHAIN; chain++))
	do  
		#Create a temprory file 
		sed '/^[#1]/d' "model${k}StanResults${chain}.csv" >> tmp.csv 
		awk -v cut="$COLCUT" -F, '{for(i=cut+1;i<=NF;i++) printf $i","; print ""}' tmp.csv > t.csv
    sed 's/.$//' t.csv >"model${k}GSC${chain}.csv" 
		rm tmp.csv
    rm t.csv
	done 

	rm columnNames.txt
	cd ..) & 
done 
 wait


#USe a python script to create matrices to be read in C. 
#Essentially .....
#Then it's a 1d array with data. First two elements will be row and column info.
cd ../
MAXPARAM=`python3 convertCMatrix.py "$NMODEL"|xargs `

#Time stamp to run gold standar modles
totalGSC=`perl -MTime::HiRes -e 'printf("%.0f\n",Time::HiRes::time()*1000)'`


##################################################
#
#							RJMCMC Crawler 
#
##################################################
#Great we have converted the gold standd chanins converted and the max paramters. 
#NOw time to write the RJMCMC Crawler (written in C)
CURDATE=$(date +"%Y-%m-%d")

#Create File 
FILE=tmpRJMCMC.c

#Final file name and name of main body for RJMCMC crawler
FINALFILE=runRJMCMC.c
SKELTON=bodypRJMCMC.c

#Initialize File 
if [ -f $FILE ]; then 
    rm $FILE 
fi

if [ -f $FINALFILE ]; then 
    rm $FINALFILE     
fi

touch "$FILE"
touch "$FINALFILE"



#Create Header
cat >> "$FILE" <<EOL
/*
DATE: $CURDATE 
NAME: runRJMCMC.c 
DESC: Bash Generated Script to Run Parallel RJMCMC Algorithm

    Run make sure to properly 
    Number of process is NMODES*CHAINS+ + CHAINS+ 1
 
    mpicc -o <name> <file>
    mpirun -np <nProc>
    Make sure to 
    export TMPDIR=/tmp

*/

/* Required Headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*MPI Header */
#include <mpi.h>
#include <string.h>

/*GNU Scientif Library Packages */
//Factorial 
//#include <gsl/gsl_sf_fact.h>
//Rando Number
#include <gsl/gsl_rng.h>



EOL


#INclude proper lenght assume burn-in is double the Sample lenght for now 
SAMPLELEN=$((NLEN*2))
#Inlcude Constants 
#cat << EOL | tee -a "$FILE"
cat >> "$FILE" <<EOL


/* Constants */
//Simulation Info
#define CHAINS $NCHAIN
#define CHAINLEN $SAMPLELEN 
#define NMODEL $NMODEL
#define COLLEN NMODEL*CHAINS
#define MAXPARAM $MAXPARAM

#define STRINGLEN (35)*sizeof(int)

//Differnt Keneral Density Estiamte (kde) implemnted
#define UNIFORM 0 
#define GAUSSIAN 1


#define PI 3.14159265358979323846

EOL


#-------------------------------------------------
#
#               WRAPPERS
#
#-------------------------------------------------
echo "Building Script..."
#-------------Proposed Models ----------
#This allows a user to jump between models 
#Assumes Uniform transition (i.e. 1/|K|)
#Here we will write the c function to jump between models 
STEP=$(echo "scale=6; 1/$NMODEL" |bc)

#Start probability
PROB="$STEP"

cat >> "$FILE" << EOL 
/* ------------------------------------------ */
/*FUNCTION: proposedModel
INPUT: current model
OUTPUT: Calculates proposed move
DESC: Get propsoed move.  */
int proposedModel(int k){
  double pU; 
  int pMove; 
  
  pU = ((double)rand()/ (double)RAND_MAX);

  if (pU < $PROB){
    pMove = 0;
  }

EOL


for ((k=1; k<"$NMODEL" - 1; k++))
do
  #echo "We were here $k times"
  PROB=$(echo "scale=6;$PROB + $STEP" |bc)
  echo "  else if (pU < $PROB){" >> "$FILE"
  echo "    pMove = $k;" >> "$FILE"
  echo "  }" >> "$FILE"
       
done 


#Last Probabilty 
cat >> "$FILE" << EOL
  else{
    pMove = $k;
  }
  return pMove;
}




EOL




#-------------Initial Models----------
#Used to pick an inital model to run 
PROB="$STEP"
cat >> "$FILE" << EOL 
/* ------------------------------------------ */
/*FUNCTION: initialModel
INPUT: 
OUTPUT: Initial Model
DESC: intialize model  */
int initModel(){
  double uInitModel; 
  int initModel;  

  uInitModel = ((double)rand()/ (double)RAND_MAX);
  

  if (uInitModel < $PROB){
    initModel = 0;
  }

EOL


for ((k=1; k<"$NMODEL" - 1; k++))
do
  #echo "We were here $k times"
  PROB=$(echo "scale=6;$PROB + $STEP" |bc)
  echo "  else if (uInitModel < $PROB){" >> "$FILE"
  echo "    initModel = $k;" >> "$FILE"
  echo "  }" >> "$FILE"
      
  #echo "$PROB"
done 


#Last Probabilty 
cat >> "$FILE" << EOL
  else{
    initModel = $k;
  }
  return initModel;
}


EOL



#-----------------------------
# 
#  
#    Include  main Script 
#
#
#
#-------------------------------
#Ok we need concat skeleton in there
#See bodypRJMCMC for more details 
cat "$FILE" "$SKELTON" >> "$FINALFILE"

rm "$FILE"

#Create directory for parallel results 
if [ -d modelSelection ]; then 
    rm -r modelSelection 
fi

mkdir modelSelection


echo "Running Crawler..."
#calculate Number of processors 
NPROC=$(($NCHAIN+1))

#For Mac we need to make sure to export tmpdir 
export TMPDIR=/tmp

#Make a directory tos reo the pictures
if [ -d  pics/ ]; then 
  rm -r pics/
fi

mkdir pics/

#Create a subdirectory for each model 
for ((i=0; i<$NMODEL; i++))
do
  mkdir pics/"model$i"
done


mpicc -o runRJMCMC runRJMCMC.c -lm -ldl -lgsl -lgslcblas  && mpirun -np "$NPROC" --oversubscribe ./runRJMCMC 

#Time stamp to run crawler
totalCrawler=`perl -MTime::HiRes -e 'printf("%.0f\n",Time::HiRes::time()*1000)'`

#Now  Let's try to change files and folders back to 
#Chagne file names first
for ((k=0; k< $NMODEL; k++))
do 
  #Grab modelName 
  MODELNAME=${MODELARRAY[k]}
  #echo "$MODELNAME"
  #Rename directories first...
  find . -type d -name "model${k}" | while read DIRECTORY; do 
    newDir="$(echo ${DIRECTORY} | sed -e "s|model${k}|${MODELNAME}|")";
    mv "$DIRECTORY" "$newDir"
  done 

  #Find Files of this nature 
  find "./goldStandardChains/${MODELNAME}/" -type f -name "model${k}*" | while read FILE; do 
    #echo "${k} ${MODELNAME}"
    newfile="$(echo ${FILE} | sed -e "s|model${k}|${MODELNAME}|")";
    #Replace names
    mv "$FILE" "$newfile"
    #echo "hey";
  done 
done

echo "Done."

#WE'll remove processed models for now....
rm -r proccessedModel/

#Should also remove redudent files in gold standard chains 
find . -type f -name "gsc*.csv" -delete


#Add  TIme infO
echo "Process,Time(s)" > $TIMEFILE 

#Add time to build models
buildTime=$((totalBuilt-totalStart)) 
gscTime=$((totalGSC-totalBuilt))
crawlerTime=$((totalCrawler-totalGSC))

echo "Build Models,$buildTime" >> $TIMEFILE 
echo "Gold Standard Chains,$gscTime" >> $TIMEFILE 
echo "RJMCMC Crawlwer,$crawlerTime" >> $TIMEFILE 


#For now keep this file as we can use it for plotting..
#find . -type f -name "*GSC*.csv" -delete

if $PLOT; then 
    python3 plotMCMCResults.py
fi


exit 
