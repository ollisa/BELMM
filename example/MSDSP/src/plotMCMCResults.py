f'''
AUTHOR: John Chavis
DATE: 2019-NOV-12
NAME: plotMCMCResults.py
DESC: This script will plots resuts from software


MIT License

Copyright (c) 2020 John T. Chavis III

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''


#Imprort packages
import numpy as np
import math

import seaborn as sns
sns.set_style('white')
#sns.set_palette("coolwarm")
#sns.set_palette('ocean')
#plt.rcParams['axes.color_cycle'][2]
import matplotlib
import matplotlib.pyplot as plt

#COllect files 
import os 
import string 

#FOr getting argument s
import sys



#Where to save pics
picDir = 'pics/' 




#Hisgotram info 
BINS = 50 

#Define colors. Thanks to Amy Cochran 

matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['xtick.labelsize'] = 10
#matplotlib.rcParams['font.family'] = "Arial"

#For colors 
lineStyleList = ['solid', 'dotted', 'dashed', 'dashdot']

brownGreen = np.array([[140,81,10], 
				[1,102,94],
				[216,179,101],
				[90,180,172],
				[246,232,195],
				[199,234,229]])/255.0

redBlueYellow = np.array([[215,48,39],
				[69, 117, 180], 
				[252, 141, 89],
				[145, 191, 219], 
				[254, 224, 144],
				[224, 243, 248]])/255.0

amyColors = np.array([[255*0.15,255*0.15,255*0.14],
					[0,112,192],
					 [0,255,255],
					 [0, 130, 65],
					 [146, 208, 80],
					 [255, 255, 0],
					 [255, 0, 102],
					 [112, 48, 160]])/255.0
#line styles
colorDic = {'brownGreen': brownGreen, 'redBlueYellow': redBlueYellow, 'amyColors':amyColors}

#COLORCHOICE = 'brownGreen'
#COLORCHOICE = 'redBlueYellow'
COLORCHOICE = 'amyColors'
colorOpt = []

for line in lineStyleList: 
	for color in colorDic[COLORCHOICE]:
		colorOpt.append([color, line])

LINEWIDTH = 0.5
ALPHA = 0.7



#Set sizes of plots
#Default font size
#Figure Sizes
HEIGHT = 20
WIDTH = 20

DEFAULTSIZE = 18
AXTITLESIZE = 45
LABELSIZE = 32
TICKSIZE = 28
LEGENDSIZE = 28
FIGTITLESIZE = 32

plt.rc('font', size=DEFAULTSIZE)

#Axes title and labe size 
plt.rc('axes', titlesize = AXTITLESIZE)
plt.rc('axes', labelsize = LABELSIZE)
plt.rc('xtick', labelsize = TICKSIZE)
plt.rc('ytick', labelsize=TICKSIZE)
plt.rc('legend', fontsize=LEGENDSIZE)
plt.rc('figure', titlesize=FIGTITLESIZE)

#Fraction of maximum to add for setting yLim 
MAXFRAC = 0.1
#print(colors)

#-------------------------------------------------- 
# 
#
#
# 				FUNCTIONS 
#--------------------------------------------------
### Should define fucntions to readability 
#NAME: readData
#INPUT: None 
#Output: Dictionary containing acceptance probabilies, data-with model distirbution, and paramters per chain 
#DESC: Read in Data from Gold Standard Chains and Model 
#  Techinically dont need the total modle count or chain count but it helps
def readData(): 
	#Collect Acceptance Probability Info and Model INde 
	dataDir = 'modelSelection/'

	#Read in model indexing
	f = open('modelIndex.txt','r')
	#Read models in 
	modelDic = {}

	#Read off header
	header = f.readline()
	modCount = 0 
	for line in f.readlines():
		line = line.strip().split(',')

		model = line[0]
		modelIndex = int(line[1])

		modelDic[modelIndex] = model
		modCount += 1


	#Sort models base on index
	modelNameList = [] 
	for k in range(modCount):
		modelNameList.append(modelDic[k])


	#Figure out the number of chains in the file  
	chainNum = []
	for file in os.listdir('{}'.format(dataDir)):
		#Get chain 
		#print(file)
		if file.endswith(".csv"):
			#print(file.split('.csv'))
			chain = int(file.split('.csv')[0][-1])
			chainNum.append(chain)
	

	chainMax = max(chainNum)

	#Great we have models and chains..
	#Now to create a dictionary to store data
	#Each  entry list will contain an array of results per chain
	dataDic = {'accept':[], 'modelIndex':[], 'modelParam':{}, 'modelName':modelNameList}

	for k in range(modCount):
		modelName = modelDic[k]
		#dataDic['modelParam']['model{}'.format(modelName)] = []
		#dataDic['modelParam']['model{}Header'.format(modelName)] = 'TBD'
		dataDic['modelParam'][modelName] = []
		dataDic['modelParam']['{}Header'.format(modelName)] = 'TBD'
		
	for name in ['accept', 'modelIndex']:		
		for chain in range(1, chainMax+1): 
			fileName = '{}{}Chain{}.csv'.format(dataDir,name, chain)

			data = []
			#Have found it faster to open line by line as opposed to a dataframe
			f = open(fileName,'r')

			#Beacuse it's a just a count of values.....for now 
			for line in f.readlines():
				line = line.strip()

				if name == 'accept':
					data.append(float(line))
				else:
					data.append(int(line))

			
			dataDic[name].append(np.array(data[:]))
			data[:] = []
			f.close()

	#Now we need to Collect the Paramter info for each model 
	#This is handy for trace plots etc, log's etc. 
	dataDir = 'goldStandardChains/'
	for k in range(modCount): 
		#modelName = 'model{}'.format(k)
		modelName = modelDic[k]
		for chain in range(1, chainMax+1): 
			fileName = '{}{}/{}GSC{}.csv'.format(dataDir, modelName,modelName,chain)

			data = []
			f = open(fileName, 'r')
			#Grab header info to figure out which variable is the log model
			#Acutally it's always the last column we should be fine. Assume float 
			header = f.readline().strip().split(',')

			for line in f.readlines():
				line = line.strip().split(',')
				#Read info 
				tmpData = [float(x) for x in line]

				data.append(tmpData[:])

			#Add this into proper place
			dataDic['modelParam'][modelName].append(np.array(data[:]))

			if dataDic['modelParam']['{}Header'.format(modelName)] == 'TBD':
				dataDic['modelParam']['{}Header'.format(modelName)] = header 

			data[:] = []
			f.close()

	#Reutrn modCount, chainMax adn data dictionary 

	return modCount, chainMax, dataDic, modelDic


#--------------------------------
#NAME: plotModelDist 
#INPUT: modCount, chainMax, dataDic
#OUTPUT: figure containing distribution of models 
#DESC: This funciton will plot distributon of model index in stationary state
def plotModelDist(modCount, dataDic): 
	#GEt length so we can figure out burn IN. Assusme burn in is 1/2 point 
	chainLength = len(dataDic['modelIndex'][0])
	burnIn = math.floor(chainLength/2.0)

	#Set up counts. Start a 0 for each modle 
	counts = np.zeros(modCount) 

	for data in dataDic['modelIndex']:
		for val in data[burnIn:]:
			counts[val] += 1

	#Convert counts to probability 
	counts = np.array(counts)
	modProbs = counts/np.sum(counts)

	x = np.arange(modCount)

	#distAx.hist(allData, 50, density= True, facecolor='blue', alpha=0.7)
	modelFig, modelAx = plt.subplots(figsize=((WIDTH,HEIGHT)))

	#print(counts)
	#print(dataList)
	barList = modelAx.bar(x, modProbs, facecolor= colorOpt[3][0], alpha=ALPHA)

	for k in range(modCount):
		barList[k].set_color(colorOpt[k][0])

	modelAx.set_title('Model Distribution')
	modelAx.set_xlabel('Model Index', fontweight='bold')
	modelAx.set_ylabel('Counts', fontweight='bold')
	plt.xticks(x,dataDic['modelName'])

	return modelFig

#--------------------------------
#NAME: plotAcceptance 
#INPUT: modCount, chainMax, dataDic
#OUTPUT: figure containing distribution of models 
#DESC: This funciton will plot distributon of model index in stationary state
def plotAcceptance(chainMax, dataDic): 
	#GEt length so we can figure out burn IN. Assusme burn in is 1/2 point 
	chainLength = len(dataDic['accept'][0])

	### PLot Acceptrances 
	acceptFig, acceptAx = plt.subplots(figsize=((WIDTH,HEIGHT)))
	
	#chainLength = chainLength/2
	for chain in range(chainMax):
		data = dataDic['accept'][chain]
		acceptAx.plot(np.arange(0,chainLength), data[:], color= colorOpt[chain][0], linestyle=colorOpt[chain][1], linewidth = LINEWIDTH, label= 'Chain {}'.format(chain+1), alpha=ALPHA)

	#fullAx.set_xlabel('Time')
	acceptAx.set_xlabel('Iteration', fontweight='bold')
	acceptAx.set_ylabel('Accept Ratio', fontweight='bold')
	acceptAx.set_ylim((0,1))
	#acceptAx.set_xlabel('Chain')
	acceptAx.legend(loc='lower right',facecolor='white', framealpha=0.3,ncol=2)
	acceptAx.set_title('Acceptance Probability')

	return acceptFig


#------------------------------
#NAME: plotParamter
#INPUT: modCount, chainMax, dataDic
#OUTPUT: figure containing Distribution of paramters and trace
#DEC: FUnction plots each models parmters distribution adn trace plot
def plotParamters(modCount, chainMax, dataDic,modelDic):
	#Generate a plot (these we are not going to return as there may be a lot)
	
	figList = False

	#print(data[:,0])
	for k in range(modCount): 
		#modelName = 'model{}'.format(k)
		modelName = modelDic[k]
		#GEt total number of partemrs (remember -1 as the last column is logModel)
		totParam = len(dataDic['modelParam']['{}Header'.format(modelName)]) - 1
		#Collect distribution adn trace for each paramter....We need the data 
		for param in range(totParam): 
			#For now we'll just plot the Paramter Distribution, Trace PLot
			paramFig, paramAx = plt.subplots(2,1, figsize=((WIDTH,HEIGHT+2)))

			parName = dataDic['modelParam']['{}Header'.format(modelName)][param]
			#For aggerating the stationary data 
			statData = []

			for chain in range(1,chainMax+1):
				#Grab chain
				#print(chain-1, param)
				data = dataDic['modelParam'][modelName][chain-1][:,param] 

				#Assuming that resulst are stationary....again we'll need to factor in 
				#burn in later
				chainLength = len(data)
				x = np.arange(0,chainLength)
				#Aggregate all the data for xthe distribution 
				statData += list(data)

				#Plot traces 
				paramAx[1].plot(x, data, color= colorOpt[chain-1][0], linestyle=colorOpt[chain-1][1], linewidth = LINEWIDTH, label= 'Chain {}'.format(chain), alpha=ALPHA)
				
			paramAx[1].set_xlabel('Iteration', fontweight='bold')
			paramAx[1].set_ylabel(parName, fontweight='bold')
			paramAx[1].set_title('Trace Plot')
			paramAx[1].legend(loc='lower right',facecolor='white', framealpha=0.3,ncol=2)
			#paramAx[1].legend(loc='lower right',facecolor='black', framealpha=0.3,ncol=2)

			#PLot Distribution 	
			paramAx[0].hist(data, BINS,  density= True, facecolor= colorOpt[k][0], alpha=ALPHA)
			paramAx[0].set_xlabel(parName, fontweight='bold')
			paramAx[0].set_ylabel('Denisty', fontweight='bold')
			paramAx[0].set_title('Parameter Distribution')

			#save fig 
			paramFig.savefig('{}{}/{}Plot.png'.format(picDir,modelName,parName))
			#plt.cla()
			#plt.clf()

			#figList.append(paramFig)
			paramFig.clear()
			plt.close(paramFig)

		#Plot Loglikelihod for the models 
		logFig, logAx = plt.subplots(figsize=((WIDTH,HEIGHT)))
		for chain in range(1,chainMax+1):
			data = dataDic['modelParam'][modelName][chain-1][:,totParam] 
			#Assuming that resulst are stationary....again we'll need to factor in 
			#burn in later
			chainLength = len(data)
			x = np.arange(0,chainLength)

			#Plot traces 
			logAx.plot(x, data, color= colorOpt[chain-1][0], linestyle=colorOpt[chain-1][1], linewidth = LINEWIDTH, label= 'Chain {}'.format(chain), alpha=ALPHA)
			
			logAx.set_xlabel('Iteration', fontweight='bold')
			logAx.set_ylabel('Log-Likelihood', fontweight='bold')
			logAx.legend(loc='lower right',facecolor='white', framealpha=0.3,ncol=2)
			logAx.set_title('{} Log-Likelihood Plot'.format(dataDic['modelName'][k]))
			#logAx.legend(loc='lower right',facecolor='black', framealpha=0.3,fontsize=12,ncol=2)

			logFig.savefig('{}{}/logLikelihood.png'.format(picDir, modelName))
       

		#figList.append(logFig)
		logFig.clear()
		plt.close(logFig)
	figList = True
	
	return figList 


#################################################
#
#
#					MAIN 
#
#
#################################################
def main():
	#Read data 
	print('Plotting Results...')
	modCount, chainMax, dataDic,modelDic = readData()

	#Plot Distirbution of Models
	modelFig = plotModelDist(modCount, dataDic)
	acceptFig = plotAcceptance(chainMax, dataDic)


	figList = plotParamters(modCount, chainMax, dataDic,modelDic)

	modelFig.savefig('{}ModelDist.png'.format(picDir))
	acceptFig.savefig('{}Accept.png'.format(picDir))

	print('Figures plotted and stored in pics/')
	#plt.show()

if __name__ == "__main__":
	main()
