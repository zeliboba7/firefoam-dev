import numpy as np
fw  = open('phiHMeanVsHeight','w')
startTime = 5.0 
nStart = 0
deltaH = 0.1
folderName = "zones/0.0011976/"

# read data from file
dataHs = np.loadtxt(folderName+"phiHs",skiprows=2)
dataHc = np.loadtxt(folderName+"phiHc",skiprows=2)
dataH  = np.loadtxt(folderName+"phiH",skiprows=2)

# transpose the array
dataHs = dataHs.transpose()
dataHc = dataHc.transpose()
dataH  = dataH.transpose()

# init arrays
time = dataHs[0]
phiHs = dataHs[1:,:]
phiHc = dataHc[1:,:]
phiH  = dataH[1:,:]

nTime = time.size
nHeight = phiHs.shape[0]

phiHsAve = np.zeros((nHeight))
phiHcAve = np.zeros((nHeight))
phiHAve  = np.zeros((nHeight))
height = np.zeros((nHeight))
timeStep = np.zeros((nTime-1))

for i in range(nTime-1):
    timeStep[i] = time[i+1]-time[i]
    if nStart==0 and time[i]>startTime:
        nStart=i

print 'nStart is ', nStart

aveTime = time[nTime-1]-time[nStart]


fw.write("height phiHs phiHc phiH\n ")
for i in range(nHeight):
    height[i]=float(i+1)*deltaH 
    phiHsAve[i] = 0.0
    phiHcAve[i] = 0.0
    phiHAve[i]  = 0.0
    for j in range(nStart,nTime-1):
	phiHsAve[i] = phiHsAve[i] + 0.5*(phiHs[i,j]+phiHs[i,j+1])*timeStep[j]
	phiHcAve[i] = phiHcAve[i] + 0.5*(phiHc[i,j]+phiHc[i,j+1])*timeStep[j]
	phiHAve[i]  = phiHAve[i]  + 0.5*(phiH[i,j] +phiH[i,j+1])*timeStep[j]

    phiHsAve[i]=phiHsAve[i]/aveTime
    phiHcAve[i]=phiHcAve[i]/aveTime
    phiHAve[i] =phiHAve[i]/aveTime
    fw.write(str(height[i])+" "+ str(phiHsAve[i])+" " + str(phiHcAve[i]) + " " + str(phiHAve[i]) + "\n")
    




