import numpy as np
fw  = open('phiHMeanVsHeight','w')
startTime = 10.0 
nStart = 0
deltaH = 0.1
#folderName = "zones/1.19998e-05/"
folderName = "zones/0.0011976/"

# read data from file
dataHs = np.loadtxt(folderName+"phiHs",skiprows=2)
dataHc = np.loadtxt(folderName+"phiHc",skiprows=2)
dataH  = np.loadtxt(folderName+"phiH",skiprows=2)
dataInlet  = np.loadtxt("patchInlet/0/faceSource.dat", skiprows=3)
dataOutlet = np.loadtxt("patchOutlet/0/faceSource.dat",skiprows=3)
dataHRR = np.loadtxt("HRR/0/cellSource.dat",skiprows=3)

# transpose the array
dataHs = dataHs.transpose()
dataHc = dataHc.transpose()
dataH  = dataH.transpose()
dataInlet   = dataInlet.transpose()
dataOutlet  = dataOutlet.transpose()
dataHRR  = dataHRR.transpose()

# init arrays
time = dataHs[0]
phiHs = dataHs[1:,:]
phiHc = dataHc[1:,:]
phiH  = dataH[1:,:]
HInlet  = dataInlet[4]
HsInlet = dataInlet[5]
HcInlet = dataInlet[6]
HOutlet  = dataOutlet[4]
HsOutlet  = dataOutlet[5]
HcOutlet  = dataOutlet[6]
HRR = dataHRR[2]

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

HRRAve = 0.0
for j in range(nStart,nTime-1):
    HRRAve = HRRAve + 0.5*(HRR[j]+HRR[j+1])*timeStep[j]
HRRAve=HRRAve/aveTime
print 'Mean HRR is ', HRRAve

# inlet
HInletAve = 0.0
HcInletAve = 0.0
HsInletAve = 0.0
for j in range(nStart,nTime-1):
    HInletAve = HInletAve + 0.5*(HInlet[j]+HInlet[j+1])*timeStep[j]
    HsInletAve = HsInletAve + 0.5*(HsInlet[j]+HsInlet[j+1])*timeStep[j]
    HcInletAve = HcInletAve + 0.5*(HcInlet[j]+HcInlet[j+1])*timeStep[j]
HInletAve=-HInletAve/aveTime
HsInletAve=-HsInletAve/aveTime
HcInletAve=-HcInletAve/aveTime
fw.write(str(0)+" "+ str(HsInletAve)+" " + str(HcInletAve) + " " + str(HInletAve) + "\n")

# outlet
HOutletAve = 0.0
HcOutletAve = 0.0
HsOutletAve = 0.0
for j in range(nStart,nTime-1):
    HOutletAve = HOutletAve + 0.5*(HOutlet[j]+HOutlet[j+1])*timeStep[j]
    HsOutletAve = HsOutletAve + 0.5*(HsOutlet[j]+HsOutlet[j+1])*timeStep[j]
    HcOutletAve = HcOutletAve + 0.5*(HcOutlet[j]+HcOutlet[j+1])*timeStep[j]
HOutletAve=HOutletAve/aveTime
HsOutletAve=HsOutletAve/aveTime
HcOutletAve=HcOutletAve/aveTime
fw.write(str(4.0)+" "+ str(HsOutletAve)+" " + str(HcOutletAve) + " " + str(HOutletAve) + "\n")


# faceZones
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
    




