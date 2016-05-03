
import matplotlib as matplotlib
matplotlib.use('Agg')
import astropy
import numpy as np
import matplotlib.pyplot as mpl
from astropy.io import fits
import pylab

#CAUTION - FIELDS REARRANGED

binz=225
boxmin=5.5 #minimum count in bins desired
minel=0.1 #cut off in terms of ellipticity
mincount=0 # minimum of the contour colour cutoff
plotwidth=30 #length of axis
cutoff=0.1
boxmin2=4

#fiddle with B settings, lines more visible, 9/11/15


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#array of input files
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

infile=["DES0025-4206-r.fits","DES0025-4249-r.fits","DES0026-4331-r.fits","DES0029-4206-r.fits","DES0029-4249-r.fits","DES0030-4331-r.fits","DES0030-4414-r.fits","DES0033-4206-r.fits", "DES0033-4249-r.fits","DES0034-4331-r.fits","DES0034-4414-r.fits","DES0037-4206-r.fits","DES0037-4249-r.fits","DES0038-4331-r.fits","DES0038-4414-r.fits","DES0039-4457-r.fits","DES0041-4249-r.fits","DES0042-4331-r.fits","DES0042-4414-r.fits","DES0043-4457-r.fits","DES0045-4331-r.fits","DES0212-0458-r.fits","DES0213-0541-r.fits","DES0215-0416-r.fits","DES0215-0458-r.fits","DES0216-0541-r.fits","DES0217-0333-r.fits","DES0217-0624-r.fits","DES0217-0707-r.fits","DES0218-0416-r.fits","DES0218-0458-r.fits","DES0219-0541-r.fits","DES0219-0624-r.fits","DES0220-0707-r.fits","DES0221-0416-r.fits","DES0221-0458-r.fits","DES0222-0541-r.fits","DES0222-0624-r.fits","DES0223-0333-r.fits","DES0223-0416-r.fits","DES0223-0707-r.fits","DES0224-0458-r.fits","DES0224-0541-r.fits","DES0225-0624-r.fits","DES0226-0333-r.fits","DES0226-0416-r.fits","DES0226-0707-r.fits","DES0227-0458-r.fits","DES0227-0541-r.fits","DES0228-0624-r.fits","DES0229-0333-r.fits","DES0229-0416-r.fits","DES0229-0707-r.fits","DES0230-0458-r.fits","DES0230-0541-r.fits","DES0231-0624-r.fits","DES0232-0333-r.fits","DES0232-0416-r.fits","DES0233-0458-r.fits","DES0233-0541-r.fits","DES0234-0624-r.fits","DES0239-0041-r.fits","DES0239-0124-r.fits","DES0242-0041-r.fits","DES0242-0124-r.fits","DES0242-0207-r.fits","DES0242+0001-r.fits","DES0245-0041-r.fits","DES0245-0124-r.fits","DES0245-0207-r.fits","DES0245+0001-r.fits","DES0248-0041-r.fits","DES0248-0124-r.fits","DES0248-0207-r.fits","DES0248+0001-r.fits","DES0248+0043-r.fits","DES0251-0041-r.fits","DES0251+0001-r.fits","DES0254+0001-r.fits","DES0251+0043-r.fits","DES0254-0041-r.fits","DES0254+0043-r.fits","DES0256+0001-r.fits","DES0325-2749-r.fits","DES0326-2832-r.fits","DES0326-2915-r.fits","DES0328-2706-r.fits","DES0328-2749-r.fits","DES0329-2832-r.fits","DES0329-2915-r.fits","DES0331-2623-r.fits","DES0331-2706-r.fits","DES0332-2749-r.fits","DES0332-2832-r.fits","DES0333-2915-r.fits","DES0333-2958-r.fits","DES0334-2623-r.fits","DES0334-2706-r.fits","DES0339-2832-r.fits","DES0339-2915-r.fits","DES0340-2958-r.fits","DES0341-2623-r.fits","DES0341-2706-r.fits","DES0341-2749-r.fits","DES0342-2832-r.fits","DES0342-2915-r.fits"]


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Reading these files into lists of data points
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RATARGET=33.6713
DECTARGET=-4.5033



print(len(infile))
j=fits.getdata(infile[0])
y=j.field('Dec')
x=np.multiply(j.field('RA'), np.cos(np.radians(y)))
e1=j.field('e1')
e2=j.field('e2')


#RATARGET=np.multiply(RATARGET, np.cos(np.radians(DECTARGET)))



for g in range(1, len(infile)):
	ja=fits.getdata(infile[g])
	yadd=ja.field('Dec')
	xadd=np.multiply(ja.field('RA'), np.cos(np.radians(yadd)))
	e1add=ja.field('e1')
	e2add=ja.field('e2')
	xs=np.append(x,xadd, axis=0)
	#ys=np.append(y,np.multiply(np.cos(yadd), xadd), axis=0)
	ys=np.append(y, yadd, axis=0)
	e1s=np.append(e1, e1add, axis=0)
	e2s=np.append(e2, e2add, axis=0)

	x=xs
	y=ys
	e1=e1s
	e2=e2s
	if g== len(infile)/2:
		print('halfway read in')
	#print(g)


print('Files read in')
if ((len(x))==(len(y))):
	if (len(x)==(len(e1))):
		if (len(x)==(len(e2))):
			if (len(x)==(len(x))):
				print('Array lengths checked: x, y, e1, e2 and gamma are all of equal length of', len(x))
			else:
				print('Error in array lengths, x =/= gamma', len(x), len(gamma))
		else:
			print('Error in array lengths, x =/= e2', len(x), len(e2))
	else:
		print('Error in array lengths, x =/= e1', len(x), len(e1))
else:
	print('Error in array lengths, x =/= y', len(x), len(y))



#def field_division(array)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Dividing data into fields, boundaries of fields input manually
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





xfieldnuma=0
xfieldnumb=0
xfieldnumc=0
xfieldnumd=0


xfielda=[0]*len(x)
yfielda=[0]*len(x)
e1fielda=[0]*len(x)
e2fielda=[0]*len(x)
gammafielda=[0]*len(x)
thetafielda=[0]*len(x)

xfieldb=[0]*len(x)
yfieldb=[0]*len(x)
e1fieldb=[0]*len(x)
e2fieldb=[0]*len(x)
gammafieldb=[0]*len(x)
thetafieldb=[0]*len(x)

xfieldc=[0]*len(x)
yfieldc=[0]*len(x)
e1fieldc=[0]*len(x)
e2fieldc=[0]*len(x)
gammafieldc=[0]*len(x)
thetafieldc=[0]*len(x)

xfieldd=[0]*len(x)
yfieldd=[0]*len(x)
e1fieldd=[0]*len(x)
e2fieldd=[0]*len(x)
gammafieldd=[0]*len(x)
thetafieldd=[0]*len(x)
print('Dividing data into 4 fields')




for o in range(0, len(x)-1):
	if y[o] <= -40:                 #selecting for Field A - ie that with RA below 12
		
		xfieldb[xfieldnumb]=x[o]
		yfieldb[xfieldnumb]=y[o]
		e1fieldb[xfieldnumb]=e1[o]
		e2fieldb[xfieldnumb]=e2[o]
		#gammafieldb[xfieldnumb]=gamma[o]
		#thetafieldb[xfieldnumb]=theta[o]
		xfieldnumb=xfieldnumb+1
	else:
		xfieldnumb=xfieldnumb


	if y[o] > -40 and y[o]<=-20:







		xfieldc[xfieldnumc]=x[o]
		yfieldc[xfieldnumc]=y[o]
		e1fieldc[xfieldnumc]=e1[o]
		e2fieldc[xfieldnumc]=e2[o]
		#gammafieldc[xfieldnumc]=gamma[o]
		#thetafieldc[xfieldnumc]=theta[o]
		xfieldnumc=xfieldnumc+1
	else:
		xfieldnumc=xfieldnumc








		


	if y[o]>=-10 and y[o]< -3:



		xfieldd[xfieldnumd]=x[o]
		yfieldd[xfieldnumd]=y[o]
		e1fieldd[xfieldnumd]=e1[o]
		e2fieldd[xfieldnumd]=e2[o]
		#gammafieldd[xfieldnumd]=gamma[o]
		#thetafieldd[xfieldnumd]=theta[o]
		xfieldnumd=xfieldnumd+1
	else:
		xfieldnumd=xfieldnumd










	if y[o]>=-300 and y[o]< 390:


		xfielda[xfieldnuma]=x[o]
		yfielda[xfieldnuma]=y[o]
		e1fielda[xfieldnuma]=e1[o]
		e2fielda[xfieldnuma]=e2[o]
		xfieldnuma=xfieldnuma+1
	else:
		xfieldnuma=xfieldnuma






























print('Data divided into fields')
print('The number of points in the fields are:')

print(xfieldnuma)
print(xfieldnumb)
print(xfieldnumc)
print(xfieldnumd)
xfieldacut=xfielda[0:xfieldnuma-1]
yfieldacut=yfielda[0:xfieldnuma-1]
e1fieldacut=e1fielda[0:xfieldnuma-1]
e2fieldacut=e2fielda[0:xfieldnuma-1]
#gammafieldacut=gammafielda[0:xfieldnuma-1]
#thetafieldacut=thetafielda[0:xfieldnuma-1]

xfieldbcut=xfieldb[0:xfieldnumb-1]
yfieldbcut=yfieldb[0:xfieldnumb-1]
e1fieldbcut=e1fieldb[0:xfieldnumb-1]
e2fieldbcut=e2fieldb[0:xfieldnumb-1]
#gammafieldbcut=gammafieldb[0:xfieldnumb-1]
#thetafieldbcut=thetafieldb[0:xfieldnumb-1]

xfieldccut=xfieldc[0:xfieldnumc-1]
yfieldccut=yfieldc[0:xfieldnumc-1]
e1fieldccut=e1fieldc[0:xfieldnumc-1]
e2fieldccut=e2fieldc[0:xfieldnumc-1]
#gammafieldccut=gammafieldc[0:xfieldnumc-1]
#thetafieldccut=thetafieldc[0:xfieldnumc-1]

xfielddcut=xfieldd[0:xfieldnumd-1]
yfielddcut=yfieldd[0:xfieldnumd-1]
e1fielddcut=e1fieldd[0:xfieldnumd-1]
e2fielddcut=e2fieldd[0:xfieldnumd-1]
#gammafielddcut=gammafieldd[0:xfieldnumd-1]
#thetafielddcut=thetafieldd[0:xfieldnumd-1]










#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Producing the arrays for use in contour/quiver plots
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Field A calculations
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



print('Number in field', xfieldnuma)
#binzb=binz
binzb1=np.round_((np.sqrt(np.divide(xfieldnuma,boxmin)))/25,0) #calculating the number of bins required to try and get boxmin
ki=np.divide(xfieldnuma, boxmin)
kj=np.sqrt(ki)
kj25=(kj)/plotwidth
kj25r=np.round_(kj25+0.5, 0)
binz=int(kj25r*plotwidth)
#print(binzb)
int(binz)
print('bins in A', binz)



print(len(xfieldacut))

bigx=np.amax(xfieldacut)
minix=np.amin(xfieldacut)
stepx=(bigx-minix)/binz


binzxst=range(binz)
binzxac=binzxst*stepx+minix

bigy=np.amax(yfieldacut)
miniy=np.amin(yfieldacut)
stepy=(bigy-miniy)/binz
binzyst=range(binz)
binzyac=binzyst*stepy+miniy

Avals=np.zeros((binz*binz))

stepxa=stepx
minixa=minix
stepya=stepy
miniya=miniy   #defining for adding to graph axes


yedgea=binzyac
xedgea=binzxac


xlocs=np.digitize(xfieldacut, binzxac)
ylocs=np.digitize(yfieldacut, binzyac)
densmapa=np.zeros((binz+1, binz+1))
e1fielda=np.zeros((binz+1, binz+1))
e2fielda=np.zeros((binz+1, binz+1))


reps=len(xfieldacut)-1

for df in range(0, reps):

	xval=xlocs[df]
	yval=ylocs[df]
	densmapa[xval, yval]=densmapa[xval, yval]+1
	e1fielda[xval, yval]=e1fielda[xval, yval]+e1fieldacut[df]
	e2fielda[xval, yval]=e2fielda[xval, yval]+e2fieldacut[df]
	
	
e1ANorm=np.zeros((binz+1, binz+1))
e2ANorm=np.zeros((binz+1, binz+1))

j=0

for d in range(0, binz+1):
	for f in range(0, binz+1):
		if boxmin <= densmapa[d,f]:
			e1ANorm[d,f]=e1fielda[d,f]/densmapa[d,f]
			e2ANorm[d,f]=e2fielda[d,f]/densmapa[d,f]
			Avals[j]=densmapa[d,f]
			j=j+1
			
		else:
			e1ANorm[d,f]=0
			e2ANorm[d,f]=0
			densmapa[d,f]=0
			


eA=np.sqrt(np.square(e1ANorm)+np.square(e2ANorm))

tot=np.sum(densmapa)
num=np.count_nonzero(densmapa)
typical=tot/num
print('Mean of nonzero bin counts in Field A', typical)

ATcut=Avals[0:j]
print('Median is hopefully', np.median(ATcut))


thetaA=np.zeros((binz+1, binz+1))

for d in range(0, binz+1):
	for f in range(0, binz+1):
		if e1ANorm[d,f] != 0:
			thetaA[d,f]=np.arctan(np.divide(e2ANorm[d,f], e1ANorm[d,f]))
	


FieldAXDir=np.multiply(eA, np.cos(thetaA))
FieldAYDir=np.multiply(eA, np.sin(thetaA))

xstretchA=np.multiply(eA, np.cos(thetaA/2))
ystretchA=np.multiply(eA, np.sin(thetaA/2))



#meda=np.nanmedian(non0a)
#print('Median bin count in Field A:', meda)

RAbin=np.divide((RATARGET-minix), stepx)
decbin=np.divide(DECTARGET-miniy, stepy)
print('RA (X) bin in A', RAbin)
print('DEC (Y) bin in A', decbin)





#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Field B Calculations
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


print('Number in field', xfieldnumb)
#binzb=binz
binzb1=np.round_((np.sqrt(np.divide(xfieldnumb,boxmin)))/25,0) #calculating the number of bins required to try and get boxmin
ki=np.divide(xfieldnumb, boxmin)
kj=np.sqrt(ki)
kj25=(kj)/plotwidth
kj25r=np.round_(kj25+0.5, 0)
binzb=int(kj25r*plotwidth)
#print(binzb)
int(binzb)
print('bins in B', binzb)


bigx=np.amax(xfieldbcut)
minix=np.amin(xfieldbcut)
stepx=(bigx-minix)/binzb
binzbxst=range(binzb)
binzbxbc=binzbxst*stepx+minix


stepxb=stepx
minixb=minix


bigy=np.amax(yfieldbcut)
miniy=np.amin(yfieldbcut)
stepy=(bigy-miniy)/binzb
binzbyst=range(binzb)
binzbybc=binzbyst*stepy+miniy



stepxb=stepx
minixb=minix
stepyb=stepy
miniyb=miniy

yedgeb=binzbybc
xedgeb=binzbxbc


#print(binzbxbc)
Bvals=np.zeros((binzb*binzb))


xlocs=np.digitize(xfieldbcut, binzbxbc)
ylocs=np.digitize(yfieldbcut, binzbybc)
densmapb=np.zeros((binzb+1, binzb+1))

e1fieldb=np.zeros((binzb+1, binzb+1))
e2fieldb=np.zeros((binzb+1, binzb+1))

reps=len(xfieldbcut)-1

for df in range(0, reps):

	xval=xlocs[df]
	yval=ylocs[df]
	densmapb[xval, yval]=densmapb[xval, yval]+1
	e1fieldb[xval, yval]=e1fieldb[xval, yval]+e1fieldbcut[df]
	e2fieldb[xval, yval]=e2fieldb[xval, yval]+e2fieldbcut[df]
 

e1BNorm=np.zeros((binzb+1, binzb+1))
e2BNorm=np.zeros((binzb+1, binzb+1))
eBNorm=np.zeros((binzb+1, binzb+1))
d=0
f=0
j=0


for d in range(0, binzb+1):
	for f in range(0, binzb+1):
		if boxmin <= densmapb[d,f]:
			e1BNorm[d,f]=e1fieldb[d,f]/densmapb[d,f]
			e2BNorm[d,f]=e2fieldb[d,f]/densmapb[d,f]
			#eBNorm[d,f]=(e1BNorm[d,f]+e2BNorm)/2
			Bvals[j]=densmapb[d,f]
			j=j+1
		else:
			e1BNorm[d,f]=0
			e2BNorm[d,f]=0
			eBNorm[d,f]=0
			densmapb[d,f]=0


eB=np.sqrt(np.square(e1BNorm)+np.square(e2BNorm))
thetaB=np.zeros((binzb+1, binzb+1))

for d in range(0, binzb+1):
	for f in range(0, binzb+1):
		if e1BNorm[d,f] != 0:
			thetaB[d,f]=np.arctan(np.divide(e1BNorm[d,f], e2BNorm[d,f]))
	


FieldBXDir=np.multiply(eB, np.cos(thetaB))
FieldBYDir=np.multiply(eB, np.sin(thetaB))


#print(eB)

nmsum=np.sum(eB, axis=0)
bnonzero=np.count_nonzero(eB)

#print(nmsum)
#print('nonz',bnonzero)


RAbin=np.divide((RATARGET-minix), stepx)
decbin=np.divide(DECTARGET-miniy, stepy)
print('RA (X) bin in B', RAbin)
print('DEC (Y) bin in B', decbin)



tot=np.sum(densmapb)
num=np.count_nonzero(densmapb)
typical=tot/num
print('Mean of nonzero bin counts in Field B', typical)

BTcut=Bvals[0:j]
print('Median is hopefully', np.median(BTcut))

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Field C Calculations
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




print('Number in field', xfieldnumc)
#binzb=binz
binzb1=np.round_((np.sqrt(np.divide(xfieldnumc,boxmin)))/25,0) #calculating the number of bins required to try and get boxmin
ki=np.divide(xfieldnumc, boxmin)
kj=np.sqrt(ki)
kj25=(kj)/plotwidth
kj25r=np.round_(kj25+0.5, 0)
binzc=int(kj25r*plotwidth)
#print(binzb)
int(binzc)
print('bins in C', binzc)




#binzc=binz

bigx=np.amax(xfieldccut)
minix=np.amin(xfieldccut)
stepx=(bigx-minix)/binzc
binzcxst=range(binzc) 
binzcxbc=binzcxst*stepx+minix

bigy=np.amax(yfieldccut)
miniy=np.amin(yfieldccut)
stepy=(bigy-miniy)/binzc
binzcyst=range(binzc)
binzcybc=binzcyst*stepy+miniy

#print(binzcxbc) 

xlocs=np.digitize(xfieldccut, binzcxbc)
ylocs=np.digitize(yfieldccut, binzcybc)
densmapc=np.zeros((binzc+1, binzc+1))
fieldcgammamap=np.zeros((binzc+1, binzc+1))
fieldcthetamap=np.zeros((binzc+1, binzc+1))


yedgec=binzcybc
xedgec=binzcxbc


stepxc=stepx
minixc=minix
stepyc=stepy
miniyc=miniy


Cvals=np.zeros((binzc*binzc))


densmapc=np.zeros((binzc+1, binzc+1))

e1fieldc=np.zeros((binzc+1, binzc+1))
e2fieldc=np.zeros((binzc+1, binzc+1))

reps=len(xfieldccut)-1

for df in range(0, reps):

	xval=xlocs[df]
	yval=ylocs[df]
	densmapc[xval, yval]=densmapc[xval, yval]+1
	e1fieldc[xval, yval]=e1fieldc[xval, yval]+e1fieldccut[df]
	e2fieldc[xval, yval]=e2fieldc[xval, yval]+e2fieldccut[df]
 

e1CNorm=np.zeros((binzc+1, binzc+1))
e2CNorm=np.zeros((binzc+1, binzc+1))
eCNorm=np.zeros((binzc+1, binzc+1))
d=0
f=0
j=0


for d in range(0, binzc+1):
	for f in range(0, binzc+1):
		if boxmin <= densmapc[d,f]:
			e1CNorm[d,f]=e1fieldc[d,f]/densmapc[d,f]
			e2CNorm[d,f]=e2fieldc[d,f]/densmapc[d,f]
			#eCNorm[d,f]=(e1CNorm[d,f]+e2CNorm)/2
			Cvals[j]=densmapc[d,f]
			j=j+1
		else:
			e1CNorm[d,f]=0
			e2CNorm[d,f]=0
			densmapc[d,f]=0


#print(e2CNorm)
#print(e1CNorm)
eC=np.sqrt(np.square(e1CNorm)+np.square(e2CNorm))
thetaC=np.zeros((binzc+1, binzc+1))

for d in range(0, binzc+1):
	for f in range(0, binzc+1):
		if e1CNorm[d,f] != 0:
			thetaC[d,f]=np.arctan(np.divide(e1CNorm[d,f], e2CNorm[d,f]))
	


FieldCXDir=np.multiply(eC, np.cos(thetaC))
FieldCYDir=np.multiply(eC, np.sin(thetaC))



RAbin=np.divide((RATARGET-minix), stepx)
decbin=np.divide(DECTARGET-miniy, stepy)
print('RA (X) bin in C', RAbin)
print('DEC (Y) bin in C', decbin)


tot=np.sum(densmapc)
num=np.count_nonzero(densmapc)
typical=tot/num
print('Mean of nonzero bin counts in Field C', typical)



CTcut=Cvals[0:j]
print('Median is hopefully', np.median(CTcut))



#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Field D Calculations
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




print('Number in field', xfieldnumd)
#binzb=binz
binzb1=np.round_((np.sqrt(np.divide(xfieldnumb,boxmin)))/25,0) #calculating the number of bins required to try and get boxmin
ki=np.divide(xfieldnumd, boxmin)
kj=np.sqrt(ki)
kj25=(kj)/plotwidth
kj25r=np.round_(kj25+0.5, 0)
binzd=int(kj25r*plotwidth)
#print(binzb)
int(binzd)
print('bins in D', binzd)




binzd=binz

bigx=np.amax(xfielddcut)
minix=np.amin(xfielddcut)
stepx=(bigx-minix)/binzd
binzdxst=range(binzd)
binzdxbc=binzdxst*stepx+minix

bigy=np.amax(yfielddcut)
miniy=np.amin(yfielddcut)
stepy=(bigy-miniy)/binzd
binzdyst=range(binzd)
binzdybc=binzdyst*stepy+miniy

yedged=binzdybc
xedged=binzdxbc

stepxd=stepx
minixd=minix
stepyd=stepy
miniyd=miniy


xlocs=np.digitize(xfielddcut, binzdxbc)
ylocs=np.digitize(yfielddcut, binzdybc)
densmapd=np.zeros((binzd+1, binzd+1))
fielddgammamap=np.zeros((binzd+1, binzd+1))
fielddthetamap=np.zeros((binzd+1, binzd+1))

densmapd=np.zeros((binzd+1, binzd+1))

e1fieldd=np.zeros((binzd+1, binzd+1))
e2fieldd=np.zeros((binzd+1, binzd+1))

reps=len(xfielddcut)-1

for df in range(0, reps):

	xval=xlocs[df]
	yval=ylocs[df]
	densmapd[xval, yval]=densmapd[xval, yval]+1
	e1fieldd[xval, yval]=e1fieldd[xval, yval]+e1fielddcut[df]
	e2fieldd[xval, yval]=e2fieldd[xval, yval]+e2fielddcut[df]
 

e1DNorm=np.zeros((binzd+1, binzd+1))
e2DNorm=np.zeros((binzd+1, binzd+1))
eDNorm=np.zeros((binzd+1, binzd+1))
d=0
f=0
j=0
Dvals=np.zeros((binzd*binzd))

for d in range(0, binzd+1):
	for f in range(0, binzd+1):
		if boxmin <= densmapd[d,f]:
			e1DNorm[d,f]=e1fieldd[d,f]/densmapd[d,f]
			e2DNorm[d,f]=e2fieldd[d,f]/densmapd[d,f]
			#eDNorm[d,f]=(e1DNorm[d,f]+e2DNorm)/2
			Dvals[j]=densmapd[d,f]
			j=j+1
		else:
			e1DNorm[d,f]=0
			e2DNorm[d,f]=0
			densmapd[d,f]=0



eD=np.sqrt(np.square(e1DNorm)+np.square(e2DNorm))
thetaD=np.zeros((binzd+1, binzd+1))

for d in range(0, binzd+1):
	for f in range(0, binzd+1):
		if e1DNorm[d,f] != 0:
			thetaD[d,f]=np.arctan(np.divide(e1DNorm[d,f], e2DNorm[d,f]))
	


FieldDXDir=np.multiply(eD, np.cos(thetaD))
FieldDYDir=np.multiply(eD, np.sin(thetaD))

RAbin=np.divide((RATARGET-minix), stepx)
decbin=np.divide(DECTARGET-miniy, stepy)


#xloctarget=np.digitize(	RATARGET, binzdxbc)
#yloctarget=np.digitize( DECTARGET, binzdybc)

print('RA (X) bin in D', RAbin)
print('DEC (Y) bin in D', decbin)

#print(xloctarget)
#print(yloctarget)


DTcut=Dvals[0:j]
print('Median is hopefully', np.median(DTcut))


tot=np.sum(densmapd)
num=np.count_nonzero(densmapd)
typical=tot/num
print('Mean of nonzero bin counts in Field D', typical)





#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Producing mass maps
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Field A
Rx=FieldAXDir
Ry=FieldAYDir


#Take a field of shear vectors in x and y and convert that into kappa and psi maps



#Define Rx as the x vectors of shear, Ry as the y vectors

Rxshift=np.fft.fftshift(Rx)
Ryshift=np.fft.fftshift(Ry)   #shifts vectors into edges for FFT

Rxsft=np.fft.fft2(Rxshift)
Rysft=np.fft.fft2(Ryshift)     # FFT of the vectors

gamma1t=np.fft.ifftshift(Rxsft)
gamma2t=np.fft.ifftshift(Rysft)   # converting FFT into gammas, by inverse shifting

gammat=gamma1t+np.multiply(gamma2t, 1j)  # gamma tilda, by definition

xsize=gammat.shape[1]
ysize=gammat.shape[0] #required length of axes


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
psit=np.zeros((ysize, xsize), dtype=np.complex)
kappat=np.zeros((ysize, xsize), dtype=np.complex) #define arrays to be filled by for loop, sizes of dimensions

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

i=0

xsizeb=int(xsize)
ysizeb=int(ysize)

if xsizeb%2:
	xcorr=0
else:
	xcorr=0.5

if ysizeb%2:
	ycorr=0
else:
	ycorr=0.5

for l in range(0,xsize):
	for i in range(0, ysize):
		xi1=(i-xsize/2)+xcorr
		xi2=(l-ysize/2)+ycorr
		if xi1 == 0 and xi2 ==0: 
			psit[l,i]=0 # otherwise is infinite
		
					
		else:
			psit[l,i]=np.divide(gammat[l,i]*(xi1**2-xi2**2-2j*xi1*xi2),2*(np.pi**2)*(xi1**2+xi2**2)**2) # calculating psi tilda, according to eqns [15/12/15]
		
		if xi1!=0 or xi2!=0:
			kappat[l,i]=np.divide((xi1**2-xi2**2-2j*(xi1*xi2))*gammat[l,i], (xi1**2+xi2**2))   #calculate kappa tilda according to eqns [15/12/15]
	
		else:
			kappat[l,i]=0 #is infinite at xi1=0, xi2=0  

	
 


psiprep=np.fft.fftshift(psit) #shifting psit for inverse FFT
kappaprep=np.fft.fftshift(kappat) #shifting kappat for inverse FFT

psireal=np.fft.ifft2(psiprep)  #inverse FFT of psit
kappareal=np.fft.ifft2(kappaprep)



psi=np.fft.ifftshift(psireal) # inverse shift of psi, should be done!
kappa=np.fft.ifftshift(kappareal) # shifting kappa back, should be it!

mpl.contourf(np.real(kappa), vmin=0.1)
mpl.colorbar()
pylab.savefig('DGraphs/MassMaps/NoNoise/KappaMap 010 boxmin4.png')
mpl.close()

mpl.contourf(np.real(kappa), vmin=0.2)
mpl.colorbar()
pylab.savefig('DGraphs/MassMaps/NoNoise/KappaMap 020 boxmin4.png')
mpl.close()

mpl.contourf(np.real(kappa), vmin=0.05)
mpl.colorbar()
pylab.savefig('DGraphs/MassMaps/NoNoise/KappaMap 005 boxmin4.png')
mpl.close()

mpl.contourf(np.real(kappa), vmin=0.15)
mpl.colorbar()
pylab.savefig('DGraphs/MassMaps/NoNoise/KappaMap 015 boxmin4.png')
mpl.close()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Graph commands
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# experimental plots
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



decplaces=4

width=25

# move to top eventually
#also need to alter binzb,c,d etc to accodmmodate for varying field size = 225 only works for one case.kinellm8

for xlimi in range(width, binz+width, width):
	for ylimi in range(width,binz+width, width):
		suf=0
		maggam=(np.square(FieldAXDir[xlimi-width:xlimi,ylimi-width:ylimi]), np.square(FieldAYDir[xlimi-width:xlimi,ylimi-width:ylimi]))#+np.square(FieldAYDir))
		cutoffm=np.amax(maggam)#np.median(maggam)
		if cutoff<= cutoffm:
			suf='Bright'
			vm=boxmin2 
		else:
			suf='Quiet'
			vm=0
		print(cutoffm)
		mpl.contourf(densmapa, cmap=mpl.cm.gray, vmin=vm)
		mpl.quiver(FieldAXDir, FieldAYDir, headwidth=0, headaxislength=50, pivot='middle', color='y')
		
		
		if cutoff<= cutoffm:
			mpl.ylim((ylimi-width, ylimi))
			ysmall=np.round_(miniya+((ylimi-width)*stepya), decplaces)
			ybig=np.round_(minixa+(ylimi*stepya), decplaces)
			xbig=np.round_(minixa+((xlimi-width)*stepxa), decplaces)
			xsmall=np.round_(minixa+(xlimi*stepxa), decplaces)
			minixar=np.round_(minixa, decplaces)
			miniyar=np.round_(miniya, decplaces)
			stepxar=np.round_(stepxa, decplaces)
			stepyar=np.round_(stepya, decplaces)
			mpl.xlabel('Xbin, becomes RAcosdec by '+str(minixar)+' plus bin by '+str(stepxar)+'\n [this case between range of RA='+str(xsmall)+'-'+str(xbig)+']')
			mpl.ylabel('Ybin, ie Dec by '+str(miniyar)+'+ bins x '+str(stepyar)+'')
			mpl.xlim((xlimi-width, xlimi))
			mpl.title('Field A section')
			pylab.savefig('DGraphs/FieldASectionat X'+str(xlimi)+' and Y'+str(ylimi)+'f'+str(suf)+'.png') #
			print('Produced Quiverf, Field D, xlim= '+str(xlimi)+' and ylim='+str(ylimi)+'')




		mpl.contourf(np.real(kappa), cmap=mpl.cm.gray, vmin=0.1)
		mpl.colorbar()
		mpl.quiver(xstretchA, ystretchA, headwidth=0, headaxislength=100, pivot='middle', color='c')
		mpl.xlim((xlimi-width, xlimi))
		mpl.ylim((ylimi-width, ylimi))
		pylab.savefig('DGraphs/MassMaps/NoNoise/'+str(xlimi)+' '+str(ylimi)+' Shear on Kappa Cyan 010 cut boxmin4 .png')
		mpl.close()
 




#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Graphs with noise
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


AX=np.zeros(((binz+1)**2))
AY=np.zeros(((binz+1)**2))
j=0
for d in range(0,binz+1):
	for f in range(0, binz+1):
		AX[j]=FieldAXDir[d,f]
		AY[j]=FieldAYDir[d,f]
		j=j+1
print(np.var(AX))
print(np.var(AY))

AXvar=np.sqrt(np.absolute(np.var(AX, axis=0)))
AYvar=np.sqrt(np.absolute(np.var(AY, axis=0)))
good=0
bad=0

for d in range(0, binz+1):
	for f in range(0, binz+1):
		if FieldAXDir[d,f]==0:
			FieldAXDir[d,f]=np.random.normal(0.0, scale=AXvar)
		if FieldAYDir[d,f]==0:
			FieldAYDir[d,f]=np.random.normal(0.0, scale=AYvar)
		if ((np.real(kappa[d,f]))/(np.multiply((np.imag(kappa[d,f])),np.conj(np.imag(kappa[d,f])))))>10:
			good=good+1
		else:
			bad=bad+1





print('The good number is', good)
print('The bad number is', bad)
print('The kappa matrix', kappa)


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Producing mass maps
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Field A
Rx=FieldAXDir
Ry=FieldAYDir


#Take a field of shear vectors in x and y and convert that into kappa and psi maps



#Define Rx as the x vectors of shear, Ry as the y vectors

Rxshift=np.fft.fftshift(Rx)
Ryshift=np.fft.fftshift(Ry)   #shifts vectors into edges for FFT

Rxsft=np.fft.fft2(Rxshift)
Rysft=np.fft.fft2(Ryshift)     # FFT of the vectors

gamma1t=np.fft.ifftshift(Rxsft)
gamma2t=np.fft.ifftshift(Rysft)   # converting FFT into gammas, by inverse shifting

gammat=gamma1t+np.multiply(gamma2t, 1j)  # gamma tilda, by definition

xsize=gammat.shape[1]
ysize=gammat.shape[0] #required length of axes


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
psit=np.zeros((ysize, xsize), dtype=np.complex)
kappat=np.zeros((ysize, xsize), dtype=np.complex) #define arrays to be filled by for loop, sizes of dimensions

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

i=0

xsizeb=int(xsize)
ysizeb=int(ysize)

if xsizeb%2:
	xcorr=0
else:
	xcorr=0.5

if ysizeb%2:
	ycorr=0
else:
	ycorr=0.5

for l in range(0,xsize):
	for i in range(0, ysize):
		xi1=(i-xsize/2)+xcorr
		xi2=(l-ysize/2)+ycorr
		if xi1 == 0 and xi2 ==0: 
			psit[l,i]=0 # otherwise is infinite
		
					
		else:
			psit[l,i]=np.divide(gammat[l,i]*(xi1**2-xi2**2-2j*xi1*xi2),2*(np.pi**2)*(xi1**2+xi2**2)**2) # calculating psi tilda, according to eqns [15/12/15]
		
		if xi1!=0 or xi2!=0:
			kappat[l,i]=np.divide((xi1**2-xi2**2-2j*(xi1*xi2))*gammat[l,i], (xi1**2+xi2**2))   #calculate kappa tilda according to eqns [15/12/15]
	
		else:
			kappat[l,i]=0 #is infinite at xi1=0, xi2=0  

	
 


psiprep=np.fft.fftshift(psit) #shifting psit for inverse FFT
kappaprep=np.fft.fftshift(kappat) #shifting kappat for inverse FFT

psireal=np.fft.ifft2(psiprep)  #inverse FFT of psit
kappareal=np.fft.ifft2(kappaprep)



psi=np.fft.ifftshift(psireal) # inverse shift of psi, should be done!
kappa=np.fft.ifftshift(kappareal) # shifting kappa back, should be it!

mpl.contourf(np.real(kappa), vmin=0)
mpl.colorbar()
pylab.savefig('DGraphs/MassMaps/Noisy/KappaMap A noise bm4.png')
mpl.close()

mpl.imshow(np.real(kappa))
pylab.savefig('DGraphs/MassMaps/Noisy/imshow kappa A with noise bm4.png')
mpl.close()

mpl.contourf(np.real(kappa), vmin=0.1)
mpl.colorbar()
pylab.savefig('DGraphs/MassMaps/Noisy/KappaMap 010 bm4.png')
mpl.close()

mpl.contourf(np.real(kappa), vmin=0.2)
mpl.colorbar()
pylab.savefig('DGraphs/MassMaps/Noisy/KappaMap 020 bm4.png')
mpl.close()

mpl.contourf(np.real(kappa), vmin=0.05)
mpl.colorbar()
pylab.savefig('DGraphs/MassMaps/Noisy/KappaMap 005.png')
mpl.close()

mpl.contourf(np.real(kappa), vmin=0.15)
mpl.colorbar()
pylab.savefig('DGraphs/MassMaps/Noisy/KappaMap 015 bm4.png')
mpl.close()


for xlimi in range(width, binz+width, width):
	for ylimi in range(width,binz+width, width):
		suf=0
		maggam=(np.square(FieldAXDir[xlimi-width:xlimi,ylimi-width:ylimi]), np.square(FieldAYDir[xlimi-width:xlimi,ylimi-width:ylimi]))#+np.square(FieldAYDir))
		cutoffm=np.amax(maggam)#np.median(maggam)
		if cutoff<= cutoffm:
			suf='Bright'
			vm=boxmin2 
		else:
			suf='Quiet'
			vm=0
		print(cutoffm)
		mpl.contourf(densmapa, cmap=mpl.cm.gray, vmin=vm)
		mpl.quiver(FieldAXDir, FieldAYDir, headwidth=0, headaxislength=50, pivot='middle', color='y')
		
		
		if cutoff<= cutoffm:
			mpl.ylim((ylimi-width, ylimi))
			ysmall=np.round_(miniya+((ylimi-width)*stepya), decplaces)
			ybig=np.round_(minixa+(ylimi*stepya), decplaces)
			xbig=np.round_(minixa+((xlimi-width)*stepxa), decplaces)
			xsmall=np.round_(minixa+(xlimi*stepxa), decplaces)
			minixar=np.round_(minixa, decplaces)
			miniyar=np.round_(miniya, decplaces)
			stepxar=np.round_(stepxa, decplaces)
			stepyar=np.round_(stepya, decplaces)
			mpl.xlabel('Xbin, becomes RAcosdec by '+str(minixar)+' plus bin by '+str(stepxar)+'\n [this case between range of RA='+str(xsmall)+'-'+str(xbig)+']')
			mpl.ylabel('Ybin, ie Dec by '+str(miniyar)+'+ bins x '+str(stepyar)+'')
			mpl.xlim((xlimi-width, xlimi))
			mpl.title('Field A section')
			pylab.savefig('DGraphs/MassMaps/Noisy/FieldASectionat X'+str(xlimi)+' and Y'+str(ylimi)+'f'+str(suf)+'bm4 .png') #
			print('Produced Quiverf, Field D, xlim= '+str(xlimi)+' and ylim='+str(ylimi)+'')




		mpl.contourf(np.real(kappa), cmap=mpl.cm.gray, vmin=0.1)
		mpl.colorbar()
		mpl.quiver(xstretchA, ystretchA, headwidth=0, headaxislength=100, pivot='middle', color='c')
		mpl.xlim((xlimi-width, xlimi))
		mpl.ylim((ylimi-width, ylimi))
		pylab.savefig('DGraphs/MassMaps/Noisy/'+str(xlimi)+' '+str(ylimi)+' Shear on Kappa Cyan 010 cut bm4 .png')
		mpl.close()


print('The kappa matrix', kappa)
