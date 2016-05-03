

import matplotlib as matplotlib
matplotlib.use('Agg')
import astropy
import numpy as np
import matplotlib.pyplot as mpl
from astropy.io import fits
import pylab
import healpy as hp
#import scipy.stats as sp
#import Benlib as ben


#========================================================================================================================================================================================================
# 							~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#========================================================================================================================================================================================================



def make_daughter(binzx, binzy, mother, gausswidth, gausswidthi):
	mother=np.asanyarray(mother, dtype=np.complex)
	daughter=mother.copy()
	daughterxorigin=(mother.shape[1]/2)
	daughteryorigin=(mother.shape[0]/2)
	meandiff=0
	for d in range(int(daughteryorigin-binzy), int(daughteryorigin+binzy)):
		for f in range(int(daughterxorigin-binzx),int(daughterxorigin+binzx)):
			alt=(np.random.normal(0.0, scale=gausswidth)+(np.random.normal(0.0, scale=gausswidthi))*1j)
			daughter[d,f]=daughter[d,f]+alt

	return daughter



def chi_squared(true, test, error_matrix):
	chi=0
	true1=np.asanyarray(true)
	true1=list((true1.copy()).reshape(-1,))
	test1=list(np.asanyarray(test.copy()).reshape(-1,))
	print(len(true1))
	print(len(test1))

	error_matrix1=list(np.asanyarray(error_matrix.copy()).reshape(-1,))
	print(len(error_matrix1))
	for d in range(0, len(true1)-1):
			dd=int(d)
		#	print(dd)
			chi=chi+((np.absolute(true1[dd]-test1[dd])**2)/(error_matrix1[dd]))
	#float(chi)
	#chi=chi/((true.shape[0])*(true.shape[1]))
	return chi 


def gamma_build(gammatilda):

#input is a matrix of some gamma tilda values, code will then reconstruct the observed shears for this map 

	shifted=np.fft.fftshift(gammatilda)
	reverset=np.fft.ifft2(shifted)
	out=np.fft.ifftshift(reverset) #real space gamma values

	return out

def gammat_build(gamma):
	shifted=np.fft.fftshift(gamma)
	transformed=np.fft.fft2(shifted)
	out=np.fft.ifftshift(transformed)

	return out


def best_daughter(mother, trueReal, error, xwidth, ywidth, loops, gausswidth, gausswidthi): 

	chilist=[0]*loops
	xorigin=trueReal.shape[1]/2
	yorigin=trueReal.shape[0]/2
	bestdaughter=mother.copy()

	daughter1=make_daughter(xwidth, ywidth, mother.copy(), gausswidth, gausswidthi)
	daughter1real=gamma_build(daughter1)
	bestchi=float(chi_squared(trueReal, daughter1real, error))

	
	
	for g in range(0, int(loops-1)):
		daughter2=make_daughter(xwidth, ywidth, mother, gausswidth, gausswidthi)
		daughter2[yorigin, xorigin]=0
		daughter2real=gamma_build(daughter2)

		newchi=float(chi_squared(trueReal, daughter2real, error))
		chilist[g]=newchi
		if newchi<bestchi:
			bestchi=newchi
			bestdaughter=daughter2


	return bestdaughter, bestchi

def real_gamma_difference(realtrue, testtilda, xwidth, ywidth):


	
	realtrue=np.asanyarray(realtrue)
	testilda=np.asanyarray(testtilda)

	xorigin=realtrue.shape[1]/2
	yorigin=realtrue.shape[0]/2

	diff=testtilda[yorigin-ywidth:yorigin+ywidth,xorigin-xwidth:xorigin+xwidth]-realtrue[yorigin-ywidth:yorigin+ywidth,xorigin-xwidth:xorigin+xwidth]
	diffline=np.real(np.reshape(diff, (1,-1)))
	magdiff=np.average(np.absolute(diffline))
	difflinei=np.imag(np.reshape(diff, (1,-1)))
	magdiffi=np.average(np.absolute(difflinei))

	return magdiff, magdiffi

def tilda_gamma_difference(truetilda, testtilda, xwidth, ywidth):


	
	realtrue=np.asanyarray(realtrue)
	testilda=np.asanyarray(testtilda)

	xorigin=realtrue.shape[1]/2
	yorigin=realtrue.shape[0]/2

	diff=testtilda[yorigin-ywidth:yorigin+ywidth,xorigin-xwidth:xorigin+xwidth]-realtrue[yorigin-ywidth:yorigin+ywidth,xorigin-xwidth:xorigin+xwidth]
	diffline=np.reshape(diff, (1,-1))
	magdiff=np.average(np.absolute(diffline))

	return magdiff


def chi_shuffle(true, mother, chideal, error, xwidth, ywidth):

#mother is the initial input matrix (hyp gammat), true being the target matrix (real gamma) which you want to minimise the difference in chi2 from. 
#xwidth is the x domain within which values are being varied, ywidth similarly. Defines a gausswidth itself although don't know a good reason to pick one yet.



	true=np.asanyarray(true)
	mother=np.asanyarray(mother)
	error=np.asanyarray(error)
	bestgt=mother


	xorigin=true.shape[1]/2
	yorigin=true.shape[0]/2

	motherr=gamma_build(mother.copy())


	kids=20
	bc=chideal+1
	generations=0
	gausswidth=20
	kl=[0]
	#cg=np.array(0)
	k=0

	truet=gammat_build(true)
	gausswidth=1*np.sqrt(real_gamma_difference(mother, truet, xwidth, ywidth)[0])
	gausswidthi=1*np.sqrt(real_gamma_difference(mother, truet, xwidth, ywidth)[1])
	bc1=chideal+1
	l=0
	chim=float(chi_squared(true, motherr, error))
	print('#C2beat', chim)
	print('initial gauss', gausswidth)
	cg=[chim]
	chix=chim
	jh=0
	motherest=mother.copy()
	z=0
	gooduns=[z]
	cgood=[chim]




	while chideal < bc1:
		bet=best_daughter(motherest, true, error, xwidth, ywidth, kids, gausswidth, gausswidthi)
		bd=bet[0]
		bc=bet[1]
		#print(bc, gausswidth, gausswidthi)
		#print('Mother generation', generations, mother[yorigin-ywidth:yorigin+ywidth,xorigin-xwidth:xorigin+xwidth])
		cg=np.append(cg, bc)
		kl=np.append(kl, k)
		mpl.scatter(kl,cg )
		mpl.xlabel('Iterations')
		mpl.ylabel('Chi2')
		mpl.savefig('cgk2 '+str(xwidth)+' by '+str(ywidth)+'.png')
		mpl.close()
		k=k+1
		chik=bc
		chidiff=chik-chix

		if chidiff <0.01:
			l=l+1
			#print(l)
		if chidiff >= 0.01:
			l=0
		if l>5000:
			bc1=chideal
		if bc < chim:
			print('fitter than mum, old was', chim, 'new is', bc)
			chim=bc
			generations=generations+1
			cgood=np.append(cgood, bc)
			gooduns=np.append(gooduns, generations)
			motherest=bd
			gausswidth=gausswidth*0.9995
			gausswidthi=gausswidthi*0.9995
			#gausswidth=np.sqrt(real_gamma_difference(bd, truet, xwidth, ywidth)[0])
			#gausswidthi=np.sqrt(real_gamma_difference(bd, truet, xwidth, ywidth)[1])

			latestmass=map_mass(gamma_build(motherest.copy()))
			mpl.contourf(np.real(latestmass))
			mpl.colorbar()
			mpl.title('Forward Fitting attempt at modelling mass, Chi2='+str(chim/dof)+'.png')
			mpl.savefig('FFmass '+str(xwidth)+' by '+str(ywidth)+'.png')
			mpl.close()
			chix=bc
			bc1=bc

			
			mpl.scatter(gooduns, cgood)
			mpl.title('Progression of best Chi estimate')
			mpl.xlabel('Generations')
			mpl.ylabel('Chi2')
			mpl.savefig('Best Chi evolutions '+str(xwidth)+' by '+str(ywidth)+' .png')
			mpl.close()
		jh=jh+1

	print('Total number of generations', generations, 'Final Chi2', chim )

#code returns the matrix and chi2 of the best fit 
	return motherest, chim, bestgt

def map_mass(gamma):
#Code to perform Kaiser-Squires analysis on a field of gamma values

	gamma=np.asanyarray(gamma)

	gamma=np.fft.fftshift(gamma)
	gamma=np.fft.fft2(gamma)
	gammat=np.fft.ifftshift(gamma)
	
	xsize=gammat.shape[1]
	ysize=gammat.shape[0]

	psit=np.zeros((gammat.shape[0], gammat.shape[1]), dtype=np.complex)
	kappat=np.zeros((gammat.shape[0], gammat.shape[1]), dtype=np.complex)

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


	kappa=np.fft.fftshift(kappat)
	kappa=np.fft.ifft2(kappa)
	kappa=np.fft.ifftshift(kappa)


	return kappa


def rad_2(true, test, xrad, yrad):
	true=true.copy()
	test=test.copy()

	truex0=true.shape[1]/2
	truey0=true.shape[0]/2


	radmax=int((xrad**2+yrad**2)**0.5)

	distances=range(0, radmax)
	differences=np.zeros((radmax+2), dtype=np.complex)
	diff=0
	count=0

	for radius in range(0, radmax):
		for x in range(int(truex0-xrad), int(truex0+xrad)):
			for y in range(int(truey0-yrad), int(truey0+yrad)):
				rad=int(((x-truex0)**2+(y-truey0)**2)**0.5)
				if rad==radius:
					diff=diff+(true[y,x]-test[y,x])
					count=count+1			
		differences[radius]=diff/count
		diff=0	
		count=0	
	
	return differences

	
def make_daughter_smart(binzx, binzy, mother, differences):
	mother=mother.copy()
	differences=differences.copy()
	daughter=mother.copy()
	daughterxorigin=(mother.shape[1]/2)
	daughteryorigin=(mother.shape[0]/2)
	meandiff=0
	for d in range(int(daughteryorigin-binzy), int(daughteryorigin+binzy)):
		for f in range(int(daughterxorigin-binzx),int(daughterxorigin+binzx)):
			dist=int(((daughteryorigin-d)**2+(daughterxorigin-f)**2)**0.5)
			realgauss=np.absolute(np.real(differences[dist]))
			imag=np.absolute(np.imag(differences[dist]))
			if realgauss!=0:
				alt=np.random.normal(0.0, scale=realgauss)
				daughter[d,f]=daughter[d,f]+alt
			if imag!=0:
				alt=(np.random.normal(0.0, scale=imag))*1j
				daughter[d,f]=daughter[d,f]+alt
			
	return daughter

def chi_shuffle_smart(true, mother, chideal, error, xwidth, ywidth):

#mother is the initial input matrix (hyp gammat), true being the target matrix (real gamma) which you want to minimise the difference in chi2 from. 
#xwidth is the x domain within which values are being varied, ywidth similarly. Defines a gausswidth itself although don't know a good reason to pick one yet.



	true=np.asanyarray(true)
	mother=np.asanyarray(mother)
	error=np.asanyarray(error)
	bestgt=mother


	xorigin=true.shape[1]/2
	yorigin=true.shape[0]/2

	motherr=gamma_build(mother.copy())


	kids=20
	bc=chideal+1
	generations=0
	gausswidth=20
	kl=[0]
	#cg=np.array(0)
	k=0

	truet=gammat_build(true)
	differences=rad_2(truet, mother, xwidth, ywidth)
	

	bc1=chideal+1
	l=0
	chim=float(chi_squared(true, motherr, error))
	print('#C2beat', chim)
	print('initial gauss', differences)
	cg=[chim]
	chix=chim

	jh=0
	motherest=mother.copy()
	differences=rad_2(truet, mother, xwidth, ywidth)
	print(differences)
	z=0
	gooduns=[z]
	cgood=[chim]
	dof=true.shape[0]*true.shape[1]


	while chideal < bc1:
		bet=best_daughter_smart(motherest, true, error, xwidth, ywidth, kids, differences*0.5)
		bd=bet[0]
		bc=bet[1]
		#print(bc, gausswidth, gausswidthi)
		#print('Mother generation', generations, mother[yorigin-ywidth:yorigin+ywidth,xorigin-xwidth:xorigin+xwidth])
		cg=np.append(cg, bc)
		kl=np.append(kl, k)
		mpl.scatter(kl,cg )
		mpl.xlabel('Iterations')
		mpl.ylabel('Chi2')
		mpl.savefig('cgk2 '+str(xwidth)+' by '+str(ywidth)+' Smart.png')
		mpl.close()
		k=k+1
		chik=bc
		chidiff=chik-chix
#		print(chidiff)


		if chidiff >0.01:
			l=l+1
			print(l)
		if chidiff <= 0.01:
			l=0
		if l>20:
			bc1=chideal
		if dof > chim:
			bc1=chideal
		if bc < chim:
			print('fitter than mum, old was', chim/dof, 'new is', bc/dof)
			chim=bc
			generations=generations+1
			cgood=np.append(cgood, bc)
			motherest=bd.copy()
			latestmass=map_mass(gamma_build(motherest.copy()))
			z=z+1
			gooduns=np.append(gooduns, generations)

			mpl.contourf(np.real(latestmass))
			mpl.colorbar()
			mpl.title('Forward Fitting attempt at modelling mass, Chi2='+str(chim/dof)+'.png')
			mpl.savefig('FFmass '+str(xwidth)+' by '+str(ywidth)+' Smart.png')
			mpl.close()
			chix=bc
			bc1=bc
			differences=rad_2(truet, bd, xwidth, ywidth)

			mpl.scatter(gooduns, cgood)
			mpl.title('Progression of best Chi estimate')
			mpl.xlabel('Generations')
			mpl.ylabel('Chi2')
			mpl.savefig('Best Chi evolutions '+str(xwidth)+' by '+str(ywidth)+' Smart.png')
			mpl.close()
		jh=jh+1

	print('Total number of generations', generations, 'Final Chi2', chim/dof )

#code returns the matrix and chi2 of the best fit 
	return motherest, chim, bestgt




def best_daughter_smart(mother, trueReal, error, xwidth, ywidth, loops, differences): 

	chilist=[0]*loops
	xorigin=trueReal.shape[1]/2
	yorigin=trueReal.shape[0]/2
	bestdaughter=mother.copy()

	daughter1=make_daughter_smart(xwidth, ywidth, mother.copy(), differences)
	daughter1real=gamma_build(daughter1)
	bestchi=float(chi_squared(trueReal, daughter1real, error))

	
	
	for g in range(0, int(loops-1)):
		daughter2=make_daughter_smart(xwidth, ywidth, mother, differences)
		daughter2[yorigin, xorigin]=0
		daughter2real=gamma_build(daughter2)

		newchi=float(chi_squared(trueReal, daughter2real, error))
		chilist[g]=newchi
		if newchi<bestchi:
			bestchi=newchi
			bestdaughter=daughter2


	return bestdaughter, bestchi

def fake_shear_normal(xpos, ypos, xsig, ysig, gammafake, mass):

	gammafake=gammafake.copy()
	cons=mass

	xo=(gammafake.shape[1]/2)+xpos
	yo=(gammafake.shape[0]/2)+ypos

	for x in range(0, gammafake.shape[1]):
		for y in range(0, gammafake.shape[0]):
			phi=np.arctan2((y-yo),(x-xo))
			gammafake[y,x]=gammafake[y,x]+cons*np.exp(-0.5*((x-xo)**2/xsig**2)-((y-yo)**2/ysig**2))*np.exp(1j*2*phi)

	return gammafake


def mass_comp(gammafake, mother, error,  motherrealkick, motherimagkick, xw, yw, spawns):

	mothere=mother.copy()
	gammafake=gammafake.copy()
	gammafakefft=gammat_build(gammafake)
	mass_composite=np.zeros((gammafake.shape[0], gammafake.shape[1]))
	chi2=[0]*spawns
	num=[0]*spawns
	mass_array=np.empty((gammafake.shape[0]*gammafake.shape[1], spawns), dtype=np.complex)

	mothere1=make_daughter(xw, yw, gammafakefft.copy(), motherrealkick, motherimagkick)
	mxo=mothere1.shape[1]/2
	myo=mothere1.shape[0]/2
	mothere=np.zeros(mothere.shape, dtype=np.complex)
	mothere[myo-yw:myo+yw, mxo-xw:mxo+xw]=mothere1[myo-yw:myo+yw, mxo-xw:mxo+xw]
	for i in range(0, spawns):
		f=chi_shuffle_smart(gammafake, mothere, 1, error, xw, yw)
		chim=f[1]

		latestmass=map_mass(gamma_build(f[0].copy()))
		massvec=np.reshape(latestmass.copy(), (-1,1))
		for h in range(0, int(massvec.shape[0])):
			mass_array[h,i]=massvec[h,0]
		#mass_array[i]=latestmass
		mpl.contourf(np.real(latestmass))
		mpl.colorbar()
		mpl.title('Forward Fitting attempt at modelling mass, Chi2='+str(chim)+' Iter='+str(i)+'.png')
		mpl.savefig('FFmass '+str(xw)+' by '+str(yw)+' Iter='+str(i)+'.png')
		mpl.close()
		mass_composite=mass_composite+latestmass

		mpl.contourf(np.real(mass_composite))
		mpl.colorbar()
		mpl.title('Forward Fitting attempt at modelling mass, Iter='+str(i)+'.png')
		mpl.savefig('FFmass '+str(xw)+' by '+str(yw)+'  Total Iter='+str(i)+'.png')
		mpl.close()

		rat=np.absolute(np.random.randn(1))
		motherrealkick=2*rat*motherrealkick+0.01
		rat=np.absolute(np.random.randn(1))
		motherimagkick=2*rat*motherimagkick+0.01

		chi2[i]=chim
		num[i]=i

#		ctest=covariance(latestmass, mass_composite)
#		mpl.imshow(np.real(ctest))
#		mpl.title('Latest mass on y covariance with total')
#		mpl.savefig('Covarianceimshow'+str(i)+'.png')
#		mpl.close()
		ctest=covariance_pixel(latestmass, mass_composite)


		mothere1=make_daughter(xw, yw, gammafakefft.copy(), motherrealkick, motherimagkick)
		mothere=np.zeros(mothere1.shape, dtype=np.complex)
		mothere[myo-yw:myo+yw, mxo-xw:mxo+xw]=mothere1[myo-yw:myo+yw, mxo-xw:mxo+xw]
		print('Mass map made, number=', i)



	mpl.scatter(num, chi2)
	mpl.xlabel('Attempt')
	mpl.ylabel('Chi2')
	mpl.savefig('Chi2iters.png')
	mpl.close()




	return mass_composite, mass_array


def covariance(matA, matB):


	first=matA.copy()
	second=matB.copy()

	first=np.reshape(first, (-1,1))
	second=np.reshape(second, (-1,1))

	if first.shape != second.shape:
		print('Error: Arrays of differing shape, refilling array with zeros')
		if first.shape[0]<second.shape[0]:
			first1=np.zeros(second.shape)
			first1[0:first.shape[0]]=first
			first=first1
		if second.shape[0]<first.shape[0]:
			first1=np.zeros(first.shape)
			first1[0:second.shape[0]]=first
			second=first1

	fave=np.average(first)
	save=np.average(second)
	cov=np.empty((1, first.shape[0] ))

	for t in range(0, first.shape[0]):
		for j in range(0, second.shape[0]):
			cov[t,j]=(first[t]-fave)*(second[j]-save)


	return cov

def covariance_pixel(matA, matB):

	first=matA.copy()
	second=matB.copy()

	first=np.reshape(first, (-1,1))
	second=np.reshape(second, (-1,1))

	if first.shape != second.shape:
		print('Error: Arrays of differing shape, refilling array with zeros')
		if first.shape[0]<second.shape[0]:
			first1=np.zeros(second.shape)
			first1[0:first.shape[0]]=first
			first=first1
		if second.shape[0]<first.shape[0]:
			first1=np.zeros(first.shape)
			first1[0:second.shape[0]]=first
			second=first1

	fave=np.average(first)
	save=np.average(second)

	firstagain=matA.copy()
	secondagain=matB.copy()
	cov=np.empty(firstagain.shape, dtype=np.complex)

	for t in range(0, int(firstagain.shape[0])):
		for g in range(0, int(secondagain.shape[1])):
			fave=(firstagain[t,g]+secondagain[t,g])*0.5
			save=fave
			cov[t,g]=(firstagain[t,g]-fave)*(secondagain[t,g]-save)

	mpl.imshow(np.real(cov))
	mpl.title('Covariance Matrix')
	mpl.colorbar()
	mpl.savefig('Covariance.png')
	mpl.close()

	return cov


def ben_men(mass_vec, spawns):

	mass_vec=np.asarray(mass_vec)
	cov=np.empty((mass_vec.shape[0],1))
	for i in range(0, mass_vec.shape[0]):
		pixave=np.sum(mass_vec[i,:])/spawns
		pixc=mass_vec[i,:]-pixave
		pixcov=np.prod(pixc)
		cov[i]=pixcov

	square=np.sqrt(mass_vec.shape[0])
	print(square)
	cov=np.reshape(cov, (square, square))
	return cov

def cov_indie(mass_vec, spawns):
			
	mass_vec=np.asarray(mass_vec)
	cov=np.zeros((mass_vec.shape[0],1), dtype=np.complex)
	colaves=[0]*spawns
	totl=cov.shape[0]
	for l in range(0, spawns):
		colaves[l]=np.sum(mass_vec[:,l])/totl
	for i in range(0, mass_vec.shape[0]-1):
		
		pixc=mass_vec[i,:]-colaves
		pixcov=np.prod(pixc)
		cov[i]=pixcov

	square=np.sqrt(mass_vec.shape[0])
	cov=np.reshape(cov, (square, square))
	return cov






def make_daughter_sp( mother, differences1, nside) :
	jh1=mother.copy()
	jh2=np.zeros(jh1.shape, dtype=np.complex) #empty array for offspring
	differences=np.asanyarray(differences1)
	#print(differences.shape)
	for g in range(0, jh1.shape[0]):
		diff=differences[g] #array of differences per lm
#		print(diff.shape)
		if (np.absolute(np.real(diff)) > 0):
			jh2[g]=jh2[g]+jh1[g]+np.random.normal(0.0, np.absolute(np.real(diff)))
		
		if (np.absolute(np.imag(diff))) >0:
			jh2[g]=jh2[g]+jh1[g]+1j*np.random.normal(0.0, np.absolute(np.imag(diff)))#add random dependent upon typical diff
#	jh2list=jh2[0].reshape(-1)#output as list of lm
#	print(jh2list.shape)

#	jh2map=hp.alm2map(jh2list, nside=nside) #function outputs map in pixel space
	return jh2

def best_daughter_sp(mother, true, error, differences, loops, nside):
	#takes inputs of mother (guessalm), true(real map of shear), error on true, differences in alm space of values between true and mother, iterations to output the best guess
	mother2=np.asanyarray(mother)
	malm=mother2.copy() #mother input as pixel space map of three arrays - a scalar (ie galaxy density), gamma1, gamma2
	true2=np.asanyarray(true)
#	[galdens, E_lm_m, B_lm_m]=hp.map2alm(mother2.copy(), pol=True)
#	malm=E_lm.copy()
	[galdensr, E_lm, B_lm]=hp.map2alm(true2.copy(), pol=True)
	fudge=error.copy()
	correct=true[1]+1j*true[2]

	
	daughter1=make_daughter_sp(malm.copy(), differences, nside=nside)
	bestdaughter=daughter1.copy()
	d1in=[galdensr, daughter1, B_lm]
	[num, d1shear1, d1shear2]=hp.alm2map(d1in, nside=nside, pol=True)
	chilist=[0]*loops
	chilist[0]=chi_squared(d1shear1+1j*d1shear2, correct, fudge)
#	print(chilist.shape)
	print(chilist[0].shape)
	bestchi=chilist[0]
	for g in range(0, loops):
		daughter2=make_daughter_sp(malm.copy(), differences, nside=nside)
		d2in=[galdensr, daughter2, B_lm]
		[num, d2shear1, d2shear2]=hp.alm2map(d2in, nside=nside, pol=True)
		chiagain=chi_squared(d2shear1+1j*d2shear2, correct, fudge)
		chilist[g]=chiagain
		print(bestchi.shape)
		if chilist[g] < bestchi:
			bestchi=chiagain
			bestdaughter=daughter2.copy()
			print('woo')
		print(g)
	dout=[galdensr, bestdaughter, B_lm]

	hp.visufunc.mollview(map_mass_lm(bestdaughter.copy(), (3*nside)-1, nside))
	mpl.savefig('fittedmass.png')
	mpl.close()
	[gals, bestmapg1, g2]=hp.alm2map(dout, nside=nside, pol=True)
	return bestmapg1+(1j*g2), bestchi


def map_mass_lm(E_lm, L, nside):
#courtesy of Boris Leistedt's code
	lfac=np.zeros(E_lm.shape)
	for el in range(L):
		for em in range(L):
			lfac[lm_indices(el, em, L)]=el

	phiE_lm=-2*E_lm/np.sqrt((lfac+2)*(lfac+1)*(lfac)*(lfac-1))
	phiE_lm[np.isnan(phiE_lm)]=0
	kappa_lm=-lfac*(lfac+1)*(phiE_lm/2)
	kappa_lm[np.isnan(kappa_lm)]=0
	kappa_map=hp.alm2map(kappa_lm, nside=nside, lmax=3*nside-1, pol=False)

	return kappa_map


def lm_indices(l, m, L):
	indices=m*(2*L-1-m)/(2+l)
	return indices


def chi_shuffle_sp(mother, true, error, differences, kids, nside, chideal):

	true=np.asanyarray(true[1]+true[2]*1j)
	mother=np.asanyarray(mother[1]+true[2]*1j)
	error=np.asanyarray(error)
	bestgt=mother


	bc=chideal+1
	generations=0
	gausswidth=20
	kl=[0]
	#cg=np.array(0)
	k=0
	loops=30

	chi2beat=best_daughter_sp(mother.copy(), true.copy(), error.copy(), differences, 2, nside)
	bc1=chi2beat[1]
	chim=chi2beat[1]


	while chideal < bc1:
		bet=best_daughter_sp(mother.copy(), true.copy(), error.copy(), differences, loops, nside)
		bd=bet[0]
		bc=bet[1]
		#print(bc, gausswidth, gausswidthi)
		#print('Mother generation', generations, mother[yorigin-ywidth:yorigin+ywidth,xorigin-xwidth:xorigin+xwidth])
		cg=np.append(cg, bc)
		kl=np.append(kl, k)
		mpl.scatter(kl,cg )
		mpl.xlabel('Iterations')
		mpl.ylabel('Chi2')
		mpl.savefig('cgk2 '+str(xwidth)+' by '+str(ywidth)+' Smart.png')
		mpl.close()
		k=k+1
		chik=bc
		chidiff=chik-chix
#		print(chidiff)


		if chidiff >0.01:
			l=l+1
			print(l)
		if chidiff <= 0.01:
			l=0
		if l>20:
			bc1=chideal
		if dof > chim:
			bc1=chideal
		if bc < chim:
			print('fitter than mum, old was', chim, 'new is', bc)
			chim=bc
			generations=generations+1
			cgood=np.append(cgood, bc)
			motherest=bd.copy()
			z=z+1
			chix=bc
			bc1=bc
			differences=rad_2(truet, bd, xwidth, ywidth)

		jh=jh+1

	print('Total number of generations', generations, 'Final Chi2', chim/dof )

#code returns the matrix and chi2 of the best fit 
	return motherest, chim, bestgt



