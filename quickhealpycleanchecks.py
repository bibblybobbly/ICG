
import matplotlib as matplotlib
matplotlib.use('Agg')
import astropy
import numpy as np
import matplotlib.pyplot as mpl
from astropy.io import fits
import pylab
import healpy as hp
#import scipy.stats as sp
import Benlibv2TEST as ben

#========================= func n ting=================================
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

def pixelise_sp(emptypix, gammaphake):

#input arrays in order to average over pixels in alm space for healpy. emptypix is an array filled with pixel numbers of the data, gammaphake is an array filled with values for those pixels

	pixcatalog=emptypix.copy().reshape((-1,1))

	densmap=np.zeros(emptypix.shape)
	outsheartot=np.zeros(gammaphake.shape, dtype=np.complex)
	donelist=[0]




	for i in range(0, pixcatalog.shape[0]):
		pix=pixcatalog[i]
		finder=any(pixcatalog[i]==donelist)
		if finder==False:
			for j in range(0, emptypix.shape[0]):
				for k in range(0, emptypix.shape[1]):
					thisone=emptypix[j,k]
					if thisone ==pix:
						outsheartot[j,k]=outsheartot[j,k]+gammaphake[j,k]
						densmap[j,k]=densmap[j,k]+1
			donelist=np.append(donelist,pixcatalog[i])

	return(outsheartot, densmap)




#==========================
#make mass blob
#=========================


xw=8
yw=8
spawns=8

gammafake=np.zeros((200,200), dtype=np.complex)

gammafake=-ben.fake_shear_normal(0,0, 500, 500, gammafake, 1)

mpl.contourf(np.real(gammafake))
mpl.savefig('ooooooflat.png')
mpl.close()

mpl.contourf(np.absolute(np.imag(gammafake.copy())))
mpl.savefig('gamma.png')
mpl.close()

mpl.contourf(np.absolute(np.real(gammafake.copy())))
mpl.savefig('gamma1.png')
mpl.close()



mpl.contourf(np.real(ben.map_mass(gammafake.copy())))
mpl.savefig('a.png')
mpl.close()

DecRange=np.zeros((200,200))
RARange=np.zeros((200,200))

RAMin=0+np.pi/20#-np.pi#/2#-np.pi+np.pi/20
DecMin=(np.pi/4)-np.pi/20

decgrad=np.pi/2000
ragrad=np.pi/2000

for x in range(0, DecRange.shape[0]):
	for y in range(0, RARange.shape[1]):
		DecRange[y,x]=DecMin+decgrad*y
		RARange[y,x]=RAMin+ragrad*x

centre=100
r=1
#decdelta=r*(DecRange)**2
#radelta=r*(np.sin(DecRange))**2*(RARange)**2

RACentre=RARange[20,20]#RAMin+(RARange.shape[0]*ragrad)
DecCentre=DecRange[20,20]#DecMin+(DecRange.shape[0]*decgrad)

RAdiffs=-(RARange-RACentre)
Decdiffs=-(DecRange-DecCentre)

decdelta=r*(Decdiffs)**2
radelta=r*(np.sin(DecRange))**2*(RAdiffs)**2
#radelta=r*(RAdiffs)**2
print(' all positive, around pi/2', DecRange)
print(np.sin(DecRange))
print(Decdiffs)




RACentre2=RARange[65,65]#RAMin+(RARange.shape[0]*ragrad)
DecCentre2=DecRange[65,65]#DecMin+(DecRange.shape[0]*decgrad)


RAdiffs2=-(RARange-RACentre2)
Decdiffs2=-(DecRange-DecCentre2)

decdelta2=r*(Decdiffs2)**2
radelta2=r*(np.sin(DecRange))**2*(RAdiffs2)**2

phi=np.arctan2(Decdiffs, np.multiply(np.sin(DecRange),RAdiffs))
phiflat=np.arctan2(Decdiffs, np.multiply(np.sin(DecRange),-RAdiffs))

phi2=np.arctan2(Decdiffs2, np.multiply(np.sin(DecRange),RAdiffs2))
phiflat2=np.arctan2(Decdiffs2, np.multiply(np.sin(DecRange),-RAdiffs2))

print(phi)
radii=np.sqrt(decdelta+radelta)
radii2=np.sqrt(decdelta2+radelta2)
#phi=np.arctan2(DecRange-100*decgrad-np.pi/2, RARange-100*ragrad-np.pi/2)
#phi=np.arctan2(Decdiffs, RAdiffs*np.absolute(np.sin(DecRange)))#np.sin(DecRange)*RAdiffs)
sigx=0.2

gammaphake=np.zeros((200,200), dtype=np.complex)
gammaflat=np.zeros((200,200), dtype=np.complex)

for d in range(0, gammaphake.shape[0]):
	for f in range(0, gammaphake.shape[1]):
		gammaphake[d,f]=(np.exp(-radii[d,f]**2/sigx**2)*np.exp(-2*1j*phi[d,f]))+(np.exp(-radii2[d,f]**2/(sigx+0.1)**2)*np.exp(-2*1j*phi2[d,f]))
		gammaflat[d,f]=(np.exp(-radii[d,f]**2/sigx**2)*np.exp(-2*1j*(phiflat[d,f])))


#gamm2t=np.imag(gammaphake.copy())
#gam1t=-np.real(gammaphake.copy())
print(gammaphake)
#gammaphake=gam1t+1j*gamm2t
flat=ben.map_mass(gammaflat.copy())
emptypix=hp.pixelfunc.ang2pix(64, DecRange, RARange)
print('EMPTY',emptypix.shape)
kappa=np.zeros((hp.nside2npix(64)), dtype=np.complex)
kappa[emptypix]=flat

pixcatalog=emptypix.copy().reshape((-1,1))

densmap=np.zeros(emptypix.shape)
outsheartot=np.zeros(gammaphake.shape, dtype=np.complex)
donelist=[0]


pixcatalog=emptypix.copy().reshape((-1,1))

densmap=np.zeros(emptypix.shape)
outsheartot=np.zeros(gammaphake.shape, dtype=np.complex)
donelist=[0]




for i in range(0, pixcatalog.shape[0]):
	pix=pixcatalog[i]
	finder=any(pixcatalog[i]==donelist)
	if finder==False:
		for j in range(0, emptypix.shape[0]):
			for k in range(0, emptypix.shape[1]):
				thisone=emptypix[d,f]
				if thisone ==pix:
					outsheartot[j,k]=outsheartot[j,k]+gammaphake[j,k]
					densmap[j,k]=densmap[j,k]+1
					print('copy!!!')
		donelist=np.append(donelist,pixcatalog[i])





outshear=outsheartot/densmap
newcategory=np.zeros(kappa.shape, dtype=np.complex)
newcategory[emptypix]=outshear

hp.visufunc.mollview(np.imag(newcategory))
mpl.savefig('filtered.png')
mpl.close()

hp.visufunc.mollview(kappa)
mpl.savefig('Masstransfer.png')
mpl.close()

mpl.contourf(np.real(flat))
mpl.colorbar()
mpl.savefig('Flatmass.png')
mpl.close()


mpl.contourf(flat)
mpl.savefig('shutuppaul.png')
mpl.close()

emptypix=hp.pixelfunc.ang2pix(64, DecRange, RARange)

shear_healpix_map=np.zeros((hp.nside2npix(64)), dtype=np.complex)
shear_healpix_map[emptypix]=gammaphake#gammafake

hp.visufunc.mollview((np.real(shear_healpix_map.copy())))
mpl.savefig('gamma1.png')
gamma2=np.real(np.imag(shear_healpix_map.copy()))
print(gamma2.shape)
print('this is gamma2', np.where(gamma2<0))
#print(gamma2[22594])
hp.visufunc.mollview(gamma2)
mpl.savefig('gamma2.png')

gamma1=np.real(shear_healpix_map.copy())

print(gamma1)

gamma1_lm=hp.sphtfunc.map2alm(gamma1, pol=False)
gamma2_lm=hp.sphtfunc.map2alm(gamma2, pol=False)

#flat_lm=hp.sphtfunc.map2alm(hp.
U=shear_healpix_map.copy()*0
T=np.absolute(shear_healpix_map.copy())#shear_healpix_map.copy()*0

T=((np.absolute(shear_healpix_map.copy()))*10)
T=kappa
tot=[T, -gamma1, -gamma2]


[T2, gamma_lm, U2]=hp.sphtfunc.map2alm(tot, pol=True)
#gamma_lm=gamma1_lm+1j*gamma2_lm
#gamma_lm=hp.map2alm(gamma1)#+1j*gamma2)

mass=map_mass_lm(gamma_lm, 191, 64)
hp.visufunc.mollview(mass)
mpl.savefig('massoutg1.png')
mpl.close()

bmode=map_mass_lm(U2, 191, 64)
hp.visufunc.mollview(bmode)
mpl.savefig('massoutb.png')
mpl.close()

hp.visufunc.mollview(hp.sphtfunc.alm2map(U2, nside=64, pol=False), 191, 64)
mpl.savefig('U2.png')
mpl.close()

hp.visufunc.mollview(hp.sphtfunc.alm2map(gamma_lm, nside=64, pol=False), 191, 64)
mpl.savefig('gamma_lm.png')
mpl.close()

hp.visufunc.mollview(hp.sphtfunc.alm2map(T2, nside=64, pol=False), 191, 64)
mpl.savefig('t2.png')
mpl.close()


mother=np.zeros(gamma_lm.shape, dtype=np.complex)
#mother=[mass, motherin, np.zeros(gamma_lm.shape, dtype=np.complex)]
true=[mass, np.real(shear_healpix_map.copy()), np.imag(shear_healpix_map.copy())]
error=np.ones(shear_healpix_map.shape)*0.02
differences=mother.copy()-gamma_lm.copy()
print('gammalm', gamma_lm.shape)

testb=ben.chi_shuffle_smart(mother, true, error, differences, 25, 64, 1)
testa=ben.best_daughter_sp(mother, true, error, differences, 250, 64)

hp.visufunc.mollview(testa, 191, 64)
mpl.savefig('fitted.png')
mpl.close()

#mymass=map_mass_lm(testa, 191, 64)
#hp.visufunc(mymas




