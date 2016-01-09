import numpy as np
import matplotlib.pyplot as mpl
import matplotlib as matplot




Rx=  #ellipticity measurements in x direction
Ry=  #ellipticity measurements in y direction

Rxshift=np.fft.fftshift(Rx)
Ryshift=np.fft.fftshift(Ry)   #shifts vectors into edges for FFT

Rxsft=np.fft.fft2(Rxshift)
Rysft=np.fft.fft2(Ryshift)     # FFT of the vectors

gamma1t=np.fft.ifftshift(Rxsft)
gamma2t=np.fft.ifftshift(Rysft)   # converting FFT into gammas, by inverse shifting

gammat=gamma1t+np.multiply(gamma2t, 1j)  # gamma tilda, by definition

ysize=gammat.shape[0]
xsize=gammat.shape[1] #finds length of the axes 


#~~~~~~~~~~~~~~~~~~~~~Transformed arrays~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
psit=np.zeros((ysize, xsize), dtype=np.complex)
kappat=np.zeros((ysize, xsize), dtype=np.complex) #define arrays to be filled by for loop, sizes of dimensions

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

i=0
for l in range(0,xsize):
	for i in range(0, ysize):
		xi1=(i-xsize/2)+0.5
		xi2=(l-ysize/2)+0.5 #coordinates in Fourier domain, 0,0 being the bottom left hand coordinate of the gamma[ysize/2, xsize/2] coordinate. 
		if xi1 == 0 or xi2 ==0: 
			psit[l,i]=0 # otherwise is infinite
		
					
		else:
			psit[l,i]=np.divide(gammat[l,i]*(xi1**2-xi2**2-2j*xi1*xi2),2*(np.pi**2)*(xi1**2*xi2**2)**2) # calculating psi tilda, according to eqns [15/12/15]
		
		kappat[l,i]=(xi1**2-xi2**2-2j*(xi1*xi2))*gammat[l,i]   #calculate kappa tilda according to eqns [15/12/15]
	l=0


psiprep=np.fft.fftshift(psit) #shifting psit for inverse FFT
kappaprep=np.fft.fftshift(kappat) #shifting kappat for inverse FFT



psireal=np.fft.ifft2(psiprep)  #inverse FFT of psit
kappareal=np.fft.ifft2(kappaprep)

psi=np.fft.ifftshift(psireal) # inverse shift of psi, should be done!
kappa=np.fft.ifftshift(kappareal) # shifting kappa back, should be it!
