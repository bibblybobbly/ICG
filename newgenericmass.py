Rx=
Ry=



#NEED TO DEFINE SIZES

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

#gammats=np.fft.fftshift(gammat)
#gammas=np.fft.ifft2



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
psit=np.zeros((ysize, xsize), dtype=np.complex)
kappat=np.zeros((ysize, xsize), dtype=np.complex) #define arrays to be filled by for loop, sizes of dimensions

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

i=0
for l in range(0,xsize):
	for i in range(0, ysize):
		xi1=(i-xsize/2)
		xi2=(l-ysize/2)
		if xi1 == 0 or xi2 ==0: 
			psit[l,i]=0 # otherwise is infinite
		
					
		else:
			psit[l,i]=np.divide(gammat[l,i]*(xi1**2-xi2**2-2j*xi1*xi2),2*(np.pi**2)*(xi1**2+xi2**2)**2) # calculating psi tilda, according to eqns [15/12/15]
		
		if xi1!=0 or xi2!=0:
			kappat[l,i]=np.divide((xi1**2-xi2**2-2j*(xi1*xi2))*gammat[l,i], (xi1**2+xi2**2))   #calculate kappa tilda according to eqns [15/12/15]
	
		else:
			kappat[l,i]=1 #is infinite at xi1=0, xi2=0  

	
 

psigr=np.real(psit)
print('kappat', kappat)
psiprep=np.fft.fftshift(psit) #shifting psit for inverse FFT
kappaprep=np.fft.fftshift(kappat) #shifting kappat for inverse FFT

psireal=np.fft.ifft2(psiprep)  #inverse FFT of psit
kappareal=np.fft.ifft2(kappaprep)
print('psir', psireal)
print('kappar', kappareal)

psi=np.fft.ifftshift(psireal) # inverse shift of psi, should be done!
kappa=np.fft.ifftshift(kappareal) # shifting kappa back, should be it!

print('psi', psi)
print('kappa', kappa)

psi=np.real(psi)
kappa=np.real(kappa)
#Produce plots of both kappa and psi




