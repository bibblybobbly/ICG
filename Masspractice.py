import numpy as np
import matplotlib.pyplot as mpl
import matplotlib as matplot



xsize=100
ysize=100

xpos=87
ypos=78 #location of centre

sigmax=5
sigmay=5
sigma=10
fac=10000
kons= 5*10**1 # constant in front of r^-2
R2=np.zeros((xsize, ysize))
theta=np.zeros((xsize, ysize))
R=np.zeros((xsize, ysize))


for i in range(0,xsize):#filling xsize x ysize array with gamma values 
	for j in range(0, ysize):
		
		
		if i-xpos ==0:
			if j-ypos<0:
				theta[j,i]=0
				
			if j-ypos >0:
				theta[j,i]=0
				
			if j-ypos ==0:
				theta[j,i]=np.pi/2
			
		else:
			theta[j,i]=(np.pi/2)+np.arctan((j-ypos+0.5)/(i-xpos+0.5))
		
			
		
		if i-xpos !=0 or j-ypos !=0:
			R[j,i]=kons/np.sqrt((i-xpos+0.5)**2+(j-ypos+0.5)**2)

print(np.max(theta), np.min(theta))

Rx=np.multiply(R, np.cos(2*theta)) #vectors in shear space
Ry=np.multiply(R, np.sin(2*theta))

mpl.quiver(Rx, Ry)

mpl.show()

mpl.imshow(Rx)

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

#gammats=np.fft.fftshift(gammat)
#gammas=np.fft.ifft2

mpl.figure(55)
mpl.imshow(np.real(gammat))
print('49,49', gammat[49,49])
print('49,51', gammat[49,51])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kappa2=np.zeros((ysize, xsize))
psit=np.zeros((ysize, xsize), dtype=np.complex)
kappat=np.zeros((ysize, xsize), dtype=np.complex) #define arrays to be filled by for loop, sizes of dimensions

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

i=0
for l in range(0,xsize):
	for i in range(0, ysize):
		xi1=(i-xsize/2)
		xi2=(l-ysize/2)
		if xi1 == 0 and xi2 ==0: 
			psit[l,i]=0 # otherwise is infinite
		
					
		else:
			psit[l,i]=np.divide(gammat[l,i]*(xi1**2-xi2**2-2j*xi1*xi2),2*(np.pi**2)*(xi1**2+xi2**2)**2) # calculating psi tilda, according to eqns [15/12/15]
		
		if xi1!=0 or xi2!=0:
			kappat[l,i]=np.divide((xi1**2-xi2**2-2j*(xi1*xi2))*gammat[l,i], (xi1**2+xi2**2))   #calculate kappa tilda according to eqns [15/12/15]
	
		else:
			kappat[l,i]=0 #is infinite at xi1=0, xi2=0  

		kappa2[l,i]=(2*(np.pi)**2)*(xi1**2+xi2**2)*psit[l,i]
 

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
 

mpl.figure(5)
mpl.subplot(1,2,1)
mpl.imshow(psi)
mpl.title('Real component of Psi results')
mpl.colorbar()

mpl.subplot(1,2,2)
mpl.imshow(kappa)
mpl.title('Real component of Kappa results')
mpl.colorbar()
mpl.show()



