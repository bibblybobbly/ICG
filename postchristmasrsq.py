import numpy as np
import matplotlib.pyplot as mpl
import matplotlib as matplot

#omgi'medittingandgitting

xsize=100
ysize=100

xpos=50
ypos=50 #location of centre

sigmax=5
sigmay=5
sigma=10
fac=10000

R2=np.zeros((xsize, ysize))
theta=np.zeros((xsize, ysize))
R=np.zeros((xsize, ysize))


for i in range(0,xsize):#filling xsize x ysize array with gamma values 
	for j in range(0, ysize):
		#R[j,i]=np.exp(-((i-xpos)**2+(j-ypos)**2)/(sigma**2))
		
		if i-xpos ==0:
			if j-ypos<0:
				theta[j,i]=-(np.pi)/2
				#R[j,i]=1/((j-ypos)**2+(i-xpos)**2)
			if j-ypos >0:
				theta[j,i]=(np.pi)/2
				#R[j,i]=1/((j-ypos)**2+(i-xpos)**2)
			if j-ypos ==0:
				theta[j,i]=np.pi/2
			theta[j,i]=theta[j,i]+np.pi/2
		else:
			theta[j,i]=(np.pi/2)+np.arctan((j-ypos)/(i-xpos))
			#R[j,i]=1/((j-ypos)**2+(i-xpos)**2)
		
			
		#print(j,i)
		if i-xpos !=0 or j-ypos !=0:
			R[j,i]=1/((i-xpos)**2+(j-ypos)**2)
		if i-xpos <0:
			#theta[j,i]=theta[j,i]
			R[j,i]=-R[j,i]	
#R[49,50]=0.00002
#R[50,49]=0.00002
#R[51,50]=0.00002
#R[50,51]=0.00002
R[50,50]=0
print(np.max(theta), np.min(theta))

Rx=np.multiply(R, np.cos(theta))
Ry=np.multiply(R, np.sin(theta))

mpl.quiver(Rx, Ry)

mpl.show()

for i in range(0, xsize):
	for j in range(0, ysize):
		if i<50:
			Rx[j,i]=-Rx[j,i]
			Ry[j,i]=-Ry[j,i]


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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
psit=np.zeros((ysize, xsize), dtype=np.complex)
kappat=np.zeros((ysize, xsize), dtype=np.complex) #define arrays to be filled by for loop, sizes of dimensions

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

i=0
for l in range(0,xsize):
	for i in range(0, ysize):
		xi1=i-xsize/2
		xi2=l-ysize/2
		if xi1 == 0: 
			psit[l,i]=0
		else:
			psit[l,i]=np.divide(gammat[l,i]*(xi1**2-xi2**2-2j*xi1*xi2),2*(np.pi**2)*(xi1**2*xi2**2)**2)
		if xi2==0:
			psit[l,i]=0 # otherwise is infinite
					
		else:
			psit[l,i]=np.divide(gammat[l,i]*(xi1**2-xi2**2-2j*xi1*xi2),2*(np.pi**2)*(xi1**2*xi2**2)**2) # calculating psi tilda, according to eqns [15/12/15]

		kappat[l,i]=(xi1**2-xi2**2-2j*(xi1*xi2))*gammat[l,i]   #calculate kappa tilda according to eqns [15/12/15]
	l=0
 

psigr=np.real(psit)

psiprep=np.fft.fftshift(psit) #shifting psit for inverse FFT
kappaprep=np.fft.fftshift(kappat) #shifting kappat for inverse FFT

psireal=np.fft.ifft(psiprep)  #inverse FFT of psit
kappareal=np.fft.ifft(kappaprep)
print('psir', psireal)
print('kappar', kappareal)

psi=np.fft.ifftshift(psireal) # inverse shift of psi, should be done!
kappa=np.fft.ifftshift(kappareal) # shifting kappa back, should be it!


#Produce plots of both kappa and psi

psideforeal=np.real(np.multiply(psi, np.conj(psi)))
kappadeforeal=np.real(np.multiply(kappa, np.conj(kappa)))


mpl.figure(5)
mpl.subplot(1,2,1)
mpl.imshow(psigr)
mpl.title('Psi results')
mpl.colorbar()

mpl.subplot(1,2,2)
mpl.imshow(kappadeforeal)
mpl.title('Kappa results')
mpl.colorbar()
mpl.show()

mpl.figure(6)
mpl.imshow(psideforeal)
mpl.title('Psi')
mpl.show()


