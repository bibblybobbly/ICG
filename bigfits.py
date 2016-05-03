import numpy as np
from astropy.io import fits


infile=["DES0025-4206-r.fits","DES0025-4249-r.fits","DES0026-4331-r.fits","DES0029-4206-r.fits","DES0029-4249-r.fits","DES0030-4331-r.fits","DES0030-4414-r.fits","DES0033-4206-r.fits", "DES0033-4249-r.fits","DES0034-4331-r.fits","DES0034-4414-r.fits","DES0037-4206-r.fits","DES0037-4249-r.fits","DES0038-4331-r.fits","DES0038-4414-r.fits","DES0039-4457-r.fits","DES0041-4249-r.fits","DES0042-4331-r.fits","DES0042-4414-r.fits","DES0043-4457-r.fits","DES0045-4331-r.fits","DES0212-0458-r.fits","DES0213-0541-r.fits","DES0215-0416-r.fits","DES0215-0458-r.fits","DES0216-0541-r.fits","DES0217-0333-r.fits","DES0217-0624-r.fits","DES0217-0707-r.fits","DES0218-0416-r.fits","DES0218-0458-r.fits","DES0219-0541-r.fits","DES0219-0624-r.fits","DES0220-0707-r.fits","DES0221-0416-r.fits","DES0221-0458-r.fits","DES0222-0541-r.fits","DES0222-0624-r.fits","DES0223-0333-r.fits","DES0223-0416-r.fits","DES0223-0707-r.fits","DES0224-0458-r.fits","DES0224-0541-r.fits","DES0225-0624-r.fits","DES0226-0333-r.fits","DES0226-0416-r.fits","DES0226-0707-r.fits","DES0227-0458-r.fits","DES0227-0541-r.fits","DES0228-0624-r.fits","DES0229-0333-r.fits","DES0229-0416-r.fits","DES0229-0707-r.fits","DES0230-0458-r.fits","DES0230-0541-r.fits","DES0231-0624-r.fits","DES0232-0333-r.fits","DES0232-0416-r.fits","DES0233-0458-r.fits","DES0233-0541-r.fits","DES0234-0624-r.fits","DES0239-0041-r.fits","DES0239-0124-r.fits","DES0242-0041-r.fits","DES0242-0124-r.fits","DES0242-0207-r.fits","DES0242+0001-r.fits","DES0245-0041-r.fits","DES0245-0124-r.fits","DES0245-0207-r.fits","DES0245+0001-r.fits","DES0248-0041-r.fits","DES0248-0124-r.fits","DES0248-0207-r.fits","DES0248+0001-r.fits","DES0248+0043-r.fits","DES0251-0041-r.fits","DES0251+0001-r.fits","DES0254+0001-r.fits","DES0251+0043-r.fits","DES0254-0041-r.fits","DES0254+0043-r.fits","DES0256+0001-r.fits","DES0325-2749-r.fits","DES0326-2832-r.fits","DES0326-2915-r.fits","DES0328-2706-r.fits","DES0328-2749-r.fits","DES0329-2832-r.fits","DES0329-2915-r.fits","DES0331-2623-r.fits","DES0331-2706-r.fits","DES0332-2749-r.fits","DES0332-2832-r.fits","DES0333-2915-r.fits","DES0333-2958-r.fits","DES0334-2623-r.fits","DES0334-2706-r.fits","DES0339-2832-r.fits","DES0339-2915-r.fits","DES0340-2958-r.fits","DES0341-2623-r.fits","DES0341-2706-r.fits","DES0341-2749-r.fits","DES0342-2832-r.fits","DES0342-2915-r.fits"]



j=fits.getdata(infile[0])
x=j.field('Dec')

for g in range(1, len(infile)):
	ja=fits.getdata(infile[g])
	xadd=ja.field('Dec')
	x=np.append(x, xadd, axis=0)


	
fits.writeto('xpoints.fits', x)
