##
# coordinates of the pixel in the 4alin.dat idl save file ii=10 & jj=25
import numpy as np

## Noisy CLE profiles
si1=np.zeros(30, dtype=np.float64)
si1[:]=[7.24176, 7.17481, 7.21602, 7.18671, 7.25229, 7.23346, 7.30228, 7.45673, 7.76848, 8.26106, 9.11638, 10.1656, 11.5156, 12.6786, 13.6026, 13.7896, 13.3566, 12.4210, 11.0947, 9.89914, 8.89439, 8.06619, 7.65103, 7.39594, 7.29253, 7.26661, 7.19432, 7.18970, 7.16826, 7.17181]

sq1=np.zeros(30, dtype=np.float64)
sq1[:]=[0.0968344,0.124722,0.148907,0.152362,0.158669,0.140984,0.146622,0.143637,0.199675,0.246524,0.361636,0.473835,0.546483,0.761311,0.835028,0.851032,0.795890,0.711007,0.584043,0.459563,0.313614,0.187505,0.176667,0.161869,0.118569,0.138428,0.0969734,0.0862800,0.117751,0.128217]

su1=np.zeros(30, dtype=np.float64)
su1[:]=[-0.044457875,0.019408356,-0.014364263,-0.0099311592,0.014877109,-0.029197676,0.0082678441,0.040854529,0.051509216,0.11864027,0.20947164,0.38501257,0.50194567,0.61440998,0.75892055,0.79546148,0.70428574,0.55108470,0.47895029,0.32626709,0.19210990,0.095756948,0.053563900,0.029628001,-0.014606290,-0.010244418,-0.0076554865,-0.012123742,-0.0050786901,-0.039388664]

sv1=np.zeros(30, dtype=np.float64)
sv1[:]=[-0.018210001,0.025144797,0.019841678,0.0028447581,-0.0049878494,-0.025167329,0.016405258,-0.033873651,0.054860491,-0.00024628281,-0.0036819680,-0.011776405,0.0074700243,0.055548724,-0.058151014,0.050382517,-0.021110022,0.033754978,0.013080808,0.0021015003,-0.010309308,-0.025866412,0.0028378896,-0.047096245,0.0052202945,-0.015673835,0.019479753,0.0052426248,0.011365019,-0.0097224098]

## pure CLE output
si2=np.zeros(30, dtype=np.float64)
si2[:]=[3.4425000e-07,0.0019310000,0.035985000,0.13800000,0.31334999,0.57800001,0.95849997,1.4810001,2.1565001,2.9710002,3.8785000,4.7995000,5.6300001,6.2650003,6.6050000,6.6050000,6.2650003,5.6300001,4.7995000,3.8785000,2.9710002,2.1565001,1.4810001,0.95849997,0.57800001,0.31334999,0.13800000,0.035985000,0.0019310000,3.4425000e-07]
sq2=np.zeros(30, dtype=np.float64)
sq2[:]=[3.8010000e-08,0.00021319999,0.0039734999,0.015234999,0.034600001,0.063800000,0.10585000,0.16350000,0.23809999,0.32804999,0.42824998,0.52999997,0.62150002,0.69150001,0.72950000,0.72950000,0.69150001,0.62150002,0.52999997,0.42824998,0.32804999,0.23809999,0.16350000,0.10585000,0.063800000,0.034600001,0.015234999,0.0039734999,0.00021319999,3.8010000e-08]

su2=np.zeros(30, dtype=np.float64)
su2[:]=[4.0305000e-08,0.00022610000,0.0042130002,0.016154999,0.036685001,0.067649998,0.11220000,0.17340001,0.25244999,0.34784999,0.45410001,0.56200004,0.65899998,0.73350000,0.77349997,0.77349997,0.73350000,0.65899998,0.56200004,0.45410001,0.34784999,0.25244999,0.17340001,0.11220000,0.067649998,0.036685001,0.016154999,0.0042130002,0.00022610000,4.0305000e-08]

sv2=np.zeros(30, dtype=np.float64)
sv2[:]=[-1.1190000e-10,-4.3730000e-07,-6.5249997e-06,-2.1554999e-05,-4.3460001e-05,-7.1650000e-05,-0.00010585001,-0.00014399999,-0.00018160001,-0.00021160000,-0.00022599999,-0.00021750000,-0.00018225001,-0.00012165000,-4.2765001e-05,4.2765001e-05,0.00012165000,0.00018225001,0.00021750000,0.00022599999,0.00021160000,0.00018160001,0.00014399999,0.00010585001,7.1650000e-05,4.3460001e-05,2.1554999e-05,6.5249997e-06,4.3730000e-07,1.1190000e-10]

## alternative definition for the wavelength array
wvl=np.array([-3.81800,-2.66100,-2.13000,-1.83600,-1.63000,-1.45700,-1.29800,-1.14300,-0.990000,-0.837000,-0.685000,-0.533000,-0.380000,-0.228000,-0.0760000,0.0760000,0.228000,0.380000,0.533000,0.685000,0.837000,0.990000,1.14300,1.29800,1.45800,1.63000,1.83600,2.13100,2.66200,3.82100], dtype=np.float)

wvl=np.zeros(30, dtype=np.float64)
wvl[:]=[-3.81800,-2.66100,-2.13000,-1.83600,-1.63000,-1.45700,-1.29800,-1.14300,-0.990000,-0.837000,-0.685000,-0.533000,-0.380000,-0.228000,-0.0760000,0.0760000,0.228000,0.380000,0.533000,0.685000,0.837000,0.990000,1.14300,1.29800,1.45800,1.63000,1.83600,2.13100,2.66200,3.82100]

np.savetxt('cle_data.csv', (si1,sq1,su1,sv1,si2,sq2,su2,sv2,wvl), delimiter=',')


# define a level-2 header structure.
from astropy.io import fits

hdu = fits.PrimaryHDU(np.zeros(wvl.shape[0]))
hdu.header.set('Naxis3',wvl.shape[0])           ## assuming the wavelength bin size is naxis3 in the header
hdu.header.set('wdt_i',0)                       ## Instrumental width
hdu.header.set('l_speed',2.9979E+5)             ## speed of light [Km s^-1]  
hdu.header.set('kb',1.3806488E-23*1.e-6)        ## Boltzman constant [km^2 kg s^-2 K^-1] # Mostly SI converted to Km^2
hdu.header.set('e_mass',9.10938356e-31 )        ## Electron mass SI [Kg]
hdu.header.set('e_charge',1.602176634e-19)      ## Electron charge SI [C]
hdu.header.set('ion_temp',6.22)                 ## Ion temperature SI [K]
hdu.header.set('ion_mass',55.845*1.672621E-27)  ## Ion mass SI [Kg] ## for ion XIII this needs to be computed for all 
hdu.header.set('Line_w',10746.8)                ## Ion referential wavelength [A]

## s is an array of 5 rows each being in order: stokes I Q, U and V + the fifth array which records the wavelength position for the observation.
s=np.zeros([hdu.header['Naxis3'],5],dtype=np.float64)
s[:,0]=si2[:]
s[:,1]=sq2[:]
s[:,2]=su2[:]
s[:,3]=sv2[:]
s[:,4]=wvl[:]

class CDF_params:
	line_obs_center   = None  ## Observed line center wavelength
	shift_lmd         = None  ## Shifts in wavelength
	shift_vel         = None  ## Shifts in velocity 
	width             = None  ## Observed line widht
	intensity         = None  ## stokes I intensity in line core avg 2 pixels.
	background        = None  ## background counts in all 4 stokes profiles
	width_th          = None  ## Thermal broadening/width
	width_nth         = None  ## Non-thermal widths
	si                = None  ## Stokes I integrated
	sq                = None  ## Stokes Q integrated
	su                = None  ## Stokes U integrated	
	sv                = None  ## Stokes V integrated
	pl                = None  ## los degree of polarization 
	pv                = None  ## full vector degree of polarization
	xi                = None  ## azimuth angle
	blos_m            = None  ## Magnetograph Blos
	blos_p            = None  ## Plowman, 2014 Blos
	mask              = None  ## mask for fit consistency check
	tmp               = None
	cdf               = None
	s1                = None

	##define a function that computes the spectroscopic properties from one pixel.
	def __init__(self,s):
		import numpy as np
		import scipy.stats as sps
		## keep a record of the original stokes I
		s1=np.zeros([hdu.header['Naxis3'],5],dtype=np.float64)
		s1[:,:]=s[:,:]
		cdf=np.zeros(hdu.header['Naxis3'],dtype=np.float64)
		intensity=np.zeros(4,dtype=np.float64)    ## array to record the peak intensity levels in stokes I,Q,U,V, in that order.
		background=np.zeros(4,dtype=np.float64)    ## array to record the background levels in stokes I,Q,U,V, in that order.
		mask=np.zeros(4,dtype=np.int)       ## mask to test the statistical significance of the data: mask[0]=1 -> no signal above background in stokes I; mask[1]=1 -> distribution is not normal (skewed); mask[2]=1 line center (observed) not found; mask[3]=1=2 one or two fwhm edges were not found
		# compute the rough cdf distribution of the stokes I component to get a measure of the statistical noise.
		for i in range (0,hdu.header['Naxis3']):
			cdf[i]=s1[0:i+1,0].sum()         ## need to check if should start from 0 or not.

		cdf[:]=cdf[:]/cdf[-1]     ##norm the cdf so its easier to use and record the noise levels of the line along the statistical 0.05 and 0.95 of a normal distribution.  
		l0=np.where(np.abs(cdf-0.05) == np.abs(cdf-0.05).min() )[0][0] 
		r0=np.where(np.abs(cdf-0.95) == np.abs(cdf-0.95).min() )[0][-1]
		## remove the background noise from the profile. the 4 positions denote the IQUV measurements, in order.
		background[0]=(s1[0:l0+1,0].mean()+s1[r0:,0].mean())/2.
		background[1]=(s1[0:l0+1,1].mean()+s1[r0:,1].mean())/2.
		background[2]=(s1[0:l0+1,2].mean()+s1[r0:,2].mean())/2.
		background[3]=(s1[0:l0+1,3].mean()+s1[r0:,3].mean())/2.

		for i in range(0,4):
			s1[:,i]=s1[:,i]-background[i]
		##now fix the possible negative values resulted from the averaging
		for i in range (0,hdu.header['Naxis3']):
			if s1[i,0] < 0:
				s1[i,0]=0.001
		## now recompute the cdf for the noise-reduced data
		for i in range (0,hdu.header['Naxis3']):
			cdf[i]=s1[0:i+1,0].sum()         ## need to check if should start from 0 or not.
		cdf[:]=cdf[:]/cdf[-1]     ##norm the cdf so its easier to use and trecord the noise levels of the line along the statistical 0.05 and 0.95 of a normal distribution.  
		## check if there is reliable signal to fit and analyze  (1) we can use cdf to fit a line and in case it does fit well, there is no (reliable) stokes profile to recover!
		##the 1.4e-5 rest in correlation corresponds to a gaussion with a peak in intensity of <5% compared to the averge signal across the spectral range, or 8% higher from the minimum value. 
		co1=sps.pearsonr(s1[:,4],cdf)[0]
		if 1.-co1 < 1.4e-5:
			print(" No signal in pixel...")
			mask[0]=1

		## compute the center of the distribution, the line width, and the doppler shifts in wavelength and velocity.
	    ##interpolate all to a continuous domain.	## use a temp variable to store all the data
		tmp=np.zeros(9,dtype=np.float64)  ## a normal distribution 0.5 is the average, sigma=34.13; FWHM=2*sqrt(2*alog(2))*sigma, order is left fwhm [0] ,center [1],rightfwhm [2], left fwhm position [3], center position [4], right position [5] ; [6],[7],[8] just record the LEFT array index where positions were recorded
		tmp[0:3]=[0.5-2*np.sqrt(2*np.log(2))*34.13/2./100, 0.5, 0.5+2*np.sqrt(2*np.log(2))*34.13/2./100 ]
		tmp[0:3]=[0.16, 0.5, 0.84]
		for j in range (0,3):
			for i in range (0,hdu.header['Naxis3']):   
				if cdf[i] < tmp[j] < cdf[i+1]:                       ## find between which bins the distribution centre is.
					k1=(tmp[j]-cdf[i])/(cdf[i+1]-cdf[i])             ## find the normed difference to theoretical centre from the left bin.
					tmp[j+3]=s1[i,4]+k1*(s1[i+1,4]-s1[i,4])
					tmp[j+6]=i
					break
				elif i == hdu.header['Naxis3']-1:
					if j == 1:
						print("Line center not found")
						mask[2]=1
					if j == 0:
						print("FWHM left margin not found")
						mask[3]+=1
					if j == 2:
						print("FWHM right margin not found")
						mask[3]+=1
		line_obs_center=tmp[4]+hdu.header['line_w'] 
		shift_lmd=line_obs_center-hdu.header['line_w']         ## and the associated center position shift
		shift_vel=shift_lmd*3e5/hdu.header['Line_w']           ## shift converted to velocities; km*s^-1
		width=tmp[5]-tmp[3]
		## record the intensity of the central wavelength for each profile. Should be ~0 for V
		intensity[0]=(s[np.int(tmp[7]),0]+s[np.int(tmp[7])+1,0])/2.
		intensity[1]=(s[np.int(tmp[7]),1]+s[np.int(tmp[7])+1,1])/2.
		intensity[2]=(s[np.int(tmp[7]),2]+s[np.int(tmp[7])+1,2])/2.
		intensity[3]=(s[np.int(tmp[7]),3]+s[np.int(tmp[7])+1,3])/2.		
		## check if there is reliable signal to fit and analyze  (2) if the distribution is skewed, the inner part of it should not fit a line.
		co2=sps.pearsonr(s1[np.int(tmp[6]):np.int(tmp[8])+1,4],cdf[np.int(tmp[6]):np.int(tmp[8])+1])[0]	
		if 1.-co2 > 5e-3:
			print("Emission does not follow a normal distribution")
			mask[1]=1
		width_th=np.sqrt(4.*np.log(2)*hdu.header['kb']*(10.**hdu.header['ion_temp'])/hdu.header['ion_mass']) 
		width_nth=np.sqrt( ((((width**2)*(hdu.header['l_speed']**2))-((hdu.header['wdt_i']**2)*(hdu.header['line_w']**2)))/(4*np.log(2)*(hdu.header['line_w']**2)))-(hdu.header['kb']*(10.**hdu.header['ion_temp'])/hdu.header['ion_mass']))
		##Integrated Stokes profiles
		si=s1[np.int(tmp[6]):np.int(tmp[8]),0].sum()
		sq=s1[np.int(tmp[6]):np.int(tmp[8]),1].sum()
		su=s1[np.int(tmp[6]):np.int(tmp[8]),2].sum()
		sv=np.absolute(s1[np.int(tmp[6]):np.int(tmp[8]),3]).sum()

		## compute the polarization quantities
		pl=np.sqrt(sq**2+su**2)/si        ## for los degree of polarization
		pv=np.sqrt(sq**2+su**2+sv**2)/si  ## technically for full vector degree of polarization...
		xi=np.arctan(su/sq)/2.              ## the azimuth angle, (vertical to horizontal (e.g. N-E) projection)

		##compute the magnetic fields
		blos_p=(-8*np.pi*sv*hdu.header['e_mass']*hdu.header['l_speed'])*hdu.header['line_w']*width_th/(3*(si+pl)*hdu.header['e_charge'])
		blos_m= 10e3*sv/(-4.6686e-7*(hdu.header['line_w']**2)*1.5*10e3*((s1[np.int(tmp[6]):np.int(tmp[8]),0]-s1[np.int(tmp[6])+1:np.int(tmp[8])+1,0]).mean())/(s1[1,4]-s1[0,4]))  ## 1.5 is the lande g factor
		## put the data in a nice class so it mimics IDL structures
		self.line_obs_center = line_obs_center
		self.shift_lmd       = shift_lmd 
		self.shift_vel       = shift_vel
		self.width           = width
		self.intensity       = intensity
		self.background      = background
		self.width_th        = width_th
		self.width_nth       = width_nth
		self.si              = si  
		self.sq              = sq
		self.su              = su
		self.sv              = sv
		self.pl              = pl
		self.pv              = pv 
		self.xi              = xi
		self.blos_m          = blos_m
		self.blos_p          = blos_p
		self.mask            = mask
		self.tmp             = tmp
		self.cdf             = cdf
		self.s1              = s1

pix=CDF_params(s) ##just call the class for each pixel of interest and record the output in arrays.




##end
#for i in range(0,3):
#	for j in range(0,10):
#		print(i,j)
#		if j < bla[i] < j+1:
#			print("bingo",bla[i], i,j )  
#			break