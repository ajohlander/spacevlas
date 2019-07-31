import numpy as np

def get_energy_spectrogram(f, v, dv, Vedges, Nbins=None):   

    print("inside get_energy spectrogram function")

    nE = Vedges.size-1;


    # array of energy and velocity values in data (long vectors, most memory usage here)
    # everything is non-relativistic
    Var = np.sqrt(np.sum(np.power(v,2),1))
    #Ear = .5*M*np.power(Var,2)/qe # in [eV]
    #print(Ear.shape)

    # if not in input, Nbins can be approximated
    if Nbins==None:
      print("Vedges = ",str(Vedges))
      EbinVolume = 4*np.pi/3*(np.power(Vedges[1:],3)-np.power(Vedges[0:-1],3))
      print("EbinVolume = "+str(EbinVolume))
      # number of bins per energy channel (no rounding needed)
      Nbins = EbinVolume/np.power(dv,3);


    # initialize spectrogram array
    fSphMean = np.zeros(nE)
    #print(fSphMean.shape)
    
    # loop through energy bins
    for iE in range(0,nE):
        #print(iE)
        #vx, vy, vz = sph2cart(phi,th,Var[iE])
      
        #print(f[np.bitwise_and(Var>Vedges[iE], Var<Vedges[iE+1])].shape)
        fSphMean[iE] = np.sum(f[np.bitwise_and(Var>Vedges[iE], Var<Vedges[iE+1])])/Nbins[iE]

        # # find bins in vlasiator velocity space where the points belong
        # idvx = np.digitize(vx,vxe)

    return fSphMean


def get_Nbins(vg,Vedges):
   #print(vg[0])
   nv = vg.size
   #print(nv)
   
   Vmat = np.zeros([nv**3,3])

   # construct matrix of velocity norms
   Vmat = np.linalg.norm(np.meshgrid(vg,vg,vg),axis=0)

   # number of energy channels
   nE = Vedges.size-1
   # loop to get number of bins in each channel
   Nbins = np.zeros(nE)
   for iE in range(0,nE):
      #print(iE)
      Nbins[iE] = np.sum(np.bitwise_and(Vmat>Vedges[iE], Vmat<Vedges[iE+1]))

   return Nbins*8


def get_reduced_dist(f,v,dv,Vedges):
  
  print(Rxyz)
  nvg = Vedges.size-1
  #dv = np.median(np.diff(v))
  dvg = np.diff(Vedges)

  # really strange behavior of meshgrid (flips order of x and y)
  #vmy,vmx,vmz = np.meshgrid(vc,vc,vc)

  x0 = Rxyz[:,0]
  x1 = Rxyz[:,1]
  x2 = Rxyz[:,2]

  #vc0 = vmx*x0[0]+vmy*x0[1]+vmz*x0[2]
  vc0 = v.dot(x0)
  vc1 = v.dot(x1)
  vc2 = v.dot(x2)
  #vc1 = vmx*x1[0]+vmy*x1[1]+vmz*x1[2]
  #vc2 = vmx*x2[0]+vmy*x2[1]+vmz*x2[2]

  #print(vc0.shape)


  Fg = np.zeros((3,nvg))
  for ivg in range(nvg):
    #print(vc0>vge[ivg])
    #print(np.bitwise_and(vc0>vge[ivg],vc0<=vge[ivg+1]))
    Fg[0,ivg] = np.sum(f[np.bitwise_and(vc0>Vedges[ivg],vc0<=Vedges[ivg+1])])*dv**3/dvg[ivg]
    Fg[1,ivg] = np.sum(f[np.bitwise_and(vc1>Vedges[ivg],vc1<=Vedges[ivg+1])])*dv**3/dvg[ivg]
    Fg[2,ivg] = np.sum(f[np.bitwise_and(vc2>Vedges[ivg],vc2<=Vedges[ivg+1])])*dv**3/dvg[ivg]

  return Fg




def get_reduced_dist2(f,v,dv,Vedges,Rxyz=np.array([[1,0,0],[0,1,0],[0,0,1]])):
	
	print(Rxyz)
	nvg = Vedges.size-1
	#dv = np.median(np.diff(v))
	dvg = np.diff(Vedges)

	# really strange behavior of meshgrid (flips order of x and y)
	#vmy,vmx,vmz = np.meshgrid(vc,vc,vc)

	x0 = Rxyz[:,0]
	x1 = Rxyz[:,1]
	x2 = Rxyz[:,2]

	#vc0 = vmx*x0[0]+vmy*x0[1]+vmz*x0[2]
	vc0 = v.dot(x0)
	vc1 = v.dot(x1)
	vc2 = v.dot(x2)
	#vc1 = vmx*x1[0]+vmy*x1[1]+vmz*x1[2]
	#vc2 = vmx*x2[0]+vmy*x2[1]+vmz*x2[2]

	#print(vc0.shape)


	Fg = np.zeros((3,nvg))
	for ivg in range(nvg):
		#print(vc0>vge[ivg])
		#print(np.bitwise_and(vc0>vge[ivg],vc0<=vge[ivg+1]))
		Fg[0,ivg] = np.sum(f[np.bitwise_and(vc0>Vedges[ivg],vc0<=Vedges[ivg+1])])*dv**3/dvg[ivg]
		Fg[1,ivg] = np.sum(f[np.bitwise_and(vc1>Vedges[ivg],vc1<=Vedges[ivg+1])])*dv**3/dvg[ivg]
		Fg[2,ivg] = np.sum(f[np.bitwise_and(vc2>Vedges[ivg],vc2<=Vedges[ivg+1])])*dv**3/dvg[ivg]

	return Fg
