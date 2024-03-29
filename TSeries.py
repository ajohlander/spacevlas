import numpy as np
import spaceunits as su 
import vlasdist as vd
import pytools as pt

def get_parameter(vlsvReader, cid, parstr, pop=None, Eedges=None, Nbins=None):
	# return (in order):
	# B [nt]
	# V [km/s]
	# n [cm^-3]
	# f (spherical mean) [s^3 m^-6]
	# E (energy) [eV]

	u = su.Units

	## Get "regular" parameters
	if parstr=='B' or parstr=='V' or parstr=='rho':
		pval = vlsvReader.read_variable(parstr,cellids=cid)
		
		if parstr=='B':
			return pval*1e9
		if parstr=='V':
			return pval*1e-3
		if parstr=='rho':
			return pval*1e-6

	elif parstr=='ft':
		## Handle ion populations
		# TODO
		M = u.mp

		## Get distribution function
		# cell ids where some psd is stored   
		velcells = vlsvReader.read_velocity_cells(cid, pop=pop)
		# velocity of those cells (nVx3 np array) [m/s]
		V = vlsvReader.get_velocity_cell_coordinates(velcells.keys(),pop=pop)
		# actual phase space density [s^3 m^-6]
		f = np.asarray(zip(*velcells.items())[1])

		# get limit of velocity bin EDGES in [m/s]
		[vxmin, vymin, vzmin, vxmax, vymax, vzmax] = vlsvReader.get_velocity_mesh_extent(pop=pop)
		# number of velspace cell blocks (every block is 4x4x4 cells)                                                                                                                                                  
		[vxsize, vysize, vzsize] = vlsvReader.get_velocity_mesh_size(pop=pop)
		# number of velocity space cells
		nV=vxsize*4
		# delta V [m/s]
		dv = (vxmax-vxmin)/nV
		# velocity cell CENTERS [m/s]
		vc = np.linspace(vxmin+dv/2,vxmax-dv/2,nV)


		## set energy channels for energy spectrogram if not in input
		if Eedges==None:
			# number of energy channels
			nE = 32.
			# energy limits for the energy spectrogram in [eV]
			Emin = 1
			E1 = 30.
			Emax = 50.*1e3
			# logarithmic energy table with one large energy bin from 0
			Eedges = np.append([0],np.logspace(np.log10(E1),np.log10(Emax),nE))

		# energy diffs
		dE = np.diff(Eedges)
		# energy centers
		Ec = Eedges[0:-1]+dE/2
		# velocity edges
		Vedges = np.sqrt(2*Eedges*u.e/M)
		# get Nbins
		if Nbins==None:
			Nbins = vd.get_Nbins(vc[nV/2:-1],Vedges)


		## get spherical mean of the distribution function
		ftval = vd.get_energy_spectrogram(f,V,dv,Vedges,Nbins=Nbins)
		return ftval, Ec, Eedges, Nbins

	else:
		print('unkown parameter')
		return None


def get_TSeries(filePath, runCode, tid, SC_coord):
	# filePath is NOT the full path
	# should have a list of outputs in input
	# TSeries should be Objects of a class!!

	# get full file path
	fullFilePath = filePath+runCode+"/"
	


	## other things
	# time resolution
	nt = tid.size
	dt = 0.5;
	timeArray = np.array(tid*dt);


	# number of points
	nsc = SC_coord.size/3

	# initiate arrays 
	#Bdata = np.zeros([nt,3])
	#Vdata = np.zeros([nt,3])
	#ndata = np.zeros(nt)
	#ftdata = np.zeros([nt,32]) # by definition

	# data is stored in dictionaries
	Bdata = {}
	Vdata = {}
	ndata = {}
	ftdata = {}

	for isc in range(nsc):
		Bdata['Bdata'+str(isc)] = np.zeros([nt,3])
		Vdata['Vdata'+str(isc)] = np.zeros([nt,3])
		ndata['ndata'+str(isc)] = np.zeros(nt)
		ftdata['ftdata'+str(isc)] = np.zeros([nt,32])

	print(ftdata)


	# First time recalculate Nbins
	Nbins = None

	# Get parameters for each time step
	startTimeInd = tid[0]
	for it in tid:
		# get FULL individual file name
		fileName = filePath+runCode+"/bulk/bulk."+str(tid[it-startTimeInd]).rjust(7,"0")+".vlsv"
		
		# get vlsvfile object
		vlsvReader = pt.vlsvfile.VlsvReader(fileName)

		# loop through spacecraft position
		for isc in range(nsc):

			# get cell id of desired cell
			if nsc == 1:
				SC_coordTemp = SC_coord
			else:
				SC_coordTemp = SC_coord[isc,:]

			print(SC_coordTemp)

			cid = vlsvReader.get_cellid(SC_coordTemp)

			# get parameters
			Bval = get_parameter(vlsvReader,cid,'B')
			Vval = get_parameter(vlsvReader,cid,'V')
			nval = get_parameter(vlsvReader,cid,'rho')
			ftval, Eval, Eedges, Nbins = get_parameter(vlsvReader,cid,'ft',pop="avgs",Nbins=Nbins)

			# put in arrays
			Bdata['Bdata'+str(isc)][it-startTimeInd,:] = Bval
			Vdata['Vdata'+str(isc)][it-startTimeInd,:] = Vval
			ndata['ndata'+str(isc)][it-startTimeInd] = nval
			ftdata['ftdata'+str(isc)][it-startTimeInd,:] = ftval

	# construct TSeries objects
	# they are also in dictionaries if more than one sc
	if nsc > 1:
		B = {}
		V = {}
		n = {}
		ft = {}

		for isc in range(nsc):
			B[str(isc)] = vector(timeArray,Bdata["Bdata"+str(isc)],name="B",unit='nT',position=SC_coord[isc,:])
			V[str(isc)] = vector(timeArray,Vdata["Vdata"+str(isc)],name="V",unit='km/s',position=SC_coord[isc,:])
			n[str(isc)] = scalar(timeArray,ndata["ndata"+str(isc)],unit='cm^-3',position=SC_coord[isc,:])
			ft[str(isc)] = spec(timeArray,ftdata["ftdata"+str(isc)],Eval,unit='s^3/m^6',position=SC_coord[isc,:],ancillary={'Eedges':Eedges})
	else:
		B = vector(timeArray,Bdata['Bdata0'],name="B",unit='nT',position=SC_coord)
		V = vector(timeArray,Vdata['Vdata0'],name="V",unit='km/s',position=SC_coord)
		n = scalar(timeArray,ndata['ndata0'],name="n",unit='cm^-3',position=SC_coord)
		ft = spec(timeArray,ftdata['ftdata0'],Eval,name="ft",unit='s^3/m^6',position=SC_coord)

	return B,V,n,ft





def scalar(time, data, unit=None, name="",position=None):

	dataShape = np.array(data.shape)
	
	# number of time steps
	nt = time.size

	if dataShape[0] != nt: 
		print('size mismatch, figure out how to throw exception')

	if dataShape.shape[0] == 1:
		return TSeries(time,data,'scalar',name=name,unit=unit,position=position)
	else:
		return None



def vector(time, data, unit=None, name="",position=None):

	dataShape = np.array(data.shape)
	
	# number of time steps
	nt = time.size

	if dataShape[0] != nt: 
		print('size mismatch, figure out how to throw exception')

	if dataShape.shape[0] == 2 and dataShape[1] == 3:
		return TSeries(time,data,'vector',name=name,unit=unit,position=position)
	else:
		return None



def tensor(time, data, unit=None, name="",position=None):

	# todo
	return None


def spec(time, data, depend, unit=None, dependUnit=None, name="",position=None,ancillary={}):

	dataShape = np.array(data.shape)
	
	# number of time steps
	nt = time.size

	if dataShape[0] != nt: 
		print('size mismatch, figure out how to throw exception')

	# add more checks here
	if dataShape.shape[0] == 2 and dataShape[1] > 3:
		return TSeries(time,data,'spec',name=name,unit=unit,depend=depend,dependUnit=dependUnit,position=position,ancillary=ancillary)
	else:
		print('error in spec:'+str(dataShape.shape[0])+','+str(dataShape[1]))
		print(time)
		print(data)




class TSeries(object):


	def __init__(self, time, data, dataType, unit=None, depend=None, dependUnit=None, name="",position=None,ancillary={}):
		# class will not complain if user misbehaves

		# number of time steps
		nt = time.size

		# determine type
		dataShape = np.array(data.shape)

		if dataShape[0] != nt: 
			print('size mismatch, figure out how to throw exception')


		self.time = time
		self.data = data
		self.type = dataType
		self.unit = unit
		self.depend = depend
		self.name = name
		self.position = position
		self.ancillary = ancillary

		# todo
		self.tensorOrder = None		
		self.cid = None


	def __add__(self,obj):
		# todo
		return None


	def abs(self,obj):
		# todo
		return None


	# and so on


	## Plot function(s)
	def plot(self,h):

		## split into cases

		if self.type == "scalar":
			h.plot(self.time,self.data,"-",linewidth=2,color="k")
		
		if self.type == "vector":
			h.plot(self.time,self.data[:,0],"-",linewidth=2,color="k")
			h.plot(self.time,self.data[:,1],"-",linewidth=2,color="b")
			h.plot(self.time,self.data[:,2],"-",linewidth=2,color="r")
		if self.type == "spec":
			cdata = self.data.transpose()
			cdata[cdata==0] = np.nan
			h.pcolor(self.time,self.depend,np.log10(cdata))


		# should also print units, must look for "^"s
		h.set_ylabel(self.name+" ["+"]")






	# def __init__(self,time,data):
		
	# 	# number of time steps
	# 	nt = time.size

	# 	# determine type
	# 	dataShape = np.array(data.shape)

	# 	if dataShape[0] != nt: 
	# 		print('size mismatch, figure out how to throw exception')

	# 	if dataShape.shape[0] == 1: # scalar (wow)
	# 		self.tensorOrder = 0
	# 		self.type = 'scalar'
	# 	elif dataShape[1] == 3:
	# 		self.tensorOrder = 1
	# 		self.type = 'vector'
	# 	elif dataShape[1] == 3 & dataShape[2] == 3:
	# 		self.tensorOrder = 2
	# 		self.type = 'tensor'
	# 	elif dataShape.shape[0] == 2 & dataShape[1] > 3:
	# 		self.tensorOrder = None
	# 		self.type = 'spec'
	# 	else:
	# 		print('unknown data shape, figure out how to throw exception')



		# self.time = time
		# self.data = data

		# self.position = None
		# self.cid = None













