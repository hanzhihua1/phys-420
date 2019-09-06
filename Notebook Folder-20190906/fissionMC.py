import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import odeint
import matplotlib.pyplot as plt          


class Neutron:
    '''Class to create and track neutron through a geometry.  We start a neutron at a position x, and our code tracks the
    neutron to the first interaction or until it leaves our volume'''
    #set up class variables- that is variables that are shared by all neutrons
    stopAfter=100000  # number of neutrons before stopping
    maxTime=500 #ns  time after which to stop
    timeBin=1 # ns  we group neutrons into bins of this width
    numberOfBins=int(maxTime/timeBin)
    mass=939.57e6 #mc**2, eV for neutron (we use eV to match NNDC)
    c=30 # cm/ns (speed of light)
    maxScatters=1000 #number of elastic scatters before we stop
    output=True
    list=[]
    for i in range(0,numberOfBins):
        list.append([])  # set up an empty list in each bin
    nFirst=0  #number of neutrons that needed to come from somewhere else(ie. not a chain reaction)
    nChain=0  # number of chain reaction neutrons

    
    def Reset(self):
        '''Resets the classes, so we can run the MC several times'''
        Neutron.nFirst=0  #number of neutrons that needed to come from somewhere else(ie. not a chain reaction)
        Neutron.nChain=0  # number of chain reaction neutrons
        for i in range(0,Neutron.numberOfBins):
             del Neutron.list[i][:] # set up an empty list in each bin
 
    def TrackNeutrons(self,nNeutrons):
        Neutron.stopAfter=nNeutrons
 #       import pdb; pdb.set_trace()
        #now we loop through the list.  If any bin is empty, we create a neutron in it at time binNumber*timeBin.  Creating a neutron
        #tracks the neutron- until it fissions or dies; the list entry corresponding to that entry is incremented to point at the parent neutron
        n=Neutron()
        n.FirstNeutron(0)
        for i in range(0,self.numberOfBins):
#            if len(self.list[i])==0 and Neutron.nChain==0:
#                n=Neutron()
#                n.FirstNeutron(i*self.timeBin)
            j=0
            while j<len(self.list[i]) and Neutron.nChain<nNeutrons:  #note: we want to do this even for the first neutron- because it could generate neutrons in this bin
                mother,e=self.list[i][j]
                n=Neutron()
                n.GenerateDaughter(mother,e)  #calculates # of daughters,times daughters die and increments the list
                self.list[i][j]=n  #replace mother with daughter 
                if Neutron.output:
                    print('daughter',n)
                j+=1
            if Neutron.nChain>=nNeutrons:
                break
    def SetEnergy(self,energy):
        self.energy=energy
        self.velocity=self.Velocity(energy)
    
    def Velocity(self,energy):
        return np.sqrt(2*self.energy/Neutron.mass)*Neutron.c
            
    def FirstNeutron(self,time):
        self.Reset()
        self.time0=time  #time of birth
        self.time=time #time during this part of the track
        self.position,self.direction,energy=geo.FirstNeutron()
        self.SetEnergy(energy)
        self.elastics=0 #number of elastic scatters
        Neutron.firstNeutron=self 
        Neutron.firstNeutronStartPosition=self.position
        self.TrackNeutron()
        if Neutron.output:
            print('mother:',self)
        
    def TrackNeutron(self):
        elastic=True
        while elastic:
            t, elastic,interactionOutput=geo.NextInteraction(self)
            self.position=self.position+t*self.velocity*self.direction
            self.time=self.time+t
            if elastic:  #for elastic scatters we just reset energy and direction and continue
                energy,ctheta=interactionOutput
                self.SetEnergy(energy)
                self.rotate(ctheta)
                self.elastics=self.elastics+1
                elastic = self.elastics<Neutron.maxScatters
            else:
                if len(interactionOutput)==0:
                    self.nDaughters='escape'
                else:    
                    self.nDaughters=interactionOutput
                    ibin=int(self.time/self.timeBin)
                    for j in self.nDaughters:
                        if j>0: # don't store 0 energy daughters (from 0 neutron fissions)
                            if ibin<self.numberOfBins :
                                self.list[ibin].append((self,j))  #we don't create the daughters, we store the mother and the daughter energy in a list.
        
        Neutron.nFirst=Neutron.nFirst+1
        
    def GenerateDaughter(self,mother,e):
        self.time0=mother.time
        self.time=self.time0
        self.position=mother.position  #starting position is end position of mother
        self.direction=geo.GetDirection()
        self.SetEnergy(e)
        Neutron.nChain=Neutron.nChain+1
        self.elastics=0 
        self.TrackNeutron()
        
    k=np.array([0,0,1])  #class variables that should never change
    i=np.array([1,0,0])

    def rotate(self,cosomega):
        '''Rotate the vector v around a random direction by an angle omega'''
        a1=np.sqrt((1+cosomega)/2)  #use half angle formula
        st=np.sqrt(1-a1**2)
        axis=np.cross(self.k,self.direction)
        if np.linalg.norm(axis)<=0.1:
            axis=np.cross(self.i,self.direction)
        axis=axis/np.linalg.norm(axis)  #unit vector
        b1,c1,d1=st*axis
    
        phi=np.pi*np.random.rand(1) #note that we left out the 2, since a is cos(phi/2)
        a2=np.cos(phi)
        st2=np.sin(phi)
        axis=self.direction/np.linalg.norm(self.direction)
        b2,c2,d2=st2*axis
    
        a=a1*a2-b1*b2-c1*c2-d1*d2
        b=a1*b2+b1*a2-c1*d2+d1*c2
        c=a1*c2+c1*a2-d1*b2+b1*d2
        d=a1*d2+d1*a2-b1*c2+c1*b2
    
        rotation=np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d),2*(b*d+a*c)],
                      [2*(b*c+a*d),a*a+c*c-b*b-d*d,2*(c*d-a*b)],
                      [2*(b*d-a*c),2*(c*d+a*b),a*a+d*d-b*b-c*c]]).reshape(3,3)
        return np.dot(rotation,self.direction)

    def __str__(self):
        if type(self.nDaughters) is list:
            d=np.array_str(np.array(self.nDaughters)*1e-6,precision=2)
        else:
            d=self.nDaughters
        return ('t=%4.1f,t0=%4.1f,e=%5.2f MeV'%(self.time,self.time0,self.energy*1e-6)+' '+d)
    
    
        
        
        
    
class Geometry:
    '''Geometry defines shapes and uses materials.  We start with a simple sphere of radius r_sphere'''
    
    throwN=1000  #number of items to throw at once
    position=throwN
    direction=throwN
     
    def ThrowPositions(self):  #throws positions within a sphere
        self.position=0
        radii=self.radius*np.random.rand(self.throwN)**(1/3)
        costheta=np.random.rand(self.throwN)*2-1
        sintheta=np.sqrt(1-costheta**2)
        phi=np.random.rand(self.throwN)*2*np.pi
        self.positions=np.array(radii*[sintheta*np.cos(phi),sintheta*np.sin(phi),costheta]).transpose()

    def ThrowDirections(self):  #throws isotropic directions
        self.direction=0
        costheta=np.random.rand(self.throwN)*2-1
        sintheta=np.sqrt(1-costheta**2)
        phi=np.random.rand(self.throwN)*2*np.pi
        self.directions=np.array([sintheta*np.cos(phi),sintheta*np.sin(phi),costheta]).transpose()
        
    def ThrowEnergies(self):  #throws energies according to fission energy distribution
        self.eIndex=0        
        self.energies=self.mat.throwFissionE(np.random.rand(self.throwN))

        
    def GetPosition(self):
        if self.position>=self.throwN:
            self.ThrowPositions()
        retval=self.positions[self.position]
        self.position+=1
        return retval
    
    def GetDirection(self):       
        if self.direction>=self.throwN:
            self.ThrowDirections()
        retval=self.directions[self.direction]
        self.direction=self.direction+1
        return retval

    def GetEnergy(self):
        if self.eIndex>=self.throwN:
            self.ThrowEnergies()
        retval=self.energies[self.eIndex]
        self.eIndex=self.eIndex+1
        return retval
    
    def FirstNeutron(self):
        return (self.GetPosition(),self.GetDirection(),self.GetEnergy())
        
    def NextInteraction(self,neutron):
        '''This routine checks for the next intersection with the sphere, or where the next nuclear interaction 
        will happen
        NextInteraction determines whether the neutron leaves the sphere, fissions, elastically scatters.  
        it returns with the time to the interaction and a tuple that contains information about outgoing neutrons
        for fissions interaction contains the number of daughters
        for geometrical events or when there are zero neutrons generated, interaction is empty'''
        
        #solve for the distance
        d=np.roots(np.array([1,2*np.dot(neutron.position,neutron.direction),np.dot(neutron.position,neutron.position)-self.radius**2]))
        dpos=1e10  #distance to next interaction
        for dd in d: #check if the interaction is a geometrical border crossing  
            if dd>0 and dd<dpos:
                dpos=dd
        if dpos==1e10:
            print('Error.  Tried finding roots when outside the sphere')
            
        #now check for neutron interactions
        elastic,dist, interactionOutput=self.mat.interaction(neutron.energy)  # returns the mean distance and the energies of the outgoing neutrons
        dd=np.random.exponential(dist)
        if dpos<dd:   #neutron escapes
            elastic=False
            interactionOutput=()
        else:
            dpos=dd
        time=dpos/neutron.velocity
        return (time,elastic,interactionOutput)
    
    def Daughter(self,energy):
        return self.mat.Daughter(energy)

    
    def __init__(self,label):
        if 'sphere' in label:  #this allows a derived class to skip these lines
            print('Initializing spherical geometry')
            self.radius=10 #cm radius of sphere

        self.mat=Material('U235')
        self.eIndex=self.throwN+1  #so that we throw new positions when we start
        self.direction=self.throwN+1
        self.position=self.throwN+1
            
    

class Material:
    '''Material contains the information about cross sections, densities, etc.  It is able to return the distance to 
    the next physics interaction and the type of interactions, and energies of elastic scattered and fission neutrons
    
    To use the class one uses three commands
    variable=Material('U235)
    fractions(energy)   to get the interaction length and the fractions of various types of scatters
    interaction(energy) to get the parameters of the interaction
    
    '''
    
        
    def fractions(self,energy):
        '''Returns the interaction length corresponding to the total cross section,
        the elastic fraction, and the elastic+fission fraction'''
        total=self.totalXS(energy)
        elastic=self.elasticXS(energy)
        fission=self.fissionXS(energy)
        interactionLength=1.0/(self.interactionLambda*self.totalXS(energy))
        return(interactionLength,elastic/total,(elastic+fission)/total)
        
    def interaction(self,energy):
        '''Calculates the probability of elastic scattering or fission.  
        Energy is the energy of the neutron in this material.
        For elastic scattering returns True, interactionLength,(cosTheta, Center of Mass velocity/velocity ).  
        Here cos theta is the cosine of the scattering angle.
        For fission returns False, interactionLength,(energy of daughters)'''
        interactionLength,elasticFraction,fissionFraction=self.fractions(energy)
        rn=np.random.rand()
        #go through the different scattering mechanisms
        if rn<elasticFraction:
            elastic=True
            retval=self.elasticCosTheta(energy)
            self.cthetaHist.append([energy,np.asscalar(retval[-1])])
        elif rn<fissionFraction:  #fissions
            elastic=False
            nDaughters=np.random.poisson(self.fissionNumber(energy))
            self.nDaughtersHist.append([energy,nDaughters])
            retval=[]
            for ix in np.arange(0,nDaughters):
                #throwFissionE returns an array
                retval.append(np.asscalar(self.throwFissionE(np.random.rand())))
                self.daughterHist.append(retval[-1])
            if nDaughters==0:  #fission without daughters
                retval.append(0)  # allows us to distinguish 0 fission captures from escapes
        else:  #treat inelastic scatters like elastic scatters, except don't histogram
            elastic=True
            retval=self.elasticCosTheta(energy)
        self.distanceHist.append([energy,interactionLength])    
        return (elastic,interactionLength,retval)    
    
    def U235Initialize(self):
        #start by reading in the angular distribution legendre polynomials
        self.energyList=[]
        self.throwAngles=[]
        self.legendreList=[]
        k=0
        lfactor=np.arange(0,30)+0.5 #(2L+1)/2 normalization of legendre polynomials
        with open('angularDistU235Elastic.txt', 'r') as f:
            for line in f:
                if 'Angular' in line:
                    pass
                elif 'eV  Coefficients' in line:
                    pass
                elif '-----' in line:
                    pass
                elif len(line)==1:
                    pass
                else:
                    if not line[0:10].isspace():
                        energy=float(line[0:10]+'e'+line[10:12])
                        self.energyList.append(energy)
                        if len(self.energyList)>1:
                            leg=np.polynomial.legendre.Legendre(np.array(l)*lfactor[:len(l)])
                            self.legendreList.append(leg)
                            leg0=leg.integ()  # integral from zero
                            legendreIntegral=leg0-leg0(-1) #integral from -1
                            xv=legendreIntegral.linspace(200)
                            intLeg=interp1d(xv[1],xv[0],kind='cubic')
                            self.throwAngles.append(intLeg)
                            k=k+1
                        l=[1.0]
                    for i in range(13,79,11):
                        if not(line[i:i+9].isspace()):
                            coeff=float(line[i:i+9]+'e'+line[i+9:i+11])
                            l.append(coeff)
        #append the last values for throwAngles and legendreList
        leg=np.polynomial.legendre.Legendre(np.array(l)*lfactor[:len(l)])
        self.legendreList.append(leg)
        leg0=leg.integ()  # integral from zero
        legendreIntegral=leg0-leg0(-1) #integral from -1
        xv=legendreIntegral.linspace(200)
        intLeg=interp1d(xv[1],xv[0],kind='cubic')
        self.throwAngles.append(intLeg)
                            
        #now read in the elastic cross section
        elasticXdata=np.genfromtxt('ElasticCrossSectionU235.txt',skip_header=1,delimiter=',')
        #the following two lines elements where two subsequent entries have the same x value.  Spline can't interpolate those!
        dd=np.where(elasticXdata[1:,0]<=elasticXdata[:-1,0])[0]
        elasticXdata=np.delete(elasticXdata,dd,axis=0)
        self.elasticXS=interp1d(elasticXdata[:,0].reshape(-1),elasticXdata[:,1].reshape(-1),
                                bounds_error=False,fill_value=0.0)
         
        #total cross section
        totalXdata=np.genfromtxt('TotalCrossSectionU235.txt',delimiter=',',skip_header=1)
        self.totalXS=interp1d(totalXdata[:,0].reshape(-1),totalXdata[:,1].reshape(-1),kind='cubic',
                              bounds_error=False, fill_value=0.001)

        #fission cross section
        fissionXdata=np.genfromtxt('FissionCrossSectionU235.txt',skip_header=1,delimiter=',')
        #the following two lines elements where two subsequent entries have the same x value.  Spline can't interpolate those!
        dd=np.where(fissionXdata[1:,0]<=fissionXdata[:-1,0])[0]
        fissionXdata=np.delete(fissionXdata,dd,axis=0)
        self.fissionXS=interp1d(fissionXdata[:,0].reshape(-1),fissionXdata[:,1].reshape(-1),bounds_error=False,
                                fill_value=0.0)

        #fission Energy
        fissionEnergy=np.genfromtxt('FissionOutgoingEnergySpectraU235.txt',skip_header=5)
        self.fissionE=interp1d(fissionEnergy[:,0].reshape(-1),fissionEnergy[:,1].reshape(-1))

        #now number of fission daughters
        fissionN=np.genfromtxt('NumberNeutronsU235.txt',skip_header=1,delimiter=',')
        self.fissionNumber=interp1d(fissionN[:,0].reshape(-1),fissionN[:,1].reshape(-1))
                                                                      

        #interpolation function to thow the fission energy                                                                   `
        integrand = lambda y, x:self.fissionE(x)
        y=0
        en=np.linspace(1,2.5e7,1000)
        dd=odeint(integrand,y,en,rtol=1e-6,atol=0.1)  # default tolerance never converges.
        dd=dd/dd[-1]  #normalize to 1

        #now we need to get rid of many duplicate bins
        xlist=[]
        ylist=[]
        j=1
        i=0
#        import pdb; pdb.set_trace()
        while i+j <len(dd):
            if i==0 or dd[i+j]/dd[i] >1.002:
                xlist.append(en[i])
                ylist.append(dd[i])
                i=i+j
            else:
                j=j+1
            xlist.append(en[i])
            ylist.append(dd[i])
        xlist.append(en[-1])
        ylist.append(dd[-1])
        xl=np.array(xlist)
        yl=np.array(ylist).reshape(-1) 
        self.throwFissionE=interp1d(yl,xl) #note that it is easy to invert an table of values by interchanging y,x
        
    def U235ElasticCosTheta(self,energy):
        k=np.searchsorted(self.energyList,energy)  #find  bin number
        if k<0 or k>=len(self.energyList):
#            print('energy overflow',energy,k,len(self.energyList),len(self.throwAngles))
            self.overFlowCounter+=1
            k=len(self.energyList)-1
        if k+1<len(self.energyList) and self.energyList[k+1]-energy<energy-self.energyList[k]:  
            k=k+1
        ctheta=self.throwAngles[k](np.random.rand())
        stheta=np.sqrt(1-ctheta**2)
        vt=Neutron.mass**2+self.mass**2+2*Neutron.mass*self.mass*ctheta
        eFinal=energy*vt/(Neutron.mass+self.mass)**2
        ct=(self.mass*ctheta+Neutron.mass)/np.sqrt(vt)
        return(eFinal,ct)


    def __init__(self, name):
        '''Initialization for each isotope.  At a minimum need to define self.mass, 
        the elastic scattering cross section and angle, self.elasticXS(E) and self.elasticCosTheta
        the fission cross section and number of neutrons, self.fissionXS(E)
        the fission energy distribution self.fissionE(E), self.fissionE(E)'''
        if 'U235' in name:
            print('Starting to initialize U235')
            self.mass=235.04*931.5e6  #mass mc**2 in eV
            self.density=19.1*235/238
            self.interactionLambda=1e-24*6.02e23/235.04*self.density #cm**2 *NA/mass per mole *gm/cm**3=1/cm
            self.U235Initialize()
            self.elasticCosTheta=self.U235ElasticCosTheta
            self.overFlowCounter=0
        

#initialize lists of variables to be histogrammed
        self.cthetaHist=[]
        self.nDaughtersHist=[]
        self.daughterHist=[]
        self.distanceHist=[]
        
  
    def validate(self):
        capture=[]
        erange=np.linspace(1,1e6,1000)
        for e in erange:
            capture.append(self.fractions(e))
        capture=np.array(capture).transpose()
        plt.figure()
        plt.plot(erange,1.0/capture[0],'r',label='fraction calculation')
        plt.plot(erange,self.totalXS(erange), label='Interpolation')
        plt.plot(erange,self.totalXS(erange)*capture[0],'g',label='ratio calculation')
        plt.legend()
        plt.title('Total Cross sections')
        plt.xlabel('Energy of neutron(eV)')
        plt.ylabel('Cross Section, Barns')        
        plt.yscale('log')
        plt.figure()
        plt.hist(self.totalXS(erange)*capture[0],bins=100,range=(20.440,20.442))
        
        plt.figure()
        plt.plot(erange,capture[1]/capture[0],'r',label='fraction calculation')
        plt.plot(erange,self.elasticXS(erange), label='Interpolation')
        plt.plot(erange,self.elasticXS(erange)*capture[0]/capture[1],'g',label='ratio calculation')
        plt.legend()
        plt.title('Elastic Cross sections')
        plt.xlabel('Energy of neutron(eV)')
        plt.ylabel('Cross Section, Barns')        
        plt.yscale('log')

        plt.figure()  #for fission we need to subtract because the fraction includes elastic
        plt.plot(erange,(capture[2]-capture[1])/capture[0],'r',label='fission calculation')
        plt.plot(erange,self.fissionXS(erange), label='Interpolation')
        plt.plot(erange,self.fissionXS(erange)*capture[0]/(capture[2]-capture[1]),'g',label='ratio calculation')
        plt.legend()
        plt.title('Fission Cross sections')
        plt.xlabel('Energy of neutron(eV)')
        plt.ylabel('Cross Section, Barns')        
        plt.yscale('log')
        
        #For interaction we will look at the number of fissions compared to the number of elastic scatters, 
        #the energy distribution of the fissions, and the angular distribution of the scatters.
    
        erange2=np.linspace(1,1e6,10)
        for e in erange2:
            for i in range(0,10000):
                elastic, length, data=self.interaction(e)
 
# now the actual work!
 
#import pdb; pdb.set_trace()
               
geo=Geometry('sphere')
n=Neutron()
while Neutron.nChain==0:
    n.TrackNeutrons(10000)
    print('Neutrons in chain:',Neutron.nChain)



