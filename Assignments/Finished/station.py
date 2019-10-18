import subprocess
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
import numpy as np
import csv

class station:
    def __init__(self,n, lat, longitude, fyear, lyear, climateid,stationid):
        self.name=n
        self.latitude=float(lat)
        self.longitude=float(longitude)
        self.firstYear=int(fyear)
        self.lastYear=int(lyear)
        self.climateID=climateid
        self.stationID=int(stationid)

    def readYear(self,year):
        if self.firstYear<=year<=self.lastYear:
            str='http://climate.weather.gc.ca/climate_data/bulk_data_e.html?format=csv&stationID={0}&Year={1}&Month=12&Day=14&timeframe=2&submit=Download+Data'.format(self.stationID,year)
            print('str=',str)
            #        import pdb
            #        pdb.set_trace()
            subprocess.run(['wget','--output-document=test.out','--quiet',str])     
            i=0
            with open('test.out','r') as fdata:
                for lines in fdata:
                    if 'Date/Time' in lines:
                        break
                    else:
                        i=i+1
            with open('test.out','r') as fdata:
                for j in range(i):
                    fdata.readline()
                datareader=csv.DictReader(fdata)
                for row in datareader:
                    try:
                        t1=float(row['Max Temp (°C)'])
                        t2=float(row['Min Temp (°C)'])
                        y=int(row['Year'])
                        m=int(row['Month'])
                        d=int(row['Day'])
                        try:
                            pos=self.dates.index(dt.date(y,m,d))
                            self.maxT[pos]=t1
                            self.minT[pos]=t2
                            print('Replacing value %s %d %d %d %f %f'%(self.name, y,m,d,t1,t2))
                        except ValueError:
                            self.dates.append(dt.date(y,m,d))
                            self.maxT=np.append(self.maxT,t1)
                            self.minT=np.append(self.minT,t2)
                            print('New value %s %d %d %d %f %f'%(self.name, y,m,d,t1,t2))

                    except ValueError:
                        continue
        return

    def readData(self):
        self.dates=[]
        self.minT=[]
        self.maxT=[]
        print(self)
        for year in range(self.firstYear, self.lastYear+1):
            self.readYear(year)

    def plot(self):
        plt.figure()
        plt.title(self.name)
        plt.gca().xaxis.set_major_formatter(mdate.DateFormatter('%m/%d/%Y'))
        plt.gca().xaxis_date()
        plt.plot(self.dates,self.maxT)
        plt.plot(self.dates,self.minT)
        plt.gcf().autofmt_xdate()
        
    def FillDayOfYear(self):
#go through and make a Day Of Year field
        self.doy=[]
        self.year=[]
        for d in self.dates:
            delta= int((d-dt.date(d.year,1,1)).days)
            if delta>364:
                delta=365
            self.doy.append(delta)
            self.year.append(d.year)
