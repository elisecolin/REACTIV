#     -*- coding: utf-8 -*-
#     « Copyright (c) 2018, Elise Koeniguer (Onera) » 
#     This file is part of REACTIV .
# 
#     REACTIV is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     REACTIV is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with REACTIV in the file copying.txt.  
#     If not, see <http://www.gnu.org/licenses/gpl.txt>.
#
# ------------------------------------------------------------------------
#
#     REACTIV est un logiciel libre ; vous pouvez le redistribuer ou le
#     modifier suivant les termes de la GNU General Public License telle
#     que publiée par la Free Software Foundation ; soit la version 3 de la
#     licence, soit (à votre gré) toute version ultérieure.
# 
#     REACTIV est distribué dans l'espoir qu'il sera utile, mais SANS
#     AUCUNE GARANTIE ; sans même la garantie tacite de QUALITÉ MARCHANDE
#     ou d'ADÉQUATION à UN BUT PARTICULIER. Consultez la GNU General Public
#     License pour plus de détails.
# 
#     Vous devez avoir reçu une copie de la GNU General Public License en
#     même temps que Gefolki dans le fichier copying.txt ; si ce n'est pas
#     le cas, consultez # <http://www.gnu.org/licenses/gpl.txt



import glob 
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt

from scipy.special import gamma
from skimage import io
from skimage import color

#l= glob.glob('/c7/MEDUSE/IDF/Saclay/Sentinel/Pile_From_SNAP/SLC/PileRecalee/DESCENDING/IW1/Recalage-Elise/AmplitudesTif/*VH_amplitude.tif') 
l= glob.glob('/home/ekoenigu/Documents/AmplitudesTif/*VH_amplitude.tif') 
N=len(l)

im=io.imread(l[0])
nx=im.shape[0]
ny=im.shape[1]

threshold=5*im.mean()

m1=im;
m2=im**2;
Imax=im;

# Compute the date vector
datenum=[]
for i in range(0,N):
    (filepath, filename) = os.path.split(l[i])
    year=(int(filename[0:4]))
    month=(int(filename[4:6]))
    day=(int(filename[6:8]))
    d1 = datetime.date(year,month,day)
    datenum.append(d1.toordinal())
    
    
ds=float(max(datenum)-min(datenum))
    
# Compute Coefficient of Variation    
M1=np.zeros([nx,ny])
M2=np.zeros([nx,ny]) 
Kmax=np.zeros([nx,ny]) 
Imax=np.zeros([nx,ny])

k=0
for i in l:
    im=io.imread(i)
    M1=M1+im
    M2=M2+im**2
    Matrix_indicek=(datenum[k]-min(datenum))/ds*np.ones([nx,ny]);
    k=k+1
    Kmax=(im>Imax)*Matrix_indicek+(im<Imax)*Kmax 
    Imax=np.maximum(Imax,im)
    
M1=M1/N
M2=M2/N
R=np.sqrt(M2-M1**2)/M1

gam=R.mean()
L= ((0.991936+0.067646*gam-0.098888*gam**2 -0.048320*gam**3)/(0.001224-0.034323*gam+4.305577*gam**2-1.163498*gam**3));
    
CV=np.sqrt((L*gamma(L)**2/(gamma(L+0.5)**2))-1); # theretical mean value
num=(L*gamma(L)**4.*(4*(L**2)*gamma(L)**2-4*L*gamma(L+1/2)**2-gamma(L+1/2)**2));
den=(gamma(L+1/2)**4.*(L*gamma(L)**2-gamma(L+1/2)**2));
alpha=1/4*num/den; # theretical standard deviation value
    
R=(R-CV)/(alpha/np.sqrt(N))/10.0+0.25;
R=(R>1)*np.ones([nx,ny])+(R<1)*R;   # Cast Coefficient of Varation R max to 1.
R=(R>0)*R;
      
I=(Imax/threshold); # Cast Intensity to threshold. 
I=(I<1)*I+(I>1)*np.ones([nx,ny]);

hsv=np.zeros([nx,ny,3])
hsv[:, :, 0]=Kmax
hsv[:, :, 1]=R 
hsv[:, :, 2]=I 
    
C=color.hsv2rgb(hsv)


gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))

from datetime import date
d1 = date.fromordinal(min(datenum))
d2 = date.fromordinal(max(datenum))
str1=str(d1)
str2=str(d2)


from matplotlib import gridspec
gs = gridspec.GridSpec(2, 1,
                       width_ratios=[1],
                       height_ratios=[10, 1]
                       )

ax1 = plt.subplot(gs[0])
plt.imshow(C, aspect=3)

ax2 = plt.subplot(gs[1])
plt.imshow(gradient, aspect=10, cmap='hsv')
ax2.set_axis_off()
ax2.set_title(str1+'                        TimeLine                      '+str2)




    
    
    
    
    
    
