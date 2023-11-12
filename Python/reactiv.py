#     -*- coding: utf-8 -*-
#     « Copyright (c) 2018, Elise Colin (Onera) » 
#     Updated November 12, 2023
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
#     même temps que REACTIV dans le fichier copying.txt ; si ce n'est pas
#     le cas, consultez # <http://www.gnu.org/licenses/gpl.txt
#
#     For more explanations:
#     https://medium.com/@elisecolin/how-to-visualize-changes-in-a-radar-timeline-fb79ef526c1e

from scipy.special import gamma

def Stack2reactiv(Stack, timeaxis=2, L=None):
    """
    Compute REACTIV representation from a Stack of Images.

    Parameters
    ----------
    Stack : float ndarray
        Input stack of images.
    timeaxis : int, optional
        Designates the dimension for the temporal axis. Default is timeaxis=2.
    L : float, optional
        Equivalent Number of Looks. If data are SLC, then L=1. If not specified,
        it will be automatically estimated from the data.

    Returns
    -------
    CV : float ndarray, between 0 and 1
        Coefficient of variation for the temporal axis.
        Scaled between 0 (theoretical mean value) and 1 (theoretical mean value + theoretical standard deviation).
    K : float ndarray, between 0 and 1
        Corresponds to the temporal fraction where maximal amplitude is found during the observation period.
    Amax : float ndarray - not scaled
        Corresponds to the image of maximal amplitude values reached during the observation period.
    """

    # Temporal dimension
    Nt = np.shape(Stack)[timeaxis]

    # Compute mean and mean of squares
    M1 = np.mean(Stack, axis=timeaxis)
    M2 = np.mean(Stack**2, axis=timeaxis)

    # Find maximum amplitude and corresponding fraction of time
    Amax = np.amax(Stack, axis=timeaxis)
    Kmax = np.argmax(Stack, axis=timeaxis)
    K = Kmax / Nt

    # Compute the Coefficient of Variation
    R = np.sqrt(M2 - M1**2) / M1
    R[M1 == 0] = 0  # Remove possible Nan output when signal is zero
    R[M1 < 0] = 0

    # Theoretical estimation of mean and std of the Coefficient of Variation for the given Equivalent Number of Looks
    if L is None:
        gam = R.mean()
        L = ((0.991936 + 0.067646 * gam - 0.098888 * gam**2 - 0.048320 * gam**3) /
             (0.001224 - 0.034323 * gam + 4.305577 * gam**2 - 1.163498 * gam**3))

    Rmean = np.sqrt((L * gamma(L)**2 / (gamma(L + 0.5)**2)) - 1)  # theoretical mean value
    num = (L * gamma(L)**4 * (4 * (L**2) * gamma(L)**2 - 4 * L * gamma(L + 1/2)**2 - gamma(L + 1/2)**2))
    den = (gamma(L + 1/2)**4 * (L * gamma(L)**2 - gamma(L + 1/2)**2))
    Rstd = 1 / 4 * num / den / np.sqrt(Nt)  # theoretical standard deviation value
    Rmax = Rmean + Rstd

    # Recast Coefficient of Variation between 0 (mean value) and 1 (max value)
    CV = (R - Rmean) / (Rmax - Rmean)
    CV = (CV >= 0) * CV
    CV = CV * (CV < 1) + (CV >= 1)

    return CV, K, Amax


def reactiv_image(CV, K, A, thresh=None):
    """
    Create a REACTIV representation as an RGB image.
    Parameters
    ----------
    CV : float ndarray, between 0 and 1
        Coefficient of variation for the temporal axis. Scaled between 0 (theoretical mean value) and 1 (theoretical mean value + theoretical standard deviation).
    K : float ndarray, between 0 and 1
        Corresponds to the temporal fraction where maximal amplitude is found during the observation period.
    A : float ndarray - not scaled
        Corresponds to the image of maximal amplitude values reached during the observation period.
    thresh : float, optional
        Threshold for amplitude values. If not specified, it is set to the mean plus the standard deviation of the amplitude.
    Returns
    -------
    rgb : ndarray
        RGB representation of the REACTIV image.
    """
    # Ensure amplitude values are positive
    A = np.abs(A)
    # Set threshold to mean plus standard deviation if not specified
    if thresh is None:
        thresh = np.mean(A) + np.std(A)
    # Normalize amplitude values between 0 and 1
    A = A / thresh
    A = A * (A < 1) + (A >= 1)
    # Stack CV, K, and normalized A to create an HSV image
    result = np.dstack((K, CV, A))
    # Convert HSV to RGB
    rgb = hsv_to_rgb(result)

    return rgb

def List2reactiv(List,L=None):
    """ Product output for REACTIV representation from a List of Images
    Parameters
    ----------
    List : List of Image files
        Example. If all your images are stored in .tif files in Directory, then call
        List= glob.glob('/home/Directory/*.tif') 
        It is preferable that .tif files containes Amplitudes (and not Intensities)
    L : float
        L is Equivalent Number of Looks. If data are SLC, then L=1. 
        if you do not specify L, it will be automatically estimated from the data
    Returns
    -------
    CV : float ndarray, between 0 and 1
        Coefficient of variation for temporal axis. 
        Casted between 0, theoretical mean value, and 1, theoretical mean value+theoretical standard deviation 
    K : float ndarray, between 0 and 1
        Corresponds to the temporal fraction where maximal Amplitude is found during the observation period
    Amax : float ndarray - not casted
        Corresponds to the image of maximal amplitude values reached during the observation period
    """
    Nt=len(List)# temporal dimension
    
    im=io.imread(List[0]) # read the first image
    nx,ny=np.shape(im) 
    Imax=im;

    # Compute Coefficient of Variation    
    M1=np.zeros([nx,ny])
    M2=np.zeros([nx,ny]) 
    Kmax=np.zeros([nx,ny]) 
    Amax=np.zeros([nx,ny])

    k=0
    for count,file in enumerate(List):
        im=io.imread(file)
        M1=M1+im
        M2=M2+im**2
        k=k+1
        Kmax=(im>Imax)*count+(im<Imax)*Kmax 
        Amax=np.maximum(Amax,im)
    M1=M1/Nt
    M2=M2/Nt
    K=Kmax/Nt
            
# Compute the Coefficient of variation       
    R=np.sqrt(M2-M1**2)/M1; 
    R[M1==0]=0 # Remove possible Nan Oupput when Signal is zero
    R[M1<0]=0
    
# Theoretical estimation of mean and std of the Coefficient of variation for the given Equivalent Number of Looks     
    if L is None:
        gam=R.mean()
        L= ((0.991936+0.067646*gam-0.098888*gam**2 -0.048320*gam**3)/(0.001224-0.034323*gam+4.305577*gam**2-1.163498*gam**3));
    
    Rmean=np.sqrt((L*gamma(L)**2/(gamma(L+0.5)**2))-1); # theretical mean value
    num=(L*gamma(L)**4.*(4*(L**2)*gamma(L)**2-4*L*gamma(L+1/2)**2-gamma(L+1/2)**2));
    den=(gamma(L+1/2)**4.*(L*gamma(L)**2-gamma(L+1/2)**2));
    Rstd=1/4*num/den/np.sqrt(Nt); # theoretical standard deviation value
    Rmax=Rmean+Rstd

# recast Coefficient of Variation between 0 (mean value) and 1 (max value)    
    CV=(R-Rmean)/(Rmax-Rmean)
    CV=(CV>=0)*CV
    CV=CV*(CV<1)+(CV>=1);
    return CV,K,Amax
