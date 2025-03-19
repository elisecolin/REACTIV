from scipy.special import gamma
import numpy as np
import matplotlib.pylab as plt
from tqdm import tqdm
from matplotlib import pyplot
from matplotlib.colors import hsv_to_rgb


def Stack2SingleCV(Stack, timeaxis=2):
    """
    Compute single coefficient of variation from a Stack of Images.

    Parameters
    ----------
    Stack : float ndarray
        Input stack of images.
    timeaxis : int, optional
        Designates the dimension for the temporal axis. Default is timeaxis=2.
    Returns
    -------
    R : float ndarray, between 0 and +infty
        Coefficient of variation for the temporal axis - Not scaled        
    """
    # Temporal dimension
    Nt = np.shape(Stack)[timeaxis]
    # Compute mean and mean of squares
    M1 = np.mean(Stack, axis=timeaxis)
    M2 = np.mean(Stack**2, axis=timeaxis)
    # Compute the Coefficient of Variation
    R = np.sqrt(M2 - M1**2) / M1
    R[M1 == 0] = 0  # Remove possible Nan output when signal is zero
    R[M1 < 0] = 0

    return R



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




def fusion2polar4reactiv(CV1, K1, A1, CV2, K2, A2):
    """
    Merge two sets of REACTIV parameters.

    Args:
        CV1 (numpy.ndarray): coefficient of variation matrix for the first polarimetric channel.
        K1 (numpy.ndarray): Matrix of temporal fraction for the first polarimetric channel .
        A1 (numpy.ndarray): Maximal Amplitude reached for the first polarimetric channel.
        CV2 (numpy.ndarray): coefficient of variation matrix for the second polarimetric channel.
        K2 (numpy.ndarray): Matrix of temporal fraction for the second polarimetric channel .
        A2 (numpy.ndarray): Maximal Amplitude reached for the second polarimetric channel.

    Returns:
        tuple: A tuple containing three matrices representing the merged reactiv variables.
            - CV (numpy.ndarray): Merged coefficient of variation values.
            - K (numpy.ndarray): Merged matrix of time fractions.
            - A (numpy.ndarray): Merged amplitude values.
    """
    # Compare matrices CV1 and CV2
    condition = CV1 > CV2
    # Select elements from K1 where CV1 > CV2, and from K2 otherwise
    K = np.where(condition, K1, K2)
    # Select elements with the maximum value between CV1 and CV2
    CV = np.maximum(CV1, CV2)
    # Select elements with the maximum value between A1 and A2
    A = A1+A2

    return CV, K, A


def CV_fromListofImages(S):
    """
    Calculate the mean vector and covariance matrix from a list of time series matrices.
    Inputs:
    - S: A list of matrices representing the time series matrices of different intensity components. 
         Each matrix has dimensions (nt, nx, ny), where nt is the number of time points,
         nx and ny are spatial dimensions. The list can contain p elements, where p is any number.
         Examples: p=2 for |Ex|,|Ey|. p=4 for Sony polarimetric Intensities. 
    Outputs:
    - CVRR: Coefficient of variation RR matrix, nd array of dimensions (nx,ny)
    - CVgv: Coefficient of variatioe gv matrix, nd array of dimensions (nx,ny)
    - CVVN: Coefficient of variatio VN matrix, nd array of dimensions (nx,ny)
    - CVAZ: Coefficient of variatio AZ matrix, nd array of dimensions (nx,ny)
    """
    p = len(S)  # Number of elements
    nx, ny = S[0].shape[1], S[0].shape[2]  # Spatial dimensions

    # Calculate mean for each matrix in S
    S_means = [np.mean(Si, axis=0) for Si in S]

    # Initialize covariance matrix C
    C = np.zeros((p, p, nx, ny))

    # Calculate covariance matrix
    for i in range(p):
        for j in range(i, p):  # Use symmetry, compute only for j >= i
            m = (S[i] - S_means[i]) * (S[j] - S_means[j])
            m_mean = np.mean(m, axis=0)
            C[i, j, :, :] = m_mean
            if i != j:
                C[j, i, :, :] = m_mean  # Use symmetry to fill in the lower half
                
    #for i in range(p):
    TraceC = np.zeros((nx, ny))
    for i in range(p):
            TraceC = TraceC + C[i, i, :, :]

    # Calculate the mean vector
    MU = np.array(S_means)
    NormU2 = np.sum(MU**2,axis=0)
    
    Cp = np.transpose(C, (2, 3, 0, 1)) # Cp has dimensions (nx,ny,p,p)
    MU_reshaped = MU.transpose(1, 2, 0) # MU_reshaped has dimensions (nx,ny,p)

    
    # Compute CVAZ - PRODUCT MU C MU 
    MU_transformed = np.repeat(MU_reshaped[:, :, :, np.newaxis], p, axis=-1) # MU_transformed has dimensions (nx,ny,p,p)
    Cmu = np.matmul(Cp,MU_transformed)  # matrix multiplication by list over the nx,ny elements
    Cmu = Cmu[:,:,:,0]  # Only vector result is needed
    muCmu=(np.sum(Cmu*MU_reshaped,axis=2))
    CVAZ=np.sqrt(muCmu)/NormU2  

    # Compute CVRR - DETERMINANT
    C_2d = Cp.reshape(-1, p, p)
    # Calculate the determinant for each (p, p) sub-matrix along the first axis
    determinants = np.linalg.det(C_2d)
    # Reshape the determinant matrix to (nx, ny)
    determinant_matrix = determinants.reshape(nx,ny)    
    
    CVRR = np.sqrt(determinant_matrix**(1/p) / NormU2)
    

    # Compute CVVN - Inverse C 
    Cinv=np.linalg.inv(C_2d)
    Cinv_p=Cinv.reshape(nx,ny,p,p)
    Cinvmu = np.matmul(Cinv_p,MU_transformed)  # matrix multiplication by list over the nx,ny elements
    Cinvmu = Cinvmu[:,:,:,0]  # Only vector result is needed
    muCinvmu=(np.sum(Cinvmu*MU_reshaped,axis=2))
    CVVN=1/np.sqrt(muCinvmu)
    
    # Compute CVgv - trace C 
    
    CVgv = np.sqrt(TraceC/NormU2)
    
    return CVRR, CVgv, CVVN, CVAZ  


def eigenimages(series):
    """
    Calculate pixel-wise covariance matrices and eigenvalues for a series of temporal image data to compute
    Inputs:
    - series: A list of p elements, where each element is a nd numpy arrays of shape (Nt, nx, ny),
              representing a time series of images.
    Outputs:
    - eigenimages: A list of p numpy arrays of shape (nx, ny), where each array represents an
                   image constructed from the eigenvalues corresponding to one of the p dimensions.
    """
    # Unpack the dimensions
    p, Nt, nx, ny = len(series), series[0].shape[0], series[0].shape[1], series[0].shape[2]
    # Initialize an array to hold the eigenvalues for each pixel
    eigenvalues_per_pixel = np.zeros((nx, ny, p))
    # Iterate over each pixel position
    for x in tqdm(range(nx)):
        for y in range(ny):
            # Extract the pixel series for this location across all p series and Nt times
            pixel_series = np.array([[series[i][t, x, y] for t in range(Nt)] for i in range(p)])
            # Calculate the covariance matrix for this pixel series
            cov_matrix = np.cov(pixel_series)*(Nt-1)/Nt
            mu=np.mean(pixel_series,axis=1)
            # Calculate the eigenvalues and store them
            eigenvalues, _ = np.linalg.eig(cov_matrix)
            eigenvalues_per_pixel[x, y,:] = eigenvalues
    # Construct the eigenimages from the eigenvalues
    eigenimages = [eigenvalues_per_pixel[:, :, i] for i in range(p)]
    return eigenimages



def CV_Generalized_Limits(series):
    """
    Calculate the limits Q-> infinity or - infinity, according to sqrt(lambda_min)/mu,  sqrt(lambda_max)/mu
    Inputs:
    - series: A list of p elements, where each element is a nd numpy arrays of shape (Nt, nx, ny),
              representing a time series of images.
    Outputs:
    - limitmin_per_pixel,limitmax_per_pixel: Two numpy arrays of shape (nx, ny), where each array corresponds to the limits 
        of Generalized CV.
    """
    # Unpack the dimensions
    p, Nt, nx, ny = len(series), series[0].shape[0], series[0].shape[1], series[0].shape[2]
    limitmin_per_pixel= np.zeros((nx, ny))
    limitmax_per_pixel= np.zeros((nx, ny))
    # Iterate over each pixel position
    for x in tqdm(range(nx)):
        for y in range(ny):
            # Extract the pixel series for this location across all p series and Nt times
            pixel_series = np.array([[series[i][t, x, y] for t in range(Nt)] for i in range(p)])
            # Calculate the covariance matrix for this pixel series
            cov_matrix = np.cov(pixel_series)*(Nt-1)/Nt
            mu=np.mean(pixel_series,axis=1)
            normMU=np.linalg.norm(mu)
            # Calculate the eigenvalues and store them
            eigenvalues, _ = np.linalg.eig(cov_matrix)
            limitmin_per_pixel[x, y] = np.sqrt(np.min(eigenvalues))/normMU
            limitmax_per_pixel[x, y] = np.sqrt(np.max(eigenvalues))/normMU
    # Construct the eigenimages from the eigenvalues
    # eigenimages = [eigenvalues_per_pixel[:, :, i] for i in range(p)]
    return limitmin_per_pixel,limitmax_per_pixel



def CV_Generalized_equally(S,Q):
    """
    Contrast definition based on the Generalized Equally Weighted Mean at order Q.
    """
    # Unpack the dimensions
    p, Nt, nx, ny = len(S), S[0].shape[0], S[0].shape[1], S[0].shape[2]
    # Initialize an array to hold the eigenvalues for each pixel
    CV = np.zeros((nx, ny))
    # Iterate over each pixel position
    for x in tqdm(range(nx)):
        for y in range(ny):
            # Extract the pixel series for this location across all p series and Nt times
            pixel_series = np.array([[S[i][t, x, y] for t in range(Nt)] for i in range(p)])
            # Calculate the covariance matrix for this pixel series
            C = np.cov(pixel_series)*(Nt-1)/Nt
            CQ=np.linalg.matrix_power(C, Q) # COvariance Matrix power Q
            mu=np.mean(pixel_series,axis=1)
            normMU=np.linalg.norm(mu)
            CV[x,y]=np.sqrt((1/p*np.trace(CQ))**(1/Q))/normMU
    return CV

def CV_Generalized_Non_equally(S,Q):
    """
    Contrast definition based on the Generalized Non Equally Weighted Mean at order Q, Q<>0.
    """
    # Unpack the dimensions
    p, Nt, nx, ny = len(S), S[0].shape[0], S[0].shape[1], S[0].shape[2]
    # Initialize an array to hold the eigenvalues for each pixel
    CV = np.zeros((nx, ny))
    # Iterate over each pixel position
    for x in tqdm(range(nx)):
        for y in range(ny):
            # Extract the pixel series for this location across all p series and Nt times
            pixel_series = np.array([[S[i][t, x, y] for t in range(Nt)] for i in range(p)])
            # Calculate the covariance matrix for this pixel series
            C = np.cov(pixel_series)*(Nt-1)/Nt
            mu=np.mean(pixel_series,axis=1)
            normMU=np.linalg.norm(mu)
            
            C_q = np.linalg.matrix_power(C, Q)
            intermediate_result = np.dot(mu.T, C_q)
            result = np.dot(intermediate_result, mu)
            
            CV[x,y]=np.sqrt((result)**(1/Q))/normMU**(1+1/Q)
    return CV

def CV_Generalized_Non_equally_Zero(S):
    """
    Contrast definition based on the Generalized Non Equally Weighted Mean at order 0.
    """
    # Unpack the dimensions
    p, Nt, nx, ny = len(S), S[0].shape[0], S[0].shape[1], S[0].shape[2]
    # Initialize an array to hold the eigenvalues for each pixel
    CV = np.zeros((nx, ny))
    # Iterate over each pixel position
    for x in tqdm(range(nx)):
        for y in range(ny):
            # Extract the pixel series for this location across all p series and Nt times
            pixel_series = np.array([[S[i][t, x, y] for t in range(Nt)] for i in range(p)])
            # Calculate the covariance matrix for this pixel series
            C = np.cov(pixel_series)*(Nt-1)/Nt
            mu=np.mean(pixel_series,axis=1)
            normMU=np.linalg.norm(mu)
            eigenvalues,U= np.linalg.eig(C)
            muprime=np.dot(U.T,mu)
            result=(np.prod(eigenvalues**(np.abs(muprime)**2)))**(1/normMU**2)
            CV[x,y]=np.sqrt(result)/normMU
    return CV



