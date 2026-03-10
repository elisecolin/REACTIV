import numpy as np
from scipy.special import gamma
from matplotlib.colors import hsv_to_rgb
from joblib import Parallel, delayed
import multiprocessing

# ============================================================
# GLOBAL UTILITIES
# ============================================================

EPS = 1e-12


def clean_array(x):
    """Replace NaN and Inf."""
    return np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)


def safe_divide(a, b):
    """Division avoiding zero."""
    return a / np.maximum(b, EPS)


def safe_sqrt(x):
    """Stable square root."""
    return np.sqrt(np.maximum(x, 0))


def covariance_matrix(X):
    """
    Stable covariance.
    X shape = (p, Nt)
    """
    Nt = X.shape[1]
    mu = np.mean(X, axis=1, keepdims=True)
    Xc = X - mu
    C = (Xc @ Xc.T) / Nt
    return C


# ============================================================
# BASIC CV
# ============================================================

def Stack2SingleCV(Stack, timeaxis=2):

    M1 = np.mean(Stack, axis=timeaxis)
    M2 = np.mean(Stack**2, axis=timeaxis)

    var = np.maximum(M2 - M1**2, 0)

    R = safe_divide(np.sqrt(var), M1)

    return clean_array(R)


# ============================================================
# REACTIV
# ============================================================

def Stack2reactiv(Stack, timeaxis=2, L=None):

    Nt = Stack.shape[timeaxis]

    M1 = np.mean(Stack, axis=timeaxis)
    M2 = np.mean(Stack**2, axis=timeaxis)

    var = np.maximum(M2 - M1**2, 0)

    R = safe_divide(np.sqrt(var), M1)

    Amax = np.max(Stack, axis=timeaxis)
    Kmax = np.argmax(Stack, axis=timeaxis)

    K = Kmax / Nt

    if L is None:

        gam = np.mean(R)

        L = ((0.991936 + 0.067646*gam -0.098888*gam**2 -0.048320*gam**3) /
             (0.001224 -0.034323*gam +4.305577*gam**2 -1.163498*gam**3))

    Rmean = np.sqrt((L * gamma(L)**2 / (gamma(L+0.5)**2)) - 1)

    num = (L*gamma(L)**4*(4*(L**2)*gamma(L)**2 - 4*L*gamma(L+0.5)**2 - gamma(L+0.5)**2))
    den = (gamma(L+0.5)**4*(L*gamma(L)**2 - gamma(L+0.5)**2))

    Rstd = 0.25 * num / den / np.sqrt(Nt)

    Rmax = Rmean + Rstd

    CV = (R - Rmean) / (Rmax - Rmean)

    CV = np.clip(CV, 0, 1)

    return clean_array(CV), K, Amax


# ============================================================
# RGB IMAGE
# ============================================================

def reactiv_image(CV, K, A, thresh=None):

    A = np.abs(A)

    if thresh is None:
        thresh = np.mean(A) + np.std(A)

    A = safe_divide(A, thresh)

    A = np.clip(A, 0, 1)

    hsv = np.dstack((K, CV, A))

    rgb = hsv_to_rgb(hsv)

    return rgb


# ============================================================
# POLARIMETRIC FUSION
# ============================================================

def fusion2polar4reactiv(CV1, K1, A1, CV2, K2, A2):

    condition = CV1 > CV2

    K = np.where(condition, K1, K2)

    CV = np.maximum(CV1, CV2)

    A = A1 + A2

    return clean_array(CV), clean_array(K), clean_array(A)


# ============================================================
# GENERALIZED CV CORE
# ============================================================

def _pixel_generalized_limits(series_pixel):

    X = np.array(series_pixel)

    if not np.isfinite(X).all():
        return 0, 0

    C = covariance_matrix(X)

    mu = np.mean(X, axis=1)

    normMU = np.linalg.norm(mu)

    if normMU < EPS:
        return 0, 0

    eigvals = np.linalg.eigvalsh(C)

    eigvals = np.maximum(eigvals, 0)

    lam_min = eigvals[0]
    lam_max = eigvals[-1]

    limitmin = np.sqrt(lam_min) / normMU
    limitmax = np.sqrt(lam_max) / normMU

    return limitmin, limitmax


# ============================================================
# PARALLEL GENERALIZED LIMITS
# ============================================================

def CV_Generalized_Limits(series, n_jobs=None):

    p = len(series)

    Nt, nx, ny = series[0].shape

    if n_jobs is None:
        n_jobs = multiprocessing.cpu_count()

    pixels = []

    for x in range(nx):
        for y in range(ny):

            pixel_series = [[series[i][t, x, y] for t in range(Nt)] for i in range(p)]

            pixels.append(pixel_series)

    results = Parallel(n_jobs=n_jobs)(
        delayed(_pixel_generalized_limits)(pix) for pix in pixels
    )

    results = np.array(results)

    limitmin = results[:, 0].reshape(nx, ny)
    limitmax = results[:, 1].reshape(nx, ny)

    return limitmin, limitmax


# ============================================================
# EIGEN IMAGES
# ============================================================

def eigenimages(series):

    p = len(series)

    Nt, nx, ny = series[0].shape

    eigenvalues = np.zeros((p, nx, ny))

    for x in range(nx):
        for y in range(ny):

            X = np.array([[series[i][t, x, y] for t in range(Nt)] for i in range(p)])

            C = covariance_matrix(X)

            eig = np.linalg.eigvalsh(C)

            eigenvalues[:, x, y] = eig

    return [eigenvalues[i] for i in range(p)]
