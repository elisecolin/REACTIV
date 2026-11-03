"""REACTIV: fast, NaN-aware NumPy implementation.

Scalar SAR inputs are *linear amplitudes* (not power or dB). Generalized
polarimetric CV inputs are real-valued temporal observables, e.g. SAR
intensities |S_HH|**2, |S_HV|**2 and |S_VV|**2.

Missing data policy:
    'omit' (default): omit non-finite acquisitions at each pixel. Polarimetric
        functions omit a date from ALL channels if any channel is missing.
    'raise': reject any non-finite input to expose a broken data pipeline.
    'propagate': an invalid date invalidates the whole corresponding pixel.

All estimators require at least two *valid* dates per pixel. Undefined results
remain NaN, never silently become genuine zero measurements. For scalar CV,
zero-amplitude pixels having >=2 valid dates retain the historical CV = 0.

Population variances are used throughout (ddof=0). K gives the index of the
maximal amplitude divided by the ORIGINAL number of acquisitions; omitted
acquisitions do not change the time coordinate. Ties select the earliest date.
"""

from __future__ import annotations

import numpy as np
from scipy.special import gammaln
from matplotlib.colors import hsv_to_rgb

EPS = 1e-12


def clean_array(x):
    """Legacy utility: replace invalid values with zero for display ONLY.

    Never use this before statistical estimations: doing so biases estimates.
    """
    return np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)


def safe_divide(a, b):
    """Legacy helper; prefer np.divide(..., where=...) for estimation."""
    return np.asarray(a) / np.maximum(b, EPS)


def safe_sqrt(x):
    return np.sqrt(np.maximum(x, 0))


def _check_nan_policy(nan_policy):
    if nan_policy not in ('omit', 'raise', 'propagate'):
        raise ValueError("nan_policy must be 'omit', 'raise', or 'propagate'.")


def _as_array(x):
    """Keep NumPy masked-array masks rather than accidentally losing them."""
    if np.ma.isMaskedArray(x):
        return np.ma.filled(x.astype(np.result_type(x.dtype, np.float64)), np.nan)
    return np.asarray(x)


def _check_scalar_stack(stack, timeaxis, nan_policy):
    _check_nan_policy(nan_policy)
    x = _as_array(stack)
    if x.ndim < 2:
        raise ValueError('Stack must have at least two dimensions.')
    if not isinstance(timeaxis, (int, np.integer)) or not -x.ndim <= timeaxis < x.ndim:
        raise ValueError('timeaxis is outside the array dimensions.')
    timeaxis = int(timeaxis) % x.ndim
    if x.shape[timeaxis] < 2:
        raise ValueError('At least two acquisitions are needed.')
    if not np.issubdtype(x.dtype, np.number) or np.iscomplexobj(x):
        raise ValueError('Expected real SAR amplitudes; use np.abs(SLC).')
    x = np.asarray(x, dtype=np.float64)
    if np.any(x[np.isfinite(x)] < 0):
        raise ValueError('Negative SAR amplitudes: expected linear amplitude, not dB.')
    if nan_policy == 'raise' and not np.isfinite(x).all():
        raise ValueError('Missing/non-finite amplitude with nan_policy="raise".')
    return x, timeaxis


def _scalar_moments(x, axis, nan_policy):
    """Numerically stable 2-pass variance with valid counts per pixel."""
    good = np.isfinite(x)
    count = np.count_nonzero(good, axis=axis)
    n_safe = np.maximum(count, 1)
    mu = np.sum(np.where(good, x, 0.0), axis=axis, dtype=np.float64) / n_safe
    centered = np.where(good, x - np.expand_dims(mu, axis=axis), 0.0)
    var = np.sum(centered * centered, axis=axis, dtype=np.float64) / n_safe
    valid = count >= 2
    if nan_policy == 'propagate':
        valid &= count == x.shape[axis]
    return mu, np.sqrt(np.maximum(var, 0.0)), count, good, valid


def _check_positive_integer(n, name):
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 2:
        raise ValueError(f'{name} must be an integer >= 2.')
    return int(n)


def reactiv_normalization(Nt, L):
    """Return theoretical mean CV and asymptotic std for Nt Rayleigh/Nakagami samples.

    L is the Equivalent Number of Looks, assumed identical for all pixels.
    The standard deviation is sqrt(variance_prefactor / Nt), *not*
    variance_prefactor / sqrt(Nt). Uses log gamma for numerical stability.
    """
    Nt = _check_positive_integer(Nt, 'Nt')
    if np.ndim(L) != 0 or not np.isfinite(L) or L <= 0:
        raise ValueError('L must be finite, scalar and strictly positive.')
    L = float(L)
    if L > 500:
        # Asymptotic expressions: avoid severe cancellation for very large L.
        return 0.5 / np.sqrt(L), 1.0 / np.sqrt(8 * L * Nt)
    log_z = 2 * (gammaln(L + 0.5) - gammaln(L)) - np.log(L)
    z = np.exp(log_z)
    rmean = np.sqrt(np.expm1(-log_z))
    variance_prefactor = (4 * L - (4 * L + 1) * z) / (4 * L * z * z * (1 - z))
    if not np.isfinite(variance_prefactor) or variance_prefactor <= 0:
        raise ValueError('Unable to calculate theoretical CV variance for L.')
    return float(rmean), float(np.sqrt(variance_prefactor / Nt))


def _estimate_L(R, valid):
    """Historical empirical inversion; short or masked stacks may bias it."""
    use = valid & np.isfinite(R) & (R > 0)
    if not np.any(use):
        return None
    gam = float(np.mean(R[use]))
    den = 0.001224 - 0.034323 * gam + 4.305577 * gam**2 - 1.163498 * gam**3
    L = (0.991936 + 0.067646 * gam - 0.098888 * gam**2
         - 0.048320 * gam**3) / den
    if not np.isfinite(L) or L <= 0:
        raise ValueError('Empirical L estimation failed; specify L explicitly.')
    return float(L)


def _normalize_cv(R, count, L, sigma_scale, valid):
    if np.ndim(sigma_scale) != 0 or not np.isfinite(sigma_scale) or sigma_scale <= 0:
        raise ValueError('sigma_scale must be finite and > 0.')
    cv = np.full(R.shape, np.nan, dtype=np.float64)
    if not np.any(valid):
        return cv
    if L is None:
        L = _estimate_L(R, valid)
        if L is None:
            cv[valid] = 0.0
            return cv
    # Nt=2 determines the prefactor; per-pixel effective count differs with NaNs.
    rmean, rstd2 = reactiv_normalization(2, L)
    rstd = rstd2 * np.sqrt(2.0 / np.maximum(count, 1))
    scaled = np.clip((R - rmean) / (sigma_scale * rstd), 0.0, 1.0)
    cv[valid] = scaled[valid]
    return cv


def Stack2SingleCV(Stack, timeaxis=2, nan_policy='omit', return_counts=False):
    """Raw temporal CV (std/mean) of a real-valued amplitude stack.

    Invalid pixels (<2 valid dates) return NaN. Valid identically-zero pixels
    return zero. When return_counts=True, also return per-pixel valid date count.
    """
    x, axis = _check_scalar_stack(Stack, timeaxis, nan_policy)
    mu, sd, counts, _, valid = _scalar_moments(x, axis, nan_policy)
    r = np.full(mu.shape, np.nan, dtype=np.float64)
    r[valid] = 0.0
    np.divide(sd, mu, out=r, where=valid & (mu > 0))
    return (r, counts) if return_counts else r


def Stack2reactiv(Stack, timeaxis=2, L=None, sigma_scale=1.0,
                  nan_policy='omit', return_counts=False):
    """REACTIV (CV-normalized, index K, maximum amplitude) of SAR amplitudes.

    Invalid pixels return NaN in all three maps; RGB rendering maps them to
    black. NaNs are ignored ONLY in nan_policy='omit' mode. The effective
    sample count determines each pixel's theoretical CV standard deviation.

    If L=None, the historical empirical global L estimator is used; with
    short stacks, explicit physical L is strongly preferred.
    """
    x, axis = _check_scalar_stack(Stack, timeaxis, nan_policy)
    Nt = x.shape[axis]
    mu, sd, counts, good, valid = _scalar_moments(x, axis, nan_policy)
    R = np.zeros_like(mu)
    np.divide(sd, mu, out=R, where=mu > 0)
    CV = _normalize_cv(R, counts, L, sigma_scale, valid)
    masked = np.where(good, x, -np.inf)
    Amax = np.max(masked, axis=axis).astype(np.float64)
    K = np.argmax(masked, axis=axis).astype(np.float64) / Nt
    Amax[~valid] = np.nan
    K[~valid] = np.nan
    return (CV, K, Amax, counts) if return_counts else (CV, K, Amax)


def List2reactiv(List, L=None, sigma_scale=1.0, nan_policy='omit',
                 return_counts=False):
    """Read temporal amplitude TIFFs in file order, one image at a time.

    TIFF NoData/masked pixels are treated as missing, not as zero intensity.
    Requires rasterio only when this function is called. The results should
    agree with Stack2reactiv(..., timeaxis=0) for the same TIFF samples.
    """
    import rasterio
    _check_nan_policy(nan_policy)
    paths = list(List)
    Nt = _check_positive_integer(len(paths), 'Number of TIFF images')
    mean = M2 = count = amax = kmax = None
    for t, path in enumerate(paths):
        with rasterio.open(path) as ds:
            if ds.count != 1:
                raise ValueError(f'Expected one image band at {path}.')
            im = ds.read(1, masked=True)
        im = _as_array(im).astype(np.float64, copy=False)
        if mean is not None and im.shape != mean.shape:
            raise ValueError(f'Inconsistent TIFF shape: {path}.')
        if np.any(im[np.isfinite(im)] < 0):
            raise ValueError(f'Negative SAR amplitude in {path}.')
        good = np.isfinite(im)
        if nan_policy == 'raise' and not np.all(good):
            raise ValueError(f'Missing/non-finite values at {path}.')
        if t == 0:
            count = np.zeros(im.shape, dtype=np.int32)
            mean = np.zeros(im.shape, dtype=np.float64)
            M2 = np.zeros(im.shape, dtype=np.float64)
            amax = np.full(im.shape, -np.inf, dtype=np.float64)
            kmax = np.zeros(im.shape, dtype=np.float64)
            seen_bad = np.zeros(im.shape, dtype=bool)
        seen_bad |= ~good
        new_count = count + good.astype(np.int32)
        delta = np.where(good, im - mean, 0.0)
        new_mean = mean + delta / np.maximum(new_count, 1)
        M2 += delta * np.where(good, im - new_mean, 0.0)
        mean = new_mean
        count = new_count
        greater = good & (im > amax)
        amax[greater] = im[greater]
        kmax[greater] = t
    valid = count >= 2
    if nan_policy == 'propagate':
        valid &= ~seen_bad
    sd = np.sqrt(np.maximum(M2 / np.maximum(count, 1), 0.0))
    R = np.zeros_like(mean)
    np.divide(sd, mean, out=R, where=mean > 0)
    CV = _normalize_cv(R, count, L, sigma_scale, valid)
    kmax[~valid] = np.nan
    amax[~valid] = np.nan
    result = (CV, kmax / Nt, amax)
    return (*result, count) if return_counts else result


def reactiv_image(CV, K, A, thresh=None):
    """HSV to RGB; undefined (NaN/Inf) pixels appear black, not colored.

    Threshold statistics are computed only on pixels with all three finite
    fields. Returns floating RGB in [0,1], shape (ny,nx,3).
    """
    cv, k, a = [np.asarray(v, dtype=np.float64) for v in (CV, K, A)]
    if cv.ndim != 2 or k.shape != cv.shape or a.shape != cv.shape:
        raise ValueError('CV, K, A must have identical 2-D shapes.')
    good = np.isfinite(cv) & np.isfinite(k) & np.isfinite(a)
    result = np.zeros(cv.shape + (3,), dtype=np.float64)
    if not np.any(good):
        return result
    amp = np.abs(a)
    if thresh is None:
        values = amp[good]
        thresh = float(np.mean(values) + np.std(values))
    if np.ndim(thresh) != 0 or not np.isfinite(thresh) or thresh < 0:
        raise ValueError('thresh must be finite and nonnegative.')
    value = (np.clip(amp[good] / thresh, 0, 1)
             if thresh > 0 else np.zeros(np.count_nonzero(good)))
    hsv = np.column_stack((np.mod(k[good], 1), np.clip(cv[good], 0, 1), value))
    result[good] = hsv_to_rgb(hsv)
    return result


def fusion2polar4reactiv(CV1, K1, A1, CV2, K2, A2):
    """Fuse two REACTIV channels, respecting missing pixels.

    When just one channel is finite use that channel; when both are valid,
    select the greater CV and sum amplitude (historical behavior).
    """
    cv1, k1, a1, cv2, k2, a2 = np.broadcast_arrays(
        *[np.asarray(x, dtype=np.float64) for x in (CV1, K1, A1, CV2, K2, A2)])
    v1 = np.isfinite(cv1) & np.isfinite(k1) & np.isfinite(a1)
    v2 = np.isfinite(cv2) & np.isfinite(k2) & np.isfinite(a2)
    choice1 = v1 & (~v2 | (cv1 > cv2))
    cv = np.where(choice1, cv1, cv2)
    k = np.where(choice1, k1, k2)
    a = np.where(v1, a1, 0.0) + np.where(v2, a2, 0.0)
    invalid = ~v1 & ~v2
    cv[invalid] = np.nan
    k[invalid] = np.nan
    a[invalid] = np.nan
    return cv, k, a


def covariance_matrix(X, nan_policy='omit'):
    """Population covariance matrix of (p,Nt) real series, with joint omission."""
    _check_nan_policy(nan_policy)
    X = _as_array(X)
    if X.ndim != 2 or X.shape[1] < 2 or not np.issubdtype(X.dtype, np.number) or np.iscomplexobj(X):
        raise ValueError('X must contain p real series of Nt>=2 samples.')
    X = np.asarray(X, dtype=np.float64)
    if nan_policy == 'raise' and not np.isfinite(X).all():
        raise ValueError('Non-finite polarimetric values.')
    good = np.isfinite(X).all(axis=0)
    if (nan_policy == 'propagate' and not np.all(good)) or np.count_nonzero(good) < 2:
        return np.full((X.shape[0], X.shape[0]), np.nan)
    Y = X[:, good]
    Yc = Y - np.mean(Y, axis=1, keepdims=True)
    return Yc @ Yc.T / Y.shape[1]


def polar_valid_counts(series, chunk_pixels=32768):
    """Map of jointly valid acquisition counts for multichannel time series.

    A date is counted at a pixel only when *all* channels are finite.
    Complex-valued input is allowed here (only validity is inspected).
    """
    if not isinstance(series, (list, tuple)) or not series:
        raise ValueError('series must be a non-empty list/tuple of arrays.')
    arrays = [_as_array(a) for a in series]
    shape = arrays[0].shape
    if len(shape) != 3 or any(a.shape != shape for a in arrays):
        raise ValueError('Every channel must have shape (Nt,nx,ny).')
    if not isinstance(chunk_pixels, (int, np.integer)) or isinstance(chunk_pixels, bool) or chunk_pixels < 1:
        raise ValueError('chunk_pixels must be a positive integer.')
    Nt, nx, ny = shape
    flats = [a.reshape(Nt, -1) for a in arrays]
    counts = np.zeros(nx*ny, dtype=np.int32)
    for i in range(0, nx*ny, chunk_pixels):
        j = min(i + chunk_pixels, nx*ny)
        common = np.ones((Nt,j-i), dtype=bool)
        for a in flats:
            common &= np.isfinite(a[:,i:j])
        counts[i:j] = np.count_nonzero(common, axis=0)
    return counts.reshape(nx,ny)


def _polar_blocks(series, chunk_pixels, nan_policy):
    """Yield (slice, eigenvalues, eigenvectors, mean, norm, valid, count).

    Covariance is the ddof=0 covariance of the SAME valid dates across all
    channels. Memory scales as p * Nt * chunk_pixels, not full scene size.
    """
    _check_nan_policy(nan_policy)
    if not isinstance(series, (list, tuple)) or not series:
        raise ValueError('series must be a list/tuple of (Nt,nx,ny) arrays.')
    arrays = [_as_array(x) for x in series]
    shape = arrays[0].shape
    if len(shape) != 3 or shape[0] < 2 or any(x.shape != shape for x in arrays):
        raise ValueError('All polarimetric series must share (Nt>=2,nx,ny).')
    if any(not np.issubdtype(x.dtype, np.number) or np.iscomplexobj(x) for x in arrays):
        raise ValueError('Polarimetric generalized CV requires real observables, e.g. |S|^2.')
    if isinstance(chunk_pixels, bool) or not isinstance(chunk_pixels, (int, np.integer)) or chunk_pixels < 1:
        raise ValueError('chunk_pixels must be a positive integer.')
    Nt, nx, ny = shape
    flats = [x.reshape(Nt, -1) for x in arrays]
    p = len(flats)
    for i in range(0, nx * ny, chunk_pixels):
        j = min(i + chunk_pixels, nx * ny)
        X = np.stack([a[:, i:j] for a in flats], axis=0).astype(np.float64)
        good = np.isfinite(X).all(axis=0)  # (Nt, pixels), shared across channels
        if nan_policy == 'raise' and not np.all(good):
            raise ValueError(f'Non-finite sample in pixel block {i}:{j}.')
        count = np.count_nonzero(good, axis=0)
        valid = count >= 2
        if nan_policy == 'propagate':
            valid &= count == Nt
        n = np.maximum(count, 1)
        mu = np.sum(np.where(good[None, :, :], X, 0.0), axis=1) / n[None, :]
        Xc = np.where(good[None, :, :], X - mu[:, None, :], 0.0)
        cov = np.einsum('ptk,qtk->kpq', Xc, Xc, optimize=True) / n[:, None, None]
        eig, U = np.linalg.eigh(cov)
        eig = np.maximum(eig, 0.0)   # clip negative floating-point roundoff
        # Detect numerically zero eigenvalues at matrix precision (relative
        # tolerance). Critical for rank-deficient covariance with Q <= 0.
        spectral_tol = (32 * p * np.finfo(np.float64).eps
                        * np.max(eig, axis=1, keepdims=True))
        eig[eig <= spectral_tol] = 0.0
        norm = np.linalg.norm(mu, axis=0)
        yield i, j, eig, U, mu, norm, valid, count


def eigenimages(series, chunk_pixels=32768, nan_policy='omit'):
    """Ascending covariance eigenvalue images; undefined pixels are NaN."""
    _, nx, ny = _polar_shape(series)
    p = len(series)
    out = np.full((p, nx * ny), np.nan)
    for i, j, eig, _, _, _, valid, _ in _polar_blocks(series, chunk_pixels, nan_policy):
        out[:, i:j] = np.where(valid[None, :], eig.T, np.nan)
    return [out[k].reshape(nx, ny) for k in range(p)]


def _polar_shape(series):
    if not isinstance(series, (tuple, list)) or not series:
        raise ValueError('series must contain at least one array.')
    shp = np.shape(series[0])
    if len(shp) != 3:
        raise ValueError('series must have shape (Nt,nx,ny).')
    return shp


def CV_Generalized_Limits(series, chunk_pixels=32768, nan_policy='omit'):
    """The Q->-inf and Q->+inf limits: sqrt(lambda_min/max)/||mu||."""
    _, nx, ny = _polar_shape(series)
    lo = np.full(nx * ny, np.nan)
    hi = np.full(nx * ny, np.nan)
    for i, j, eig, _, _, norm, valid, _ in _polar_blocks(series, chunk_pixels, nan_policy):
        ok = valid & (norm > 0)
        l = np.full(j-i, np.nan)
        h = np.full(j-i, np.nan)
        np.divide(np.sqrt(eig[:, 0]), norm, out=l, where=ok)
        np.divide(np.sqrt(eig[:, -1]), norm, out=h, where=ok)
        lo[i:j], hi[i:j] = l, h
    return lo.reshape(nx, ny), hi.reshape(nx, ny)


def _check_Q(Q):
    if np.ndim(Q) != 0 or not np.isrealobj(Q) or not np.isfinite(Q):
        raise ValueError('Q must be a finite real scalar.')
    return float(Q)


def _power_mean(lam, Q, weights):
    """Generalized spectral power mean of nonnegative eigenvalues.

    lam, weights both (pixels, p); weights sum to 1 per pixel. Exact zero
    eigenvalues contribute zero for Q>=0; for Q<0, any zero eigenvalue with
    positive spectral weight makes the inverse-power mean zero.
    """
    # Ignore zero-weight components even for negative Q; avoid 0 * inf.
    active = weights > 0
    zeros_active = np.any((lam == 0) & active, axis=1)
    if Q == 0:
        log_lam = np.log(np.where(lam > 0, lam, 1.0))
        result = np.exp(np.sum(np.where(active, weights * log_lam, 0.0), axis=1))
        result[zeros_active] = 0.0
    else:
        powered = np.power(np.where((lam > 0) & active, lam, 1.0), Q)
        summed = np.sum(np.where(active, weights * powered, 0.0), axis=1)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            result = np.power(summed, 1.0/Q)
        if Q < 0:
            result[zeros_active] = 0.0
    return result


def CV_Generalized_equally(S, Q, chunk_pixels=32768, nan_policy='omit'):
    """Equal eigenvalue weights: sqrt((mean lambda**Q)**(1/Q))/||mu||.

    Q=0 uses the continuous geometric mean, avoiding the original 1/Q
    singularity. Q may be positive, negative or zero.
    """
    Q = _check_Q(Q)
    _, nx, ny = _polar_shape(S)
    out = np.full(nx * ny, np.nan)
    for i, j, eig, _, _, norm, valid, _ in _polar_blocks(S, chunk_pixels, nan_policy):
        p = eig.shape[1]
        weights = np.full_like(eig, 1/p)
        spectral = _power_mean(eig, Q, weights)
        ok = valid & (norm > 0)
        np.divide(np.sqrt(spectral), norm, out=out[i:j], where=ok)
    return out.reshape(nx, ny)


def CV_Generalized_Non_equally(S, Q, chunk_pixels=32768, nan_policy='omit'):
    """Mean-vector-weighted spectral CV of order Q.

    This is sqrt(power_mean(lambda, w, Q))/||mu||, where weights w are the
    squared projections of mean vector on covariance eigenvectors.
    Q=0 is the continuous weighted geometric mean. For rank-deficient
    covariance, negative Q uses its continuous generalized-mean limit.
    """
    Q = _check_Q(Q)
    _, nx, ny = _polar_shape(S)
    out = np.full(nx * ny, np.nan)
    for i, j, eig, U, mu, norm, valid, _ in _polar_blocks(S, chunk_pixels, nan_policy):
        projected = np.einsum('kpa,pk->ka', U, mu, optimize=True)
        weights = np.zeros_like(eig)
        np.divide(projected*projected, norm[:, None]**2, out=weights,
                  where=norm[:, None] > 0)
        # Ignore roundoff-only projections onto orthogonal eigenspaces.
        weights[weights <= 64 * np.finfo(np.float64).eps] = 0.0
        # Renormalize to remove floating-point norm errors.
        wsum = np.sum(weights, axis=1, keepdims=True)
        np.divide(weights, wsum, out=weights, where=wsum > 0)
        spectral = _power_mean(eig, Q, weights)
        ok = valid & (norm > 0)
        np.divide(np.sqrt(spectral), norm, out=out[i:j], where=ok)
    return out.reshape(nx, ny)


def CV_Generalized_Non_equally_Zero(S, chunk_pixels=32768, nan_policy='omit'):
    """Weighted generalized CV at Q=0, computed via spectral logarithms."""
    return CV_Generalized_Non_equally(S, 0, chunk_pixels, nan_policy)


def CV_fromListofImages(S, chunk_pixels=32768, nan_policy='omit'):
    """Historical (CVRR, CVgv, CVVN, CVAZ), NaN-aware and batched.

    Preserves *historical* CVgv = sqrt(trace(C))/||mu||. Note this is sqrt(p)
    TIMES CV_Generalized_equally(S, 1), whose average is trace(C)/p.
    CVRR equals equal-weight Q=0; CVVN equals weighted Q=-1 if C is
    invertible; CVAZ equals weighted Q=1. Singular cases use spectral limits
    instead of failing during inversion.
    """
    _, nx, ny = _polar_shape(S)
    output = [np.full(nx*ny, np.nan) for _ in range(4)]
    for i, j, eig, U, mu, norm, valid, _ in _polar_blocks(S, chunk_pixels, nan_policy):
        p = eig.shape[1]
        equal = np.full_like(eig, 1/p)
        proj = np.einsum('kpa,pk->ka', U, mu, optimize=True)
        weights = np.zeros_like(eig)
        np.divide(proj*proj, norm[:, None]**2, out=weights,
                  where=norm[:, None] > 0)
        weights[weights <= 64 * np.finfo(np.float64).eps] = 0.0
        s = weights.sum(axis=1, keepdims=True)
        np.divide(weights, s, out=weights, where=s > 0)
        spectral = (_power_mean(eig, 0, equal), eig.sum(axis=1),
                    _power_mean(eig, -1, weights), _power_mean(eig, 1, weights))
        ok = valid & (norm > 0)
        for dst, value in zip(output, spectral):
            np.divide(np.sqrt(value), norm, out=dst[i:j], where=ok)
    return tuple(x.reshape(nx, ny) for x in output)

