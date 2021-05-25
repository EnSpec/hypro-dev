import os
import logging

import itertools
from itertools import chain

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from ENVI import read_envi_header
from Radiometry import bad_detectors_from_hyspex_set

logger = logging.getLogger(__name__)


def neighborhood(X, idcs, steps, axes, write=False):
    ''' Pure `numpy` approach using `lib.stride_tricks.sliding_window_view`;
         Create a sliding window view with the desired shape, returning
         the subset of views that are aligned with the neighborhood.
    '''
    # Get sliding window view
    window_shape = tuple(2*s+1 for s in steps)
    S = sliding_window_view(X, window_shape, axis=axes, writeable=write)
    
    # Initialize all indices with empty slice
    indices = [slice(None)]*X.ndim
    
    # Update specified axes with adjusted index
    for axis, idx, step in zip(axes, idcs, steps):
        indices[axis] = idx-step
    
    # Get views over the specified neighborhood
    views = S[tuple(indices)]
    
    return views

def get_neighborhood_mask(steps, axes, ndim=None):
    
    # Form of mask depends only on `step`
    shape = tuple(1+2*S for S in steps)
    mask = np.zeros(shape, dtype=np.int)
    mask[steps] = 1
    
    # Sort the axes
    mask = np.moveaxis(mask, np.argsort(axes), np.arange(len(axes)))
    
    if ndim is None: ndim = max(axes)+1
    
    # Any axes which are unspecified must be added to the mask
    missing_axes = [i for i in range(ndim) if not i in axes]
    mask = np.expand_dims(mask, missing_axes)
    
    return mask


### ----------------------------------------

def matsolve__lstsq_(P, Q):
    return np.linalg.lstsq(P, Q, rcond=None)[0]


def matsolve__pinv_(P, Q):
    return np.linalg.pinv(P).dot(Q)


def polynomial_combinations(dof, order=3, inter_only=False):
    ''' Get all combinations of polynomial features of degree less than
         or equal to `order` among input variables representing one or more
         degrees of freedom, as given by `dof`.
        
        This method separates the combinatorics of polynomial fitting
         from the matrix math.
    '''
    # Get combinations algorithm & number of resulting features
    if inter_only:
        # Interactions only, i.e. sampling without replacement
        # No term may be more than first order in any variable 
        combinations = itertools.combinations
        # # Number of features is nCk
        # n_features = sum(comb(dof, i) for i in range(order+1))
    else:
        # Sampling with replacement
        # All combinations are considered
        combinations = itertools.combinations_with_replacement
        # # Number of features is (n+k-1)Ck
        # n_features = sum(comb(dof+i-1, i) for i in range(order+1))
    
    # Get combinations of input variables with degree <= `order`
    combos = chain.from_iterable(
                    np.array(list(combinations(range(dof), i))).astype(np.int64)
                    for i in range(order+1)
                   )
    
    return list(combos)


def polynomial_matrix(X, combos):
    ''' Generate a matrix of polynomial terms given a set of observations
         and a list of feature combinations to be included.
    '''
    n_samples, _ = np.shape(X)
    
    # Initialize output array
    XP = np.empty((n_samples, len(combos)), dtype=X.dtype)
    
    # Populate with polynomial terms
    for i, combo in enumerate(combos):
        XP[:, i] = X[:, combo].prod(axis=1)
    # build_poly_matrix__numba_(X, combos, XP)
    
    return XP


### ----------------------------------------

def interpolate_detectors(sensor_dict, mode='bicurvilinear', *args, **kwargs):
    ''' Interpolate radiance over bad detectors. '''
    
    if mode == 'bicurvilinear':
        interpolate = interpolate_detectors__bicurvilinear_
    
    # Initialize radiance cube from disk
    hdr = read_envi_header(os.path.splitext(sensor_dict['raw_rdn_image_file'])[0]+'.hdr')
    rdn_image = np.memmap(sensor_dict['raw_rdn_image_file'],
                          shape=(hdr['lines'], hdr['bands'], hdr['samples']),
                          offset=hdr['header offset'],
                          mode='r+', dtype='float32')#, dtype='int16')
    
    ### TODO: Read dtype from header!
    
    # Load bad detectors registry from settings file
    bad_detectors = bad_detectors_from_hyspex_set(sensor_dict['setting_file'])
    
    # Interpolate over bad detectors
    logger.info('Do interpolation.')
    wavelengths = np.array(hdr['wavelength'], dtype=np.float64)
    
    interpolate(rdn_image, bad_detectors, wavelengths, *args, **kwargs)
    
    # Commit changes to disk
    del rdn_image


def interpolate_detectors__bicurvilinear_(image, detectors, wavelengths,
                                          axes=(1,2), step=1, order=2,
                                          mode='lstsq'):
    '''
        NOTE: regardless of `step`, only the N±1 bands are aggregated
                to predict the radiance at the central detector in band N
                
        TODO:
        ----
         * Consider how grid coordinates are supplied
         * Aggregate model should be weighted
         - Maybe `step` should only be applied along the column axis (not band)
         - Handle different band interleave orderings
         - Write output in-place, or write to a separate object?
             - A separate object would be especially useful for
                 e.g. comparing with raw signal to estimate new gains & offsets
         - Can we implement at least the inner loop with `numba`,
            e.g. using `@stencil` or `@guvectorize`?
    '''
        
    # Choose matrix equation solver
    solver = {'lstsq': matsolve__lstsq_, 'pinv': matsolve__pinv_}[mode]
    
    steps = (step, step)
    
    n_var = 1
    n_obs = 2*step+1
    
    for band, column in detectors:
        
        # Get column coordinates
        cols = np.arange(column-step, column+step+1)[:,None]
        # Get wavelengths
        wvls = wavelengths[band-step:band+step+1]
        
        # Get neighborhood as a subset of the original array
        S = neighborhood(image, (band, column), steps, axes, write=True)
            
        # Build boolean mask for bad detector
        I = np.ones((n_obs,)).astype(bool)
        I[step] = False
        
        # Generate polynomial feature combinations
        combos = polynomial_combinations(n_var, order=order)
        # Convert column coordinates to polynomial features
        P = polynomial_matrix(cols, combos)
        
        for line in S:
            
            # Solve best-fit coefficients for N±1 bands
            # M0 = solver(P, line[0:step].T)
            # M2 = solver(P, line[step+1:].T)
            M0 = solver(P, line[step-1:step].T)
            M2 = solver(P, line[step+1:step+2].T)
            
            # Calculate weights for neighboring bands
            wts = wvls[[step-1, step+1]]
            wts /= wts.sum()
            
            # Get aggregate coefficients
            M1 = (np.concatenate([M0, M2], axis=1)*wts).sum(axis=1)
            
            # Predict cross-track radiance curve
            T = polynomial_matrix(cols, combos).dot(M1)
            
            # Compute mean of differences between unmasked radiance values
            #  in central band relative to the predicted curve
            R = (line[step:step+1,I]-T[I]).mean()
            
            # Predict radiance at the central detector
            X = polynomial_matrix(cols[step,None], combos).dot(M1)+R
            
            ### TODO: Check for datatype conflicts!
            
            line[~I,~I] = X
            