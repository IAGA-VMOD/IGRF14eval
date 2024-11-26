#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 22:21:03 2023

@author: Ciaran Beggan (British Geological Survey, UK)

Functions for computing metrics of Gauss coefficients and file checking. 

Contributions from D. Kerridge and E. Thebault

"""
import numpy as np
import os
from src import sha_lib as sha
r2d = np.rad2deg
d2r = np.deg2rad

def msum(coeffs, nmax=13):
    """
    Compute the spatial power spectrum at the Earth's standard radius 6371.2 km

    Parameters
    ----------
    coeffs : ndarray, shape (..., N*(N+2))
        Spherical harmonic coefficients for degree `N`.
    nmax : int, optional
        Maximum spherical degree (defaults to `N`).


    Returns
    -------
    P : array 
        Power spectrum for degrees up to degree ``nmax``

    Notes
    -----
    The spatial power spectrum for a potential field is defined as
    dP_i = sqrt(sum_1^N gi^2 + hi^2)

    References
    ----------
    Lowes, F. J. (1974). Spatial power spectrum of the main geomagnetic field, 
    and extrapolation to the core. 
    Geophysical Journal International, 36(3), 717-730.

    """

    R = np.empty(nmax)
    l1 = 1
    
    for n in range(1, nmax+1):
        R[n-1] = 0
        l2 = l1 + 2*n
        for m in range(l1,l2+1):
            R[n-1] = R[n-1] + coeffs[m-1]**2
            
        l1 = l1 + 2*n +1 

    return R

def msumnsq(coeffs, nmax=13):
    """
    Compute the spatial unsquared 
    power spectrum at the Earth's standard radius 6371.2 km

    Parameters
    ----------
    coeffs : ndarray, shape (..., N*(N+2))
        Spherical harmonic coefficients for degree `N`.
    nmax : int, optional
        Maximum spherical degree (defaults to `N`).


    Returns
    -------
    P : array 
        Power spectrum for degrees up to degree ``nmax``

    Notes
    -----
    The spatial power spectrum for a potential field is defined as
    dP_i = sqrt(sum_1^N gi^2 + hi^2)

    References
    ----------
    Lowes, F. J. (1974). Spatial power spectrum of the main geomagnetic field, 
    and extrapolation to the core. 
    Geophysical Journal International, 36(3), 717-730.

    """

    R = np.empty(nmax)
    l1 = 1
    
    for n in range(1, nmax+1):
        R[n-1] = 0
        l2 = l1 + 2*n
        for m in range(l1,l2+1):
            R[n-1] = R[n-1] + coeffs[m-1]
            
        l1 = l1 + 2*n +1 

    return R

def check_cof_file_format(files, main_sv):
    """
    Check format of IGRF candidate files
    File should have three header lines starting with #
    Line 4 onwards should have four or six columns
    The final n should be 8 or 13 depending whether it is a SV or DGRF/IGRf
    candidate model
    
    

    Parameters
    ----------
    files : list of file names
        A list of the filenames of the candidate model.
        
    main_sv: str
        'main' - a degree 13 model having 13*15 = 195 coefficients
        'sv' - a degree 8 model have 8*10 = 80 coefficients
        

    Returns
    -------
    pass_fail:
        list of Boolean pass (True) or fail (False) test of each file format

    """
    pass_or_fail = np.empty(len(files))
    
    # Walk through each file
    for l in range(0, len(files)):

        # Assume the file format is accepted
        pass_or_fail[l] = True
        count = 0

        fp = open(files[l], 'r')
        Lines = fp.readlines()

        # Test the first 3 header lines have a #
        while count < 3:   
            if Lines[count][0] != '#':
               pass_or_fail[l] = False
               print('Missing # in header in file: ' + files[l] )

            count += 1
                    

        # Test last value is either 8 for SV or 13 for main
        last_line = Lines[-1].split()       
        if main_sv == 'main':
            if not(int(last_line[0]) == 13):
                pass_or_fail[l] = False
                print('Missing degrees e.g. 13 in file: ' + files[l] + '?')

        
        if main_sv == 'sv':
            if not(int(last_line[0]) == 8):
                pass_or_fail[l] = False
                print('Missing degrees e.g. 8 in file: ' + files[l] + '?')

                
    return pass_or_fail
    

            
def load_coeffs_to_numpy(files, main_sv, model_acronym):
    """
    Load the candidate model files into a pandas dataframe
    
    
    Parameters
    ----------
    files : list of file names
        A list of the filenames of the candidate model.
        
    main_sv: str
        'main' - a degree 13 model having 13*15 = 195 coefficients
        'sv' - a degree 8 model have 8*10 = 80 coefficients
        
    model_acronym: str
        'DGRF', 'IGRF' or 'SV' for stripping off the model type from the 
        institute name

    Returns
    -------
    coeffs: numpy array
        numpy array containing colunms of Gauss coefficients 
        
    institute: str list 
        names of institute or team from the filenames
        
    degree: int
        degree of the model whether main or SV
        
    """
    
    # Create an empty numpy array for the coefficients
    if main_sv == 'main':
        degree = 13
        coeffs = np.empty( (degree*(degree+2),len(files)) )
    else:
        degree = 8 
        coeffs = np.empty( (degree*(degree+2),len(files)) )
    
    institute = []
    
    # Walk through each file
    for l in range(0, len(files)):
        
        data = np.loadtxt(files[l])
        
        
        # As a shortcut, stack the n,m values, then remove np.where(nm=0)
        nm = np.stack((data[:,0],data[:,1]),axis=1)
        nm = np.reshape(nm, (len(nm)*2,1))
        
        # Extract the coefficients only
        gh = np.stack((data[:,2],data[:,3]),axis=1)
        gh = np.reshape(gh, (len(gh)*2,1))
        
        gh = gh[np.where(nm!=0)]
        
        coeffs[:,l] = gh
        
        # Isolate the institute from the filename
        institute.append( os.path.splitext(
            os.path.basename(files[l]))[0].replace(model_acronym + '_',''))
        
    return coeffs, institute, degree
            
def shm_calculator(gh, nmax, altitude, colat, long, coord):
    
    """
    Compute the magnetic field values at lat/long/altitude for 
    geodetic or geocentric frame
    
    
    Parameters
    ----------
    gh : numpy array
        Gauss coefficients of the magnetic field with
        a (with monopole) in the usual arrangement
        
    nmax: int
        maximum degree to evaluate the Gauss coefficients to
        
    altitude: numpy array
        array of altitudes in km above the Earth's centre for geocentric
        or km above/below the standard radius in km
        
    colat: numpy array
        array of colatitude values (90-latitude)
        
    long: numpy array
        array of longitudes
    
    coord: str
        'Geodetic' or 'Geocentric' 
        convert colat and altitude to geocentric coords if geodetic
     
    Returns
    -------
    X, Y, Z: numpy array
        values of the magnetic field three orthogonal components
        
    """
    
    RREF     = 6371.2 #The reference radius assumed by the IGRF
    degree   = nmax
    phi      = long

    if (coord == 'Geodetic'):
        # Geodetic to geocentric conversion using the WGS84 spheroid
        rad, theta, sd, cd = sha.gd2gc(altitude, colat)
    else:
        rad   = altitude
        theta = colat

    # Function 'rad_powers' to create an array with values of (a/r)^(n+2) 
    # for n = 0,1, 2 ..., degree
    rpow = sha.rad_powers(degree, RREF, rad)

    # Function 'csmphi' to create arrays with cos(m*phi), sin(m*phi) 
    # for m = 0, 1, 2 ..., degree
    cmphi, smphi = sha.csmphi(degree,phi)

    # Function 'gh_phi_rad' to create arrays with terms such 
    # as [g(3,2)*cos(2*phi) + h(3,2)*sin(2*phi)]*(a/r)**5 
    ghxz, ghy = sha.gh_phi_rad(gh, degree, cmphi, smphi, rpow)

    # Function 'pnm_calc' to calculate arrays of the 
    # Associated Legendre Polynomials for n (&m) = 0,1, 2 ..., degree
    pnm, xnm, ynm, znm = sha.pxyznm_calc(degree, theta)

    # Geomagnetic field components are calculated as a dot product
    X =  np.dot(ghxz, xnm)
    Y =  np.dot(ghy,  ynm)
    Z = -np.dot(ghxz, znm)

    # Convert back to geodetic (X, Y, Z) if required
    if (coord == 'Geodetic'):
        t = X
        X = X*cd + Z*sd
        Z = Z*cd - t*sd

    return((X, Y, Z))

def xyz2dhif(x, y, z):
    """Calculate D, H, I and F from (X, Y, Z)
      
    Input parameters
    ---------------
    X: north component (nT) 
    Y: east component (nT)
    Z: vertical component (nT)
    
    Output
    ------
    A tuple: (D, H, I, F)
    D: declination (degrees)
    H: horizontal intensity (nT)
    I: inclination (degrees)
    F: total intensity (nT)
    
    Dependencies
    ------------
    numpy
         
    Revision date
    -------------
    28 April 2019

    """
    hsq = x*x + y*y
    hoz  = np.sqrt(hsq)
    eff = np.sqrt(hsq + z*z)
    dec = np.arctan2(y,x)
    inc = np.arctan2(z,hoz)
    return((r2d(dec), hoz, r2d(inc), eff))


def xyz2dhif_sv(x, y, z, xdot, ydot, zdot):
    """Calculate secular variation in D, H, I and F from (X, Y, Z) and
    (Xdot, Ydot, Zdot)
      
    Input parameters
    ---------------
    X: north component (nT), and Xdot=dX/dt 
    Y: east component (nT), and Xdot=dX/dt 
    Z: vertical component (nT), and Xdot=dX/dt 
    
    Output
    ------
    A tuple: (Ddot, Hdot, Idot, Fdot)
    Ddot: rate of change of declination (degrees/year)
    Hdot: rate of change of horizontal intensity (nT/year)
    Idot: rate of change of inclination (degrees/year)
    Fdot: rate of change of total intensity (nT/year)
    
    Dependencies
    ------------
    numpy
         
    Revision date
    -------------
    28 April 2019

    """
    h2  = x*x + y*y
    h   = np.sqrt(h2)
    f2  = h2 + z*z
    hdot = (x*xdot + y*ydot)/h
    fdot = (x*xdot + y*ydot + z*zdot)/np.sqrt(f2)
    ddot = r2d((xdot*y - ydot*x)/h2)*60
    idot = r2d((hdot*z - h*zdot)/f2)*60
    return((ddot, hdot, idot, fdot))


def ra_spectrum(nmax, nmin, avg_window, coeffs):
    """
    Compute the magnetic field azimuthal spectrum
    
    Based on Matlab code from E. Thebault (U. Nantes)
    
    Parameters
    ----------
    
    nmax: int
        maximum degree to evaluate the Gauss coefficients to
        
    nmin: int
        minimum degree to evaluate the Gauss coefficients from
        
    avg_window: float
        averaging window for the m/n width
    
    coeffs: numpy array
        Gauss coefficients in the usual arrangement
     
    Returns
    -------
    Xa: numpy array
        values of the m/n azimuthal ratio averages
    
    Ra: numpy array
            values of the azimuthal spectrum
    
        
    Dependencies
    ------------
    numpy
         
    Revision date
    -------------
    20-Dec-2023
    
    """
    
    noff = nmin**2 -1
    ncof = 0
    ar = np.zeros(((nmax+1)**2-1,1))
    pr = np.zeros(((nmax+1)**2-1,1))
    
    for n in range(nmin, nmax+1):
        for m in range(0, n+1):
            if m == 0:
                ar[ncof] = 0
                pr[ncof] = (n+1) * coeffs[noff+ncof]**2
                
            else:
                # positive side
                ncof += 1
                ar[ncof] = m/n
                pr[ncof] = (n+1) * coeffs[noff+ncof]**2
               
                # negative side
                ncof += 1
                ar[ncof] = -m/n
                pr[ncof] = (n+1) * coeffs[noff+ncof]**2
        
        ncof += 1
    
    # Use a 'mergesort' to create a correctly ordered list of occurrence of the
    # sorted values. It turns out 'quicksort' {default} and 'heapsort' do not
    # retain the order of occurence - they start at the end of the list or
    # in the middle. Not the required behaviour for azimuthal sorting
    po = pr[np.argsort(ar,axis=0, kind='mergesort')]
    po = np.squeeze(po)
    ao = np.sort(ar, axis=0)
    
    bins = np.arange(-1-avg_window, 1+avg_window, avg_window)
    nbins = len(bins)-2
    xa = np.zeros((nbins,1))
    ra = np.zeros((nbins,1))
    
    # Now, as a sad hack, we have to forcibly truncate the decimal places of 
    # the bins and ao arrays, as >= or <= does not match the matlab 
    # implementation. *sighs loudly*
    bins = np.around(bins, 4)
    ao = np.around(ao, 4)
    
    for i in range(nbins):
        bin_start = np.max((-1, bins[i]))
        bin_end = np.min((1, bins[i+2]))
        psum = 0
        for j in range(ncof):
            if (ao[j] >= bin_start) and (ao[j] <= bin_end):
                psum += po[j]
        
        xa[i] = bins[i+1]
        ra[i] = psum
    
    return xa, ra

def gh2triangle(coeffs,nmax):
    """
    Prepare a 2D matrix with the coefficients (or differences) arranged
    by degree in the rows and order by colunms moving out from the centre
        
    Parameters
    ----------
   
    coeffs: numpy array
        Gauss coefficients in the usual arrangement in a 1D array (vector)
        
    nmax: int
        maximum degree to evaluate the Gauss coefficients to
    
     
    Returns
    -------
    tri_gh: numpy array
        values of the m/n arranged ready for a triangle plot
    
        
    Dependencies
    ------------
    numpy
         
    Revision date
    -------------
    20-Dec-2023
    
    """
    
    tri_gh =  np.full((nmax, 2*nmax+1),np.nan)
    
    k = 0
    for n in range(1,nmax+1):
        for m in range(0,n+1):
            tri_gh[n-1, nmax+m] = coeffs[k]
            k += 1
            if m > 0:
                tri_gh[n-1, nmax-m] = coeffs[k]
                k += 1
                
    
    return tri_gh

def write_coefficients_to_cof(filename, coeff_array, uncertainty_array, 
                              degree, header, dp=2):
    """
    Write out the coefficient files 
    Truncate Gauss coefficients to two decimal places unless other
    

    Parameters
    ----------
    filename: str
        path of the file e.g DGRF_Median.cof
    coeff_array : float, n x 1
        numpy array of Gauss coefficients to be written out
    uncertainty_array : float, , n x 1
        numpy array of uncertainties coefficients to be written out
    degree : str, int
        Degree of the file
    header : string
        name of the model for the header e.g. DGRF Median 
    dp : int
        number of decimal places to truncate the numpy array to
        
    

    Returns
    -------
    None.

    """
    # check if degree is an int or str
    
    if degree == 'main':
        degree = 13
    elif degree == 'sv':
        degree = 8
    
    outarray_len = np.cumsum(np.arange(2,degree+2))[-1]
    
    # Generate the ordering of the g/h elements
    gh_labels = []
    for l in range(1,degree+1):
        gh_labels.append(np.stack([np.ones(l+1)*l, np.arange(0,l+1)]).T)
        
    ghl = np.vstack(gh_labels)

    # Check the uncertainty array exists and has values
    # If not, set up an n x 2 array
    if np.size(uncertainty_array) == 0:
        uncertainty_array = np.zeros([outarray_len,2])

                 
    # Format the array in to outarray_len x 2 with the h_1^0 coefficient set to 0
    ghs = np.insert(coeff_array,np.cumsum(np.arange(1,degree+1)*2-1),0)
    gho = np.round(np.reshape(ghs, [outarray_len, 2]), dp)
    
    # Write out the coefficients
    out = np.hstack([ghl, gho, uncertainty_array])
    
    with open(filename,"w") as f:
        f.writelines(['# Candidate model for IGRF-14 \n'])
        f.writelines(['# Analysis for' + header + '\n'])
        f.writelines(['# n  m         gnm         hnm       uncertainty_gnm       uncertainty_hnm \n'])
        np.savetxt(f, out, fmt='  %d  %d   %10.' +str(dp)+ 'f '
                                +  '%10.' +str(dp)+ 'f '
                                + ' %8.' +str(dp)+ 'f ' 
                                + ' %8.' +str(dp)+ 'f', delimiter=' ', 
                                                           newline='\n')
    
    f.close()
    
    
def load_shcfile(filepath):
    """
    Load shc-file and return coefficient arrays.

    Parameters
    ----------
    filepath : str
        File path to spherical harmonic coefficient shc-file.


    Returns
    -------
    time : ndarray, shape (N,)
        Array containing `N` times for each model snapshot 
    coeffs : ndarray, shape (nmax(nmax+2), N)
        Coefficients of model snapshots. Each column is a snapshot up to
        spherical degree and order `nmax`.
    parameters : dict, {'SHC', 'nmin', 'nmax', 'N', 'order', 'step'}
        Dictionary containing parameters of the model snapshots and the
        following keys: ``'SHC'`` shc-file name, `nmin` minimum degree,
        ``'nmax'`` maximum degree, ``'N'`` number of snapshot models,
        ``'order'`` piecewise polynomial order and ``'step'`` number of
        snapshots until next break point. Extract break points of the
        piecewise polynomial with ``breaks = time[::step]``.

    """
    
    read_values_flag = False

    with open(filepath, 'r') as f:

        data = np.array([])
        for line in f.readlines():

            if line[0] == '#':
                continue

            read_line = np.fromstring(line, sep=' ')
            if (read_line.size == 7) and not(read_values_flag):
                name = os.path.split(filepath)[1]  # file name string
                values = [name] + read_line.astype(int).tolist()
                read_values_flag = True
            else:
                data = np.append(data, read_line)

        # unpack parameter line
        keys = ['SHC', 'nmin', 'nmax', 'N', 'order', 'step', 'start_year', 'end_year']
        parameters = dict(zip(keys, values))
        
        time = data[:parameters['N']]
        coeffs = data[parameters['N']:].reshape((-1, parameters['N']+2))
        coeffs = np.squeeze(coeffs[:, 2:])  # discard columns with n and m


    return time, coeffs, parameters