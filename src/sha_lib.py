# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 14:46:31 2019

@author: djk
"""

import numpy as np


def gd2gc(h, gdcolat):
    """Geodetic to geocentric colatitude calculation.

    Use WGS84 values for equatorial radius and reciprocal flattening

    """
    eqrad = 6378.137
    flat  = 1/298.257223563
    plrad = eqrad*(1-flat)
    ctgd  = np.cos(np.deg2rad(gdcolat))
    stgd  = np.sin(np.deg2rad(gdcolat))
    a2    = eqrad*eqrad
    a4    = a2*a2
    b2    = plrad*plrad
    b4    = b2*b2
    c2    = ctgd*ctgd
    s2    = 1-c2
    rho   = np.sqrt(a2*s2 + b2*c2)
    rad   = np.sqrt(h*(h+2*rho) + (a4*s2+b4*c2)/rho**2)
    cd    = (h+rho)/rad
    sd    = (a2-b2)*ctgd*stgd/(rho*rad)
    cthc  = ctgd*cd - stgd*sd           # Also: sthc = stgd*cd + ctgd*sd
    thc   = np.rad2deg(np.arccos(cthc)) # arccos returns values in [0, pi]
    return((rad, thc, sd, cd))


def rad_powers(n, a, r):
    """Calculate values of (a/r)^(n+2) for n=0, 1, 2 ..., nmax"""
    arp = np.zeros(n+1)
    t0  = a/r
    arp[0] = t0*t0
    for i in range(n):
        arp[i+1] = t0*arp[i]
    return(arp)


def csmphi(m,phi):
    """Populate arrays with cos(m*phi), sin(m*phi)"""
    cmp = np.zeros(m+1)
    smp = np.zeros(m+1)
    cmp[0] = 1
    smp[0] = 0
    cp = np.cos(np.deg2rad(phi))
    sp = np.sin(np.deg2rad(phi))
    for i in range(m):
        cmp[i+1] = cmp[i]*cp - smp[i]*sp
        smp[i+1] = smp[i]*cp + cmp[i]*sp
    return((cmp,smp))


def gh_phi_rad(gh, nmax, cp, sp, rp):
    """Populate arrays with terms such as g(3,2)*cos(2*phi)*(a/r)**5 """
    rx = np.zeros(nmax*(nmax+3)//2+1)
    ry = np.zeros(nmax*(nmax+3)//2+1)
    igx=-1
    igh=-1
    for i in range(nmax+1):
        igx += 1
        igh += 1
        rx[igx]= gh[igh]*rp[i]
        for j in range(1,i+1):
            igh += 2
            igx += 1
            rx[igx] = (gh[igh-1]*cp[j] + gh[igh]*sp[j])*rp[i]
            ry[igx] = (gh[igh-1]*sp[j] - gh[igh]*cp[j])*rp[i]
    return((rx, ry))


def gh_phi(gh, nmax, cp, sp):
    """As for gh_phi_rad but without the radial function"""
    rx = np.zeros(nmax*(nmax+3)//2+1)
    igx=-1
    igh=-1
    for i in range(nmax+1):
        igx += 1
        igh += 1
        rx[igx]= gh[igh]
        for j in range(1,i+1):
            igh += 2
            igx += 1
            rx[igx] = (gh[igh-1]*cp[j] + gh[igh]*sp[j])
    return(rx)


def pnmindex(n,m):
    """Index for terms of degree=n and order=m in arrays pnm, xnm, ynm and znm"""
    return(n*(n+1)//2+m)


def gnmindex(n,m):
    if(m==0):
        igx = n*n
    else:
        igx = n*n+2*m-1
    return(igx)


def hnmindex(n,m):
    return(n*n+2*m)


def pnm_calc(nmax, th):
    """Calculate arrays of the Associated Legendre Polynomials pnm"""

    # Initialise
    nel   = nmax*(nmax+3)//2+1
    pnm   = np.zeros(nel)
    ct    = np.cos(np.deg2rad(th)); st = np.sin(np.deg2rad(th))
    pnm[0] = 1
    if(nmax==0): return(pnm)
    pnm[1] = ct
    pnm[2] = st

    for i in range(2,nmax+1): # Loop over degree
        idx0 = pnmindex(i,i)
        idx1 = pnmindex(i-1,i-1)
        t1   = np.sqrt(1-1/(2*i))
        pnm[idx0] = t1*st*pnm[idx1]       
        for j in range(i):   # Loop over order
            idx0 = pnmindex(i,j)
            idx1 = pnmindex(i-1,j)
            idx2 = pnmindex(i-2,j)
            t1 = (2*i-1)
            t2 = np.sqrt((i-1+j)*(i-1-j))
            t3 = np.sqrt((i+j)*(i-j))
            pnm[idx0] = (t1*ct*pnm[idx1] - t2*pnm[idx2])/t3           
    return(pnm)


def pxyznm_calc(nmax, th):
    """Calculate arrays of the Associated Legendre Polynomials pnm and the related 
    values xnm, ynm and znm which are needed to compute the X, Y and Z
    geomagnetic field components
    Reverse sign on znm: 16 April 2019
    """
    # Initialise
    nel   = nmax*(nmax+3)//2+1
    pnm   = np.zeros(nel); xnm = np.zeros(nel)
    ynm   = np.zeros(nel); znm = np.zeros(nel)
    ct    = np.cos(np.deg2rad(th)); st = np.sin(np.deg2rad(th))
    pnm[0] = 1
    if(nmax==0): return(pnm)
    pnm[1] =  ct; pnm[2] = st
    xnm[0] = 0; xnm[1] = -st; xnm[2] = ct
    ynm[2] = 1
    znm[0] = 1; znm[1] = 2*ct; znm[2] = 2*st #Check znm[0]
    eps    = 10**(-6)

    for i in range(2,nmax+1): # Loop over degree
        idx0 = pnmindex(i,i)
        idx1 = pnmindex(i-1,i-1)
        t1   = np.sqrt(1-1/(2*i))
        pnm[idx0] = t1*st*pnm[idx1]
        xnm[idx0] = t1*(st*xnm[idx1] + ct*pnm[idx1])
        znm[idx0] = (i+1)*pnm[idx0]
        if(np.abs(st) > eps):
            ynm[idx0] = i*pnm[idx0]/st
        else:
                ynm[idx0] = xnm[idx0]*ct

        for j in range(i):   # Loop over order
            idx0 = pnmindex(i,j)
            idx1 = pnmindex(i-1,j)
            idx2 = pnmindex(i-2,j)
            t1 = (2*i-1)
            t2 = np.sqrt((i-1+j)*(i-1-j))
            t3 = np.sqrt((i+j)*(i-j))
            pnm[idx0] = (t1*ct*pnm[idx1] - t2*pnm[idx2])/t3
            xnm[idx0] = (t1*(ct*xnm[idx1] - st*pnm[idx1]) - t2*xnm[idx2])/t3
            znm[idx0] = (i+1)*pnm[idx0]
            if(np.abs(st) > eps):
                ynm[idx0] = j*pnm[idx0]/st
            else:
                ynm[idx0] = xnm[idx0]*ct

    return((pnm, xnm, ynm, -znm))


def shm_calculator(gh, nmax, altitude, colat, long, coord):
    """Compute values of the geomagnetic field from a model gh of 
    maximum degree and order nmax for either geodetic ot geocentric coordinates.
    """
    RREF     = 6371.2
    degree   = nmax
    phi      = long

    if (coord == 'Geodetic'):
        # Geodetic to geocentric conversion
        rad, theta, sd, cd = gd2gc(altitude, colat)
    else:
        rad   = altitude
        theta = colat

    # Create array with  values of (a/r)^(n+2) for n=0,1, 2 ..., degree
    rpow = rad_powers(degree, RREF, rad)
    # Create arrays with cos(m*phi), sin(m*phi)
    cmphi, smphi = csmphi(degree,phi)
    # Create arrays with terms such as g(3,2)*cos(2*phi)*(a/r)**5 
    ghxz, ghy = gh_phi_rad(gh, degree, cmphi, smphi, rpow)
    # Calculate arrays of the Associated Legendre Polynomials
    pnm, xnm, ynm, znm = pxyznm_calc(degree, theta)
    # Calculate geomagnetic field components are calculated as a dot product
    X =  np.dot(ghxz, xnm)
    Y =  np.dot(ghy,  ynm)
    Z =  np.dot(ghxz, znm)
    # Convert back to geodetic (X, Y, Z) if required
    if (coord == 'Geodetic'):
        t = X
        X = X*cd + Z*sd
        Z = Z*cd - t*sd

    return((X, Y, Z))
