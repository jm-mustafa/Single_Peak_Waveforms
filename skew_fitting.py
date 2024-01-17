#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 20:03:54 2023

@author: jackmustafa
"""

def skewed_wave_fit(xdata,ydata,attempts=1,lin=True):
    """
    Calculates the best-fit skew-permitting waveform of an input array.

    Parameters
    ----------
    xdata : 1-D np.array-like object, length N
        Equally spaced independent data.
    ydata : 1-D np.array-like object, length N
        Dependent data.
    attempts : Positive integer, optional
        The maximum number of initial guesses that will be attempted if the optimization process fails to converge. The default is 1.
    lin : Boolean, optional
        Convert the best-fit skew parameter to a linearised value if True. The default is True.

    Returns dictionary object containing the best-fit parameters (A, phi, alpha and c) and a np.array containing the values of the best-fit waveform evaluated at each point in xdata.
    -------
    Dependent on numpy and scipy. The xdata are assumed to be cyclic but without repeated end values.

    """
    
    # Import modules (if not already imported)
    try:
        if type(np) == 'module':
            np = np
    except NameError:
        try:
            if type(numpy) == 'module':
                np = numpy
        except NameError:
            import numpy as np
    try:
        if type(curve_fit) == 'function':
            curve_fit = curve_fit
    except NameError:
        from scipy.optimize import curve_fit
    
    # Standardise input data type
    xdata_np = np.array(xdata)
    ydata_np = np.array(ydata)
    
    # Check xdata is equally spaced
    step = np.diff(xdata_np)[0]
    if all(np.diff(xdata_np) == step) == False:
        return 'xdata input are not equally spaced.'
    
    # Check equal length of xdata and ydata
    len_x = len(xdata_np)
    len_y = len(ydata_np)
    if len_x != len_y:
        return 'xdata and ydata have unequal lengths ({} and {}).'.format(len_x,len_y)
    
    # Shift x-axis to set first x-coordinate to 0, if not already 0
    xshift = xdata_np[0]
    xdata_np = xdata_np - xshift
    
    # Calculate xdata cycle period/wavelength
    T = len_x * step
    
    # Subtract mean
    c = ydata_np.mean()
    ydata_anom = ydata_np - c
    
    # Use Fourier analysis to extract the first diurnal harmonic parameters
    fft_FDH = np.sum(np.exp(-2 * np.pi * xdata_np / T * 1j) * ydata_anom) / len_x
    FDH_A = 2 * np.abs(fft_FDH)
    if FDH_A > 0:
        FDH_phi = ((np.arctan(-np.imag(fft_FDH) / np.real(fft_FDH)) + np.pi * (np.real(fft_FDH) < 0)) * T/(2 * np.pi)) % T
    else:
        FDH_phi = 0
    
    # Define internal function to linearise to skew value
    def lin_skew(alpha):
        poly = np.array([-4.73503677e-02, 0, 5.16562873e-02, 0, -2.60647935e-02, 0, -3.73032305e-01, 0, 1.39603118e+00, 0])
        scaler = np.polyval(poly,1)
        poly = poly / scaler # the polynomial is scaled (with negligible effect) to ensure that the extreme values of -1 and +1 precisely align for the original and linearised skew values
        
        # Solve for linearised skew (alpha) by finding the real root of the polynomial in the -1 to +1 range
        alpha_allroots = np.roots([*poly[:-1], -SKEW_alpha_0])
        alpha_lin = np.nan
        for j in alpha_allroots:
            if np.abs(j) <=1 and np.imag(j) == 0:
                alpha_lin = np.real(j)
        if alpha_lin == np.nan:
            return 'Failed to linearise the best-fit skew value.'
        return alpha_lin
    
    # Define internal functions to evaluate the first diurnal harmonic and skew-permitting waveform for given input parameters
    def wave_eval_FDH(x,A,phi):
        """
        Evaluates first diurnal harmonic (cosine) function with given parameters.
        
        Parameters
        ----------
        x : 1-D np.array-like object.
            Equally spaced independent data, length N.
        A : Float
            First diurnal harmonic amplitude.
        phi : Float
            First diurnal harmonic phase.
        
        Returns np.array object of length N containing the values of the first diurnal harmonic waveform evaluated at each point in x.
        -------
        """
        
        # Calculate the cycle frequency
        k = (2 * np.pi) / T
        
        # Return first diurnal harmonic waveform output
        return A * np.cos(k * (x - phi))
    
    def wave_eval_skew(x,A,phi,alpha_0):
        """
        Evaluates skew-permitting waveform function with given parameters.
        
        Parameters
        ----------
        x : 1-D np.array-like object.
            Equally spaced independent data, length N.
        A : Float
            Skew-permitting waveform amplitude.
        phi : Float
            Skew-permitting waveform phase.
        alpha_0 : Float
            Skew-permitting waveform skew. Accepts values between and including -1 and +1.
        
        Returns np.array object of length N containing the values of the skew-permitting waveform evaluated at each point in x.
        -------
        """
        
        # If alpha_0 = 0, default to first diurnal harmonic
        if alpha_0 == 0:
            return wave_eval_FDH(x, A, phi)
        
        # Calculate the cycle frequency
        k=(2*np.pi)/T
        
        # Return skew-permitting waveform function output
        return ( (A / np.arctan(alpha_0 / (1 - alpha_0 ** 2) ** 0.5)) *
                np.arctan(alpha_0 * np.sin(k * (x - phi) + np.arccos(alpha_0)) /
                          (1 - alpha_0 * np.cos(k * (x - phi) + np.arccos(alpha_0)))
                          )
                )
    
    # Set the best-fit parameter bounds
    bounds = ((0, -T, -1), (np.inf, 2*T, 1)) # allowing extra space for phi to prevent entrapment at 0 or T, return to 0-T range later
    
    # Perform best-fit optimisation
    attempts_left = attempts # a maximum of kwarg "attempts" attempts are allowed for the optimization to succeed, with the initial guess of phase altered between attempts
    phi_adj = 0
    while attempts_left > 0:
        try:
            popt_skew, pcov_skew = curve_fit(wave_eval_skew, xdata_np, ydata_anom, [FDH_A, (FDH_phi + phi_adj) % T, 0], bounds=bounds)
            break
        except RuntimeError:
            attempts_left = attempts_left - 1
            if attempts_left == 0:
                return 'Optimisation procedure failed to converge. Use kwarg "attempts" to increase the number of initial guess scenarios trialled.'
            else:
                phi_adj = phi_adj + T/attempts
    
    # Assign best-fit output
    [SKEW_A, SKEW_phi, SKEW_alpha_0] = [*popt_skew]
    SKEW_phi = (SKEW_phi % T) + xshift
    
    # Evaluate the best-fit waveform along the input xdata
    y_out = wave_eval_skew(xdata_np, SKEW_A, SKEW_phi, SKEW_alpha_0)
    
    # Calculate the linearised the skew value (or return output if retaining original skew value)
    if not lin:
        return {'A': SKEW_A, 'phi': SKEW_phi, 'alpha_0': SKEW_alpha_0, 'c': c, 'y_out': y_out+c}
    else:
        # Use skew linearisation function to calculate linearised skew
        SKEW_alpha = lin_skew(SKEW_alpha_0)
    
    # Return output with linearised skew value
    return {'A': SKEW_A, 'phi': SKEW_phi, 'alpha': SKEW_alpha, 'c': c, 'y_out': y_out+c}
