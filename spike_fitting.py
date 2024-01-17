#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 20:05:14 2023

@author: jackmustafa
"""

def spiked_wave_fit(xdata,ydata,attempts=1,lin=True):
    """
    Calculates the best-fit spike-permitting waveform of an input array.

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

    Returns dictionary object containing the best-fit parameters (A, phi, beta and c) and a np.array containing the values of the best-fit waveform evaluated at each point in xdata.
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
    
    # Define internal function to linearise to spike value
    def lin_spike(beta_orig):
        poly = np.array([1.47098869e-02, 0, -1.75130704e-01, 0, 8.91851778e-01, 0, -2.53500696e+00, 0, 4.40870944e+00, 0, 
                         -4.84510857e+00, 0, 3.38301987e+00, 0, -1.48042804e+00, 0, 5.68165374e-01, 0,  1.95859171e-02, 0])
        scaler = np.polyval(poly,np.pi/2)
        poly = poly / scaler
        power = np.arange(len(poly)-1,-1e-6,-1)
        if np.max(np.abs(beta_orig)) > np.pi/2:
            return "Non-linearised beta values outside of -pi/2 to +pi/2 range not accepted."
        beta_lin = sum(poly[i]*beta_orig**power[i] for i in range(len(poly)))
        return beta_lin
    
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
    
    def wave_eval_spike(x,A,phi,beta_0):
        """
        Evaluates spike-permitting waveform function with given parameters.
        
        Parameters
        ----------
        x : 1-D np.array-like object.
            Equally spaced independent data, length N.
        A : Float
            Spike-permitting waveform amplitude.
        phi : Float
            Spike-permitting waveform phase.
        beta_0 : Float
            Spike-permitting waveform spike. Accepts values between and including -pi/2 and +pi/2.
        
        Returns np.array object of length N containing the values of the spike-permitting waveform evaluated at each point in x.
        -------
        """
        
        # If beta_0 = 0, default to first diurnal harmonic
        if beta_0 == 0:
            return wave_eval_FDH(x, A, phi)
        
        # Calculate approximate offset
        offset = A * lin_spike(beta_0)
        
        # If beta_0 < 0, invert parameters
        if beta_0 < 0:
            A = -A
            phi = (phi - T/2)%T
            beta_0 = -beta_0
        
        # Return spike-permitting waveform function output
        return ( A * np.cos( np.pi/beta_0 * 
                              np.arctan( (2/T * ((x - phi - T/2)%T) -1) * np.tan(beta_0)) )
                    + offset)
    
    # Set the best-fit parameter bounds
    bounds = ((0, -T, -np.pi/2), (np.inf, 2*T, np.pi/2)) # allowing extra space for phi to prevent entrapment at 0 or T, return to 0-T range later
    
    # Perform best-fit optimisation
    attempts_left = attempts # a maximum of kwarg "attempts" attempts are allowed for the optimization to succeed, with the initial guess of phase altered between attempts
    phi_adj = 0
    while attempts_left > 0:
        try:
            popt_spike, pcov_spike = curve_fit(wave_eval_spike, xdata_np, ydata_anom, [FDH_A, (FDH_phi + phi_adj) % T, 0.0001], bounds=bounds)
            break
        except RuntimeError:
            attempts_left = attempts_left - 1
            if attempts_left == 0:
                return 'Optimisation procedure failed to converge. Use kwarg "attempts" to increase the number of initial guess scenarios trialled.'
            else:
                phi_adj = phi_adj + T/attempts
    
    # Assign best fit output
    [SPIKE_A, SPIKE_phi, SPIKE_beta_0] = [*popt_spike]
    SPIKE_phi = (SPIKE_phi % T) + xshift
    
    # Evaluate the best-fit waveform along the input xdata
    y_out = wave_eval_spike(xdata_np, SPIKE_A, SPIKE_phi, SPIKE_beta_0)
    
    # Calculate the linearised the spike value (or return output if retaining original spike value)
    if not lin:
        return {'A': SPIKE_A, 'phi': SPIKE_phi, 'beta_0': SPIKE_beta_0, 'c': c, 'y_out': y_out+c}
    else:
        # Use spike linearisation function to calculate linearised spike
        SPIKE_beta = lin_spike(SPIKE_beta_0)
    
    # Return output with linearised spike value
    return {'A': SPIKE_A, 'phi': SPIKE_phi, 'beta': SPIKE_beta, 'c': c, 'y_out': y_out+c}