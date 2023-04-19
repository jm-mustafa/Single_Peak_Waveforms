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
        Equally spaced independent data. First entry should be 0.
    ydata : 1-D np.array-like object, length N
        Dependent data.
    attempts : Positive integer, optional
        The maximum number of initial guesses that will be attempted if the optimization process fails to converge. The default is 1.
    lin : Boolean, optional
        Convert the best-fit skew parameter to a linearised value if True. The default is True.

    Returns dictionary object containing the best-fit parameters (A, phi, beta, g_beta0 and c) and a np.array containing the values of the best-fit waveform evaluated at each point in xdata.
    -------
    Dependent on numpy and scipy. The xdata are assumed to be cyclic but without repeated end values.

    """
    
    # Import modules
    import numpy as np
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
    
    def wave_eval_spike(x,A,phi,beta_0,g_beta0,suppress_warn=True):
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
        g_beta0 : Float
            Spike-permitting waveform vertical readjustment function. Should be the value between -1 and +1 that enforces a waveform integral of 0.
        
        Returns np.array object of length N containing the values of the spike-permitting waveform evaluated at each point in x.
        -------
        """
        
        # If beta_0 = 0, default to first diurnal harmonic
        if beta_0 == 0:
            return wave_eval_FDH(x, A, phi)
        
        # If beta_0 < 0, invert parameters
        if beta_0 < 0:
            A = -A
            phi = (phi - T/2)%T
            beta_0 = -beta_0
            g_beta0 = -g_beta0
        
        # Calculate the spike-permitting waveform function output at high resolution to test alignment of beta_0 and g_beta0
        x_hires = np.linspace(0,T,1000)
        y = A * (np.cos( np.pi/beta_0 * 
                          np.arctan( (2/T * ((x_hires - phi - T/2)%T) -1) * np.tan(beta_0)) )
                    + g_beta0)
        
        # Check whether g_beta0 approximately aligns with beta_0 (such that the waveform cycle integral is approximately zero)
        diff = np.mean(y)
        if np.abs(diff) > np.abs(A)/20 and suppress_warn == False:
            print("Warning: unmatched values of spike and the vertical readjustment function result in a mean waveform displacement from the mean of around {:.2e}.".format(diff))
        
        # Return spike-permitting waveform function output
        return A * (np.cos( np.pi/beta_0 * 
                          np.arctan( (2/T * ((x - phi - T/2)%T) -1) * np.tan(beta_0)) )
                    + g_beta0)
    
    # Set the best-fit parameter bounds
    bounds1 = ((0, -T, 0, 0), (np.inf, 2*T, np.pi/2, 1)) # allowing extra space for phi to prevent entrapment at 0 or T, return to 0-T range later
    bounds2 = ((-np.inf, -T, 0, 0), (0, 2*T, np.pi/2, 1)) # inverted bounds with negative amplitude to achieve a narrow minimum waveform
    
    # Perform best-fit optimisation
    attempts_left = attempts # a maximum of kwarg "attempts" attempts are allowed for the optimization to succeed, with the initial guess of phase altered between attempts
    phi_adj = 0
    while attempts_left > 0:
        try:
            popt_spike1, pcov_spike1 = curve_fit(wave_eval_spike, xdata_np, ydata_anom, [FDH_A, (FDH_phi + phi_adj) % T, 0.001, 0], bounds=bounds1)
            popt_spike2, pcov_spike2 = curve_fit(wave_eval_spike, xdata_np, ydata_anom, [-FDH_A, (FDH_phi + T/2 + phi_adj) % T, 0.001, 0], bounds=bounds2)
            break
        except RuntimeError:
            attempts_left = attempts_left - 1
            if attempts_left == 0:
                return 'Optimisation procedure failed to converge. Use kwarg "attempts" to increase the number of initial guess scenarios trialled.'
            else:
                phi_adj = phi_adj + T/attempts
    
    # Compare least-squares residual of each solution and assign best-fit output of best solution
    anom1 = ydata_anom - wave_eval_spike(xdata_np, *popt_spike1)
    anom2 = ydata_anom - wave_eval_spike(xdata_np, *popt_spike2)
    sqrs1 = np.sum(anom1 ** 2)
    sqrs2 = np.sum(anom2 ** 2)
    if sqrs1 <= sqrs2:
        [SPIKE_A, SPIKE_phi, SPIKE_beta_0, SPIKE_g_beta0] = [*popt_spike1]
        print('Pos')
    else:
        [SPIKE_A, SPIKE_phi, SPIKE_beta_0, SPIKE_g_beta0] = [*popt_spike2]
        SPIKE_A = -SPIKE_A
        SPIKE_phi = (SPIKE_phi - T/2) % T
        SPIKE_beta_0 = -SPIKE_beta_0
        SPIKE_g_beta0 = -SPIKE_g_beta0
        print('Neg')
    
    # Evaluate the best-fit waveform along the input xdata
    y_out = wave_eval_spike(xdata_np, SPIKE_A, SPIKE_phi, SPIKE_beta_0, SPIKE_g_beta0, suppress_warn=False)
    
    # Calculate the linearised the spike value (or return output if retaining original spike value)
    if not lin:
        return {'A': SPIKE_A, 'phi': SPIKE_phi, 'beta_0': SPIKE_beta_0, 'g_beta0': SPIKE_g_beta0, 'c': c, 'y_out': y_out+c}
    else:
        # Define the approximate relation between original and linearised spike
        poly = np.array([1.47098869e-02, 0, -1.75130704e-01, 0, 8.91851778e-01, 0, -2.53500696e+00, 0, 4.40870944e+00, 0, 
                     -4.84510857e+00, 0, 3.38301987e+00, 0, -1.48042804e+00, 0, 5.68165374e-01, 0,  1.95859171e-02, 0]) # a best-fit 19th order polynomial was calculated to relate the original skew (alpha_0) to the desired linearised value (alpha)
        scaler = np.polyval(poly,np.pi/2)
        poly = poly / scaler # the polynomial is scaled (with negligible effect) to ensure that the extreme values of -1 and +1 for the linearised spike precisely align with -pi/2 and +pi/2 for the original spike
        power = np.arange(len(poly)-1,-1e-6,-1)
        
        # Solve for linearised spike (beta) by summing the polynomial components at the original spike
        SPIKE_beta = sum(poly[i]*SPIKE_beta_0**power[i] for i in range(len(poly)))
    
    # Return output with linearised spike value
    return {'A': SPIKE_A, 'phi': SPIKE_phi, 'beta': SPIKE_beta, 'g_beta0': SPIKE_g_beta0, 'c': c, 'y_out': y_out+c}
