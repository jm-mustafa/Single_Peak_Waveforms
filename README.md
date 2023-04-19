# Single_Peak_Waveforms
 Public repo with cyclic form-fitting code

This repository contains two Python scripts: 'skew_fitting.py' and 'spike_fitting.py'.

'skew_fitting.py' : Contains the function 'skewed_wave_fit' which calculates the best-fit skew-permitting waveform for input data.

	Accepts two positional input arguments: 'xdata' and 'ydata'.
		'xdata' : a 1-D np.array-like object containing equally-spaced x-coordinates.
		'ydata' : a 1-D np.array-like object of equal length to xdata containing the corresponding y-coordinates.
	Returns a dictionary containing:
		A     : the best-fit amplitude parameter, A.
		phi   : the best-fit phase parameter, phi.
		alpha : the rescaled best-fit skew parameter, alpha.
		c     : the mean value of the input ydata.
		y_out : a 1-D np.array object containing the values of the best-fit waveform evaluated for the input xdata.

'spike_fitting.py' : Contains the function 'spiked_wave_fit' which calculates the best-fit spike-permitting waveform for input data.

	Accepts two positional input arguments: 'xdata' and 'ydata'.
		'xdata' : a 1-D np.array-like object containing equally-spaced x-coordinates.
		'ydata' : a 1-D np.array-like object of equal length to xdata containing the corresponding y-coordinates.
	Returns a dictionary containing:
		A       : the best-fit amplitude parameter, A.
		phi     : the best-fit phase parameter, phi.
		beta    : the rescaled best-fit skew parameter, beta.
		g_beta0 : the best-fit vertical readjustment function parameter, g_beta0.
		c       : the mean value of the input ydata.
		y_out   : a 1-D np.array object containing the values of the best-fit waveform evaluated for the input xdata.