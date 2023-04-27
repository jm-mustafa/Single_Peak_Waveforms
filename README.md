# Single_Peak_Waveforms
 Public repo with cyclic form-fitting code

This repository contains two Python scripts: 'skew_fitting.py' and 'spike_fitting.py'.

Script dependencies: numpy, scipy(.optimize).

NOTE (27/Apr/2023): These scripts have not been rigorously tested. Please contact me at j.mustafa@uea.ac.uk if you encounter problems.



'skew_fitting.py' : Contains the function 'skewed_wave_fit' which calculates the best-fit skew-permitting waveform for input data.

	Requires two positional input arguments: (xdata, ydata).
		'xdata' : a 1-D np.array-like object containing equally-spaced x-coordinates.
		'ydata' : a 1-D np.array-like object of equal length to xdata containing the corresponding y-coordinates.
	Accepts two keyword input arguments:
		'attempts' : positive integer specifying the maximum number of initial guesses that will be attempted if the optimization process fails to converge. The default is 1.
		'lin' : Boolean specifying whether the output skew parameter should be converted to a linearised value (recommended), or retained in its raw mathematical value (for purposes directly related to the function mathematics). The default is True.
	Returns a dictionary containing:
		A     : the best-fit amplitude parameter, A.
		phi   : the best-fit phase parameter, phi.
		alpha : the rescaled best-fit skew parameter, alpha.
		c     : the mean value of the input ydata.
		y_out : a 1-D np.array object containing the values of the best-fit waveform evaluated for the input xdata.



'spike_fitting.py' : Contains the function 'spiked_wave_fit' which calculates the best-fit spike-permitting waveform for input data.

	Requires two positional input arguments: (xdata, ydata).
		'xdata' : a 1-D np.array-like object containing equally-spaced x-coordinates.
		'ydata' : a 1-D np.array-like object of equal length to xdata containing the corresponding y-coordinates.
	Accepts two keyword input arguments:
		'attempts' : positive integer specifying the maximum number of initial guesses that will be attempted if the optimization process fails to converge. The default is 1.
		'lin' : Boolean specifying whether the output spike parameter should be converted to a linearised value (recommended), or retained in its raw mathematical value (for purposes directly related to the function mathematics). The default is True.
	Returns a dictionary containing:
		A       : the best-fit amplitude parameter, A.
		phi     : the best-fit phase parameter, phi.
		beta    : the rescaled best-fit skew parameter, beta.
		g_beta0 : the best-fit vertical readjustment function parameter, g_beta0.
		c       : the mean value of the input ydata.
		y_out   : a 1-D np.array object containing the values of the best-fit waveform evaluated for the input xdata.