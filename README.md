# Single_Peak_Waveforms
 Public repo with cyclic form-fitting code

This repository contains two Python scripts: 'skew_fitting.py' and 'spike_fitting.py'.

Script dependencies: numpy, scipy. Ensure numpy is imported (as 'np' or 'numpy') and scipy.optimize.curve_fit is imported (as 'curve_fit') for efficiency.

NOTE (10/Feb/2025): Further details on the waveforms and the form-fitting approach are found in Section 5 of Mustafa et al., 2024 (doi:10.1002/qj.4725). These scripts have been created with care, but have not been rigorously tested! Please contact me at j.m.mustafa@leeds.ac.uk if you encounter problems.



'skew_fitting.py' : Contains the function 'skewed_wave_fit' which calculates the best-fit skew-permitting waveform for input data.

	skewed_wave_fit(xdata, ydata, attempts=1, lin=True)

	Requires two positional input arguments: (xdata, ydata).
		'xdata' : a 1-D np.array-like object containing equally-spaced x-coordinates.
		'ydata' : a 1-D np.array-like object of equal length to xdata containing the corresponding y-coordinates.
	Accepts two keyword input arguments:
		'attempts' : positive integer specifying the maximum number of initial guesses that will be attempted if the optimization process fails to converge. The default is 1.
		'lin' : Boolean specifying whether the output skew parameter should be converted to a linearised value (recommended), or retained in its raw mathematical value (for purposes directly related to the function mathematics). The default is True.
	Returns a dictionary containing:
		A     : the best-fit amplitude parameter, A.
		phi   : the best-fit phase parameter, phi.
		alpha : the rescaled best-fit skew parameter, alpha. (OR the original best-fit skew parameter, alpha_0, if lin=False.)
		c     : the mean value of the input ydata.
		y_out : a 1-D np.array object containing the values of the best-fit waveform evaluated for the input xdata.



'spike_fitting.py' : Contains the function 'spiked_wave_fit' which calculates the best-fit spike-permitting waveform for input data.

	spiked_wave_fit(xdata, ydata, attempts=1, lin=True)

	Requires two positional input arguments: (xdata, ydata).
		'xdata' : a 1-D np.array-like object containing equally-spaced x-coordinates.
		'ydata' : a 1-D np.array-like object of equal length to xdata containing the corresponding y-coordinates.
	Accepts two keyword input arguments:
		'attempts' : positive integer specifying the maximum number of initial guesses that will be attempted if the optimization process fails to converge. The default is 1.
		'lin' : Boolean specifying whether the output spike parameter should be converted to a linearised value (recommended), or retained in its raw mathematical value (for purposes directly related to the function mathematics). The default is True.
	Returns a dictionary containing:
		A       : the best-fit amplitude parameter, A.
		phi     : the best-fit phase parameter, phi.
		beta    : the rescaled best-fit skew parameter, beta. (OR the original best-fit spike parameter, beta_0, if lin=False.)
		c       : the mean value of the input ydata.
		y_out   : a 1-D np.array object containing the values of the best-fit waveform evaluated for the input xdata.
