import numpy as np

def parse_SPC(filename, skip_rows=5):
	dtype = [('p', float),		# pressure, mb
	         ('z', float),		# altitude, m
		 ('T', float),		# temperature, C
		 ('Td', float),		# dewpoint, C
                 ('RH', float),		# relative humidity (%)
                 ('w', float),          # mixing ratio (g/kg)
		 ('wind_dir', float),	# wind direction, degrees
		 ('wind_spd', float),	# wind speed, knots
                 ('theta', float),      # theta, K
                 ('theta-e', float),    # theta-e, K
                 ('theta-v', float),    # theta-v, K

		 ]
	data = np.genfromtxt(filename, dtype=dtype, skip_header=skip_rows, skip_footer=80)
	return data
