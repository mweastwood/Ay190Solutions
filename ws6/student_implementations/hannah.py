import sys
sys.path.insert(0, "..")
import numpy as np
import matplotlib.pyplot as pl

# Calculates the discrete Fourier transform of a vector x using matrix
# multiplication.
def dft(x):
  # Generate the matrix
  w = np.exp(-2 * np.pi * np.complex(0,1) / len(x))
  w_matrix = np.array([np.array([w ** (j * k) for j in range(len(x))]) \
          for k in range(len(x))])
  # Return its product with the input vector
  return np.dot(w_matrix, x)

# Uncomment to test that the DFT algorithm matches numpy's FFT results.
#test = np.random.random(20) 
#np.testing.assert_allclose(dft(test), np.fft.fft(test))
