import numpy as np
import pyPDWT

v=np.arange(18)

w=pyPDWT.Wavelet1D(3,False,"Daub4")
w.Initialize(v)
w.forward()
w.backward()

