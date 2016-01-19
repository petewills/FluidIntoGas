__author__ = 'peterwills'
import math as math

q = 500.0 / (3600.0 * 24.0)                                       # flow rate m3/sec
mu = 1.0 * 9.0 / 10000.0                                          # viscosity of water in Pa-s at 20C
mu = mu * 1.4                                                     # to explain data
k = 9.869 / 1000.0 / 1000.0 / 1000.0 / 1000.0 / 10.0              # permeability in m2. equal to 1 D.
k *= 0.1                                                          # to explain data
mu_4k = mu / (4.0*k)

h = 10.0                                     # reservoir thickness(m)
phi = 0.28                                   # porosity
cr = 1.0 / 1000.0 / 1000.0 / 1000.0          # compressibility reservoir field example(Pa-1)
cr *= 10.0                                 # This can be big if there is gas

D = phi * mu_4k * cr                         # phi * mu * C / ( 4.8 * k )        units are s / m2
p0 = q * mu_4k / (3.1415 * h)                # Pressure multiplies exponential integral in Pa
pi = 4.0*1000.0 * 1000.0                     # background reservoir pressure(Pa) prior to injection

nu = k / (phi * mu * cr)                     # fluid diffusivity constant. multiplies spatial term.
qnorm = q * mu / (2.0 * math.pi * k * h)     # injection provides boundary condition at r=0.

fault = 'YES'                                # uses mirror well to produce a N-S strike fault at well

gsat = 0.2                                  # Initial gas saturation
tankrad = 250.0                             # radius of stimulated reservoir
