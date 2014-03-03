import mpmath
import numpy
import matplotlib.pyplot as plt

dn = mpmath.ellipfun('dn')
data = []
xpts = numpy.linspace(-100,0,100)
for x in xpts:
	a = dn(x, 1)
	print x,a
	data.append(a)

print xpts


plt.plot(xpts, data)
plt.show()