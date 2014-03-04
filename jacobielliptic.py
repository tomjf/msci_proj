import mpmath
import numpy
import matplotlib.pyplot as plt

dn = mpmath.ellipfun('dn')
sn = mpmath.ellipfun('sn')
data1, data2, data3 = [],[],[]
grads= []

xmin,xmax,xsteps = -5,0.0,100.0
xpts = numpy.linspace(xmin,xmax,xsteps)
h = ((xmax-xmin)/xsteps)
adashdash = [0,0,0]
addova = [None, None, None]
addova_analytical = []
for i in range(int(xsteps)):
	x = xpts[i]
	a = dn(x,1)
	b = sn(x,1)
	data1.append(10*a)
	data2.append(abs(1/x))
	data3.append(10*(b+1))
	if i > 2:
		seconderiv = (data1[i] - 2*data1[i-1] + data1[i-2])/(h*h)
		adashdash.append(seconderiv)
		addova.append(adashdash[i]/a)

plt.plot(xpts, data1, label = "dn fnc")
plt.plot(xpts, adashdash, color = 'm', label = "a''")
plt.plot(xpts, addova, color = 'c', label = "a''/a")
plt.plot(xpts, data2 , color ='r', label = "1/n")
plt.plot(xpts, data3, color = 'g', label = "sn fnc")
plt.legend(loc=3)
plt.show()