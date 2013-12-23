from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy
import math
global k

#--------------------------------------------------------------------------------------------------------------------------------
## Define Parameter Class so don't have to pass multiple variables between functions and can have preset values ie for inflation
# class ParamsType(object):
# 	w = -1
# 	kspace = numpy.linspace(1,10,100)
# 	n_a = -1000/(numpy.amin(kspace))
# 	n_b = -1/(1000*numpy.amax(kspace))
# 	coeff = 2/((3*w)+1)
# 	C = coeff*(coeff-1)
#-----------------------------------------------
def deriv(y,t):
	a = -((k*k)-(C/(t*t)))
	b = 0
	return numpy.array([ y[1], a*y[0]+b*y[1] ])
#-----------------------------------------------------------------------------------------------

def calc_gradient():
	coefficients = numpy.polyfit(data[:,0], data[:,3], 1)
	polynomial = numpy.poly1d(coefficients)
	return polynomial
#-----------------------------------------------------------------------------------------------------

def RK4(k, n_inside, n_outside, num_steps):
	times = numpy.linspace(n_inside,n_outside,num_steps)
	h = times[1]-times[0]

	data = numpy.zeros((len(times),4))
	xn = 1/(math.sqrt(2*k))*math.cos(k*n_inside)
	yn = -k/(math.sqrt(2*k))*math.sin(k*n_outside)

	# RK4 method of integrating ODE where 2nd order mukhanov-sasaki equation has been re-written as 2 1st order coupled ODEs using: x = v, y = dv/dn
	# --> ODE_1: dx/dn = y
	# --> ODE_2: dy/dn = -[(k^2)+(2/n^2)]x

	for idx, time in enumerate(times):
		n = times[idx]
		if idx == 0:
			data[idx,0], data[idx,1], data[idx,2] = n, xn, yn
		else:
			xk1 = h*yn
			yk1 = -h*((k*k)+(2/n*n))*xn

			xk2 = h*(yn+0.5*yk1)
			yk2 = -h*((k*k)+(2/(n+0.5*h)*(n+0.5*h)))*(xn+0.5*xk1)

			xk3 = h*(yn+0.5*yk2)
			yk3 = -h*((k*k)+(2/(n+0.5*h)*(n+0.5*h)))*(xn+0.5*xk2)

			xk4 = -h*(yn+0.5*yk3)
			yk4 = -h*((k*k)+(2/(n+h)*(n+h)))*(xn+xk3)

			xadd = xk1 +2*xk2 + 2*xk3 + xk4

			xn = xn + (1.0/6.0)*(xk1 + (2*xk2) + (2*xk3) + xk4)
			yn = yn + (1.0/6.0)*(yk1 + (2*yk2) + (2*yk3) + yk4)

			data[idx,0], data[idx,1], data[idx,2] = n, xn, yn

			# Trying to plot analytical solution for de-Sitter space ('http://www.damtp.cam.ac.uk/user/db275/Cosmology/Chapter5.pdf' eqn 5.1.23 but when alpha=beta=1)
			data[idx,3] = (2/math.sqrt(2*k))*(math.cos(k*n)-(1/(k*n))*math.sin(k*n))
	return data
#-----------------------------------------------------------------------------------------------------
k = 1


kspace = numpy.linspace(1,10,100)
n_inside = -1000/(numpy.amin(kspace))
n_outside = -1/(1000*numpy.amax(kspace))
num_steps = 100000
times = numpy.linspace(n_inside,n_outside,num_steps)

w = -1
coeff = 2/((3*w)+1)
C = coeff*(coeff-1)

# odeint method for comparison with RK4 and analytical methods
xinit = numpy.array([1/(math.sqrt(2*k))*math.cos(k*n_inside), -k/(math.sqrt(2*k))*math.sin(k*n_inside)])
odemethod = odeint(deriv,xinit,times) 

data = RK4(k,n_inside,n_outside,num_steps)

plt.plot(data[:,0], data[:,1], color = 'r')
# plt.scatter(data[:,0], data[:,1], color = 'r')
plt.plot(data[:,0], odemethod[:,0], color = 'c')
# plt.scatter(data[:,0], odemethod[:,0], color = 'c')
plt.plot(data[:,0], data[:,3], color = 'm')
plt.axvline(x=-1/k, ymin=-99, ymax=999, color='m')
plt.axvline(x=-1)
# plt.axvline(x=-1/freeze_time,ymin=-999,ymax=999,color='m')
plt.axhline(y=0,xmin=-9999,xmax=9999,color='k')
plt.xlim(-1.5,0)												# Comment in to see how it behaves outside horizon (should go v ~ a but doesnt change until the final point regardless of step size)
# plt.ylim(-1,1)
plt.show()
