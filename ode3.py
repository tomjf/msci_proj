from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy
import math
global k
import os
base_path = os.path.dirname(os.path.abspath(__file__))
#--------------------------------------------------------------------------------------------------------------------------------
# Define Parameter Class so don't have to pass multiple variables between functions and can have preset values ie for inflation
class ParamsType(object):
	k = 1
#--------------------------------------------------------------------------------------------------------------------------------
def deriv(y,t):
	a = -((k*k)-(C/(t*t)))
	return numpy.array([ y[1], a*y[0] ])
#--------------------------------------------------------------------------------------------------------------------------------
def calc_gradient():
	coefficients = numpy.polyfit(data[:,0], data[:,3], 1)
	polynomial = numpy.poly1d(coefficients)
	return polynomial
#--------------------------------------------------------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------------------------------------------------------
def plotgraphs(data, odemethod):
	# rk4 = plt.plot(data[:,0], data[:,1], color = 'r')
	# plt.scatter(data[:,0], data[:,1], color = 'r')
	ode_scipy = plt.plot(data[:,0], odemethod[:,0], color = 'c')
	# plt.scatter(data[:,0], odemethod[:,0], color = 'c')
	# analytical = plt.plot(data[:,0], data[:,3], color = 'm')
	plt.axvline(x=-1/k, ymin=-99, ymax=999, color='m')
	plt.axvline(x=-1)
	# plt.axvline(x=-1/freeze_time,ymin=-999,ymax=999,color='m')
	plt.axhline(y=0,xmin=-9999,xmax=9999,color='k')
	plt.xlim(-0.1,0)												# Comment in to see how it behaves outside horizon (should go v ~ a but doesnt change until the final point regardless of step size)
	plt.ylim(-1,2)
	# plt.legend([rk4[0], ode_scipy[0], analytical[0]], ["RK4 Method", "Odeint RK4 Method", "Analytical Solution"])
	plt.show()
#--------------------------------------------------------------------------------------------------------------------------------
def kspectrum(kspace, n_a, n_b, nt):
	for i in range(0,len(kspace)):
		k = kspace[i]
		times = numpy.linspace(n_a,n_b,nt) 
		# odeint method for comparison with RK4 and analytical methods
		xinit = numpy.array([1/(math.sqrt(2*k))*math.cos(k*n_inside), -k/(math.sqrt(2*k))*math.sin(k*n_inside)])
		odemethod = odeint(deriv,xinit,times) 
		data = RK4(k,n_inside,n_outside,num_steps)
		plotgraphs(data, odemethod)
		
		# data[i,3] = math.log(k*k*k*data[i,2]*data[i,2]*time[999]*time[999])
#--------------------------------------------------------------------------------------------------------------------------------
params = ParamsType()

w_list, ns_list = [], []

w_llim, w_ulim, w_steps = -0.05,0.05,10
wspace = numpy.linspace(w_llim,w_ulim,w_steps)
for w_val in wspace:
	w = w_val
	print w
	if -0.34 < w < -0.32:
		print '@@@@@@@@@@@@@@@@@ in the condition'
		w_list.append(w)
		ns_list.append(False)
	else:
		kspace = numpy.linspace(1,100,10)
		n_inside = -100/(numpy.amin(kspace))
		n_outside = -1/(1000*numpy.amax(kspace))
		if w > 1/3:
			n_inside, n_outside = -n_outside, -n_inside
		num_steps = 100000
		times = numpy.linspace(n_inside,n_outside,num_steps)

		# Use w to work out: {numerator in M-S eqn, a dependence on eta}
		# Also (in general) the power eta (conformal time) is raised to in the equation for a is calculated (a_power)
		a_power = 2/((3*w)+1)
		C = a_power*(a_power-1)

		''' Comment this bit back in to iterate for just one k value and see the point where the solution blows up outside the horizon
		k=1
		xinit = numpy.array([1/(math.sqrt(2*k))*math.cos(k*n_inside), -k/(math.sqrt(2*k))*math.sin(k*n_inside)])
		odemethod = odeint(deriv,xinit,times) 
		data = RK4(k,n_inside,n_outside,num_steps)
		plotgraphs(data, odemethod) '''

		data = numpy.zeros((len(kspace),4))

		for i in range(0,len(kspace)):
			k = kspace[i]
			time = numpy.linspace(n_inside,n_outside,num_steps) 
			# time = numpy.logspace(n_a,n_b,num=100) 
			# print time
			xinit = numpy.array([1/(math.sqrt(2*k))*math.cos(k*n_inside), -k/(math.sqrt(2*k))*math.sin(k*n_inside)])
			x = odeint(deriv,xinit,time) 
			yinit = numpy.array([-1/(math.sqrt(2*k))*math.sin(k*n_inside), -k/(math.sqrt(2*k))*math.cos(k*n_inside)])
			y = odeint(deriv,yinit,time)
			x_len = int(len(x))-3
			data[i,0] = math.log(k)
			data[i,1] = math.sqrt(x[5,0]*x[5,0] + y[5,0]*y[5,0]) 
			data[i,2] = math.sqrt(x[x_len,0]*x[x_len,0] + y[x_len,0]*y[x_len,0]) 
			data[i,3] = math.log(k*k*k*data[i,2]*data[i,2]*math.pow(abs(time[x_len]),a_power)*math.pow(abs(time[x_len]),a_power))		

		#--------------------------------------------------------------------------------------------------------------------------------
		coefficients = numpy.polyfit(data[:,0], data[:,3], 1)
		polynomial = numpy.poly1d(coefficients)
		ys = polynomial(data[:,0])
		# print polynomial
		w_list.append(w)
		ns_list.append(polynomial[1])

# plt.scatter(w_list, ns_list)
# plt.show()

# add 1 to each element to convert from ns-1 to ns
ns_list = [x+1 for x in ns_list]

# find ns values of interest closte to 0.96
ns_scale_invariant, w_scale_invariant = [], []
for a,b in enumerate(ns_list):
	if 0.9 < b < 1.06:
		ns_scale_invariant.append(b)
		w_scale_invariant.append(w_list[a])

# Plot w vs ns for the region of w being studied
plt.plot(w_list, ns_list)
plt.xlabel('w')
plt.ylabel('ns')
plt.axhline(y=0.96,xmin=-9999,xmax=9999)
if len(ns_scale_invariant) > 0:
	plt.scatter(w_scale_invariant, ns_scale_invariant, color = 'r')
plt.savefig('~Results' '{' + str(w_llim) + ',' + str(w_ulim) + '}' + ',' + str(w_steps) + '_' + 'steps' + '_' + 'graph.png')
plt.show()

# save data to text file so don't have to run again
text_data = numpy.vstack((w_list,ns_list))
text_data = numpy.reshape(text_data,[len(w_list),2])
numpy.savetxt('~Results' '{' + str(w_llim) + ',' + str(w_ulim) + '}' + ',' + str(w_steps) + '_' + 'steps' + '_' + 'data.txt', text_data)	

	# f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
	# ax1.set_title('inside horizon')
	# ax1.plot(data[:,0],data[:,1])
	# ax2.set_title('outside horizon')
	# ax2.plot(data[:,0],data[:,2], color = 'r')
	# ax3.set_title('outside horizon powe spec')
	# ax3.plot(data[:,0],data[:,3], color = 'g')
	# ax3.plot(data[:,0],ys)
	# # plt.xlim(1,2.5)
	# plt.show()


# a_pow = calc_a_dependence(1/3)
# print a_pow

# test_w_vs_a(-10,10)