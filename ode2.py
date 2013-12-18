from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy
import math
global k
#-----------------------------------------------
def deriv(y,t):
	print t
	# print C, '@@@@.....@@@@@@'
	a = -((k*k)-(C/(t*t)))
	b = 0
	return numpy.array([ y[1], a*y[0]+b*y[1] ])
#-----------------------------------------------
w = raw_input("w = ... ")
w = float(w)
coeff = 2/((3*w)+1)
C = coeff*(coeff-1)

kspace = numpy.linspace(1,10,100)
n_a = -1000/(numpy.amin(kspace))
n_b = -1/(1000*numpy.amax(kspace))

if w > 1/3:
	n_a, n_b = -n_b, -n_a


data = numpy.zeros((len(kspace),4))
nt = 100000
#------------------------------------------------
for i in range(0,len(kspace)):
	# print i
	k = kspace[i]
	time = numpy.linspace(n_a,n_b,nt) 
	# time = numpy.logspace(n_a,n_b,num=100) 
	# print time
	xinit = numpy.array([1/(math.sqrt(2*k))*math.cos(k*n_a), -k/(math.sqrt(2*k))*math.sin(k*n_a)])
	x = odeint(deriv,xinit,time) 
	yinit = numpy.array([-1/(math.sqrt(2*k))*math.sin(k*n_a), -k/(math.sqrt(2*k))*math.cos(k*n_a)])
	y = odeint(deriv,yinit,time)
	print x
	data[i,0] = math.log(k)
	data[i,1] = math.sqrt(x[5,0]*x[5,0] + y[5,0]*y[5,0]) 
	data[i,2] = math.sqrt(x[999,0]*x[999,0] + y[999,0]*y[999,0]) 
	

	plt.plot(time, x[:,1], color = 'r')
	plt.plot(time, x[:,0], color = 'c')
	plt.axvline(x=-1/k, ymin=-99, ymax=999, color='m')
	# plt.axvline(x=-1/freeze_time,ymin=-999,ymax=999,color='m')
	plt.axhline(y=0,xmin=-9999,xmax=9999,color='k')
	# plt.xlim(-1,0)
	plt.ylim(-1,1)
	plt.show()
	
	data[i,3] = math.log(k*k*k*data[i,2]*data[i,2]*time[999]*time[999])
#-----------------------------------------------------

coefficients = numpy.polyfit(data[:,0], data[:,3], 1)
polynomial = numpy.poly1d(coefficients)
ys = polynomial(data[:,0])
print polynomial

f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
ax1.set_title('inside horizon')
ax1.plot(data[:,0],data[:,1])
ax2.set_title('outside horizon')
ax2.plot(data[:,0],data[:,2], color = 'r')
ax3.set_title('outside horizon powe spec')
ax3.plot(data[:,0],data[:,3], color = 'g')
ax3.plot(data[:,0],ys)
# plt.xlim(1,2.5)
plt.show()