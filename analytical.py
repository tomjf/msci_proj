import matplotlib.pyplot as plt
import numpy

w_vals = numpy.linspace(-2, 2,1000)
analytical1, analytical2, real = [],[],[]
for w in w_vals:
	if -0.33333<w<1.0:
		real.append(4+((3*(w-1))/(1+(3*w))))
	else:
		real.append(4-((3*(w-1))/(1+(3*w))))

	# analytical1.append(4-((3*(w-1))/(1+(3*w))))
	# analytical2.append(4+((3*(w-1))/(1+(3*w))))
plt.ylim(-80,10)
plt.plot(w_vals,real)
plt.axhline(y=0.96,xmin=-9999,xmax=9999)
# plt.plot(w_vals,analytical2, color='g')
# plt.scatter(w_vals,real, color='y')
# plt.axvline(x=-1.0/3.0,ymin=-9999,ymax=9999, color='r')
# plt.axvline(x=1.0,ymin=-9999,ymax=9999, color='m')
plt.show()

