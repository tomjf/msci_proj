import matplotlib.pyplot as plt
import numpy

x = numpy.logspace(0.01,3,num=5,base=10)
y = [1,2,3,4,5]
# y=[10,20,30,40,50,60,70,80,90]

a = 1.0
b = 'sn'

plt.scatter(x,y)
plt.title("Jacobi Elliptic" + "__"+ b + "__" + "m = " + str(a))
plt.xlabel('Amplitude')
plt.ylabel('ns')
plt.axhline(y=0.96,color='r')
plt.show()
