
import os
import numpy as np
from math import sin, cos

pi = np.pi
nn = 64

x = np.zeros(nn)
y = np.zeros(nn)

for k in range(nn):
    x[k] = cos(2*pi*k/nn)
    y[k] = sin(2*pi*k/nn)

fid = open("example.poly", "w")

fid.write("{0} {1} {2} {3}\n".format(nn, 2, 0, 1))
for k in range(nn):
    fid.write("{0} {1} {2} 1\n".format(k+1, x[k], y[k]))

fid.write("{0} 1\n".format(nn))
for k in range(nn):
    fid.write("{0} {1} {2} 1\n".format(k+1, k+1, (k+1)%nn + 1))

fid.write("0\n")
fid.close()

area = 0.0005
os.system("triangle -pqnea" + str(area) + " example.poly")
