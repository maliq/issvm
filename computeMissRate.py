from sklearn.datasets import load_svmlight_file
import numpy as np
import sys

def sign(x):
	if x >= 0:
		return 1
	return -1

def hinge(y1,y2):
	if (1.0-(y1*y2)) < 0.0:
		return 0.0
	return 1.0-(y1*y2)

predictions = np.loadtxt(sys.argv[1])
data_real = load_svmlight_file(sys.argv[2])
y = data_real[1]

mistakes = 0
hinge_loss = 0.0


for idx in range(1,len(y)):
	if(sign(predictions[idx])!=y[idx]):
		mistakes+=1
	hinge_loss = hinge_loss + hinge(predictions[idx],y[idx])

hinge_loss = hinge_loss/float(len(y))
miss_loss = float(mistakes)/float(len(y))
print "Hinge Loss=%f"%hinge_loss
print "Miss Loss=%f"%miss_loss
