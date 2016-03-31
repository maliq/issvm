from sklearn.datasets import load_svmlight_file
import pylab as plt
import sys
import numpy as np 

def sign(x):
	if x >= 0:
		return 1
	return -1

def hinge(y1,y2):
	if (1.0-(y1*y2)) < 0.0:
		return 0.0
	return 1.0-(y1*y2)

def comp_loss(pred_name,true_name):
	predictions = plt.loadtxt(pred_name)
	data_real = load_svmlight_file(true_name)

	text = sys.argv[4]
	output_f = open(sys.argv[5],'a')

	y = data_real[1]

	mistakes = 0
	hinge_loss = 0.0


	for idx in range(1,len(y)):
		if(sign(predictions[idx])!=y[idx]):
			mistakes+=1
		hinge_loss = hinge_loss + hinge(predictions[idx],y[idx])

	hinge_loss = hinge_loss/float(len(y))
	miss_loss = float(mistakes)/float(len(y))
	return miss_loss

miss_train = comp_loss(sys.argv[1],sys.argv[2])
miss_test = comp_loss(sys.argv[3],sys.argv[4])
model_data = plt.loadtxt(sys.argv[5],dtype={'names': ('bias', 'weights', 'example'),
          'formats': (np.float, np.float, '|S1024')}, skiprows=0) #delimiter=','
supp_size = len(model_data)

text = sys.argv[6]
output_f = open(sys.argv[7],'a')

output_f.write(text+" ")
output_f.write(str(miss_test)+" ")
output_f.write(str(miss_train)+" ")
output_f.write(str(supp_size)+"\n")

