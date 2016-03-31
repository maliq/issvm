from sklearn.datasets import load_svmlight_file
import pylab as plt
import sys
import numpy as np

model_data = plt.loadtxt(sys.argv[1],dtype={'names': ('bias', 'weights', 'example'),
          'formats': (np.float, np.float, '|S1024')}, skiprows=0) #delimiter=','

print len(model_data)