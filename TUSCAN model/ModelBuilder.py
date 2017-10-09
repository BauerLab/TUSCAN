#!/usr/bin/python

#USAGE: ./this.program.py -c/r matrix.txt features.txt
#choose if doing classification or regression
#Given a feature table, this program will build a model and store it in a joblib file

import sys
import FeatureMatrix
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.externals import joblib
from numpy import genfromtxt, savetxt

#get important features for feature selection
l = []
f = open(sys.argv[3])
for line in f:
	feat = line.split()
	for b in feat:
		l.append(b.strip('"'))
f.close()

#grabs full list of features from matrix
lines = []
with open(sys.argv[2], 'r') as f:
	features = f.readline()
features = features.split()
num_features = len(features)
features = features[1:]

#gets index of important features
a = [features.index(i) for i in l]

data = genfromtxt(sys.argv[2], dtype = 'f8', skip_header = 1, usecols = range(1,num_features))
if sys.argv[1] == '-r':
	target = data[:, features.index('Activity')]
elif sys.argv[1] == '-c':
	target = data[:, features.index('Class')]
train = data[:, a]

#Create and train the randomForest
print('Building Random Forest')
if sys.argv[1] == '-r':
	rf = RandomForestRegressor(n_estimators=1001, max_features=2)
	print('Fitting Random Forest Regressor')
	rf.fit(train, target)
	with open('rfModelregressor.joblib', 'wb') as f:
		joblib.dump(rf, f)
elif sys.argv[1] == '-c':
	rf = RandomForestClassifier(n_estimators=1001, max_features=2)
	print('Fitting Random Forest Classifier')
	rf.fit(train, target)
	with open('rfModelclassifier.joblib', 'wb') as f:
		joblib.dump(rf, f)
else:
	print('Please specify -r or -c for regressor of classifier (first argument)')
