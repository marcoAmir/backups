#!/usr/bin/env python

# Here I will run outlier-detection approach using One-class SVM, in order to detect transcripts that deviate from normal transcripts
						
import numpy as np
import random
import sys

import os

from sklearn import svm
from sklearn import preprocessing
from sklearn.covariance import EllipticEnvelope

#import matplotlib.pyplot as plt
#import matplotlib.font_manager
#from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

if len(sys.argv)!=2:
	sys.exit('usage:\n\ttranscriptOutlierDetection.py querySpecies\n')
querySpecies = sys.argv[1]

if not os.path.exists('output'):
	os.makedirs('output')

# specify which features to include - the values indicate the column numbers to read from (0 based)
# to exclude a feature replace with a negative number
featureDict = {}
featureDict['stop']         = 4		# col4  - stop codons
featureDict['fs']           = 5		# col5  - frameshift indels
featureDict['delExons']     = 6		# col6  - deleted Exons
featureDict['totalExons']   = 7		# col7  - total Exons
featureDict['fracDelExons'] = 8		# col8  - fraction deleted exons
featureDict['delBases']     = 9		# col9  - deleted coding bases
featureDict['totalBases']   = 10	# col10 - total coding bases
featureDict['fracDelBases'] = 11	# col11 - fraction deleted coding bases
featureDict['percentID']    = -1	# col12 - percent Identity
featureDict['KaKs']         = 13	# col13 - Ka/Ks
columns = []
for key in featureDict:
	if (featureDict[key] > 0): 
		columns.append(featureDict[key])
columns.sort()

# define two outlier detection tools to be compared
#classifiers = {
#	"One-Class_SVM": svm.OneClassSVM(nu=0.02, kernel="rbf", gamma=0.1),
#	"robust_covariance_estimator": EllipticEnvelope(contamination=.1)}
classifiers = {
	"One-Class_SVM": svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.1)}

# load data files
# There are three categories: intact, others, curated with 0, 1, 2 respectively at column 16 (or 15, zero based...)
data = np.genfromtxt('/cluster/u/amirma/geneLoss/intactGenes/Classifiers/featuresBySpecies/' + querySpecies + '.allFeatures')
TxIntact  = data[data[:,15]==0][:,columns]
TxOthers  = data[data[:,15]==1][:,columns]
TxCurated = data[data[:,15]==2][:,columns]

# Split the intact transcript (TxIntact) at random to 80% & 20% training and test set, respectively:
u = range(0, len(TxIntact))
random.shuffle(u)
trainSize = int(round( len(TxIntact) * 0.80 ))
TxIntact_train = TxIntact[u[0:trainSize],:]
TxIntact_test  = TxIntact[u[trainSize:len(TxIntact)],:]

# Standartization of the column:
TxIntact_train = preprocessing.scale(TxIntact_train)
TxIntact_test  = preprocessing.scale(TxIntact_test)
TxOthers       = preprocessing.scale(TxOthers)
TxCurated      = preprocessing.scale(TxCurated)

# Now load the data for all the transcripts align to the species and write down what you predict
#allData = np.genfromtxt("Mouse.allTx.allFeatures")
#TxAll = allData[:,columns]
#TxAll = preprocessing.scale(TxAll)
TxAll  = np.empty([0, len(columns)])
geneTxs = []
f = open('/cluster/u/amirma/geneLoss/intactGenes/Classifiers/featuresBySpecies/' + querySpecies + '.allTx.allFeatures')
for line in f:
	geneTxs.append('\t'.join(line.split('\t')[:4]))
	a = np.array([float(line.split('\t')[i]) for i in columns])
	TxAll = np.append(TxAll, a.reshape(1, len(a)), axis = 0)
f.close()
TxAll = preprocessing.scale(TxAll)
	

#'\t'.join(content.split('\t')[:4])

logFile = open('output/classifiers.' + querySpecies + '.log', 'a')
# fit the model - based on the canonical transcripts of the intact genes
for i, (clf_name, clf) in enumerate(classifiers.items()):
	logFile.write('\n##########################################################################')
	print "\n\n\t...Learning with " + clf_name + "...\n"
	clf.fit(TxIntact_train)
	y_pred_TxIntact_train = clf.predict(TxIntact_train)				# should return 1 for regular observation, -1 for novelty
	y_pred_TxIntact_test  = clf.predict(TxIntact_test)
	y_pred_others         = clf.predict(TxOthers)
	y_pred_curated        = clf.predict(TxCurated)
	n_novels_train        = y_pred_TxIntact_train[y_pred_TxIntact_train == -1].size	# how many of the training fall outside the boundary
	n_novels_test         = y_pred_TxIntact_test[y_pred_TxIntact_test == -1].size	# how many of the test fall outside the boundary
	n_novels_others       = y_pred_others[y_pred_others == -1].size			# how many of the others fall outside the boundary
	n_novels_curated      = y_pred_curated[y_pred_curated == -1].size		# how many of the curated fall outside the boundary
	# Report results - for train and test in canonical transcripts
	logFile.write("\nOutelier Detection with " + clf_name + "\n\t" + str(clf) + "\n\n")
	logFile.write("%d of %d training examples (%.3g percent) fall outside learned rule\n" %(n_novels_train, len(TxIntact_train), 100*float(n_novels_train)/len(TxIntact_train)))
	logFile.write("%d of %d test cases (%.3g percent) fall outside learned rule\n" %(n_novels_test, len(TxIntact_test), 100*float(n_novels_test)/len(TxIntact_test)))
	logFile.write("%d of %d others (%.3g percent) fall outside learned rule\n" %(n_novels_others, len(TxOthers), 100*float(n_novels_others)/len(TxOthers)))
	logFile.write("%d of %d curated cases (%.3g percent) fall outside learned rule\n" %(n_novels_curated, len(TxCurated), 100*float(n_novels_curated)/len(TxCurated)))
	# now predict for all the transcripts
	y_TxAll = clf.predict(TxAll)
	n_all = y_TxAll[y_TxAll == -1].size
	logFile.write("\nAll Transcripts:\n%d of all %d transcripts (%.3g percent) fall outside learned rule\n\n" %(n_all, len(TxAll), 100*float(n_all)/len(TxAll)))
	allPreds = open('output/TxPredictions.' + querySpecies + '.' + clf_name + '.gamma_' + str(clf.gamma) + '.nu_' + str(clf.nu), 'w')
	for j, currTx in enumerate(geneTxs):
		allPreds.write(currTx + '\t' + str(int(y_TxAll[j]<0)) + '\n')
	allPreds.close()

logFile.close()



