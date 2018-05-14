#!/usr/bin/env python

import numpy as np
import random

from sklearn import svm
from sklearn import preprocessing

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

# load data files
# There are three categories: intact, others, curated with 0, 1, 2 respectively at column 16 (or 15, zero based...)
data = np.genfromtxt("Rat.allFeatures")
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


# fit the model - based on the canonical transcripts of the intact genes
gamma_val = 1.0
nu_val = 0.15
clf = svm.OneClassSVM(nu=nu_val, kernel="rbf", gamma=gamma_val)		# tune gamma for the regularization
clf.fit(TxIntact_train)
y_pred_TxIntact_train = clf.predict(TxIntact_train)				# should return 1 for regular observation, -1 for novelty
y_pred_TxIntact_test  = clf.predict(TxIntact_test)
y_pred_others         = clf.predict(TxOthers)
y_pred_curated        = clf.predict(TxCurated)
n_novels_train        = y_pred_TxIntact_train[y_pred_TxIntact_train == -1].size	# how many of the training fall outside the boundary
n_novels_test         = y_pred_TxIntact_test[y_pred_TxIntact_test == -1].size	# how many of the test fall outside the boundary
n_novels_others       = y_pred_others[y_pred_others == -1].size			# how many of the others fall outside the boundary
n_novels_curated      = y_pred_curated[y_pred_curated == -1].size		# how many of the others fall outside the boundary

# Report results
print("\nOneClassSVM with: nu=" + str(nu_val) + ", kernel=rbf, gamma=" + str(gamma_val))
print("%d of %d training examples (%.3g percent) fall outside learned rule" %(n_novels_train, len(TxIntact_train), 100*float(n_novels_train)/len(TxIntact_train)))
print("%d of %d test cases (%.3g percent) fall outside learned rule" %(n_novels_test, len(TxIntact_test), 100*float(n_novels_test)/len(TxIntact_test)))
print("%d of %d others (%.3g percent) fall outside learned rule" %(n_novels_others, len(TxOthers), 100*float(n_novels_others)/len(TxOthers)))
print("%d of %d curated cases (%.3g percent) fall outside learned rule\n" %(n_novels_curated, len(TxCurated), 100*float(n_novels_curated)/len(TxCurated)))

# Now load the data for all the transcripts align to the species and write down what you predict
allData = np.genfromtxt("Rat.allTx.allFeatures")
TxAll = allData[:,columns]
TxAll = preprocessing.scale(TxAll)
y_TxAll = clf.predict(TxAll)
n_all      = y_TxAll[y_TxAll == -1].size

out = open('predictions.txt', 'w')
for i in y_TxAll:
	out.write("%d\n" %i)
out.close()
print("All transcripts\n\t%d of %d all transcripts (%.3g percent) fall outside learned rule\n" %(n_all, len(TxAll), 100*float(n_all)/len(TxAll)))
	
# Write output for grid search:
out = open('output.tsv', 'a')
out.write("%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\t%.3g\n" %(gamma_val, nu_val, 100*float(n_novels_train)/len(TxIntact_train), 100*float(n_novels_test)/len(TxIntact_test), 100*float(n_novels_others)/len(TxOthers), 100*float(n_novels_curated)/len(TxCurated), 100*float(n_all)/len(TxAll)))
out.close()
