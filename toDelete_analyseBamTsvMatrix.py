import numpy as np
import umap
from sklearn.cluster import DBSCAN
from collections import Counter

def run(matrix, sample,region):
    results=[]
    variableColumns=[]
    for i in range(matrix.shape[1]):#iterate over columns which are A/C/T/G positins
        if len(np.unique(matrix[:,i]))>1:
            variableColumns.append(i)
    #variableColumns=np.asarray(variableColumns)
    for dist in range(5,30,5):
        print(dist)
        X_embedded = umap.UMAP(n_neighbors=dist, min_dist=0.01, metric='hamming').fit_transform(matrix[:,variableColumns])

        for eps in range(6,20,2):
            clustering = DBSCAN(eps=eps/10, min_samples=5).fit(X_embedded).labels_
            counts=Counter(clustering)
            for item in counts:
                if counts[item]/len(clustering)>0.05:
                    results.append([sample,region,counts[item]/len(clustering)])
    return results