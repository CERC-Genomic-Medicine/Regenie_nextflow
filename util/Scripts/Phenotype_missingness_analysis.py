#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
import seaborn as sns
from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Clustering Phenotype based on missingness similarity with a threshold limit')
parser.add_argument('-P','--phenotype', metavar='FILE', dest='Phenotype', required=True, type=str, help='Phenotype file.')
parser.add_argument('-T', '--threshold', dest='threshold', type=float,required=False, default=0.15, help='Maximum jaccard distance (i.e. (present-common)/(all-present) within cluster')
parser.add_argument('-m', '--maxColumn', dest='max_col', type=int,required=False, default=50, help='Maximum number of phenotype per cluster')
parser.add_argument('-o', '--output', dest='output', type=str,required=False, default='output', help='output name for the figure generated')


def calculate_jaccard_distance(df):
    # Create a boolean matrix where True indicates presence (non-missing values)
    bool_matrix = df.notna().astype(int)
# Compute the pairwise Jaccard distances
    jaccard_distances = pdist(bool_matrix.T, metric='jaccard')
    fig = plt.figure(figsize=(25, 25))
    dist = squareform(jaccard_distances)
    sns.heatmap(dist, cmap="mako")
    plt.savefig(f'heatmap_{args.output}.png')
    return jaccard_distances

def calculate_unshared_nas(df, cluster):
    # Calculate the percentage of rows with unshared NaNs in a cluster
    bool_matrix = df[cluster].dropna(how='all').isna()
    shared_absence = bool_matrix.sum(axis=0).max()
    #print(bool_matrix.iloc[1,:])
    return shared_absence / len(bool_matrix.iloc[:,0]) * 100

def hierarchical_clustering(df, max_unshared, max_cols):
    # Calculate the Jaccard distance between columns
    jaccard_distances = calculate_jaccard_distance(df)
    #print(jaccard_distances)
    # Perform hierarchical clustering using the 'complete' linkage method
    linkage_matrix = linkage(jaccard_distances, method='complete')
    fig = plt.figure(figsize=(25, 10))
    dn = dendrogram(linkage_matrix)
    plt.savefig(f'Dendogram_{args.output}.png')
    # Calculate a threshold for clustering based on maximum unshared missingness
    threshold = max_unshared
    
    # Form flat clusters from the hierarchical clustering defined by the linkage matrix
    clusters = fcluster(linkage_matrix, t=threshold, criterion='distance')
    
    # Group columns based on their cluster assignments
    cluster_dict = {}
    for idx, cluster_id in enumerate(clusters):
        if cluster_id not in cluster_dict:
            cluster_dict[cluster_id] = []
        cluster_dict[cluster_id].append(df.columns[idx])
    
    # Limit the number of columns per cluster to the max_cols
    final_clusters = []
    unshared_nas_stats = []
    for cols in cluster_dict.values():
        for i in range(0, len(cols), max_cols):
            cluster = cols[i:i + max_cols]
            final_clusters.append(cluster)
            unshared_nas_stats.append(calculate_unshared_nas(df, cluster))
    
    return final_clusters, unshared_nas_stats
if __name__ == '__main__':
    args = parser.parse_args()
    df = pd.read_csv(args.Phenotype, sep='\t', header=0, index_col='IID')
    df.pop('FID')
    clusters, unshared_nas_stats = hierarchical_clustering(df, max_unshared=args.threshold, max_cols=args.max_col)
    #Output clusters and their maximum unshared NaNs percentage
    for i, (cluster, unshared_na) in enumerate(zip(clusters, unshared_nas_stats)):
        clustered=','.join(cluster)# usage for phenocol
        print(f"{clustered}")
