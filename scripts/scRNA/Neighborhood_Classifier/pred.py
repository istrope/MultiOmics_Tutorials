import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
import scanpy as sc
import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from scipy.spatial.distance import cdist


def generate_clusters(adata, resolution=0.3):
    clusters = np.unique(adata.obs.cell_type)
    cluster = []
    ind = []
    for clus in clusters:
        ndata = adata[adata.obs.cell_type == clus]
        sc.pp.neighbors(adata,use_rep='X_scANVI')
        if clus == 'lumhr':
            sc.tl.leiden(ndata,resolution = 0.7)
        else:
            sc.tl.leiden(ndata,resolution = resolution)
        ndata.obs['cluster'] = ndata.obs.cell_type.astype(str) + '_' + ndata.obs.leiden.astype(str)
        cluster.extend(ndata.obs['cluster'].values.tolist())
        ind.extend(ndata.obs.index)

    df = pd.DataFrame({'cluster':cluster})
    df.index = ind
    adata.obs = adata.obs.merge(df,left_index=True,right_index=True,how='left')
    v = adata.obs.cluster.value_counts()
    adata = adata[adata.obs.cluster.isin(v.index[v.gt(100)])]
    return adata

def compute_convex_hulls_and_centroids(adata, embedding):
    cluster_hulls = {}
    centroids = {}
    coordinates = adata.obsm[embedding]
    for cluster in adata.obs['leiden'].unique():
        points_in_cluster = coordinates[adata.obs['leiden'] == cluster]
        hull = ConvexHull(points_in_cluster)
        cluster_hulls[cluster] = hull
        centroids[cluster] = points_in_cluster.mean(axis=0)
    return cluster_hulls, centroids

# Function to check if a point is in a convex hull
def point_in_hull(point, hull):
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull.points[hull.vertices])
    return hull.find_simplex(point) >= 0

# Compute convex hulls for each cluster in the latent space

def predict_clusters_with_hulls_and_scores(new_data, cluster_hulls, centroids, cluster_radii):

    new_data_latent = new_data.obsm['X_scANVI']
    # Initialize lists to store results
    predicted_clusters = []
    likelihood_scores = []
    
    # Predict clusters and calculate likelihood scores
    for point in new_data_latent:
        assigned_cluster = 'Unassigned'
        point_scores = {}
        for cluster, hull in cluster_hulls.items():
            if point_in_hull(point, hull):
                assigned_cluster = cluster
                centroid = centroids[cluster]
                radius = cluster_radii[cluster]
                distance = np.linalg.norm(point - centroid)
                score = 1 - (distance / radius) if distance <= radius else 0
                point_scores[cluster] = score
        predicted_clusters.append(assigned_cluster)
        likelihood_scores.append(point_scores)
    
    return predicted_clusters, pd.DataFrame(likelihood_scores)



def compute_cluster_radii(adata, centroids,embedding):
    cluster_radii = {}
    coordinates = adata.obsm[embedding]
    for cluster in adata.obs['cluster'].unique():
        points_in_cluster = coordinates[adata.obs['cluster'] == cluster]
        centroid = centroids[int(cluster)]
        distances = cdist(points_in_cluster, [centroid], 'euclidean')
        cluster_radii[cluster] = distances.max()
    return cluster_radii

#RANDOM FOREST CLASSIFIER
def random_forest(adata,embedding):
    # Train Random Forest Classifier
    X_train = adata.obsm[embedding]
    y_train = adata.obs['cluster']

    # Encode the cluster labels to ensure they are numeric
    label_encoder = LabelEncoder()
    y_train_encoded = label_encoder.fit_transform(y_train)

    random_forest = RandomForestClassifier(n_estimators=100, random_state=42)
    calibrated_rf = CalibratedClassifierCV(random_forest, method='isotonic', cv=5)
    calibrated_rf.fit(X_train,y_train_encoded)
    return calibrated_rf, label_encoder

# Function to predict clusters for new data points and get probabilities
def predict_clusters_with_proba(new_data, calibrated_rf, label_encoder, threshold=0.8):

    new_adata_latent = new_data.obsm['X_scANVI']
    # Predict probabilities using Random Forest
    predicted_proba = calibrated_rf.predict_proba(new_adata_latent)
    
    # Determine the cluster with the highest probability and check the threshold
    predicted_clusters = []
    for proba in predicted_proba:
        max_proba = np.max(proba)
        if max_proba >= threshold:
            cluster = label_encoder.inverse_transform([np.argmax(proba)])[0]
        else:
            cluster = 'Unassigned'
        predicted_clusters.append(cluster)
    return predicted_clusters, predicted_proba

def main():
    adata = sc.read('merged.h5ad')
    navin = adata[adata.obs.project == 'navin']
    fariba = adata[adata.obs.project == 'fariba']
    sc.pp.neighbors(navin,use_rep='X_scANVI')
    navin = generate_clusters(navin)

    cluster_counts = navin.obs['cluster'].value_counts()
    clusters_to_keep = cluster_counts[cluster_counts > 100].index
    navin = navin[navin.obs['cluster'].isin(clusters_to_keep)]
    #cluster_hulls,centroids = compute_convex_hulls_and_centroids(navin, embedding="X_scVI")
    #cluster_radii = compute_cluster_radii(navin, centroids, embedding="X_scVI")

    #predict_clusters, distance_scores = predict_clusters_with_hulls_and_scores(fariba,cluster_hulls,centroids,cluster_radii)

    #predict_clusters.to_csv('hull_pred_clusters.csv')
    #distance_scores.to_csv('hull_distance_scores.csv')

    rf, label_encoder = random_forest(navin,'X_scANVI')
    predicted_clusters, prob = predict_clusters_with_proba(fariba,rf,label_encoder)
    predicted_clusters.extend(navin.obs.cluster.values.tolist())
    ind = []
    ind.extend(fariba.obs.index.astype(str).tolist())
    ind.extend(navin.obs.index.astype(str).tolist())
    pred_clusters = pd.DataFrame(predicted_clusters, columns=["Predicted Cluster"])
    pred_clusters.index = ind
    navin.obs['cluster'].to_csv('navin_clusters.csv')
    pred_clusters.to_csv("predicted_clusters.csv")
    prob = pd.DataFrame(prob)
    prob.to_csv('predicted_probabilities.csv')

if __name__ == '__main__':
    main()
