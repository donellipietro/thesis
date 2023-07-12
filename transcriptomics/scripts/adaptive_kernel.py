import numpy as np
from scipy.linalg import sqrtm

def adaptive_kernel(X, V, x_ref, v_ref, N_kn_loc, N_kn_val, max_iter, toll=1e-6):
    # Initialization
    iter = 0
    convergence = False
    cov = np.identity(2)  # identity matrix

    # Difference between each spot and the reference one
    d = X - x_ref

    # L2 distance for the gene expression
    dist_val = np.sqrt(np.sum((V - v_ref[:,np.newaxis]) ** 2, axis=0))

    # Iterative algorithm
    while iter < max_iter and not convergence:
        # Inverse of the covariance matrix
        inv_cov = np.linalg.inv(cov)

        # Distance computation
        dist_loc = np.sum(d @ inv_cov * d, axis=1)

        # Identify the closest N_kn_loc points to the reference spot in terms of location
        order_loc = np.argsort(dist_loc)
        closest_points_loc = X[order_loc[:N_kn_loc]]  # location of the closest points
        closest_points_val = dist_val[order_loc[:N_kn_loc]]  # (distance) value of the closest points

        # Among the spots just selected using the spatial criteria, identify the
        # closest N_kn_val spots to the reference spot in terms of values
        order_values = np.argsort(closest_points_val)
        interesting_points_loc = closest_points_loc[order_values[:N_kn_val]]

        # Save the previous covariance matrix
        cov_prev = cov

        # Compute the covariance matrix of the cloud of interesting points
        cov = np.cov(interesting_points_loc.T)

        # Convergence check
        if np.linalg.norm(cov - cov_prev) < toll:
            convergence = True

        iter += 1

    if not convergence:
        cov = np.identity(2)

    # Kernel computation
    inv_cov = np.linalg.inv(cov)
    dist_loc = np.sum(d @ inv_cov * d, axis=1)
    ker = np.exp(-0.5 * dist_loc)

    return ker

# Data
X = np.genfromtxt("../data/DLPFC/DLPFC_xy_coords.csv", delimiter=",", skip_header=True)
V = np.genfromtxt("../data/DLPFC/DLPFC_gene_by_spot_mat_normalized_zero_mean_1stddev.csv", delimiter=",", skip_header=True)[:, 1:]
Clustering = np.genfromtxt("../data/DLPFC/DLPFC_Ground_Truth_Labels_Integers.csv", delimiter=",")

# Algorithm test
ref_point_index = np.random.randint(0, X.shape[0])
x_ref = X[ref_point_index]
v_ref = V[:,ref_point_index]

# Hyperparameters
N_kn_loc = 25
N_kn_val = 10

# Stopping criteria
toll = 1e-6
max_iter = 10

# Execution
result = adaptive_kernel(X, V, x_ref, v_ref, N_kn_loc, N_kn_val, max_iter, toll)