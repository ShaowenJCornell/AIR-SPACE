# Chase Holdener (ch2228@cornell.edu)
# March 2024

# IN-PLACE AND GRID-BASED GAUSSIAN SMOOTHING FUNCTIONS (NON-PARALLELIZED)

import pandas as pd
import numpy as np
import math
import anndata as ad
import time

from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool, cpu_count


## FUNCTIONS FOR CONSTRUCTING GRID
# function to make square or hexagonal grid
def make_grid2(max_X, min_X, max_Y, min_Y, grid_size, is_hex_grid=False):
    x_vals = np.arange(min_X, max_X + grid_size, grid_size)
    y_vals = np.arange(min_Y, max_Y + grid_size, grid_size)

    if is_hex_grid:
        grid_points = []
        # Shift every other row horizontally by grid_size * 0.5
        for i, row in enumerate(y_vals):
            if i % 2 != 0:
                # shifted row
                grid_points.extend([(x + grid_size * 0.5, y) for x, y in zip(x_vals[:-1], [row] * (len(x_vals) - 1))])
            else:
                # non-shifted row
                grid_points.extend([(x, y) for x, y in zip(x_vals[:-1], [row] * (len(x_vals) - 1))])
            
    else:
        grid_points = [(x, y) for x in x_vals for y in y_vals]

    return grid_points

# impose grid onto anndata structure and filter for grid points within target distance to anndata point
def fit_grid_to_adata(anndata, grid, dist):
    fit_grid = []
    for i in range(len(grid)):
        x_center, y_center = grid[i] 
        grid_point_neighbors = anndata.obsm['spatial'][(anndata.obsm['spatial'][:,0] <= x_center + dist) &
                                       (anndata.obsm['spatial'][:,0] >= x_center - dist) &
                                       (anndata.obsm['spatial'][:,1] <= y_center + dist) &
                                       (anndata.obsm['spatial'][:,1] >= y_center - dist)].copy()
        if len(grid_point_neighbors) != 0:
            fit_grid.append(grid[i])
    
    return np.array(fit_grid)


# VARIOUS GAUSSIAN KERNELS
# for a given point, (x,y), finds its height in a type of gaussian kernel centered at (0,0) with std. dev. s 
def gaussian_kernel_height(x, y, s):
    # Calculate the height of the Gaussian kernel at point (x, y)
    height = math.exp(-((x**2 + y**2) / (2 * s**2)))
    return height

def x_partial_gaussian_kernel_height(x, y, s):
    # Calculate the height of the X partial Gaussian kernel at point (x, y)
    height = ((-x)/s**2) * math.exp(-((x**2 + y**2) / (2 * s**2)))
    return height

def y_partial_gaussian_kernel_height(x, y, s):
    # Calculate the height of the Y partial Gaussian kernel at point (x, y)
    height = ((-y)/s**2) * math.exp(-((x**2 + y**2) / (2 * s**2)))
    return height


## FUNCTIONS FOR APPLYING A SMOOTHING KERNEL ACROSS THE DATA
# apply some gaussian kernel to each point in fit_grid, using gene expr. data from adata
def apply_gauss_kernel_grid(adata, fit_grid, gauss_kernel_function, gaussian_sd, kernel_radius):
    
    # preallocate the .X array with len(fit_grid) empty rows to fill in with smoothed points
    X_sm = np.empty((len(fit_grid), adata.n_vars))

    # Extract x and y columns from adata.X
    x_values = adata.obs['x'].values
    y_values = adata.obs['y'].values
    
    # for each point in the grid
    for i in range(len(fit_grid)):
        
        # Create boolean masks for extracting neighboring beads
        x_mask = (x_values >= (fit_grid[i, 0] - kernel_radius)) & (x_values <= (fit_grid[i, 0] + kernel_radius))
        y_mask = (y_values >= (fit_grid[i, 1] - kernel_radius)) & (y_values <= (fit_grid[i, 1] + kernel_radius))

        # Apply boolean masks to filter rows in adata.X
        temp_X = adata.X[x_mask & y_mask]

        # Extract corresponding x, y coordinates
        temp_x_values = x_values[x_mask & y_mask]
        temp_y_values = y_values[x_mask & y_mask]

        # define kernel
        kernel = np.array([gauss_kernel_function(temp_x_values[j] - fit_grid[i,0], 
                                       temp_y_values[j] - fit_grid[i,1],
                                       gaussian_sd) for j in range(len(temp_y_values))])

        
        # ONLY FOR gauss_kernel_function = x_partial or y_partial: divide positive side of kernel by sum of positive
        #  kernel weights and divide negative side of kernel by abs(sum of negative kernel weights) to normalize each side 
        #  of the kernel to have the same volume
        if gauss_kernel_function in (x_partial_gaussian_kernel_height, y_partial_gaussian_kernel_height):
            
            # Extract positive and negative values
            pos_values = kernel[kernel > 0]
            neg_values = kernel[kernel < 0]

            # If bead count on positive or negative half of kernel is less than threshold,
            # set partial deriv value to 0 for all genes as kernel is on very edge of tissue!
            PARTIAL_BEAD_THRESHOLD = 50
            if (len(pos_values) < PARTIAL_BEAD_THRESHOLD) or (len(neg_values) < PARTIAL_BEAD_THRESHOLD):
                X_sm[i:i+1, :] = np.zeros(adata.n_vars)
                continue
            
            # Calculate sums
            pos_sum = np.sum(pos_values)
            neg_sum = np.sum(neg_values)
            
            # Normalize values by their sum
            pos_values /= pos_sum
            neg_values /= np.abs(neg_sum)
            
            # Update the kernel with normalized values
            kernel[kernel > 0] = pos_values
            kernel[kernel < 0] = neg_values

        
        # 1) reshape and transpose gauss_kernel array and temp_adata.X
        # 2) apply dot prod of 'gauss_kernel' to temp_adata to get smoothed genes for bead of interest [1 x n_genes] matrix!
        temp_boi_sm = (temp_X.T).dot(kernel.reshape(-1, 1)).T
    
        # ONLY FOR gauss_kernel_function = gaussian_kernel_height: you have to divide by the sum of gauss_kernel weights 
        #  to normalize volume under gaussian kernel (differing #s of points in each kernel for bead i)
        if (gauss_kernel_function == gaussian_kernel_height):  

            # If bead count under kernel is less than threshold, set smoothed value to
            # 0 for all genes as kernel center is on very edge of tissue!
            GAUSS_BEAD_THRESHOLD = 50
            if (len(kernel) < GAUSS_BEAD_THRESHOLD):
                X_sm[i:i+1, :] = np.zeros(adata.n_vars)
                continue

            temp_boi_sm = temp_boi_sm / np.sum(kernel) 

        # add smoothed bead to .X array we want to construct
        X_sm[i:i+1, :] = temp_boi_sm

    return X_sm


# apply some gaussian kernel to each point in original adata
def apply_gauss_kernel_inplace(adata, gauss_kernel_function, gaussian_sd, kernel_radius):
    
    # preallocate the .X array with len(fit_grid) empty rows to fill in with smoothed points
    X_sm = np.empty((len(adata), adata.n_vars))

    # Extract x and y columns from adata.X
    x_values = adata.obs['x'].values
    y_values = adata.obs['y'].values
    
    # for each point in the grid
    for i in range(len(adata)):

        # Create boolean masks for extracting neighboring beads
        x_mask = (x_values >= (adata.obsm['spatial'][i, 0] - kernel_radius)) & (x_values <= (adata.obsm['spatial'][i, 0] + kernel_radius))
        y_mask = (y_values >= (adata.obsm['spatial'][i, 1] - kernel_radius)) & (y_values <= (adata.obsm['spatial'][i, 1] + kernel_radius))

        # Apply boolean masks to filter rows in adata.X
        temp_X = adata.X[x_mask & y_mask]

        # Extract corresponding x, y coordinates
        temp_x_values = x_values[x_mask & y_mask]
        temp_y_values = y_values[x_mask & y_mask]

        # define kernel
        kernel = np.array([gauss_kernel_function(temp_x_values[j] - adata.obsm['spatial'][i,0], 
                                       temp_y_values[j] - adata.obsm['spatial'][i,1],
                                       gaussian_sd) for j in range(len(temp_y_values))])

        # ONLY FOR gauss_kernel_function = x_partial or y_partial: divide positive side of kernel by sum of positive
        #  kernel weights and divide negative side of kernel by abs(sum of negative kernel weights) to normalize each side 
        #  of the kernel to have the same volume
        if gauss_kernel_function in (x_partial_gaussian_kernel_height, y_partial_gaussian_kernel_height):
            
            # Extract positive and negative values
            pos_values = kernel[kernel > 0]
            neg_values = kernel[kernel < 0]

            # If bead count on positive or negative half of kernel is less than threshold,
            # set partial deriv value to 0 for all genes as kernel is on very edge of tissue!
            PARTIAL_BEAD_THRESHOLD = 50
            if (len(pos_values) < PARTIAL_BEAD_THRESHOLD) or (len(neg_values) < PARTIAL_BEAD_THRESHOLD):
                X_sm[i:i+1, :] = np.zeros(adata.n_vars)
                continue
            
            # Calculate sums
            pos_sum = np.sum(pos_values)
            neg_sum = np.sum(neg_values)
            
            # Normalize values by their sum
            pos_values /= pos_sum
            neg_values /= np.abs(neg_sum)
            
            # Update the kernel with normalized values
            kernel[kernel > 0] = pos_values
            kernel[kernel < 0] = neg_values

        
        # 1) reshape and transpose gauss_kernel array and temp_adata.X
        # 2) apply dot prod of 'gauss_kernel' to temp_adata to get smoothed genes for bead of interest [1 x n_genes] matrix!
        temp_boi_sm = (temp_X.T).dot(kernel.reshape(-1, 1)).T
    
        # ONLY FOR gauss_kernel_function = gaussian_kernel_height: you have to divide by the sum of gauss_kernel weights 
        #  to normalize volume under gaussian kernel (differing #s of points in each kernel for bead i)
        if (gauss_kernel_function == gaussian_kernel_height):  

            # If bead count under kernel is less than threshold, set smoothed value to
            # 0 for all genes as kernel center is on very edge of tissue!
            GAUSS_BEAD_THRESHOLD = 50
            if (len(kernel) < GAUSS_BEAD_THRESHOLD):
                X_sm[i:i+1, :] = np.zeros(adata.n_vars)
                continue
            
            temp_boi_sm = temp_boi_sm / np.sum(kernel) 
        
        # add smoothed bead to .X array we want to construct
        X_sm[i:i+1, :] = temp_boi_sm

    return X_sm


## CONVERT LAYERS OF SMOOTHED SURFACES TO AN ANNDATA STRUCTURE
# converts X_sm, X_gradient_mag, and X_gradient_dir arrays to an anndata structure, using ref_adata and fit_adata_grid to get info
def X_smoothed_w_gradient_to_anndata_grid(X_sm, X_gradient_mag, X_gradient_dir, ref_adata, fit_adata_grid):
    
    # construct anndata for grid_smoothed
    grid_sm_grad_adata = ad.AnnData(X_sm)

    # add layers for X_smoothed, X_gradient_mag, X_gradient_dir
    grid_sm_grad_adata.layers['X_grad_mag'] = X_gradient_mag
    grid_sm_grad_adata.layers['X_grad_dir'] = X_gradient_dir

    # add spatial coordinates from fit_adata_grid
    grid_sm_grad_adata.obsm['spatial'] = fit_adata_grid
    grid_sm_grad_adata.obs['x'] = fit_adata_grid[:,0]
    grid_sm_grad_adata.obs['y'] = fit_adata_grid[:,1]

    # add all .var fields for genes from ref_adata
    grid_sm_grad_adata.var = ref_adata.var.copy()

    ## CHOOSE/EDIT: here you should add any custom anndata.obs fields you want on your grid data points
    # re-add reovirus stats to grid points
    grid_sm_grad_adata.obs['total_counts'] = np.sum(grid_sm_grad_adata.X, axis=1)

    return grid_sm_grad_adata


# converts the smoothed_X output array to an anndata structure, using ref_adata and fit_adata_grid to get info
def X_smoothed_w_gradient_to_anndata_inplace(X_sm, X_gradient_mag, X_gradient_dir, ref_adata):
    
    # add layers for X_smoothed, X_gradient_mag, X_gradient_dir to original adata
    ref_adata.layers['X_sm'] = X_sm
    ref_adata.layers['X_grad_mag'] = X_gradient_mag
    ref_adata.layers['X_grad_dir'] = X_gradient_dir

    return ref_adata

## RUNS ALL CODE TOGETHER TO COMPLETE SMOOTHING 
# This runs all the functions above to correctly output an AnnData structure of smoothed gene expression
def run_grid_based_smoothing_and_gradient_analysis(adata, IS_HEX_GRID, GRID_PTS_ACROSS, GRID_SIZE, GRID_PT_MIN_DIST, GAUSSIAN_SD, KERNEL_RADIUS):
    
    ## Make grid of points and apply the grid to adata coordinates
    
    # find edges of grid
    max_X, min_X = adata.obs['x'].max(), adata.obs['x'].min()
    max_Y, min_Y = adata.obs['y'].max(), adata.obs['y'].min()
    
    # report ranges of grid
    print("Variation in X-direction: " + str(max_X - min_X))
    print("Variation in Y-direction: " + str(max_Y - min_Y))

    if GRID_SIZE > 0:
        grid_size = GRID_SIZE
    else:
        # Generate grid_size using NUM_PTS_ACROSS
        grid_size = (max_X - min_X) / GRID_PTS_ACROSS

    # Determine DIST_THRESHOLD
    dist_threshold = GRID_PT_MIN_DIST
    
    # generate grid
    adata_grid = make_grid2(max_X, min_X, max_Y, min_Y, grid_size, IS_HEX_GRID)
    
    # fit grid
    fit_adata_grid = fit_grid_to_adata(adata, adata_grid, dist_threshold)

    # report number of grid points
    print("Number of points in grid: " + str(len(fit_adata_grid)))

    
    ## Run grid-based guassian and partial gaussian smoothing (takes 7 minutes with ~24000 grid points)
    
    # run apply_gauss_kernel_grid function to get smoothed gaussian
    X_grid_sm = apply_gauss_kernel_grid(adata, 
                                  fit_adata_grid, 
                                  gaussian_kernel_height, 
                                  GAUSSIAN_SD, 
                                  KERNEL_RADIUS)
    
    # run apply_gauss_kernel_grid function to get x partial gaussian derivative
    X_xpartial_grid_sm = apply_gauss_kernel_grid(adata, 
                                           fit_adata_grid, 
                                           x_partial_gaussian_kernel_height, 
                                           GAUSSIAN_SD, 
                                           KERNEL_RADIUS)
    
    # run apply_gauss_kernel_grid function to get y partial gaussian derivative
    X_ypartial_grid_sm = apply_gauss_kernel_grid(adata, 
                                           fit_adata_grid, 
                                           y_partial_gaussian_kernel_height, 
                                           GAUSSIAN_SD, 
                                           KERNEL_RADIUS)
    
    # find the magnitude of the gradient at each point in the grid
    X_gradient_magnitude_grid = np.sqrt(X_xpartial_grid_sm**2 + X_ypartial_grid_sm**2)
    
    # find the direction of the gradient at each point in the grid
    X_gradient_direction_grid = np.arctan2(X_ypartial_grid_sm, X_xpartial_grid_sm)
    
    # convert several X arrays to new anndata structure with several layers
    grid_sm_grad_adata = X_smoothed_w_gradient_to_anndata_grid(X_grid_sm, 
                                                          X_gradient_magnitude_grid, 
                                                          X_gradient_direction_grid, 
                                                          adata, 
                                                          fit_adata_grid)
    return grid_sm_grad_adata



# This runs all the functions above to correctly output an AnnData structure of smoothed gene expression
def run_inplace_smoothing_and_gradient_analysis(adata, GAUSSIAN_SD, KERNEL_RADIUS):
    # WARNING - THIS IS LIKELY SLOW (TODO: Make parallelized version)
    
    ## Run in-place guassian and partial gaussian smoothing (takes 7 minutes with ~24000 points)
    
    # run apply_gauss_kernel_inplace function to get smoothed gaussian
    X_inplace_sm = apply_gauss_kernel_inplace(adata, 
                                  gaussian_kernel_height, 
                                  GAUSSIAN_SD, 
                                  KERNEL_RADIUS)
    
    # run apply_gauss_kernel_inplace function to get x partial gaussian derivative
    X_xpartial_inplace_sm = apply_gauss_kernel_inplace(adata, 
                                           x_partial_gaussian_kernel_height, 
                                           GAUSSIAN_SD, 
                                           KERNEL_RADIUS)
    
    # run apply_gauss_kernel_inplace function to get y partial gaussian derivative
    X_ypartial_inplace_sm = apply_gauss_kernel_inplace(adata, 
                                           y_partial_gaussian_kernel_height, 
                                           GAUSSIAN_SD, 
                                           KERNEL_RADIUS)
    
    # find the magnitude of the gradient at each point in the inplace
    X_gradient_magnitude_inplace = np.sqrt(X_xpartial_inplace_sm**2 + X_ypartial_inplace_sm**2)
    
    # find the direction of the gradient at each point in the inplace
    X_gradient_direction_inplace = np.arctan2(X_ypartial_inplace_sm, X_xpartial_inplace_sm)

    # add the smoothed surfaces to the original adata as layers
    adata = add_inplace_X_smoothed_w_gradient_to_anndata(X_inplace_sm,
                                                         X_gradient_magnitude_inplace,
                                                         X_gradient_direction_inplace, 
                                                         adata)
    return adata



# TODO: IN_PLACE AND GRID-BASED GAUSSIAN SMOOTHING FUNCTIONS (PARALLELIZED)


# # This runs all the functions above to correctly output an AnnData structure of smoothed gene expression
# def P_run_inplace_smoothing_and_gradient_analysis(adata, GAUSSIAN_SD, KERNEL_RADIUS, NUM_PROCESSES):
    
#     # Function to apply Gaussian kernel in place
#     def apply_gauss_kernel_inplace_wrapper(args):
#         return apply_gauss_kernel_inplace(*args)

    
#     # PARALLELIZE gaussian smoothing
#     with ThreadPoolExecutor(max_workers=NUM_PROCESSES) as executor:
#         # Submit the computation to the executor to run in parallel
#         future = executor.submit(apply_gauss_kernel_inplace_wrapper, 
#                                  (adata, gaussian_kernel_height, GAUSSIAN_SD, KERNEL_RADIUS))
        
#         # Retrieve the result from the future
#         X_sm_inplace = future.result()

#     # PARALLELIZE partial X gaussian smoothing
#     with ThreadPoolExecutor(max_workers=NUM_PROCESSES) as executor:
#         # Submit the computation to the executor to run in parallel
#         future = executor.submit(apply_gauss_kernel_inplace_wrapper, 
#                                  (adata, x_partial_gaussian_kernel_height, GAUSSIAN_SD, KERNEL_RADIUS))
        
#         # Retrieve the result from the future
#         X_xpartial_inplace_sm = future.result()

#     # PARALLELIZE partial Y gaussian smoothing
#     with ThreadPoolExecutor(max_workers=NUM_PROCESSES) as executor:
#         # Submit the computation to the executor to run in parallel
#         future = executor.submit(apply_gauss_kernel_inplace_wrapper, 
#                                  (adata, y_partial_gaussian_kernel_height, GAUSSIAN_SD, KERNEL_RADIUS))
        
#         # Retrieve the result from the future
#         X_ypartial_inplace_sm = future.result()

    
#     # find the magnitude of the gradient at each point in the grid
#     X_gradient_magnitudes_inplace = np.sqrt(X_xpartial_inplace_sm**2 + X_ypartial_inplace_sm**2)
    
#     # find the direction of the gradient at each point in the grid
#     X_gradient_direction_inplace = np.arctan2(X_ypartial_inplace_sm, X_xpartial_inplace_sm)

#     # add the smoothed surfaces to the original adata as layers
#     adata = add_inplace_X_smoothed_w_gradient_to_anndata(X_sm_inplace,
#                                                          X_gradient_magnitude_inplace,
#                                                          X_gradient_direction_inplace, 
#                                                          adata)
#     return adata
