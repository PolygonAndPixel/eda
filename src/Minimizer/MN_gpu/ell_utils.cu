/* @brief Contains most functions for MultiNest from utils1.F90 from
 * MultiNest v3.8
 *
 * References:
 * 'MultiNest: an efficient and robust Bayesian inference tool for cosmology and particle physics'
 * F. Feroz, M.P. Hobson, M. Bridges. Sep 2008. 14 pp.
 * Published in Mon.Not.Roy.Astron.Soc. 398 (2009) 1601-1614
 * DOI: 10.1111/j.1365-2966.2009.14548.x
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 */

/* @brief: Evolve an ellipsoid where a point has been rejected/inserted.
 *
 * \param inserted          Point inserted or rejected
 * \param n_pt              Number of points after rejection/before insertion
 * \param n_dims            Dimensionality
 * \param new_pt            Point to be rejected/inserted
 * \param pts               Points after rejection/before insertion
 * \param mean              Centroid of the ellipsoid
 * \param eval              Eigenvalues of the ellipsoid
 * \param inv_cov           Inverse covariance matrix of the ellipsoid
 * \param k_fac             Point enlargement before (on out: after) 
                            rejection/insertion
 * \param eff               Volume enlargement before (on out: after) 
                            rejection/insertion
 * \param vol               Volume before (on out: after) 
                            rejection/insertion
 * \param target_vol        Target volume
 */

// Overload some functions according to the type.
// Needed for the template value_t
inline double forward_max(double a, double b) { return fmax(a,b);}
inline float forward_max(float a, float b) { return fmaxf(a,b);}
// inline cublasStatus_t forward_mm(
//     cublasHandle_t handle,
//     cublasOperation_t transa, cublasOperation_t transb,
//     int m, int n, int k,
//     const float           *alpha,
//     const float           *A, int lda,
//     const float           *B, int ldb,
//     const float           *beta,
//     float           *C, int ldc) 
// {
//     return cublasSgemm(cublasHandle_t handle,
//         cublasOperation_t transa, cublasOperation_t transb,
//         int m, int n, int k,
//         const float           *alpha,
//         const float           *A, int lda,
//         const float           *B, int ldb,
//         const float           *beta,
//         float           *C, int ldc);
// }
// inline cublasStatus_t forward_mm(
//     cublasHandle_t handle,
//     cublasOperation_t transa, cublasOperation_t transb,
//     int m, int n, int k,
//     const double          *alpha,
//     const double          *A, int lda,
//     const double          *B, int ldb,
//     const double          *beta,
//     double          *C, int ldc)
// {
//     cublasStatus_t cublasDgemm(cublasHandle_t handle,
//         cublasOperation_t transa, cublasOperation_t transb,
//         int m, int n, int k,
//         const double          *alpha,
//         const double          *A, int lda,
//         const double          *B, int ldb,
//         const double          *beta,
//         double          *C, int ldc);
// }

#define CUDART_PI_F 3.141592654f
// See https://docs.nvidia.com/cuda/cublas/index.html#device-api
#define IDX2F(i,j,ld) ((((j)-1)*(ld))+((i)-1))

template<typename value_t>
__device__
void scale_factor(
    const int & n_pt,
    const int & n_dims,
    const value_t * pt,
    const value_t & mean,
    const value_t * inv_cov,
    value_t & k_max)
{
    // TODO: Move that to a global helper array
    value_t temp_p[n_dims*n_pt];
    value_t pt_m[n_dims];
    // cublasHandle_t cnpHandle;
    // cublasStatus_t status = cublasCreate(&cnpHandle);
    k_max = 0.0;
    value_t scale_fac = 0;
    for(int k=0; k<n_pt; k++)
    {
        for(int i=0; i<n_dims; i++)
            pt_m[i] = temp_p[i + k*n_dims];
            scale_fac = mat_mul(mat_mul(pt_m, inv_cov), transpose(pt_m));
            if(scale_fac > k_max) k_max = scale_fac;
    }
    // cublasDestroy(cnpHandle);
}


template<typename value_t>
__device__
value_t volume_of_ell(
    int n_dims,
    value_t * eval,
    value_t k_fac)
{
    value_t prod = 1.0;
    for(int i=0; i<n_dims; i++) prod *= eval[i];
    value_t vol = sqrtf(prod * (k_fac**n_dims));

    if(n_dims%2 == 0)
    {
        for(int i=2; i<=n_dims; i+=2)
            vol *= 2.0*CUDART_PI_F/(value_t) i;
    } else 
    {
        vol *= 2;
        for(int i=3; i<=n_dims; i+=2)
            vol *= 2.0*CUDART_PI_F/(value_t) i;
    }
    return vol;
}



template<typename value_t>
__device__
void evolve_ell(
    bool inserted,
    int n_pt,
    int n_dims,
    value_t * new_pt,
    value_t * pts,
    value_t * mean,
    value_t * eval,
    value_t * inv_cov,
    value_t k_fac,
    value_t eff,
    value_t vol,
    value_t target_vol)
{
    // This check should be a bit more like if vol < epsilon
    if(target_vol == 0.0 || eval[n_dims] == 0)
    {
        k_fac = 0;
        eff = 1;
        vol = 0;
        return;
    }
    if(inserted)
    {
        // sanity says, if n_pt == 0, we should abort
        if(n_pt == 0) 
        {
            fprintf(stderr, 
                "Can not insert point in an ellipsoid with no points\n");
            exit(0);
        }
        // Calculate the scale factor
        value_t new_k_fac;
        scale_factor(1, n_dims, new_pt, mean, inv_cov, new_k_fac);
        if(new_k_fac > k_fac)
        {
            eff *= k_fac/new_k_fac;
            k_fac = new_k_fac;
            if(eff < 1.0)
            {
                eff = 1.0;
                vol = volume_of_ell(n_dims, eval, k_fac);
            }
        }
        // If the point has been inserted to a cluster with n_pt > 0,
        // then k_fac remains the same.
        // Scale eff if target vol is bigger
        if(target_vol > vol)
        {
            eff *= (target_vol/vol) ** (2.0/(value_t) n_dims);
            vol = target_vol;
        } else 
        {
            // Target volume is smaller, so calculate the point volume
            // and scale eff only if required
            vol = volume_of_ell(n_dims, eval, k_fac);
            // Scale enlargement factor
            if(vol < target_vol)
            {
                eff = forward_max(1.0, 
                    (target_vol/vol) ** (2.0/(value_t) n_dims));
                vol = target_vol;
            } else 
            {
                eff = 1.0;
            }
        }
    } else 
    {
        if(n_pt == 0)
        {
            vol = k_fac = 0.0;
            eff = 1.0;
            return;            
        }
        // If point has been rejected, check if it was at the boundary
        scale_factor(1, n_dims, new_pt, mean, inv_cov, new_k_fac);
        // In case of boundary point, find new k_fac
        if(new_k_fac == k_fac)
        {
            scale_factor(n_pt, n_dims, pts, mean, inv_cov, k_fac);
            if(k_fac > new_k_fac*eff && n_pt > 1)
            {
                fprintf(stderr, 
                    "Problem within evolve_ell: %f, %f, %f, %f\n",
                    k_fac, new_k_fac, eff, n_pt);
                exit(0);
            }
            eff = 1.0;
            vol = ell_vol(n_dims, eval, k_fac);
        }
        // Scale enlargement factor
        if(vol < target_vol)
        {
            eff = (target_vol/vol) ** (2.0/(value_t) n_dims);
            vol = target_vol;
        }
    }
}

template<typename value_t>
__device__
value_t calc_ell(
    value_t * S,
    value_t * E,
    int n_dims)
{

}

/** Calculate mean, covariance matrix, inverse covariance matrix, eigenvalues,
 *  eigenvectors, determinante of the covariance matrix and enlargement of
 *  a given point set.
 *  One block calculates on one set.
 *
 */
template<typename value_t>
__device__
void calc_ell_prop(
    int n_pt,
    int n_dims,
    value_t * pts,
    value_t target_vol,
    value_t * mean,
    value_t * cov,
    value_t * invcov,
    value_t * t_mat,
    value_t * eigen_vec,
    value_t * eigen_val,
    value_t & det_cov,
    value_t & k_fac,
    value_t & eff,
    value_t & vol,
    value_t * cache)
{
    auto tid = threadIdx.x;
    sum_reduce_2D(pts, n_dims, mean, 1.0/n_pt, cache);

}
