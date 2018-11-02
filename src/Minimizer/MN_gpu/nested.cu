/* @brief Contains most functions for MultiNest from nested.F90 from
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
#include <cooperative_groups.h>
namespace cg = cooperative_groups;

/* Generate some livepoints at the beginning of the minimization process.
 * This function does not resume from a file.
 *
 * \param p                 Points in unit hyperspace
 * \param l          
 * \param n_live            The number of livepoints to use       
 * \param func              Indication for the actual likelihood function
 * \param n_accepted        Number of accepted points
 * \param n_dims            Number of dims
 */
template<typedef value_t>
__device__
void gen_initial_live(
    value_t * p,
    value_t * p_tmp,
    value_t * l,
    int n_live,
    int func,
    int n_accepted,
    int n_dims)
{
    auto g = this_thread_block();
    int id = threadIdx.x + blockDim.x*blockIdx.x;
    int bid = blockDim.x;
    int tid = g.thread_rank();
    int nend = 0;
    int n_gen = n_live;
    int pt_per_thread = SDIV(n_gen, blockDim.x*gridDim.x);
    int k = 0;
    int n_generated = 0;
    int j = 0;
    int p_idx = 0;

    for(int generated=0; generated<n_live; generated++)
    {
        if(tid == 0)
            l[generated] = log_zero;
        g.sync();
        while(l[generated] <= log_zero)
        {   
            // Sample 32 points
            if(tid < n_dims*32)
                get_random(n_dims, p_tmp, bid); // start points
            // this should be another kernel where several 
            // blocks might calculate this
            for(p_idx=0; p_idx<32; p_idx++)
            {
                g.sync();
                loglike(p_tmp+p_idx*n_dims, n_dims, tot_par, 
                        l+generated, bid, g); 
                g.sync();
                if(l[generated] > log_zero) {
                    if(tid == 0)
                        for(int i=0; i<n_dims; i++)
                            p[i + n_dims*generated] = p_tmp[i + p_idx*n_dims];
                    generated++;
                    if(generated >= n_live) break;
            }
            if(generated >= n_live) break;
        }
    }
}

template<typedef value_t>
__device__
void get_random(
    int n_dims,
    value_t * x,
    int id)
{
    x[id] = ranmarns(id);
}
