// TODO: mat_mul, transpose

/** Warp reduce as seen in 
 * https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
 *
 */

#include "Minimizer/MN_gpu/linear_algebra.cuh"

template<
    typename value_t,
    uint block_size>
__device__
void warp_reduce(volatile value_t * sdata, uint tid, uint n_dims)
{
    if(block_size >= 64) sdata[tid] += sdata[tid + 32 + 32%n_dims];
    if(block_size >= 32) sdata[tid] += sdata[tid + 16 + 16%n_dims];
    if(block_size >= 16) sdata[tid] += sdata[tid + 8 + 8%n_dims];
    if(block_size >= 8) sdata[tid] += sdata[tid + 4 + 4%n_dims];
    if(block_size >= 4) sdata[tid] += sdata[tid + 2 + 2%n_dims];
    if(block_size >= 2) sdata[tid] += sdata[tid + 1];
}

/** Calculate the sum over all values along all dimensions separately. 
 *  Might divide by the value given in norm. Uses one block and shared memory.
 *
 * \param pts           The point set for the reduction
 * \param n_dims        The dimensionality of the points
 * \param reduced       On out: the reduction for each dimension
 * \param norm          Factor to multiply the reduced values with
 */
template<
    typename value_t,
    uint block_size>
__device__
void sum_reduce_2D(
    value_t * pts,
    uint n_dims,
    uint n_pts,
    value_t * reduced,
    value_t norm,
    value_t * cache)
{
    uint tid = threadIdx.x;
    value_t * sdata = cache;
    sdata[tid] = 0;
    uint stride = SDIV(block_size, n_dims)*n_dims;
    for(uint i=tid; i<n_pts-stride; i+=stride*2)
    {
        sdata[tid] += pts[i] + pts[i + stride];
    }
    // In case stride is bigger than block_size, we have to add the missing
    // points
    if(stride > block_size)
        for(uint i=tid + stride-block_size; i<n_pts-stride; i+=2*stride)
        {   
            if(tid > stride-block_size)
                sdata[tid] += pts[i] + pts[i + stride];
        }
    __syncthreads();

    // Reduce with all threads until we are small enough for warp size
    // We need to carefully calculate which values to add
    if(block_size >= 512) 
    {
        if(tid < 256)
            sdata[tid] += sdata[tid + 256 + 256%n_dims];
        if(256%n_dims != 0)
            if(tid > 256%n_dims && tid < n_dims)
                sdata[tid] += sdata[tid + 256 - n_dims];
        __syncthreads();
    }
    if(block_size >= 256) 
    {
        if(tid < 128)
            sdata[tid] += sdata[tid + 128 + 128%n_dims];
        if(128%n_dims != 0)
            if(tid > 128%n_dims && tid < n_dims)
                sdata[tid] += sdata[tid + 128 - n_dims];
        __syncthreads();
    }
    if(block_size >= 128) 
    {
        if(tid < 64)
            sdata[tid] += sdata[tid + 64 + 64%n_dims];
        if(64%n_dims != 0)
            if(tid > 64%n_dims && tid < n_dims)
                sdata[tid] += sdata[tid + 64 - n_dims];
        __syncthreads();
    }
    if(tid < 32) warp_reduce<value_t, block_size>(sdata, tid, n_dims);
    if(tid < n_dims) reduced[tid] = sdata[tid] * norm;
}