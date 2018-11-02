/* @brief Contains some utilities to construct random numbers and seed
 * those properly.
 *
 * References:
 * 'MultiNest: an efficient and robust Bayesian inference tool for cosmology and particle physics'
 * F. Feroz, M.P. Hobson, M. Bridges. Sep 2008. 14 pp.
 * Published in Mon.Not.Roy.Astron.Soc. 398 (2009) 1601-1614
 * DOI: 10.1111/j.1365-2966.2009.14548.x
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 */

/* Initialize random numbers according to the number of threads.
 *
 * 
 */
template<typename value_t>
__device__ init_random_ns(
    int i = 42)
{
    int kl, ij;
    int id = threadIdx.x + blockDim.x*blockIdx.x;

    if(i >= 0)
    {
        kl = 9373;
        ij = (i+(k-1))*45;
    } else
    {
        // I could add something that uses system clock
    }
    rmarinns(ij, kl, id);
}




/* Taken from utils.f90 from MultiNest v3.8
 * This is the initialization routine for the randomNS number generator 
 * ranmarns()
 * The seed variables can have values between:    0 <= IJ <= 31328
 *                                                0 <= KL <= 30081
 * The randomNS number sequences created by these two seeds are of sufficient 
 * length to complete an entire calculation with. For example, if sveral 
 * different groups are working on different parts of the same calculation,
 * each group could be assigned its own IJ seed. This would leave each group
 * with 30000 choices for the second seed. That is to say, this randomNS 
 * number generator can create 900 million different subsequences -- with 
 * each subsequence having a length of approximately 10^30.
 *
 * Use IJ = 1802 & KL = 9373 to test the randomNS number generator. The
 * subroutine ranmarns should be used to generate 20000 randomNS numbers.
 * Then display the next six randomNS numbers generated multiplied by 4096*4096
 * If the randomNS number generator is working properly, the randomNS numbers
 * should be:
 *           6533892.0  14220222.0  7275067.0
 *           6172232.0  8354498.0   10633180.0
 *
 * \param ij        seed variable between 0 and 31328
 * \param kl        seed variable between 0 and 30081
 * \param id        thread id
 */
template<typename value_t>
__device__ void rmarinns(
    int ij,
    int kl,
    int id)
{
    if(ij < 0) ij += 31328;
    if(ij > 31328) ij = ij%31328;
    if(kl < 0) kl += 30081;
    if(kl > 30081) kl = kl%30081;

    // if(ij < 0 || ij > 31328 || kl < 0 || kl > 30081)
    // {
    //     // Abort...
    // }

    int i = ((ij/177)%177)+2;
    int j = (ij%177)+2;
    int k = ((kl/169)%178)+1;
    int l = kl%169;

    for(int ii=1; ii<98; ii+)
    {
        value_t s = 0.0;
        value_t t = 0.5;
        for(int jj=1; jj<25; jj+)
        {
            int m = (((i*j)%179)*k)%179;
            i = j;
            j = k;
            k = m;
            l = (53*l+1)/169;
            if((l*m)%64 >= 32) s += t;
            t *= 0.5;
        }
        u[id + ii*blockDim.x] = s;
    }
    c[id] = 362436.0 / 16777216.0;
    cd[id] = 7654321.0 / 16777216.0;
    cm[id] = 16777213.0 / 16777216.0;
    i97[id] = 97;
    j97[id] = 33;
}

/* Geroge Marsaglia's random number generator (Florida State University)
 * Report: FSU-SCRI-87-50
 * Modified by F. James to produce an array of pseudorandom numbers.
 */
template<typename value_t>
value_t ranmarns(int idg)
{
    
}