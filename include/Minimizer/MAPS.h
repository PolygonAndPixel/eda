#ifndef MAPS_H_INCLUDED
#define MAPS_H_INCLUDED

#include "MinimizerResult.h"
#include "Minimizer.h"
#include "helper/abbreviations.h"

#include <random>
#include <lapacke.h>
#include <cblas.h>
#include <fstream>

class MAPS : public Minimizer {
public:

    MAPS(double tolerance, int max_iter, int min_iter,
         int max_points=0, int n_start_points=1000,
         int size_sub_pop=100, int max_sub_pops=9,
         int n_selected=500, int n_sub_selected=25, double size_factor=1.5,
         int seed=1025, bool dump_points=false);

    virtual std::unique_ptr<Minimizer> clone() const override {
        return std::make_unique<MAPS>(*this);
    }

    std::string get_name();

    /// core method: minimizer a given function with given initial conditions
    MinimizerResult Minimize(TestFunctions test_func, v_d lower_bounds,
                             v_d upper_bounds);

    // Check if a population is premature and store it to the discarded ones
    // if true
    bool check_premature(m_d pop, uint32_t idx, uint32_t ndims,
                         double epsilon=1e-3);

    // Construct a histogram in 1D sub-space. The direction is the direction
    // to go in this subspace. This can be the Cartesian coordinates given
    // in the problem or transformed axes. Here we use PCA
    v_i make_histogram(m_d pop, m_d direction, uint32_t dim, m_d &freq_pop);
    // Calculate the product pairwise of all rows with b
    v_d prod(m_d a, v_d b);
    // Find a higher bin to determine fast changings in the likelihood region
    // Return a vector with each bin that was higher than another one within
    // some threshold
    v_i find_higher_bin(v_i freq);
    // Calculate the sub populations around the higher bins
    std::vector<m_d> confirm_bins(v_i higher_bins, v_i freq, m_d freq_pop,
                                  uint32_t ndims);
    // Calculate the final sub populations
    std::vector<m_d> iterative_observation(m_d sub_pop, m_d direction,
                                           uint32_t dim, uint32_t ndims);

    // Given a truncated population, calculate the directions
    // to search for via PCA and call iterative_observation
    std::vector<m_d> maintaining(m_d offspring, uint32_t ndims);
    // Ignore premature promising areas aka delete populations that are
    // attracted to such areas
    std::vector<m_d> processing(std::vector<m_d> estimated_sub_pops,
        uint32_t ndims);

    void pca(m_d in, m_d & eigen_v, m_d & cov, v_d & eigen_values,
             uint32_t ndims, bool real_cov=false);

    // Check by using Mahalanobis distance of parameters + llh with identity matrix as
    // covariance matrix. Might change that later to some covariance matrix
    // which leads to the question, which one?
    // TODO: Change epsilon to (h-l)/100 which only makes sense if all
    // parameters are in the same range though...
    bool is_similar(v_d a, v_d b, m_d & cov, double epsilon=1e-3);
    bool is_similar(m_d A, m_d B, m_d & cov, uint32_t ndims, double epsilon=1e-3);

    v_d get_center(m_d pop, bool ignore_last_col=true);
    m_d get_cov(m_d pop, uint32_t ndims, bool ignore_last_col=true,
                bool real_cov=false);

    void execute_maps(uint32_t ndims);

    // Also sorts pop
    m_d truncatedly_select(m_d &pop, uint32_t n, uint32_t ndims);

    /// Transform point from hypercube to physical space
    v_d to_physics(v_d cube, uint32_t ndims);

    // Sample with given population.
    m_d evolve_population(m_d pop, uint32_t ndims);


private:

    m_d discarded_pops_;
    uint32_t n_start_points_, max_sub_pops_;
    uint32_t size_sub_pop_, n_selected_, n_sub_selected_;
    m_d cov_;
    std::vector<std::pair<uint32_t, double>> premature_list_;
    double size_factor_;
    uint32_t n_init_samples;
};

#endif
