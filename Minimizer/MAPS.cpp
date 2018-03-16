/**
 * @brief EA algorithm "Maintaining and Processing Sub-models (MAPS)"
 * Calculates the maximum likelihood of a given function.
 * It assumes that promising areas, aka areas with an optimum, have large
 * changes in the likelihood. This may be a pitfal for edgy regions.
 *
 * References:
 * Improving Estimation of Distribution Algorithm on
 * Multimodal Problems by Detecting Promising Area
 * by Peng Yang, Ke Tang and Xiaofen Lu
 * Submitted to IEEE Transactions on Cybernetics (2014)
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 */

#include "MAPS.h"

/** Constructor and destructor **/
MAPS::MAPS(
    double tolerance,
    uint32_t max_iter,
    uint32_t min_iter,
    uint32_t max_points,
    uint32_t n_start_points,
    uint32_t size_sub_pop,
    uint32_t max_sub_pops,
    uint32_t n_selected,
    uint32_t n_sub_selected,
    uint32_t seed,
    bool dump_points) : Minimizer(tolerance, max_iter, min_iter, 
                                  max_points, seed, dump_points)
{
    n_start_points = n_start_points_;
    size_sub_pop = size_sub_pop_;
    max_sub_pops = max_sub_pops_;
    n_selected = n_selected_;
    n_sub_selected = n_sub_selected_;
}

/** Check if a population is similar to a discarded one and return true.
 *  Also stores the mean of the discarded population in discarded_pops_.
 *
 *  \param pop      The population to check
 *  \param real_cov Calculate the real covariance matrix if true. Else just
 *                  use a tenth of the identity matrix.
 *
 *  \return         True if population is similar to a discarded one
 * */
bool MAPS::check_premature(
    m_d pop,
    bool real_cov) {

    v_d mean = get_center(pop);
    m_d cov = get_cov(pop, real_cov);
    for(auto discarded_means=discarded_pops_.begin();
        discarded_means<discarded_pops_.end(); discarded_means++) {

        if(is_similar(mean, *discarded_means, cov)) {
            discarded_pops_.push_back(mean);
            return true;
        }
    }
    return false;
}



/** Construct a histogram in 1D sub-space. The direction is the direction
 *  to go in this subspace. This can be the Cartesian coordinates given
 *  in the problem or transformed axes. Here we use PCA
 *
 *  \param pop              The population aka the active samples
 *  \param direction        The axes along we want to explore the likelihood
 *                          May be Cartesian or better eigenvectors
 *  \param dim              The current dimension in the parameter above we
 *                          are interested in
 *  \param freq_pop         The population sorted in its bins. The row is the
 *                          bin and the vectors are concatenated in each row.
 *
 *  \return                 A histogram
 * */
v_i MAPS::make_histogram(
    m_d pop,
    m_d direction,
    uint32_t dim,
    m_d &freq_pop) {

    // Project all the points to the new coordinate in the given direction
    v_d projected_pop = prod(pop, direction[dim]);
    double max_value = projected_pop[0];
    double min_value = projected_pop[0];
    for(uint32_t d=1; d<projected_pop.size(); d++) {
        if(projected_pop[d] > max_value) max_value = projected_pop[d];
        if(projected_pop[d] < min_value) min_value = projected_pop[d];
    }
    uint32_t n_bins = SDIV(projected_pop.size(), 5);
    v_i freq(n_bins);
    freq_pop.resize(n_bins);
    uint32_t width = SDIV(max_value-min_value, n_bins);
    auto not_projected = pop.begin();
    for(auto p=projected_pop.begin(); p<projected_pop.end();
        p++, not_projected++) {

        uint32_t idx = SDIV(*p-min_value, width);
        freq[idx]++;
        freq_pop[idx].insert(freq_pop[idx].end(), not_projected->begin(),
            not_projected->end());
    }
    return freq;
}

/** Project all vectors in a to the vector b.
 *
 *  \param a        A vector of vector, e.g. a population.
 *  \param b        A vector to project to.
 *
 *  \return         The projected values.
 * */
v_d prod(
    m_d a,
    v_d b) {

    v_d projected(a.size());
    uint32_t i=0;
    for(auto m_iter=a.begin(); m_iter<a.end(); ++m_iter, ++i) {
        double sum = 0;
        for(auto a_iter=m_iter->begin(), b_iter=b.begin();
            (a_iter<m_iter->end() && b_iter<b.end()); ++a_iter, ++b_iter) {

            sum += (*a_iter) * (*b_iter);
        }
        projected[i] = sum;
    }
    return projected;
}

/** Find a higher bin to determine fast changings in the likelihood region
 *  Return a vector with each bin that was higher than another one within
 *  some threshold
 *
 *  \param freq     The histogram which had been calculated along one dimension
 *
 *  \return
 * */
v_i MAPS::find_higher_bin(
    v_i freq) {

    bool reset_flag = false;
    uint32_t num = 0;
    uint32_t large = 1;
    v_i higher_bins(1);
    for(uint32_t i=1; i<freq.size(); i++) {
        if(freq[i] > freq[large]) large = i;
        if(EULER_CONST * freq[i] < freq[large]) {
            num++;
            // Line above is also possible
            // higher_bins.insert_element(higher_bins.size(), large);
            higher_bins.push_back(large);
            reset_flag = true;
        }
        if(reset_flag && freq[i] > freq[i-1]) {
            large = i;
            reset_flag = false;
        }
    }
    return higher_bins;
}

/** Calculate the sub populations around the higher bins
 *
 *  \param higher_bins      The bins which are bigger than a bin before within
 *                          a threshold
 *  \param freq             The histogram which had been calculated along one
 *                          dimension
 *
 *  \return
 * */
std::vector<m_d> MAPS::confirm_bins(
    v_i higher_bins,
    v_i freq,
    m_d freq_pop) {

    std::vector<m_d> sub_pops;
    uint32_t n_bins = higher_bins.size();
    for(uint32_t i=1; i<n_bins; i++) {
        uint32_t left_frontier = n_bins;
        uint32_t right_frontier = n_bins;
        for(uint32_t j=higher_bins[i]; j>0; j--) {
            if(EULER_CONST * freq[j] < freq[higher_bins[i]]) {
                left_frontier = j;
            }
        }
        for(uint32_t j=higher_bins[i]; j<n_bins; j++) {
            right_frontier = j;
        }
        // Form a sub population if and only if a right and left frontier have
        // been found
        if(left_frontier < n_bins && right_frontier < n_bins) {
            m_d tmp_pop;
            for(uint32_t j=left_frontier; j<=right_frontier; j++) {
                tmp_pop.push_back(freq_pop[j]);
            }
            sub_pops.push_back(tmp_pop);
        }
    }
    return sub_pops;
}

/** Calculate the final sub populations
 *
 *  \param sub_pops      Sub population calculated around "higher" bins in one
                         dimension.
 *  \param direction     The axes along we want to explore the likelihood
 *                       May be Cartesian or better eigenvectors
 *  \param dim           The current dimension in the parameter above we
 *                       are interested in
 *  \param ndims         Dimensionallity of parameter space in terms of
 *                       free parameter for minimization
 *
 *  \return
 * */
std::vector<m_d> MAPS::iterative_observation(
    m_d sub_pop,
    m_d directions,
    uint32_t dim,
    uint32_t ndims) {

    m_d freq_pop;
    std::vector<m_d> fine_sub_pop;
    v_i freq = make_histogram(sub_pop, directions, dim, freq_pop);
    v_i higher_bins = find_higher_bin(freq);
    std::vector<m_d> new_sub_pops = confirm_bins(higher_bins, freq, freq_pop);
    if(dim < ndims) {
        for(auto sub_pop_p=new_sub_pops.begin();
            sub_pop_p<new_sub_pops.end(); sub_pop_p++) {

            fine_sub_pop = iterative_observation(*sub_pop_p,
                directions, dim+1, ndims);
        }
    }
    return fine_sub_pop;
}

/** Given a truncated population, calculate the directions to search
 *  for via PCA and call iterative_observation
 *
 *  \param offspring    A truncated population
 *
 *  \return             All the accepted sub populations with the llh as
 *                      extra dimension for each point
 * */
std::vector<m_d> MAPS::maintaining(
    m_d offspring,
    uint32_t ndims) {

    m_d eigen_v;
    v_d eigen_values(ndims);
    m_d cov(ndims, v_d(ndims));
    pca(offspring, eigen_v, cov, eigen_values, ndims);

    double eigen_threshold = 0.85; // We want to take the eigenvectors that
                                   // make up 85%, see paper page 4
    double eigen_sum = 0;

    for(auto eigen=eigen_values.begin(); eigen<eigen_values.end();
        eigen++) {

        eigen_sum += abs(*eigen);
    }
    m_d directions;
    eigen_sum *= eigen_threshold;
    double current_sum = 0;
    auto eigen_vec=eigen_v.begin();
    for(auto eigen=eigen_values.begin(); eigen<eigen_values.end();
        eigen++, eigen_vec++) {

        current_sum += abs(*eigen);
        if(current_sum > eigen_sum) break;
        directions.push_back(*eigen_vec);
    }
    std::vector<m_d> fine_sub_pops = iterative_observation(offspring,
        directions, 1, ndims);
    return fine_sub_pops;

}

/** Ignore premature promising areas aka delete populations that are
 *  attracted to such areas. Saves discarded populations to discarded_points.
 *
 *  \param estimated_sub_pops       All truncated populations
 *  \param ndims                    Dimensionallity of parameter space in terms
 *                                  of free parameter for minimization
 *
 *  \return                         All populations that are not discarded
 * */
std::vector<m_d> MAPS::processing(
    std::vector<m_d> estimated_sub_pops,
    uint32_t ndims) {

    boost::uniform_real<> uni_dist(0,1);
    boost::uniform_01<boost::minstd_rand> uf(intgen);
    for(auto pop=estimated_sub_pops.begin();
        pop<estimated_sub_pops.end(); pop++) {

        m_d tmp_pop;
        // Sample size_sub_pop_ many points
        for(uint32_t i=0; i<n_sub_selected_; i++) {
            double rnd = uf();
            uint32_t idx = round(rnd * (pop->size()));
            tmp_pop.push_back((*pop)[idx]);
        }
        // Truncatedly select n_sub_selected_ many points for this
        // population
        m_d selected_pop = truncatedly_select(tmp_pop, n_sub_selected_, ndims);

        // Update the parameters with the selected individuals
        *pop = evolve_population(selected_pop, ndims);
    }

    // Calculate the means of each population first to avoid further calculations
    m_d centers(estimated_sub_pops.size(), v_d(ndims));
    auto c_iter = centers.begin();
    for(auto pop=estimated_sub_pops.begin(); pop<estimated_sub_pops.end(); 
        pop++, c_iter++) {

        c_iter->push_back(get_center(*pop));
    }
    std::vector<m_d> final_selected_pops;
    uint32_t idx_current_pop = 0;
    for(auto pop=estimated_sub_pops.begin();
        pop<estimated_sub_pops.end(); pop++) {

        // Check if this population is premature
        // Need to store the value of the best fit and how long it didn't
        // change. If it didn't change within 1e-4 for then generations,
        // it is considered premature.
        if(check_premature(*pop)) {
            idx_current_pop++;
            continue;
        }

        // Check if population is similar to any other one
        uint32_t compare_idx = 0;
        bool break_up = false;
        for(auto c_iter=centers.begin(); c_iter<centers.end();
            c_iter++, compare_idx++) {

            if(compare_idx == idx_current_pop) continue;
            if(is_similar(*c_iter, centers[idx_current_pop], cov_)) {
                // Check which one is better
                if((*pop)[ndims] < (*c_iter)[ndims]) {
                    break_up = true;
                    break;

                }
            }
        }
        if(break_up) {
            idx_current_pop++;
            continue;
        }

        // Check if population is similar with any discarded one
        for(auto dis_iter=discarded_pops_.begin();
            dis_iter<discarded_pops_.end(); dis_iter++) {

            if(is_similar(*dis_iter, centers[idx_current_pop], cov_)) {
                break_up = true;
                break;
            }
        }
        if(break_up) {
            idx_current_pop++;
            continue;
        }

        // Add current population to final one
        final_selected_pops.push_back(*pop);
        idx_current_pop++;
    }

    return final_selected_pops;
}

/** Function to map from the unit hypercube to Theta in the physical space.
 *
 * \param cube              Hypercube coordinates of point that
 *                          shall be evaluated.
 * \param ndims             Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 *
 *  \return                 Physical coordinates
 * */
v_d MAPS::to_physics(
    v_d cube,
    uint32_t ndims) {

    v_d theta(ndims);

    for (int i=0; i<ndims; i++) {
        theta(i) = (this->lower_bnds[i]
            + (this->upper_bnds[i] - this->lower_bnds[i])
            * cube[i]);
    }
    return theta;
}

/** Calculate the principal component analysis of a given matrix.
 *  Calculates the covariance matrix first.
 *  See http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html#ga442c43fca5493590f8f26cf42fed4044
 *  for the used PCA.
 *
 *  \param in               The input matrix
 *  \param eigen_v          The eigenvectors on return
 *  \param cov              The covariance matrix on return
 *  \param eigen_values     The eigenvalues in ascending order
 *  \param ndims            Dimensionallity of parameter space in terms of
 *                          free parameter for minimization
 *  \param real_cov         Calculate the real covariance matrix if true, else
 *                          calculate a diagonal matrix
 *                          ((h-l)/max_sub_pops_) * I
 *                          with I the identity matrix, h and l the boundaries
 *                          of the parameter space in the current population
 *                          in the current dimension.
 *
 *  \return                 = 0: sucessful
 *                          < 0: if return = -i, the i-th argument had an
 *                               illegal value in the covariance matrix
 *                          > 0: if return = i, the pca failed to converge. i
 *                               diagonal elements of an intermediate
 *                               tridiagonal form  did not converge to zero.
 * */
int MAPS::pca(
    m_d in,
    m_d eigen_v,
    m_d cov,
    v_d eigen_values,
    uint32_t ndims,
    bool real_cov) {

    in.resize(in.size(), ndims);
    if(real_cov) {
        // Adjust data s.t. it is centered around 0
        v_d means(ndims);
        for(auto v=in.begin(); v<in.end(); v++) {
            uint32_t i = 0;
            for(auto e=v.begin(); e<v.end(); e++, i++) {
                means(i) += *e;
            }
        }
        for(auto e=means.begin(); e<means.end(); e++) {
            *e /= in.size();
        }

        for(auto v=in.begin(); v<in.end(); v++) {
            *v -= means;
        }
        // Easier to read but is it also faster than plain for loops?
        cov = 1/(in.size()-1) * outer_prod(trans<m_d>(in), in);

        // // Calculate the symmetric covariance matrix
        // for(uint32_t row=0; row<ndims; row++) {
        //     for(uint32_t col=0; col<ndims; col++) {
        //         if(row <= col) {
        //             double cov_v = 0.0;
        //             for(uint32_t i=0; i<in.size(); i++) {
        //                 cov_v += in(i, row)*in(i, col);
        //             }
        //             cov(row, col) = cov(row, col) = cov_v/in.size();
        //         }
        //     }
        // }
    } else {
        if(cov.size() != ndims || cov[0].size() != ndims)
            cov.resize(ndims, ndims, false);
        for(uint32_t i=0; i<ndims; i++) {
            double high = in[0][i];
            double low = in[0][i];
            for(uint32_t j=1; j<in.size(); j++) {
                if(in[j][i] > high) high = in[i][j];
                if(in[j][i] < low) low = in[i][j];
            }
            const double constant = SDIV((high-low), max_sub_pops_);
            cov[i, i] = constant;
         }
    }
    // Calculate the eigenvectors v^-1 C v = D
    // matrix_layout, jobz, uplo, n, a, lda, w
    // Resize if someone forgot to allocate memory.
    if(eigen_v.size() != ndims || eigen_v[0].size() != ndims)
        eigen_v.resize(cov.size(), cov[0].size(), false);
    if(eigen_values.size != ndims) eigen_values.resize(ndims, false);
    std::copy(cov.begin(), cov.end(), eigen_v.begin());
    int info = LAPACKE_dsyev(
        LAPACK_ROW_MAJOR,               // The matrix layout
        'V',                            // Calculate the eigenvectors too
        'L',                            // Handle the input matrix as lower triangular matrix
        eigen_v.size()*in[0].size(),    // Order (size) of the matrix
        &eigen_v[0],                    // Input matrix and the eigenvectors on output
        eigen_v.size(),                // Leading order dimension
        &eigen_values[0]);              // The eigenvalues on output

    return info;
}

/** Check if two vectors are similar by calculating the Mahalanobis distance
 *  of their means in parameter space.
 *
 *  \param a            Vector a
 *  \param b            Vector b
 *  \param cov          The covariance matrix to use.
 *  \param epsilon      Precision at which both are similar
 *
 *  \return             True if both are similar
 * */
bool MAPS::is_similar(
    v_d a,
    v_d b,
    m_d cov,
    double epsilon) {

    // sqrt((a-b)T cov^-1 (a-b))
    v_d a_b(a.size());
    for(uint32_t i=0; i<a.size(); i++) {
        a_b = a[i] - b[i];
    }
    double distance = inner_prod(prod<v_d>(prod<v_d>(trans<v_d>(a_b),
        inverse<m_d>(cov)) ,a_b)))
    return (distance < epsilon);

}

/** Check if two populations are similar by calculating the Mahalanobis distance
 *  of their means in parameter space.
 *
 *  \param A            Matrix of a population
 *  \param B            Matrix of a population
 *  \param cov          The covariance matrix to use.
 *  \param ndims        Dimensionallity of parameter space in terms of
 *                      free parameter for minimization
 *  \param epsilon      Precision at which both are similar
 *
 *  \return             True if both are similar
 * */
bool MAPS::is_similar(
    m_d A,
    m_d B,
    m_d cov,
    uint32_t ndims,
    double epsilon) {

    v_d a = get_center(A);
    v_d b = get_center(B);
    return is_similar(a, b, cov, epsilon, ndims);
}

/** Calculates the center of a population.
 *
 *  \param pop                  The matrix containing the population. Each row
 *                              is one individual of this populatio
 *  \param ignore_last_col      Ignore the last column in the population
 *
 *  \return                     The center of the population
 * */
v_d MAPS::get_center(
    m_d pop,
    bool ignore_last_col) {

    v_d centre(pop[0].size());
    for(auto v=pop.begin(); v<pop.end(); v++) {
        centre += *v;
    }
    centre /= pop.size();
    return centre;
}

/** Calculates either the real covariance matrix or a fake one,
 *  (h-l)/10 * I
 *  with h the higher bound, l the lower bound of the parameter space and I
 *  the parameter space.
 *
 *  \param pop                  The population to calculate the covariance
 *                              for
 *  \param ndims                Dimensionallity of parameter space in terms of
 *                              free parameter for minimization
 *  \param ignore_last_col      Ignore the last column in the population
 *  \param real_cov             Calculate the real covariance matrix if true,
 *                              else calculate a diagonal matrix
 *                              ((h-l)/max_sub_pops_) * I
 *                              with I the identity matrix, h and l the
 *                              boundaries of the parameter space in the
 *                              current population in the current dimension.
 *
 *  \return                     The covariance matrix
 * */
m_d MAPS::get_cov(
    m_d pop,
    uint32_t ndims,
    bool ignore_last_col,
    bool real_cov) {

    if(ignore_last_col) pop.resize(pop.size(), ndims, true);
    m_d cov(ndims, v_d(ndims));
    if(real_cov) {
        // Adjust data s.t. it is centered around 0
        v_d means(ndims);
        for(auto v=pop.begin(); v<pop.end(); v++) {
            uint32_t i = 0;
            for(m_d::iterator2 e=v.begin(); e<v.end(); e++, i++) {
                means[i] += *e;
            }
        }
        for(auto e=means.begin(); e<means.end(); e++) {
            *e /= pop.size();
        }

        for(auto v=pop.begin(); v<pop.end(); v++) {
            *v -= means;
        }
        // Easier to read but is it also faster than plain for loops?
        cov = 1/(pop.size()-1) * outer_prod(trans<m_d>(pop), pop);

        // // Calculate the symmetric covariance matrix
        // for(uint32_t row=0; row<ndims; row++) {
        //     for(uint32_t col=0; col<ndims; col++) {
        //         if(row <= col) {
        //             double cov_v = 0.0;
        //             for(uint32_t i=0; i<in.size(); i++) {
        //                 cov_v += in(i, row)*in(i, col);
        //             }
        //             cov(row, col) = cov(row, col) = cov_v/in.size();
        //         }
        //     }
        // }
    } else {
        // for(uint32_t i=0; i<ndims; i++) {
        //     double high = pop[0, i];
        //     double low = pop[0, i];
        //     for(uint32_t j=1; j<in.size(); j++) {
        //         if(pop[j, i] > high) high = pop[i, j];
        //         if(pop[j, i] < low) low = pop[i, j];
        //     }
        //     const double constant = SDIV((high-low), max_sub_pops_);
        //     cov[i, i] = constant;
        //  }
         return cov_;
    }
    return cov;
}

/** Execute the MAPS algorithm and store the current best and worst fit.
 *
 *  \param ndims        Dimensionallity of parameter space in terms of
 *                      free parameter for minimization
 * */
void MAPS::execute_maps(
    uint32_t ndims) {

    uint32_t n_iter = 0;
    // This is just preliminarily
    identity_matrix<double, ndims> identity;
    boost::uniform_real<> uni_dist(0,1);
    boost::uniform_01<boost::minstd_rand> uf(intgen);
    while(true) {
        m_d init_samples(n_start_points_, v_d(ndims+1));
        for(uint32 i=0; i<n_start_points_; i++) {
            for(uint32_t j=0; j<ndims; j++) {
                init_samples(i, j) = uf();
            }
            init_samples(i, ndims) = get_llh(to_physics(init_samples[i],
                ndims));
        }

        m_d offspring = truncatedly_select(init_samples, n_selected_, ndims);

        // sub_pops is SS from the paper, page 5
        vector<m_d> sub_pops = maintaining(offspring);

        // Sort points in each sub_pop descending of its fitness
        for(auto pop=sub_pops.begin(); pop<sub_pops.end(); pop++) {

            std::sort(pop->begin(), pop->end(),
                boost::bind(&v_d[ndims], _1) <
                boost::bind(&v_d[ndims], _2));
        }
        // estimated_pops is ES from the paper, page 5
        std::vector<m_d> estimated_pops(max_sub_pops_, v_d(1, v_d(ndims+1)));
        // Pick point from the first order on
        // if estimated population is not full and point != any of the
        // points in the estimated population and in the premature populations
        // then add point to estimated population
        int n_estimated_models = 0;
        int offset = 0;
        for(auto pop=sub_pops.begin();
            pop<sub_pops.end(); pop++) {

            v_d best_fit_individual = (*pop)[0];
            // Check for similarity in discarded points and estimated_pop
            bool add_this = true;
            for(auto es=estimated_pops.begin();
                es<estimated_pops.end(); es++) {

                if(is_similar(*es, best_fit_individual, identity)) {
                    add_this = false;
                    break;
                }
            }
            if(!add_this) continue;
            for(auto ds=discarded_points.begin();
                ds<discarded_points.end(); ds++) {

                if(is_similar(*ds, best_fit_individual, identity, ndims)) {
                    add_this = false;
                    break;
                }
            }
            if(!add_this) continue;

            std::copy(pop->begin(), pop->end(), estimated_pops.begin()+offset);
            n_estimated_models++;
            offset += sizeof(*pop);
            // Check if estimated_pop is "full" although I am not sure
            // what the maximal size should be... Hence I just leave it be.
            // MAPS usually uses "less than 10" (page 6, Table 1)
            if(n_estimated_models == max_sub_pops_) break;
        }
        while(true) {
            estimated_pops = processing(estimated_pops);
            // Get best and worst llh from the each population and check if all
            // reach the stopping criterion (best and worst fit in this
            // population are smaller than the tolerance) or if the maximum
            // iterations is reached.
            n_iter++;
            if(n_iter == max_iter_ and min_iter_ < n_iter) return;
            lh_bestFit_ = estimated_pops[0][0][ndims];
            for(uint32_t d=0; d<ndims; d++) {
                params_best_fit[d] = estimated_pops[0][0][d];
            }
            lh_worstFit_ = lh_bestFit_;
            for(auto e_iter=estimated_pops.begin();
                e_iter<estimated_pops.end(); e_iter++) {

                for(m_d::iterator pop_iter=e_iter->begin(); pop_iter<e_iter->end();
                    pop_iter++) {

                    if(lh_bestFit_ > (*pop_iter)[ndims]) {
                        for(uint32_t d=0; d<ndims; d++) {
                            params_best_fit[d] = (*pop_iter)[d];
                        }
                        lh_bestFit_ = (*pop_iter)[ndims];
                    }
                    if(lh_worstFit_ < (*pop_iter)[ndims]) {
                        lh_worstFit_ = (*pop_iter)[ndims];
                    }
                }
            }

            if(fabs(lh_bestFit_ - lh_worstFit_) < precision_criterion_) {
                return;
            }
            // If n_estimated_models == 0 (the models are the null set), start over
            // Else restart processing.
            if(estimated_pops.empty()) break;
        }
    }
    
}

/** Truncatedly selects n elements from pop and stores them in pop_out. On
 *  return, pop is sorted in a descending order.
 *
 *  \param pop          The input population to take samples from.
 *  \param n            The number of samples to take. Needs to be smaller
 *                      than the number of rows in pop.
 *  \param ndims        Dimensionallity of parameter space in terms of
 *                      free parameter for minimization
 * 
 *  \return             The selected samples.
 * */
m_d MAPS::truncatedly_select(
    m_d &pop,
    uint32_t n
    uint32_t ndims) {

    // Sort points descending of its fitness
    std::sort(pop.begin(), pop.end(),
        boost::bind(&v_d(ndims), _1) <
        boost::bind(&v_d(ndims), _2));

    m_d pop_out(n, v_d(n));
    // We truncatedly select R \leq N of them
    for(uint32_t i=0; i<n; i++) {
        pop_out[i] = pop[i];
    }
    return pop_out;
}

/** Sample new points with the given set of points. Sample only as many
 *  points as there are individuals in the population.
 *  This can be Metropolis Hastings, chordial sampling, nested sampling with
 *  the population as active points, etc...
 *
 *  \param pop      The input population to take samples from.
 *  \param ndims    Dimensionallity of parameter space in terms of
 *                  free parameter for minimization.
 *
 *  \return         A new set of sampled points.
 * */
m_d MAPS::evolve_population(
    m_d pop,
    uint32_t ndims) {

    m_d new_pop(pop.size(), pop[0].size());

    // Metropolis Hastings using a normal distribution
    boost::uniform_real<> uni_dist(0,1);
    boost::uniform_01<boost::minstd_rand> uf(intgen);
    double variance_2 = 0.4;
    double multiplier = 1.0/(std::sqrt(variance_2*M_PI));
    variance_2 = -1.0/variance_2;
    uint32_t p=0;
    for(auto pop_iter=pop.begin();
        pop_iter<pop.end(); pop_iter++, p++) {

        v_d new_p(ndims+1);
        for(uint32_t i=0; i<ndims; i++) {
            double value = uf();
            value = value - (*pop_iter)[i];
            new_p[i] = multiplier * std::exp(value*value*variance_2);
        }
        v_d new_p_phys = to_physics(new_p, ndims);
        new_p[ndims] = get_llh(new_p_phys);
        if(new_p[ndims] > (*pop_iter)[ndims]) {
            new_pop[p] = new_p;
            result.lh_efficiency += 1;
            continue;
        }
        if(new_p[ndims] / (*pop_iter)[ndims] > uf()) {
            new_pop[p] = new_p;
            result.lh_efficiency += 1;
            continue;
        }
        new_pop[p] = *pop_iter;
    }
    return new_pop;
    /////// End of Metropolis Hastings

}

/** Required Minimize() function for every minimizer. Sets the bounds.
 *
 *  \param test_func        The function which shall be minimized
 *  \param lower_bounds     The lower bounds for each dimension
 *  \param upper_bounds     The upper bounds for each dimension
 * 
 *  \return                 The result of the minimization
 * */
MinimizerResult
MAPS::Minimize(
    TestFunctions test_func,
    v_d lower_bounds,
    v_d upper_bounds ) {
    
    reset_calls();
    results.clear();
    
    upper_bnds = upper_bounds;
    lower_bnds = lower_bounds;
    test_func_ = &test_func;
    file_name_ = test_func_->get_name();
    // Build the identity covariance matrix
    cov_ = m_d(test_func_->get_ndims(), v_d(test_func_->get_ndims(), 0));
    for(uint32_t row=0; row<test_func_->get_ndims(); row++) {
        cov_[row][row] = 1;
    }

    execute_maps(test_func_->get_ndims());
    
    result.minimizer_name = "MAPS";
    result.best_fit = lh_bestFit_;
    result.params_best_fit = params_best_fit;
    result.lh_efficiency /= result.n_lh_calls;
    
    return result;
}
