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
const bool DEBUG (false);

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
    n_start_points_ = n_start_points;
    size_sub_pop_ = size_sub_pop;
    max_sub_pops_ = max_sub_pops;
    n_selected_ = n_selected;
    n_sub_selected_ = n_sub_selected;
    size_factor_ = EULER_CONST;
    size_factor_ = 1.1;
}

/** Check if a population is premature.
 *  If the fit of the best individual didn't change within 1e-4 for ten
 *  generations, it is considered premature (aka doesn't contain a global
 *  optimum). If it is premature, it is recorded to discarded_pop_.
 *
 *  \param pop      The population to check
 *  \param idx      The index of this population in premature_list_.
 *  \param ndims    Dimensionallity of parameter space in terms of
 *                  free parameter for minimization
 *  \param epsilon  Precision at which two fits are similar.
 *
 *  \return         True if population is premature
 * */
bool MAPS::check_premature(
    m_d pop,
    uint32_t idx,
    uint32_t ndims,
    double epsilon) {

    double best_fit = pop[0][ndims];
    for(auto pop_iter=pop.begin(); pop_iter<pop.end(); ++pop_iter) {
        if(best_fit > (*pop_iter)[ndims]) {
            best_fit = (*pop_iter)[ndims];
        }
    }
    // Check first if this population had been stored before.
    if(premature_list_.size() > idx+1) {
        if(abs(premature_list_[idx].second - best_fit) < epsilon) {
            if(premature_list_[idx].first == 10) {
                v_d mean = get_center(pop);
                discarded_pops_.push_back(mean);
                return true;
            } else {
                premature_list_[idx].first++;
                return false;
            }
        } else {
            premature_list_[idx].first = 1;
            premature_list_[idx].second = best_fit;
        }
    } else {
        premature_list_.emplace_back(1, best_fit);
    }
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
 *                          bin and the vectors are concatenated in each row
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
    double width = (max_value-min_value) / (double)n_bins;
    auto not_projected = pop.begin();
    for(auto p=projected_pop.begin(); p<projected_pop.end();
        p++, not_projected++) {

        uint32_t idx = (*p-min_value)/width;
        if(idx==n_bins) idx--;
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
v_d MAPS::prod(
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

        if(size_factor_ * freq[i] < freq[large]) {
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
 *  \param freq_pop         The population sorted in its bins. The row is the
 *                          bin and the vectors are concatenated in each row
 *
 *  \return                 The populations within the bins that are around
 *                          higher bins
 * */
std::vector<m_d> MAPS::confirm_bins(
    v_i higher_bins,
    v_i freq,
    m_d freq_pop,
    uint32_t ndims) {

    std::vector<m_d> sub_pops;
    uint32_t n_bins = higher_bins.size();
    for(uint32_t i=1; i<n_bins; i++) {
        uint32_t left_frontier = n_bins;
        uint32_t right_frontier = n_bins;
        for(uint32_t j=higher_bins[i]; j>0; j--) {
            if(size_factor_ * freq[j] < freq[higher_bins[i]]) {
                left_frontier = j;
            }
        }
        for(uint32_t j=higher_bins[i]; j<n_bins; j++) {
            if(size_factor_ * freq[j] < freq[higher_bins[i]]) {
                right_frontier = j;
            }
        }
        // Form a sub population if and only if a right and left frontier have
        // been found
        if(left_frontier < n_bins && right_frontier < n_bins) {
            m_d tmp_pop;
            for(uint32_t j=left_frontier; j<=right_frontier; j++) {
                uint32_t d = 0;
                v_d v(ndims+1);
                for(auto val: freq_pop[j]) {
                    v[d] = val;
                    // We also store the likelihood at the end of each vector
                    // hence the ndims+1
                    d = (d+1)%(ndims+1);
                    if(d==0) tmp_pop.push_back(v);
                }
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
    std::vector<m_d> new_sub_pops = confirm_bins(higher_bins, freq, freq_pop,
                                                 ndims);
    if(dim < ndims-1) {
        for(auto &sub_pop_p : new_sub_pops) {
            std::vector<m_d> fine_sub_pop_tmp = iterative_observation(sub_pop_p,
                directions, dim+1, ndims);
            for(auto &pop: fine_sub_pop_tmp) {
                fine_sub_pop.push_back(pop);
            }
        }
        return fine_sub_pop;
    }
    return new_sub_pops;
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
        directions, 0, ndims);
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

    for(auto &pop: estimated_sub_pops) {
        m_d tmp_pop;
        // Sample size_sub_pop_ many points
        for(uint32_t i=0; i<n_sub_selected_; i++) {
            double rnd = uf(intgen);
            uint32_t idx = round(rnd * (pop.size()-1));
            tmp_pop.push_back(pop[idx]);
        }
        // Truncatedly select n_sub_selected_ many points for this
        // population
        m_d selected_pop = truncatedly_select(tmp_pop, n_sub_selected_, ndims);

        // Update the parameters with the selected individuals
        pop = evolve_population(selected_pop, ndims);
    }

    // Calculate the means of each population first to avoid further calculations
    m_d centers(estimated_sub_pops.size(), v_d(ndims));
    for(auto const &pop: estimated_sub_pops) {
       centers.push_back(get_center(pop));
    }
    std::vector<m_d> final_selected_pops;
    uint32_t idx_current_pop = 0;
    uint32_t idx_premature = 0;

    std::vector<std::pair<uint32_t, double>> premature_list_tmp(premature_list_.size());
    std::copy(premature_list_.begin(), premature_list_.end(), premature_list_tmp.begin());

    for(auto pop=estimated_sub_pops.begin();
        pop<estimated_sub_pops.end(); ++pop, ++idx_premature, ++idx_current_pop) {

        if(DEBUG) {
            std::cout << "Checking another population" << std::endl;
        }
        // Check if this population is premature
        // Need to store the value of the best fit and how long it didn't
        // change. If it didn't change within 1e-4 for then generations,
        // it is considered premature.
        if(check_premature(*pop, idx_premature, ndims, precision_criterion_)) {
            if(DEBUG)
                std::cout << "This population is premature" << std::endl;
            // Rearrange the list of best fits and generations.
            premature_list_.erase(premature_list_.begin()+idx_premature);
            // Need to reduce the idx because we just erased an element.
            idx_premature--;
            continue;
        }

        // Check if population is similar to any other one
        uint32_t compare_idx = 0;
        bool break_up = false;
        for(auto c_iter=centers.begin(); c_iter<centers.end();
            ++c_iter, ++compare_idx) {

            if(compare_idx == idx_current_pop) continue;

            if(is_similar(*c_iter, centers[idx_current_pop], cov_,
                precision_criterion_)) {

                // Check which one is better
                if(premature_list_[idx_premature].second < premature_list_tmp[compare_idx].second) {
                    if(DEBUG) {
                        std::cout << "This population is similar to another one and inferior" << std::endl;
                    }
                    break_up = true;
                    break;

                }
            }
        }
        if(break_up) {
            // Rearrange the list of best fits and generations.
            premature_list_.erase(premature_list_.begin()+idx_premature);
            // Need to reduce the idx because we just erased an element.
            idx_premature--;
            continue;
        }

        // Check if population is similar with any discarded one
        // for(auto dis_iter=discarded_pops_.begin();
        //     dis_iter<discarded_pops_.end(); dis_iter++) {
        for(auto const &dis: discarded_pops_) {
            if(is_similar(dis, centers[idx_current_pop], cov_,
                precision_criterion_)) {

                if(DEBUG) {
                    std::cout << "This population is similar to a discarded one" << std::endl;
                }
                break_up = true;
                break;
            }
        }
        if(break_up) {
            // Rearrange the list of best fits and generations.
            premature_list_.erase(premature_list_.begin()+idx_premature);
            // Need to reduce the idx because we just erased an element.
            idx_premature--;
            continue;
        }
        if(DEBUG)
            std::cout << "Pushing a population\n";
        // Add current population to final one
        final_selected_pops.push_back(*pop);
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
        theta[i] = (this->lower_bnds[i]
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
void MAPS::pca(
    m_d in,
    m_d & eigen_v,
    m_d & cov,
    v_d & eigen_values,
    uint32_t ndims,
    bool real_cov) {

    cov = get_cov(in, ndims, true, real_cov);

    // Calculate the eigenvectors v^-1 C v = D
    // matrix_layout, jobz, uplo, n, a, lda, w
    // Resize if someone forgot to allocate memory.
    if(eigen_v.size() != ndims || eigen_v[0].size() != ndims) {
        for(auto &v: eigen_v) {
            v.resize(ndims);
        }
    }
    if(eigen_values.size() != ndims) {
        eigen_values.resize(ndims);
    }
    eigen_v = cov;
    char triang = 'U';
    char eigen = 'V';
    lapack_int n = eigen_v.size();
    lapack_int lda = eigen_v.size();
    v_d flat_eigen_v(eigen_v.size() * eigen_v[0].size());

    int info = LAPACKE_dsyev(
        LAPACK_ROW_MAJOR,               // The matrix layout
        eigen,                          // Calculate the eigenvectors too
        triang,                         // Store the output matrix as upper triangular matrix
        n,                              // Order (sqrt(size)) of the matrix
        &flat_eigen_v[0],               // Input matrix and the eigenvectors on output
        lda,                            // Leading order dimension
        &eigen_values[0]);              // The eigenvalues on output

    uint32_t v_size = eigen_v[0].size();
    for(uint32_t row=0; row<eigen_v.size(); ++row) {
        for(uint32_t col=0; col<v_size; ++col) {
            eigen_v[row][col] = flat_eigen_v[row*v_size + col];
        }
    }
}

/** Check if two vectors are similar by calculating the Mahalanobis distance
 *  of their means in parameter space.
 *                  distance = sqrt((a-b)^T cov^-1 (a-b))
 *  Calculate a-b first. Factorize cov, s.t. it is an upper (or lower)
 *  triangular matrix. Cholesky decomposition using forward substitution:
 *                  cov = L D L*
 *  Solve for y in
 *                  L y = (a-b)
 *                  distance = y^T*y
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
    m_d & cov,
    double epsilon) {

    v_d a_b(a.size());
    for(uint32_t i=0; i<a.size(); i++) {
        a_b[i] = a[i] - b[i];
    }
    // If we use the identity matrix, we won't need any fancy inverse matrix.
    // Just use Euclidean distance.
    if(&cov == &cov_) {
        double distance = 0.0;
        for(auto const &val : a_b) {
            distance += val*val;
        }
        return (sqrt(distance) < epsilon);
    }
    int info = 0;
    int lda = cov.size();
    int nrhs = 1;
    char triang = 'U';
    double * L = new double[cov.size()*cov.size()];
    std::copy(&(cov[0][0]), &(cov[0][0])+(cov.size()*cov.size()), L);
    double * y = new double[a_b.size()];
    std::copy(&(a_b[0]), &(a_b[0])+a_b.size(), y);
    // http://www.netlib.org/lapack/explore-html/d1/d7a/group__double_p_ocomputational_ga167aa0166c4ce726385f65e4ab05e7c1.html
    // cov * x = a_b
    LAPACK_dpotrs(&triang,        // Store the upper triangular part
            &lda,           // The order of cov (N)
            &nrhs,          // The number of columns on the right hand side
            L,            // Array of size (lda,N)
            &lda,           // Leading dimension of cov (lda)
            y,            // The right hand side. on exit the solution
            &nrhs,          // Leading dimension of a_b (ldb)
            &info);

    if(info < 0) {
        std::cout << "Error in MAPS::is_similar: Factorization of ";
        std::cout << "covariance matrix failed!. Illegal value at ";
        std::cout << abs(info) << std::endl;
        return false;
    } else if(info > 0) {
        std::cout << "Error in MAPS::is_similar: Factorization of ";
        std::cout << "covariance matrix failed!. leading minor of order ";
        std::cout << abs(info) << " is not positive definite!" << std::endl;
        return false;
    }
    double distance = 0;
    for(uint32_t i=0; i<a_b.size(); i++) {
        distance += y[i]*y[i];
    }

    return (sqrt(distance) < epsilon);
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
    m_d & cov,
    uint32_t ndims,
    double epsilon) {

    v_d a = get_center(A);
    v_d b = get_center(B);
    return is_similar(a, b, cov, precision_criterion_);
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

    v_d centre(pop[0].size(), 0.0);
    for(auto const &v : pop) {
        for(uint32_t i=0; i<v.size(); i++) {
            centre[i] += v[i];
        }
    }
    for(auto &c : centre) c /= pop.size();

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

    if(ignore_last_col) {
        for(auto &v : pop) {
            v.resize(ndims);
        }
    }
    m_d cov(ndims, v_d(ndims));
    if(real_cov) {
        // Adjust data s.t. it is centered around 0
        v_d means(ndims);
        for(auto const &v : pop) {
            uint32_t i = 0;
            for(auto const &e : v) {
                means[i] += e;
            }
        }
        for(auto &e : means) {
            e /= pop.size();
        }

        for(auto &v : pop) {
            for(uint32_t i=0; i<v.size(); i++) {
                v[i] -= means[i];
            }
        }

        // Calculate the symmetric covariance matrix
        for(uint32_t row=0; row<ndims; row++) {
            for(uint32_t col=0; col<ndims; col++) {
                if(row <= col) {
                    double cov_v = 0.0;
                    for(uint32_t i=0; i<pop.size(); i++) {
                        cov_v += pop[i][row]*pop[i][col];
                    }
                    cov[row][col] = cov[row][col] = cov_v/pop.size();
                }
            }
        }
    } else {
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
    premature_list_.clear();
    while(true) {
        m_d init_samples(n_start_points_, v_d(ndims+1));
        for(uint32_t i=0; i<n_start_points_; i++) {
            for(uint32_t j=0; j<ndims; j++) {
                init_samples[i][j] = uf(intgen);
            }
            init_samples[i][ndims] = get_llh(to_physics(init_samples[i],
                ndims));
        }

        m_d offspring = truncatedly_select(init_samples, n_selected_, ndims);

        // sub_pops is SS from the paper, page 5
        std::vector<m_d> sub_pops = maintaining(offspring, ndims);

        // Sort points in each sub_pop descending of its fitness
        for(auto &pop : sub_pops) {
            std::sort(pop.begin(), pop.end(),
                [](const v_d & a, const v_d & b) {
                    return a[a.size()-1] < b[b.size()-1];});
        }
        // estimated_pops is ES from the paper, page 5
        std::vector<m_d> estimated_pops;
        // Pick point from the first order on
        // if estimated population is not full and mean != any of the
        // means of the estimated populations and in the premature populations
        // then add point to estimated population
        uint32_t n_estimated_models = 0;
        // what is offset for?
        uint32_t offset = 0;
        // Precalculate the means of the populations to compare their
        // similarity later on
        m_d means;
        for(auto const &pop : sub_pops) {
            means.push_back(get_center(pop));
        }
        uint32_t means_idx = 0;
        v_i means_est_idx;
//         std::cout << "Amount of sub_pops: " << sub_pops.size() << std::endl;
        // Build the list of best individuals of each population
        for(auto const &pop : sub_pops) {
            // Check for similarity in discarded points and estimated_pop
            bool add_this = true;
            for(auto es : means_est_idx) {
                if(is_similar(means[es], means[means_idx], cov_,
                    precision_criterion_)) {

                    add_this = false;
                        break;
                }
            }

            if(!add_this) {means_idx++; continue;}
            for(auto const &ds : discarded_pops_) {
                if(is_similar(ds, means[means_idx], cov_,
                    precision_criterion_)) {

                    add_this = false;
                    break;
                }
            }
            if(!add_this) {means_idx++; continue;}

            // formerly offset
            estimated_pops.push_back(pop);
            n_estimated_models++;
            means_est_idx.push_back(means_idx);
            means_idx++;
            offset += sizeof(pop);
            // Add the best individual to the premature list to see if
            // the population is going to be premature
            premature_list_.emplace_back(1, pop[0][ndims]);

            // Check if estimated_pop is "full" although I am not sure
            // what the maximal size should be... Hence I just leave it be.
            // MAPS usually uses "less than 10" (page 6, Table 1)
            if(n_estimated_models == max_sub_pops_) break;
        }
        while(true) {
            if(DEBUG)
                std::cout << "starting processing. " << estimated_pops[0][0].size() << std::endl;
            estimated_pops = processing(estimated_pops, ndims);

            // Get best and worst llh from the each population and check if all
            // reach the stopping criterion (best and worst fit in this
            // population are smaller than the tolerance) or if the maximum
            // iterations is reached.
            // Start over if no population has been found.
            n_iter++;
            if(n_iter == max_iter_ and min_iter_ < n_iter) return;
            if(estimated_pops.size() == 0) {
                break;
            }
            if(DEBUG) {
                std::cout << "Now1 " << estimated_pops.size() << std::endl;
                std::cout << "Now2 " << estimated_pops[0].size() << std::endl;
                std::cout << "Now3 " << estimated_pops[0][0].size() << std::endl;
            }
            lh_bestFit_ = estimated_pops[0][0][ndims];
            for(uint32_t d=0; d<ndims; d++) {
                params_best_fit[d] = estimated_pops[0][0][d];
            }
            lh_worstFit_ = lh_bestFit_;
            for(auto const &e : estimated_pops) {
                for(auto const &pop : e) {
                    if(lh_bestFit_ > pop[ndims]) {
                        for(uint32_t d=0; d<ndims; d++) {
                            params_best_fit[d] = pop[d];
                        }
                        lh_bestFit_ = pop[ndims];
                    }
                    if(lh_worstFit_ < pop[ndims]) {
                        lh_worstFit_ = pop[ndims];
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
    uint32_t n,
    uint32_t ndims) {

    // Sort points descending of its fitness
    std::sort(pop.begin(), pop.end(),
        [](const v_d & a, const v_d & b) {
            return a[a.size()-1] < b[b.size()-1];});

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

    m_d new_pop(pop.size(), v_d(pop[0].size()));

    // Metropolis Hastings using a normal distribution
    // For each point: Sample a new point and take the new one if it satisfies
    // the Metropolis criterion, else take the old point
    double variance_2 = 0.5;
    double multiplier = 1.0/(std::sqrt(variance_2*M_PI));
    variance_2 = -1.0/variance_2;
    uint32_t p=0;
    for(auto pop_iter=pop.begin(); pop_iter<pop.end(); pop_iter++, p++) {
        v_d new_p(ndims+1);
        for(uint32_t i=0; i<ndims; i++) {
            double value = uf(intgen);
            value = value - (*pop_iter)[i];
            new_p[i] = multiplier * std::exp(value*value*variance_2);
        }
        v_d new_p_phys = to_physics(new_p, ndims);
        new_p[ndims] = get_llh(new_p_phys);
        if(new_p[ndims] > (*pop_iter)[ndims]) {
            new_pop[p] = new_p;
            result.lh_efficiency += 1;
            if(dump_points_) {
                std::ofstream ofile((base_dir_+file_name_).c_str(),
                    std::ofstream::out  | std::ofstream::app);
                for(auto &p: new_p_phys)
                    ofile << p << "\t";
                ofile << new_p[ndims] << "\t" << std::endl;
                ofile.close();
            }
            continue;
        }
        if(new_p[ndims] / (*pop_iter)[ndims] > uf(intgen)) {
            new_pop[p] = new_p;
            result.lh_efficiency += 1;
            if(dump_points_) {
                std::ofstream ofile((base_dir_+file_name_).c_str(),
                    std::ofstream::out  | std::ofstream::app);
                for(auto &p: new_p_phys)
                    ofile << p << "\t";
                ofile << new_p[ndims] << "\t" << std::endl;
                ofile.close();
            }
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
    if(DEBUG)
        std::setbuf(stdout, NULL);

    upper_bnds = upper_bounds;
    lower_bnds = lower_bounds;
    test_func_ = &test_func;
    file_name_ = test_func_->get_name();
    params_best_fit.resize(test_func_->get_ndims());
    // Build the identity covariance matrix
    cov_ = m_d(test_func_->get_ndims(), v_d(test_func_->get_ndims(), 0));
    for(uint32_t row=0; row<test_func_->get_ndims(); row++) {
        cov_[row][row] = 1;
    }
    if(dump_points_) {
        std::ofstream ofile((base_dir_+file_name_).c_str(),
            std::ofstream::out  | std::ofstream::app);

        for(int j=0; j<test_func_->get_ndims(); j++)
            ofile << "Param" << j << "\t";
        ofile << std::endl;
        ofile.close();
    }
    execute_maps(test_func_->get_ndims());
    result.function_name = file_name_;
    result.minimizer_name = "MAPS";
    result.best_fit = lh_bestFit_;
    result.params_best_fit = params_best_fit;
    result.lh_efficiency /= result.n_lh_calls;

    return result;
}
