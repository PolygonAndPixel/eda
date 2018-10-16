/* Run and compare different algorithms to estimate a distribution or to
 * find a maximum/minimum for functions, that take lots of time to evaluate
 * (not for the test functions proposed here though) and that are not smooth
 * (i.e. not everywhere differentiable), that have degeneracies or are not
 * convex. In short: find a maximum/minimum for ill-posed functions and
 * as a bonus provide an estimation of the density (e.g. for error calculations).
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 * */

#include "helper/read_xml.h"

#include <iostream>

void print_result(
    MinimizerResult &result) {

    std::cout << std::endl << "Best fit for " << result.function_name
        << " with " << result.minimizer_name << " is "
        << -result.best_fit << std::endl;
    std::cout << "It took " << result.n_lh_calls << " likelihood evaluations"
        << std::endl;
    std::cout << "It had an efficiency of " << result.lh_efficiency
        << std::endl;
    std::cout << "Found parameters are ";
    for(auto &v: result.params_best_fit) std::cout << v << ", ";
    std::cout << std::endl << std::endl;
}

std::vector<MinimizerResult> run_tests(
    Minimizer &minimizer,
    std::string &func_path) {

    m_d lower_bounds;
    m_d upper_bounds;
    std::vector<TestFunctions> test_funcs;
    load_likelihood_xml(func_path, test_funcs, lower_bounds, upper_bounds);

    std::vector<MinimizerResult> results;
    for(index_t i=0; i<test_funcs.size(); ++i) {
        std::cout << "Running " << test_funcs[i].get_name()
        << " with " << minimizer.get_name() << std::endl;
        MinimizerResult result = minimizer.Minimize(test_funcs[i],
            lower_bounds[i], upper_bounds[i]);
        print_result(result);
        results.push_back(result);
    }
    return results;
}
//
// /* Can be used to get the best estimated time. Not sure why I would need it.
//  *
//  */
// value_t best_time_seed(
//     Track track,
//     PulseMap t_map) {
//
//     ESource casc = track.get_source(0);
//     v_d shell_times;
//     v_d shell_charges;
//     v_d shell_distances;
//
//     DOM dom;
//     v_d times;
//     v_d charges;
//     while(t_map.get_next(dom, times, charges)) {
//         if(times.empty()) continue;
//         value_t d = dist(dom.get_pos(), casc.get_pos());
//         shell_times.push_back(times[0] - d/0.3);
//         shell_distances.push_back(d);
//         shell_charges.push_back(charges[9]);
//     }
//
//     value_t best_time = shell_times[0];
//     value_t best_val = std::numeric_limits<value_t>::min();
//
//     for(value_t &s_time: shell_times) {
//         v_d gamma_singles(shell_charges.size());
//         value_t gamma_sum = 0;
//         for(index_t i=0; i<shell_charges.size(); i++) {
//             gamma_singles[i] = shell_charges[i] * log(
//                 hit_time(shell_distances[i], shell_times[i] - s_time)
//                 * hit_charge(shell_distances[i]) + 1e-3);
//
//             gamma_sum += gamma_singles[i];
//         }
//         if(gamma_sum > best_val) {
//             best_val = gamma_sum;
//             best_time = s_time;
//         }
//     }
//     return best_time;
// }

index_t main(index_t argc, char* argv[]) {

    std::string instr;
    std::string testf;
	index_t n_runs = 1;
    if(argc > 2) {
        instr = argv[1];
        testf = argv[2];
    } else {
        std::cout << "Please enter an argument which minimizer to use\n";
        std::cout << "The first argument should point to your configuration"
            << " file in ../../xml/Minimizer/\n";
        std::cout << "The second argument should point to your likelihood"
            << " configuration in ../../xml/likelihood/\n";
		std::cout << "An optional third argument can be used for the number "
			<< "of times you want to execute the algorithm\n";
        return 1;
    }
	if(argc > 3) {
		n_runs = atoi(argv[3]);
	}
    // We could basically load different minimizers by looping over
    // configuration files and execute every one of them.
    std::vector<std::unique_ptr<Minimizer>> minimizers;
    load_minimizer_xml(instr, minimizers);

	for(index_t i=0; i<n_runs; ++i) {
	    for(auto &minimizer: minimizers) {
	        // Do whatever you feel like with the result.
	        std::vector<MinimizerResult> results = run_tests(*minimizer, testf);
			for(auto &r: results) {
				std::string filename = r.function_name + "_" + r.minimizer_name;
				std::ofstream ofile(filename.c_str(),
					std::ofstream::out  | std::ofstream::app);
				ofile << r.lh_efficiency << "\t" << r.n_lh_calls << "\t" << r.best_fit
					<< "\t" << "0\t";
				for(auto &p: r.params_best_fit) {
					ofile << p << "\t";
				}
				ofile << "\n";
				ofile.close();
			}
	    }
	}

    return 0;
}
