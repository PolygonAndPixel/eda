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
        << result.best_fit << std::endl;
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
    for(uint32_t i=0; i<test_funcs.size(); ++i) {
        std::cout << "Running " << test_funcs[i].get_name()
        << " with " << minimizer.get_name() << std::endl;
        MinimizerResult result = minimizer.Minimize(test_funcs[i],
            lower_bounds[i], upper_bounds[i]);
        print_result(result);
        results.push_back(result);
    }
    return results;
}

int main(int argc, char* argv[]) {

    std::string instr;
    std::string testf;
    if(argc > 2) {
        instr = argv[1];
        testf = argv[2];
    } else {
        std::cout << "Please enter an argument which minimizer to use\n";
        std::cout << "The first argument should point to your configuration"
            << " file in ../../xml/Minimizer/\n";
        std::cout << "The second argument should point to your likelihood"
            << " configuration in ../../xml/likelihood/\n";
        return 1;
    }
    // We could basically load different minimizers by looping over
    // configuration files and execute every one of them.
    std::vector<std::unique_ptr<Minimizer>> minimizers;
    load_minimizer_xml(instr, minimizers);
    for(auto &minimizer: minimizers) {
        // Do whatever you feel like with the result.
        std::vector<MinimizerResult> results = run_tests(*minimizer, testf);
    }

    return 0;
}
