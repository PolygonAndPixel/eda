#include <boost/cstdint.hpp>
#include <iostream>

#include "helper/abbreviations.h"
#include "Minimizer/SampleSpace.h"
#include "Minimizer/MAPS.h"
#include "likelihood/TestFunctions.h"
#include "Minimizer/MinimizerResult.h"


int main(int argc, char* argv[]) {

    double tolerance = 1e-4;
    uint32_t max_iter = 1<<12;
    uint32_t min_iter = 1<<10;
    uint32_t ndims = 2;
    uint32_t n_start_points = 1000;
    uint32_t size_sub_pop = n_start_points/10;
    uint32_t n_selected = n_start_points/2;

    MAPS sampler(tolerance, max_iter, min_iter, 0, n_start_points,
                 size_sub_pop, 9, n_selected);
    std::string path = "output/MAPS/";
    sampler.set_output(path);
    v_d lower_bounds(ndims);
    v_d upper_bounds(ndims);
    TestFunctions test_func;
    /*
    // Townsend function
    test_func.set_func("townsend", ndims);
    lower_bounds[0] = -2.25;
    lower_bounds[1] = -2.5;
    upper_bounds[0] = 2.5;
    upper_bounds[1] = 1.75;
    MinimizerResult result = sampler.Minimize(test_func, lower_bounds,
        upper_bounds);
    std::cout << std::endl << "Best fit for " << result.function_name
        << " with " << result.minimizer_name << " is "
        << result.best_fit << std::endl;
    std::cout << "It took " << result.n_lh_calls << " likelihood evaluations"
        << std::endl;
    std::cout << "It had an efficiency of " << result.lh_efficiency
        << std::endl;

    // Try the eggholder function
    test_func.set_func("eggholder", ndims);
    lower_bounds[0] = -512;
    lower_bounds[1] = -512;
    upper_bounds[0] = 512;
    upper_bounds[1] = 512;
    result = sampler.Minimize(test_func, lower_bounds,
        upper_bounds);
    std::cout << std::endl << "Best fit for " << result.function_name
        << " with " << result.minimizer_name << " is "
        << result.best_fit << std::endl;
    std::cout << "It took " << result.n_lh_calls << " likelihood evaluations"
        << std::endl;
    std::cout << "It had an efficiency of " << result.lh_efficiency
        << std::endl;

    // Rosenbrock function
    test_func.set_func("rosenbrock", ndims);
    lower_bounds[0] = -2;
    lower_bounds[1] = -1;
    upper_bounds[0] = 2;
    upper_bounds[1] = 3;
    result = sampler.Minimize(test_func, lower_bounds,
        upper_bounds);
    std::cout << std::endl << "Best fit for " << result.function_name
        << " with " << result.minimizer_name << " is "
        << result.best_fit << std::endl;
    std::cout << "It took " << result.n_lh_calls << " likelihood evaluations"
        << std::endl;
    std::cout << "It had an efficiency of " << result.lh_efficiency
        << std::endl;
    */
    // Himmelblau's function
    test_func.set_func("himmelblau", ndims);
    lower_bounds[0] = -6;
    lower_bounds[1] = -6;
    upper_bounds[0] = 6;
    upper_bounds[1] = 6;
    MinimizerResult result = sampler.Minimize(test_func, lower_bounds,
        upper_bounds);
    std::cout << std::endl << "Best fit for " << result.function_name
        << " with " << result.minimizer_name << " is "
        << result.best_fit << std::endl;
    std::cout << "It took " << result.n_lh_calls << " likelihood evaluations"
        << std::endl;
    std::cout << "It had an efficiency of " << result.lh_efficiency
        << std::endl;

    // Gaussian shell
    test_func.set_func("gaussian_shell", ndims);
    result = sampler.Minimize(test_func, lower_bounds,
        upper_bounds);
    std::cout << std::endl << "Best fit for " << result.function_name
        << " with " << result.minimizer_name << " is "
        << result.best_fit << std::endl;
    std::cout << "It took " << result.n_lh_calls << " likelihood evaluations"
        << std::endl;
    std::cout << "It had an efficiency of " << result.lh_efficiency
        << std::endl;
    /*
    uint32_t n_iterations = 1000000;
    uint32_t max_points = n_iterations/10;
    uint32_t seed = 1025;
    uint32_t ndims = 2;

    std::string path = "output/SampleSpace/";
    SampleSpace sampler(n_iterations, max_points, seed, true);
    sampler.set_output(path);
    v_d lower_bounds(ndims);
    v_d upper_bounds(ndims);

    // Townsend function
    TestFunctions test_func("townsend", ndims);
    lower_bounds[0] = -2.25;
    lower_bounds[1] = -2.5;
    upper_bounds[0] = 2.5;
    upper_bounds[1] = 1.75;
    MinimizerResult result = sampler.Minimize(test_func, lower_bounds,
        upper_bounds);
    std::cout << std::endl << "Best fit for " << result.function_name
        << " with " << result.minimizer_name << " is "
        << result.best_fit << std::endl;

    // Try the eggholder function
    test_func.set_func("eggholder", ndims);
    lower_bounds[0] = -512;
    lower_bounds[1] = -512;
    upper_bounds[0] = 512;
    upper_bounds[1] = 512;
    result = sampler.Minimize(test_func, lower_bounds,
        upper_bounds);
    std::cout << std::endl << "Best fit for " << result.function_name
        << " with " << result.minimizer_name << " is "
        << result.best_fit << std::endl;

    // Rosenbrock function
    test_func.set_func("rosenbrock", ndims);
    lower_bounds[0] = -2;
    lower_bounds[1] = -1;
    upper_bounds[0] = 2;
    upper_bounds[1] = 3;
    result = sampler.Minimize(test_func, lower_bounds,
        upper_bounds);
    std::cout << std::endl << "Best fit for " << result.function_name
        << " with " << result.minimizer_name << " is "
        << result.best_fit << std::endl;

    // Himmelblau's function
    test_func.set_func("himmelblau", ndims);
    lower_bounds[0] = -6;
    lower_bounds[1] = -6;
    upper_bounds[0] = 6;
    upper_bounds[1] = 6;
    result = sampler.Minimize(test_func, lower_bounds,
        upper_bounds);
    std::cout << std::endl << "Best fit for " << result.function_name
        << " with " << result.minimizer_name << " is "
        << result.best_fit << std::endl;

    // Gaussian shell
    test_func.set_func("gaussian_shell", ndims);
    result = sampler.Minimize(test_func, lower_bounds,
        upper_bounds);
    std::cout << std::endl << "Best fit for " << result.function_name
        << " with " << result.minimizer_name << " is "
        << result.best_fit << std::endl;
    */
    return 0;
}
