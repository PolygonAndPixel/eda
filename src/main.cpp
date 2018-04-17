#include <boost/cstdint.hpp>
#include <iostream>
#include <string>

#include "helper/abbreviations.h"
#include "Minimizer/SampleSpace.h"
#include "Minimizer/MAPS.h"
#include "Minimizer/PolyChord.h"
#include "Minimizer/Minimizer.h"
#include "likelihood/TestFunctions.h"
#include "Minimizer/MinimizerResult.h"

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

void run_tests(
    Minimizer &sampler) {

    uint32_t ndims = 2;
    v_d lower_bounds(ndims);
    v_d upper_bounds(ndims);
    TestFunctions test_func;
    MinimizerResult result;

    // Himmelblau's function
    test_func.set_func("himmelblau", ndims);
    lower_bounds[0] = -6;
    lower_bounds[1] = -6;
    upper_bounds[0] = 6;
    upper_bounds[1] = 6;
    std::cout << "Running Himmelblau" << std::endl;
    result = sampler.Minimize(test_func, lower_bounds, upper_bounds);
    print_result(result);

    // Townsend function
    test_func.set_func("townsend", ndims);
    lower_bounds[0] = -2.25;
    lower_bounds[1] = -2.5;
    upper_bounds[0] = 2.5;
    upper_bounds[1] = 1.75;
    std::cout << "Running Townsend" << std::endl;
    result = sampler.Minimize(test_func, lower_bounds, upper_bounds);
    print_result(result);

    // Rosenbrock function
    test_func.set_func("rosenbrock", ndims);
    lower_bounds[0] = -2;
    lower_bounds[1] = -1;
    upper_bounds[0] = 2;
    upper_bounds[1] = 3;
    std::cout << "Running Rosenbrock" << std::endl;
    result = sampler.Minimize(test_func, lower_bounds, upper_bounds);
    print_result(result);

    // Eggholder function
    test_func.set_func("eggholder", ndims);
    lower_bounds[0] = -512;
    lower_bounds[1] = -512;
    upper_bounds[0] = 512;
    upper_bounds[1] = 512;
    std::cout << "Running Eggholder" << std::endl;
    result = sampler.Minimize(test_func, lower_bounds, upper_bounds);
    print_result(result);

    // Gaussian shell
    test_func.set_func("gaussian_shell", ndims);
    std::cout << "Running Gaussian Shell" << std::endl;
    result = sampler.Minimize(test_func, lower_bounds, upper_bounds);
    print_result(result);
}

int main(int argc, char* argv[]) {

    char *instr = nullptr;
    if(argc > 1) {
        instr = argv[1];
    } else {
        std::cout << "Please enter an argument which minimizer to use\n";
        std::cout << "Possible arguments are:\n";
        std::cout << "maps     Using MAPS\n";
        std::cout << "poly     Using PolyChord\n";
        std::cout << "sample   Using SampleSpace to sample the function\n";
        std::cout << "multi    Using MultiNest\n";
        return 1;
    }

    if(strcmp(instr, "maps") == 0) {

        /////// MAPS
        double tolerance = 1e-4;
        uint32_t max_iter = 1<<12;
        uint32_t min_iter = 1<<10;
        uint32_t n_start_points = 1000;
        uint32_t size_sub_pop = n_start_points/10;
        uint32_t n_selected = n_start_points/2;

        MAPS sampler(tolerance, max_iter, min_iter, 0, n_start_points,
                     size_sub_pop, 9, n_selected);
        std::string path = "../../output/MAPS/";
        sampler.set_output(path);
        run_tests(sampler);
    }

    if(strcmp(instr, "poly") == 0) {

        /////// PolyChord
        double tolerance = 1e-4;

        PolyChord sampler(tolerance);
        std::string path = "../../output/PolyChord";
        sampler.set_output(path);
        run_tests(sampler);
    }

    if(strcmp(instr, "sample") == 0) {

        /////// SampleSpace
        // Parameters for sampling the space
        uint32_t n_iterations = 1000000;
        uint32_t max_points = n_iterations/10;
        uint32_t seed = 1025;

        SampleSpace sampler(n_iterations, max_points, seed, true);
        std::string path = "../../output/SampleSpace";
        sampler.set_output(path);
        run_tests(sampler);
    }

    return 0;
}
