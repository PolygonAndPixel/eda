#include <iostream>

#include "helper/abbreviations.h"
#include "Minimizer/SampleSpace.h"
#include "Minimizer/MAPS.h"
#include "Minimizer/PolyChord.h"
#include "Minimizer/Minimizer.h"
#include "Minimizer/MultiNest.h"
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
    Minimizer &sampler,
    std::string &func_name) {

    uint32_t ndims = 2;
    v_d lower_bounds(ndims);
    v_d upper_bounds(ndims);
    TestFunctions test_func;
    MinimizerResult result;

    if(func_name == HIMMEL || func_name == ALL) {
        // Himmelblau's function
        test_func.set_func(HIMMEL, ndims);
        lower_bounds[0] = -6;
        lower_bounds[1] = -6;
        upper_bounds[0] = 6;
        upper_bounds[1] = 6;
        std::cout << "Running " << HIMMEL << std::endl;
        result = sampler.Minimize(test_func, lower_bounds, upper_bounds);
        print_result(result);
    }

    if(func_name == TOWN || func_name == ALL) {
        // Townsend function
        test_func.set_func(TOWN, ndims);
        lower_bounds[0] = -2.25;
        lower_bounds[1] = -2.5;
        upper_bounds[0] = 2.5;
        upper_bounds[1] = 1.75;
        std::cout << "Running " << TOWN << std::endl;
        result = sampler.Minimize(test_func, lower_bounds, upper_bounds);
        print_result(result);
    }

    if(func_name == ROSEN || func_name == ALL) {
        // Rosenbrock function
        test_func.set_func(ROSEN, ndims);
        lower_bounds[0] = -2;
        lower_bounds[1] = -1;
        upper_bounds[0] = 2;
        upper_bounds[1] = 3;
        std::cout << "Running " << ROSEN << std::endl;
        result = sampler.Minimize(test_func, lower_bounds, upper_bounds);
        print_result(result);
    }

    if(func_name == EGG || func_name == ALL) {
        // Eggholder function
        test_func.set_func(EGG, ndims);
        lower_bounds[0] = -512;
        lower_bounds[1] = -512;
        upper_bounds[0] = 512;
        upper_bounds[1] = 512;
        std::cout << "Running " << EGG << std::endl;
        result = sampler.Minimize(test_func, lower_bounds, upper_bounds);
        print_result(result);
    }

    if(func_name == GAUSS || func_name == ALL) {
        // Gaussian shell
        test_func.set_func(GAUSS, ndims);
        std::cout << "Running " << GAUSS << std::endl;
        result = sampler.Minimize(test_func, lower_bounds, upper_bounds);
        print_result(result);
    }
}

int main(int argc, char* argv[]) {

    std::string instr;
    std::string testf;
    if(argc > 2) {
        instr = argv[1];
        testf = argv[2];
    } else {
        std::cout << "Please enter an argument which minimizer to use\n";
        std::cout << "Possible arguments are:\n";
        std::cout << MAPSNAME   << "     Using MAPS\n";
        std::cout << POLY       << "     Using PolyChord\n";
        std::cout << SAMPLE     << "   Using SampleSpace to sample the function\n";
        std::cout << MULTI      << "    Using MultiNest\n";
        std::cout << "also enter which test function to use\n";
        std::cout << "Possible arguments are:\n";
        std::cout << HIMMEL << std::endl;
        std::cout << TOWN << std::endl;
        std::cout << ROSEN << std::endl;
        std::cout << EGG << std::endl;
        std::cout << GAUSS << std::endl;
        std::cout << ALL << std::endl;
        return 1;
    }

    if(instr == MAPSNAME) {

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
        run_tests(sampler, testf);
    }

    if(instr == POLY) {

        /////// PolyChord
        double tolerance = 1e-4;

        PolyChord sampler(tolerance);
        std::string path = "../../output/PolyChord/";
        sampler.set_output(path);
        run_tests(sampler, testf);
    }

    if(instr == MULTI) {

        /////// MultiNest
        double tolerance = 1e-4;

        MultiNest sampler(tolerance);
        std::string path = "../../output/MultiNest/";
        sampler.set_output(path);
        run_tests(sampler, testf);
    }

    if(instr == SAMPLE) {

        /////// SampleSpace
        // Parameters for sampling the space
        uint32_t n_iterations = 1000000;
        uint32_t max_points = n_iterations/10;
        uint32_t seed = 1025;

        SampleSpace sampler(n_iterations, max_points, seed, true);
        std::string path = "../../output/SampleSpace/";
        sampler.set_output(path);
        run_tests(sampler, testf);
    }

    return 0;
}
