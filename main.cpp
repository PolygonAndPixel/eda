#include <boost/cstdint.hpp>
#include <iostream>

#include "Minimizer/SampleSpace.h"
#include "likelihood/TestFunctions.h"
#include "Minimizer/MinimizerResult.h"


int main(int argc, char* argv[]) {
    
    TestFunctions gauss_shell("gauss", 2);
    SampleSpace sampler(1000, 100, 1025);
    v_d lower_bounds(2);
    v_d upper_bounds(2);
    lower_bounds[0] = -6;
    lower_bounds[1] = -6;
    upper_bounds[0] = 6;
    upper_bounds[1] = 6;
    MinimizerResult result = sampler.Minimize(gauss_shell, lower_bounds, 
        upper_bounds);
    std::cout << "Best fit is " << result.best_fit << std::endl;
    
    return 0;
}