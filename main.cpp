#include <boost/cstdint.hpp>
#include <iostream>

#include "Minimizer/SampleSpace.h"
#include "likelihood/TestFunctions.h"
#include "Minimizer/MinimizerResult.h"


int main(int argc, char* argv[]) {
    
    TestFunctions gauss_shell("egg", 2);
    SampleSpace sampler(100000000, 1000000000, 1025);
    v_d lower_bounds(2);
    v_d upper_bounds(2);
    lower_bounds[0] = -512;
    lower_bounds[1] = -512;
    upper_bounds[0] = 512;
    upper_bounds[1] = 512;
    MinimizerResult result = sampler.Minimize(gauss_shell, lower_bounds, 
        upper_bounds);
    std::cout << std::endl << "Best fit is " << result.best_fit << std::endl;
    
    return 0;
}