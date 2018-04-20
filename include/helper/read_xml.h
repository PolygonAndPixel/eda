#ifndef INI_H_INCLUDED
#define INI_H_INCLUDED

/* Just some methods to read an xml file.
 *
 * Author: Maicon Hieronymus <mhierony@students.uni-mainz.de>
 * */
// Boost ptree and lapacke have conflicts if included in the wrong order.
 #include <boost/property_tree/ptree.hpp>
 #include <boost/property_tree/xml_parser.hpp>
 #include <boost/algorithm/string.hpp>

#include "helper/abbreviations.h"
#include "Minimizer/Minimizer.h"
#include "Minimizer/SampleSpace.h"
#include "Minimizer/MAPS.h"
#include "Minimizer/PolyChord.h"
#include "Minimizer/Minimizer.h"
#include "Minimizer/MultiNest.h"
#include "likelihood/TestFunctions.h"
#include "Minimizer/MinimizerResult.h"

#include <memory>

/** Load a configuration file for a minimizer, create it with the configuration
 *  and return the object.
 *
 *  \param xml_file     The path and name to the file
 *
 *  \return             The configured minimizer
 * */
std::unique_ptr<Minimizer> load_minimizer_xml(
    std::string &xml_file) {

    boost::property_tree::ptree pt;
    boost::property_tree::xml_parser::read_xml(xml_file, pt);
    std::string name = pt.get<std::string>("Minimizer.name");

    // Load specific values depending on the method
    if(name == MAPSNAME) {
        double tolerance    = pt.get<double>("Minimizer.tolerance");
        int max_iter        = pt.get<int>("Minimizer.max_iter");
        int min_iter        = pt.get<int>("Minimizer.min_iter");
        int max_points      = pt.get<int>("Minimizer.max_points");
        int seed            = pt.get<int>("Minimizer.seed");
        bool dump_points    = pt.get<bool>("Minimizer.dump_points");

        int n_start_points      = pt.get<int>("Minimizer.n_start_points");
        int size_sub_pop        = pt.get<int>("Minimizer.size_sub_pop");
        int max_sub_pops        = pt.get<int>("Minimizer.max_sub_pops");
        int n_selected          = pt.get<int>("Minimizer.n_selected");
        int n_sub_selected      = pt.get<int>("Minimizer.n_sub_selected");
        double size_factor      = pt.get<double>("Minimizer.size_factor");
        MAPS minimizer(tolerance, max_iter, min_iter, max_points,
            n_start_points, size_sub_pop, max_sub_pops, n_selected, size_factor,
            seed, dump_points);
        std::string path = "../../output/MAPS/";
        if(dump_points) minimizer.set_output(path);
        return minimizer.clone();
    }

    if(name == POLY) {
        double tolerance    = pt.get<double>("Minimizer.tolerance");
        int max_iter        = pt.get<int>("Minimizer.max_iter");
        int min_iter        = pt.get<int>("Minimizer.min_iter");
        int max_points      = pt.get<int>("Minimizer.max_points");
        int seed            = pt.get<int>("Minimizer.seed");
        bool dump_points    = pt.get<bool>("Minimizer.dump_points");

        int n_prior             = pt.get<int>("Minimizer.n_prior");
        int n_grade             = pt.get<int>("Minimizer.n_grade");
        std::string values      = pt.get<std::string>("Minimizer.grade_frac");
        v_d grade_frac;
        std::vector<std::string> values_split;
        std::string::size_type sz;
        boost::split(values_split, values, boost::is_any_of(" "));
        for(std::string &item: values_split)
            grade_frac.push_back(std::stod(item));
        int n_live              = pt.get<int>("Minimizer.n_live");
        int feedback            = pt.get<int>("Minimizer.feedback");
        int max_dead            = pt.get<int>("Minimizer.max_dead");
        double boost_posterior  = pt.get<double>("Minimizer.boost_posterior");
        int num_repeats         = pt.get<int>("Minimizer.num_repeats");
        bool posteriors         = pt.get<bool>("Minimizer.posteriors");
        bool equals             = pt.get<bool>("Minimizer.equals");
        bool cluster_posteriors = pt.get<bool>("Minimizer.cluster_posteriors");
        bool do_clustering      = pt.get<bool>("Minimizer.do_clustering");

        PolyChord minimizer(tolerance, max_iter, min_iter, max_points,
            n_prior, n_grade, grade_frac.data(), n_live, feedback, max_dead,
            boost_posterior, num_repeats, posteriors, equals,
            cluster_posteriors, do_clustering, seed, dump_points);

        std::string path = "../../output/PolyChord/";
        if(dump_points) minimizer.set_output(path);
        return minimizer.clone();
    }

    if(name == MULTI) {
        double tolerance    = pt.get<double>("Minimizer.tolerance");
        int max_iter        = pt.get<int>("Minimizer.max_iter");
        int max_points      = pt.get<int>("Minimizer.max_points");
        int seed            = pt.get<int>("Minimizer.seed");
        bool dump_points    = pt.get<bool>("Minimizer.dump_points");

        bool ins                = pt.get<bool>("Minimizer.ins");
        bool mode_separation    = pt.get<bool>("Minimizer.mode_separation");
        bool const_eff          = pt.get<bool>("Minimizer.const_eff");
        int n_live              = pt.get<int>("Minimizer.n_live");
        double enlargement       = pt.get<double>("Minimizer.enlargement");
        int feedback_interval   = pt.get<int>("Minimizer.feedback_interval");
        int max_modes           = pt.get<int>("Minimizer.max_modes");
        bool feedback           = pt.get<bool>("Minimizer.feedback");

        MultiNest minimizer(tolerance, max_iter, ins, mode_separation,
            const_eff, n_live, enlargement, feedback_interval, max_modes,
            feedback, seed, dump_points);

        std::string path = "../../output/MultiNest/";
        if(dump_points) minimizer.set_output(path);
        return minimizer.clone();
    }

    if(name == SAMPLE) {
        int max_iter        = pt.get<int>("Minimizer.max_iter");
        int max_points      = pt.get<int>("Minimizer.max_points");
        int seed            = pt.get<int>("Minimizer.seed");
        bool dump_points    = pt.get<bool>("Minimizer.dump_points");

        SampleSpace minimizer(max_iter, max_points, seed, dump_points);

        std::string path = "../../output/SampleSpace/";
        if(dump_points) minimizer.set_output(path);
        return minimizer.clone();
    }

    std::cerr << "No such minimizer! "
        << "Did you specify the right name in your configuration file?"
        << std::endl << "Falling back to SampleSpace." << std::endl;
    SampleSpace minimizer(10000, 1000, 1025, false);
    return minimizer.clone();
}

/** Load a configuration file for a minimizer, create it with the configuration
 *  and return the object.
 *
 *  \param xml_file         The path and name to the file
 *  \param test_funcs       A vector of TestFunctions object that shall
 *                          be initialized
 *  \param lower_bounds     A vector of vectors that will be used for the
 *                          lower bounds.
 *  \param upper_bounds     A vector of vectors that will be used for the
 *                          upper bounds
 *
 *  \return             The configured minimizer
 * */
void load_likelihood_xml(
    std::string &xml_file,
    std::vector<TestFunctions> &test_funcs,
    m_d &lower_bounds,
    m_d &upper_bounds) {

    lower_bounds.clear();
    upper_bounds.clear();
    test_funcs.clear();
    boost::property_tree::ptree pt;
    boost::property_tree::xml_parser::read_xml(xml_file, pt);
    bool first = true;
    for(boost::property_tree::ptree::value_type &v: pt.get_child("likelihoods")) {
        boost::property_tree::ptree likelihood = v.second;
        // The first entry is just about the encoding
        if(first) {first = false; continue;}
        TestFunctions test_func;
        std::string func_name = likelihood.get<std::string>("function_name");

        v_d lower;
        std::string values = likelihood.get<std::string>("lower_bounds");
        std::vector<std::string> values_split;
        boost::split(values_split, values, boost::is_any_of(" "));
        for(std::string &item: values_split) {
            lower.push_back(std::stod(item));
        }
        v_d upper;
        values = likelihood.get<std::string>("upper_bounds");
        boost::split(values_split, values, boost::is_any_of(" "));
        for(std::string &item: values_split)
            upper.push_back(std::stod(item));

        int n_dims = upper.size();
        test_func.set_func(func_name, n_dims);
        test_funcs.push_back(test_func);
        lower_bounds.push_back(lower);
        upper_bounds.push_back(upper);
    }
}

#endif
