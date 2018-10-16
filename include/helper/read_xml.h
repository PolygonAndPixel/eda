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
#include "Minimizer/DalexMinimizer.h"
#include "Minimizer/GradientDescent.h"
#include "Minimizer/MAPS.h"
#include "Minimizer/Minimizer.h"
#include "Minimizer/MinimizerResult.h"
#include "Minimizer/MultiNest.h"
#include "Minimizer/PolyChord.h"
#include "Minimizer/SampleSpace.h"
#include "Minimizer/Scan.h"
#include "likelihood/TestFunctions.h"

#include <memory>

/** Load a configuration file for a minimizer, create it with the configuration
 *  and return the object.
 *
 *  \param xml_file     The path and name to the file
 *  \param minimizers   Vector of minimizer where the new minimizers shall be
 *                      stored.
 *
 *  \return             The configured minimizer
 * */
void load_minimizer_xml(
    std::string &xml_file,
    std::vector<std::unique_ptr<Minimizer>> &minimizers) {

    boost::property_tree::ptree ptree;
    boost::property_tree::xml_parser::read_xml(xml_file, ptree);

    bool first = true;
    for(boost::property_tree::ptree::value_type &v: ptree.get_child("Minimizers")) {
        boost::property_tree::ptree pt = v.second;
        // The first entry is just about the encoding
        if(first) {first = false; continue;}

        std::string name = pt.get<std::string>("name");
        // Load specific values depending on the method
        if(name == MAPSNAME) {
            value_t tolerance    = pt.get<value_t>("tolerance");
            index_t max_iter        = pt.get<index_t>("max_iter");
            index_t min_iter        = pt.get<index_t>("min_iter");
            index_t max_points      = pt.get<index_t>("max_points");
            index_t seed            = pt.get<index_t>("seed");
            bool dump_points    = pt.get<bool>("dump_points");

            index_t n_start_points      = pt.get<index_t>("n_start_points");
            index_t size_sub_pop        = pt.get<index_t>("size_sub_pop");
            index_t max_sub_pops        = pt.get<index_t>("max_sub_pops");
            index_t n_selected          = pt.get<index_t>("n_selected");
            index_t n_sub_selected      = pt.get<index_t>("n_sub_selected");
            value_t size_factor      = pt.get<value_t>("size_factor");
            MAPS minimizer(tolerance, max_iter, min_iter, max_points,
                n_start_points, size_sub_pop, max_sub_pops, n_selected, size_factor,
                seed, dump_points);
            std::string path = "../../output/MAPS/";
            if(dump_points) minimizer.set_output(path);
            minimizers.push_back(minimizer.clone());
            continue;
        }

        if(name == POLY) {
            value_t tolerance    = pt.get<value_t>("tolerance");
            index_t max_iter        = pt.get<index_t>("max_iter");
            index_t min_iter        = pt.get<index_t>("min_iter");
            index_t max_points      = pt.get<index_t>("max_points");
            index_t seed            = pt.get<index_t>("seed");
            bool dump_points    = pt.get<bool>("dump_points");

            index_t n_prior             = pt.get<index_t>("n_prior");
            index_t n_grade             = pt.get<index_t>("n_grade");
            std::string values      = pt.get<std::string>("grade_frac");
            v_d grade_frac;
            std::vector<std::string> values_split;
            std::string::size_type sz;
            boost::split(values_split, values, boost::is_any_of(" "));
            for(std::string &item: values_split)
                grade_frac.push_back(std::stod(item));
            index_t n_live              = pt.get<index_t>("n_live");
            index_t feedback            = pt.get<index_t>("feedback");
            index_t max_dead            = pt.get<index_t>("max_dead");
            value_t boost_posterior  = pt.get<value_t>("boost_posterior");
            index_t num_repeats         = pt.get<index_t>("num_repeats");
            bool posteriors         = pt.get<bool>("posteriors");
            bool equals             = pt.get<bool>("equals");
            bool cluster_posteriors = pt.get<bool>("cluster_posteriors");
            bool do_clustering      = pt.get<bool>("do_clustering");
            PolyChord minimizer(tolerance, max_iter, min_iter, max_points,
                n_prior, n_grade, grade_frac.data(), n_live, feedback, max_dead,
                boost_posterior, num_repeats, posteriors, equals,
                cluster_posteriors, do_clustering, seed, dump_points);

            std::string path = "../../output/PolyChord/";
            if(dump_points) minimizer.set_output(path);
            minimizers.push_back(minimizer.clone());
            continue;
        }

        if(name == MULTI) {
            value_t tolerance    = pt.get<value_t>("tolerance");
            index_t max_iter        = pt.get<index_t>("max_iter");
            index_t seed            = pt.get<index_t>("seed");
            bool dump_points    = pt.get<bool>("dump_points");

            bool ins                 = pt.get<bool>("ins");
            bool mode_separation     = pt.get<bool>("mode_separation");
            bool const_eff           = pt.get<bool>("const_eff");
            index_t n_live               = pt.get<index_t>("n_live");
            value_t enlargement       = pt.get<value_t>("enlargement");
            index_t feedback_interval    = pt.get<index_t>("feedback_interval");
            index_t max_modes            = pt.get<index_t>("max_modes");
            bool feedback            = pt.get<bool>("feedback");

            MultiNest minimizer(tolerance, max_iter, ins, mode_separation,
                const_eff, n_live, enlargement, feedback_interval, max_modes,
                feedback, seed, dump_points);

            std::string path = "../../output/MultiNest/";
            if(dump_points) minimizer.set_output(path);
            minimizers.push_back(minimizer.clone());
            continue;
        }

        if(name == SAMPLE) {
            index_t max_iter        = pt.get<index_t>("max_iter");
            index_t max_points      = pt.get<index_t>("max_points");
            index_t seed            = pt.get<index_t>("seed");
            bool dump_points    = pt.get<bool>("dump_points");

            SampleSpace minimizer(max_iter, max_points, seed, dump_points);

            std::string path = "../../output/SampleSpace/";
            if(dump_points) minimizer.set_output(path);
            minimizers.push_back(minimizer.clone());
            continue;
        }

        if(name == SCAN) {
            index_t n_points_per_dim = pt.get<index_t>("n_points_per_dim");
            index_t max_points       = pt.get<index_t>("max_points");
            index_t seed;
            // Those have default values 
            try{
                 seed = pt.get<index_t>("seed");
            } catch(const boost::property_tree::ptree_bad_path &e) {
                // seed = NULL;
            }
            bool dump_points     = pt.get<bool>("dump_points");

            Scan minimizer(n_points_per_dim, max_points, 
                seed, dump_points);

            std::string path = "../../output/Scan/";
            if(dump_points) minimizer.set_output(path);
            minimizers.push_back(minimizer.clone());
            continue;
        }

        if(name == DALEX) {
            index_t max_iter        = pt.get<index_t>("max_iter");
            index_t max_points       = pt.get<index_t>("max_points");
            index_t seed;
            // Those have default values 
            try{
                 seed = pt.get<index_t>("seed");
            } catch(const boost::property_tree::ptree_bad_path &e) { }
            bool dump_points;
            try{
                 bool dump_points    = pt.get<bool>("dump_points");
            } catch(const boost::property_tree::ptree_bad_path &e) { }

            DalexMinimizer minimizer(max_iter, max_points, 
                seed, dump_points);

            std::string path = "../../output/Dalex/";
            if(dump_points) minimizer.set_output(path);
            minimizers.push_back(minimizer.clone());
            continue;
        }

        if(name == GD) {
            index_t max_iter        = pt.get<index_t>("max_iter");
            index_t min_iter        = pt.get<index_t>("min_iter");
            index_t max_points       = pt.get<index_t>("max_points");

            // Those have default values 
            value_t conv_crit;
            try{
                conv_crit = pt.get<value_t>("conv_crit");
            } catch(const boost::property_tree::ptree_bad_path &e) { }
            index_t n_gradients;
            try{
                 n_gradients = pt.get<index_t>("n_gradients");
            } catch(const boost::property_tree::ptree_bad_path &e) { }
            index_t seed;
            try{
                 seed = pt.get<index_t>("seed");
            } catch(const boost::property_tree::ptree_bad_path &e) { }
            bool dump_points;
            try{
                 bool dump_points    = pt.get<bool>("dump_points");
            } catch(const boost::property_tree::ptree_bad_path &e) { }

            value_t stepsize;
            try{
                 stepsize = pt.get<value_t>("stepsize");
            } catch(const boost::property_tree::ptree_bad_path &e) { }

            GradientDescent minimizer(max_iter, min_iter, max_points, 
                n_gradients, seed, stepsize, conv_crit, dump_points);

            std::string path = "../../output/GD/";
            if(dump_points) minimizer.set_output(path);
            minimizers.push_back(minimizer.clone());
            continue;
        }

        std::cerr << "No such minimizer! "
            << "Did you specify the right name in your configuration file?"
            << std::endl << "Falling back to SampleSpace." << std::endl;
        SampleSpace minimizer(10000, 1000, 1025, false);
        minimizers.push_back(minimizer.clone());
    }
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

        index_t n_dims = upper.size();
        test_func.set_func(func_name, n_dims);
        test_funcs.push_back(test_func);
        lower_bounds.push_back(lower);
        upper_bounds.push_back(upper);
    }
}

#endif
