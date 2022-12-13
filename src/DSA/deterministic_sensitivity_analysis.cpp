//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Thomas Saigre <saigre@math.unistra.fr>
//! @date 10 Dec 2022
//! @copyright 2022 Feel++ Consortium
//!
#include <boost/dll.hpp>
#include <boost/algorithm/string/split.hpp>

#include <feel/feelcore/table.hpp>
#include <feel/feelmor/options.hpp>
#include <feel/feelmor/crbplugin_interface.hpp>
#include <feel/feelmor/crbmodeldb.hpp>

#include <iostream>
#include <ctime>
#include <execution>

#if defined(FEELPP_HAS_MONGOCXX )
#include <bsoncxx/json.hpp>

#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>
#include <bsoncxx/builder/basic/document.hpp>
#include <bsoncxx/builder/basic/kvp.hpp>
#include <bsoncxx/builder/basic/array.hpp>
#endif

// #include <omp.h>
#include "../tqdm/tqdm.h"

typedef Feel::ParameterSpaceX::element_type element_t;
typedef std::shared_ptr<Feel::CRBPluginAPI> plugin_ptr_t;
typedef std::shared_ptr<Feel::ParameterSpaceX> parameter_space_ptr_t;


std::shared_ptr<Feel::CRBPluginAPI>
loadPlugin()
{
    using namespace Feel;

    std::string crbmodelName = Environment::expand( soption(_name="crbmodel.name") );
    CRBModelDB crbmodelDB{ crbmodelName, uuids::nil_uuid() };

    std::string attribute = soption(_name="crbmodel.attribute" );
    std::string attribute_data;
    if ( attribute == "id"  || attribute == "name")
    {
        attribute_data = Environment::expand( soption(_name=fmt::format("crbmodel.db.{}",attribute) ) );
    }
    else if ( attribute == "last_created" || attribute == "last_modified" )
    {
        std::vector<std::string> split_;
        boost::split(split_, attribute, boost::is_any_of("_"));
        attribute_data = split_[1];
    }
    else
    {
        throw std::runtime_error( "no crbmodel selection, crbmodel.db.id or crbmodel.db.last should be defined" );
    }
    auto meta = crbmodelDB.loadDBMetaData( attribute, attribute_data );
    Feel::cout << "-- crbmodelDB::dbRepository()=" << crbmodelDB.dbRepository() << std::endl;

    return crbmodelDB.loadDBPlugin( meta, soption(_name="crbmodel.db.load" ) );
}



inline Feel::AboutData makeAbout()
{
    Feel::AboutData about( "deterministic_sensitivity_analysis",
                     "DSA" ,
                     "0.1",
                     "Deterministic sensitivity analysis",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2022 Feel++ Consortium" );

    about.addAuthor( "Thomas Saigre", "developer", "saigre@math.unistra.fr", "" );
    return about;
}



std::vector<double> linspace(double start, double end, size_t num_in)
{
    std::vector<double> linspaced;
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1) 
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);
    for(int i=0; i < num-1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end);
    return linspaced;
}


#if 0
OT::Sample output(OT::Sample const& input, plugin_ptr_t const& plugin, Eigen::VectorXd &time_crb, double online_tol, int rbDim)
{
    size_t n = input.getSize();
    OT::Sample output(n, 1);
    {
        parameter_space_ptr_t Dmu = plugin->parameterSpace();
        std::vector<std::string> names = Dmu->parameterNames();
        Feel::cout << "Start to compute outputs" << std::endl;
#if 0
        auto exec_rb = [&input, &Dmu, &time_crb, rbDim, &plugin, online_tol, &output] (int i) {
            OT::Point X = input[i];
            element_t mu = Dmu->element();
            for (size_t j = 0; j < Dmu->dimension(); ++j)
            {
                mu.setParameter(j, X[j]);
            }
            Feel::CRBResults crbResult = plugin->run( mu, time_crb, online_tol, rbDim, false );
            output[i] = OT::Point({boost::get<0>( crbResult )[0]});
        };

        std::vector<int> r(n);
        std::for_each( std::execution::par, r.begin(), r.end(), exec_rb );
#else
        for (size_t i: tqdm::range(n))          // std::for_each
        {
            element_t mu = Dmu->element();
            OT::Point X = input[i];
            for (size_t j = 0; j < Dmu->dimension(); ++j)
            {
                mu.setParameter(j, X[j]);
            }
            Feel::CRBResults crbResult = plugin->run( mu, time_crb, online_tol, rbDim, false );
            output[i] = OT::Point({boost::get<0>( crbResult )[0]});
        }
#endif
        Feel::cout << "output computed" << std::endl;
    }
    return output;
}
#endif

void print_vector(std::vector<double> vec)
{
  std::cout << "size: " << vec.size() << std::endl;
  for (double d : vec)
    std::cout << d << " ";
  std::cout << std::endl;
}

void print_results_to_file(std::vector<double> param_vect, std::vector<double> output_vect, std::string filename)
{
    std::ofstream file;
    file << std::scientific;
    file.open(filename);
    file << "param,output" << std::endl;
    for (size_t i = 0; i < param_vect.size(); ++i)
    {
        file << param_vect[i] << "," << output_vect[i] << std::endl;
    }
    file.close();
}



/**
 * @brief Compute sobol indices
 *
 * @param plugin std::vector containing the plugin from load_plugin
 * @param sampling_size size of the input sample used for computation of sobol indices
 * @param computeSecondOrder boolean to compute second order sobol indices
 */
int runSensitivityAnalysis( std::vector<plugin_ptr_t> plugin, size_t sampling_size )
{
    using namespace Feel;

    bool loadFiniteElementDatabase = boption(_name="crb.load-elements-database");
    std::cout << __LINE__ << std::endl;

    Eigen::VectorXd/*typename crb_type::vectorN_type*/ time_crb;
    double online_tol = 1e-2;               //Feel::doption(Feel::_name="crb.online-tolerance");
    int rbDim = ioption(_name="rb-dim");
    bool print_rb_matrix = false;           //boption(_name="crb.print-rb-matrix");
    parameter_space_ptr_t muspace = plugin[0]->parameterSpace();

    std::vector<std::string> tableRowHeader = muspace->parameterNames();
    size_t dim = muspace->dimension();
    std::cout << __LINE__ << std::endl;

    std::string param = soption( _name="parameter.name" );
    std::map<std::string, double> param_map = {
        {"h_bl", 65},         // [W / m^2 / K]
        {"h_amb", 10},        // [W / m^2 / K]
        {"h_tau", 6},         // [W / m^2 / K]
        {"T_bl", 310},        // [K] 36.85°C
        {"T_amb", 294},       // [K] 24.85°C
        {"E", 40},            // [W / m^2]
    };

    std::vector<double> results(sampling_size);
    element_t mu = muspace->element();
    mu.setParameters(param_map);
    std::cout << __LINE__ << std::endl;

    element_t mu_min = muspace->min();
    element_t mu_max = muspace->max();

    std::cout << "mu_min = " << mu_min << std::endl;
    std::cout << "mu_max = " << mu_max << std::endl;
    std::cout << __LINE__ << std::endl;

    double min_from_option = doption( _name="parameter.min" ),
           max_from_option = doption( _name="parameter.max" );
    double min_value = (min_from_option == DBL_MAX) ? mu_min.parameterNamed(param) : min_from_option,
           max_value = (max_from_option == DBL_MIN) ? mu_max.parameterNamed(param) : max_from_option;

    if (min_value > max_value)
    {
        Feel::cout << tc::red << "Error: min value is greater than max value" << std::endl;
        return 1;
    }
    std::cout << __LINE__ << std::endl;

    Feel::cout << tc::cyan << "Running deterministic analysis for parameter " << param
        << " in [" << min_value << ", " << max_value << "] of size " << sampling_size << tc::reset << std::endl;

    std::vector<double> params_vect = linspace(min_value, max_value, sampling_size);
    std::cout << __LINE__ << std::endl;

    for (size_t i: tqdm::range(sampling_size))
    // for (size_t i = 0; i < sampling_size; ++i)
    {
        mu.setParameterNamed(param, params_vect[i]);
        Feel::CRBResults crbResult = plugin[0]->run( mu, time_crb, online_tol, rbDim, print_rb_matrix );
        results[i] = boost::get<0>( crbResult )[0];
    }
    std::cout << __LINE__ << std::endl;

    print_results_to_file(params_vect, results, "deterministic_analysis_" + param + ".csv");

    return 0;
}

int main( int argc, char** argv )
{
    using namespace Feel;
    po::options_description crbonlinerunoptions( "crb online run options" );
    crbonlinerunoptions.add_options()
        ( "plugin.dir", po::value<std::string>()->default_value( Info::libdir().string() ) , "plugin directory" )

        ( "crbmodel.name", po::value<std::string>(), "CRB online code name" )
        ( "crbmodel.attribute", po::value<std::string>()->default_value( "last_modified" ), "last_created, last_modified, id, name" )
        ( "crbmodel.db.id", po::value<std::string>(), "CRB online code id" )
        ( "crbmodel.db.name", po::value<std::string>(), "CRB online code name" )
        ( "crbmodel.db.last", po::value<std::string>()->default_value( "modified" ), "use created or modified" )
        ( "crbmodel.db.load", po::value<std::string>()->default_value( "rb" ), "load rb, fe or all (fe and rb)" )
        ( "crbmodel.db.root_directory", po::value<std::string>()->default_value( "${repository}/crbdb" ), "root directory of the CRB database " )

        ( "sampling.size", po::value<int>()->default_value( 2000 ), "size of sampling" )
        ( "sampling.type", po::value<std::string>()->default_value( "random" ), "type of sampling" )
        ( "rb-dim", po::value<int>()->default_value( -1 ), "reduced basis dimension used (-1 use the max dim)" )
        ( "output_results.save.path", po::value<std::string>(), "output_results.save.path" )


        ( "parameter.name", po::value<std::string>(), "selected parameter")
        ( "parameter.min", po::value<double>()->default_value( DBL_MAX ), "minimal value taken by the varying parameter" )
        ( "parameter.max", po::value<double>()->default_value( DBL_MIN ), "maximal value taken by the varying parameter")

        ( "query", po::value<std::string>(), "query string for mongodb DB feelpp.crbdb" )
        ( "compare", po::value<std::string>(), "compare results from query in mongodb DB feelpp.crbdb" )
        ( "list", "list registered DB in mongoDB  in feelpp.crbdb" )
        ;
    po::options_description crbonlinerunliboptions( "crb online run lib options" );
#if 1
    crbonlinerunliboptions.add(crbOptions())
        .add(crbSEROptions())
        .add(eimOptions())
        .add(podOptions())
        .add(backend_options("backend-primal"))
        .add(backend_options("backend-dual"))
        .add(backend_options("backend-l2"))
        ;
#endif

    Environment env( _argc = argc, _argv = argv,
                     _desc = crbonlinerunoptions,
                     _desc_lib = crbonlinerunliboptions.add( feel_options() ),
                     _about = makeAbout() );

    plugin_ptr_t plugin = loadPlugin();
    // runCrbOnline( { plugin } );
    int err = runSensitivityAnalysis( { plugin }, ioption(_name="sampling.size") );

    if (err == 0)
        Feel::cout << tc::green << "Done ✓" << tc::reset << std::endl;
    else
        Feel::cout << tc::red << "Error ✗" << tc::reset << std::endl;
    return err;
}
