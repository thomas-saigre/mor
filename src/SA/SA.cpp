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
//! @date 25 Jul 2022
//! @copyright 2022 Feel++ Consortium
//!
#include <boost/dll.hpp>
#include <boost/algorithm/string/split.hpp>
#include <openturns/OT.hxx>

#include <feel/feelcore/table.hpp>
#include <feel/feelmor/options.hpp>
#include <feel/feelmor/crbplugin_interface.hpp>
#include <feel/feelmor/crbmodeldb.hpp>

#include <iostream>
#include <ctime>
#include <fstream>
#include <algorithm>
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
    Feel::AboutData about( "sensitibity_analysis",
                     "SA" ,
                     "0.1",
                     "Sensitivity analysis",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2022 Feel++ Consortium" );

    about.addAuthor( "Thomas Saigre", "developer", "saigre@math.unistra.fr", "" );
    return about;
}

OT::ComposedDistribution composedFromModel(parameter_space_ptr_t Dmu )
{
    using namespace Feel;
    element_t mumin = Dmu->min();
    element_t mumax = Dmu->max();
    std::vector<std::string> names = Dmu->parameterNames();
    Feel::cout << tc::cyan << "names = " << names << tc::reset << std::endl;

    OT::Collection<OT::Distribution> marginals(Dmu->dimension());

    for ( uint16_type d=0; d<Dmu->dimension(); ++d)
    {
        OT::Distribution dist;
        if (names[d] == "h_amb")
            dist = OT::LogNormal(2., 0.6137056388801094, 8.);
        else
            dist = OT::Uniform( mumin(d), mumax(d) );
        dist.setDescription({names[d]});
        marginals[d] = dist;
        Feel::cout << tc::blue << "Distribution " << d << " (" << names[d] << ") = " << dist << tc::reset << std::endl;
    }
    
    return OT::ComposedDistribution( marginals );
}

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

/**
 * @brief Export values computed in a json file
 * 
 * @param filename path to the exported file
 * @param names names of the parameters
 * @param algo algo used to compute sobol indices
 * @param size size of the sample used for computation
 * @param firstOrder first order sobol indices
 * @param totalOrder total order sobol indices
 * @param firstOrderMin lower bound of the first order sobol indices
 * @param firstOrderMax upper bound of the first order sobol indices
 * @param totalOrderMin lower bound of the total order sobol indices
 * @param totalOrderMax upper bound of the total order sobol indices
 */
void exportValues(const std::string filename,
                  const std::vector<std::string>& names, const std::string algo, size_t size,
                  const std::vector<OT::Scalar>& firstOrder, const std::vector<OT::Scalar>& totalOrder,
                  const std::vector<OT::Scalar>& firstOrderMin, const std::vector<OT::Scalar>& firstOrderMax,
                  const std::vector<OT::Scalar>& totalOrderMin, const std::vector<OT::Scalar>& totalOrderMax)
{
    std::ofstream file;
    file.open(filename);
    size_t n = names.size();
    file << "{\n\t\"N\": " << n << ",\n";
    file << "\t\"sampling-size\": " << size << ",\n";
    file << "\t\"algo\": \"" << algo << "\",\n";
    file << "\t\"Names\": [";
    for (size_t i = 0; i < n; ++i)
    {
        file << "\"" << names[i] << "\"";
        if (i != n - 1)
            file << " ,";
    }
    file << "]," << std::endl;
    file << "\t\"FirstOrder\":\n\t{" << std::endl;
    file << "\t\t\"values\": [";
    assert(n == firstOrder.size());
    for (size_t i = 0; i < n; ++i)
    {
        file << firstOrder[i];
        if (i != n - 1)
        {
            file << ", ";
        }
    }
    file << "]," << std::endl;
    file << "\t\t\"intervals\": [";
    for (size_t i = 0; i < n; ++i)
    {
        file << "[" << firstOrderMin[i] << ", " << firstOrderMax[i] << "]";
        if (i != n - 1)
        {
            file << ", ";
        }
    }
    file << "]" << std::endl;
    file << "\t}," << std::endl;
    file << "\t\"TotalOrder\":\n\t{" << std::endl;
    file << "\t\t\"values\": [";
    assert(n == totalOrder.size());
    for (size_t i = 0; i < n; ++i)
    {
        file << totalOrder[i];
        if (i != n - 1)
        {
            file << ", ";
        }
    }
    file << "]," << std::endl;
    file << "\t\t\"intervals\": [";
    for (size_t i = 0; i < n; ++i)
    {
        file << "[" << totalOrderMin[i] << ", " << totalOrderMax[i] << "]";
        if (i != n - 1)
        {
            file << ", ";
        }
    }
    file << "]" << std::endl;
    file << "\t}" << std::endl;
    file << "}" << std::endl;
    file.close();
}


/**
 * @brief Compute sobol indices
 * 
 * @param plugin std::vector containing the plugin from load_plugin
 * @param sampling_size size of the input sample used for computation of sobol indices
 * @param computeSecondOrder boolean to compute second order sobol indices
 */
void runSensitivityAnalysis( std::vector<plugin_ptr_t> plugin, size_t sampling_size, bool computeSecondOrder=true )
{
    using namespace Feel;

    Feel::cout << "Running sensisivity analysis with a sample of size " << sampling_size << std::endl;

    bool loadFiniteElementDatabase = boption(_name="crb.load-elements-database");

    Eigen::VectorXd/*typename crb_type::vectorN_type*/ time_crb;
    double online_tol = 1e-2;               //Feel::doption(Feel::_name="crb.online-tolerance");
    int rbDim = ioption(_name="rb-dim");
    bool print_rb_matrix = false;           //boption(_name="crb.print-rb-matrix");
    parameter_space_ptr_t muspace = plugin[0]->parameterSpace();

    OT::ComposedDistribution composed_distribution = composedFromModel( muspace );
    std::vector<std::string> tableRowHeader = muspace->parameterNames();
    size_t dim = muspace->dimension();


    if ( !boption("algo.poly") )
    {
        OT::SobolIndicesExperiment sobol(composed_distribution, sampling_size, computeSecondOrder);
        tic();
        OT::Sample inputDesign = sobol.generate();
        toc("input design");
        Feel::cout << "inputDesign generated" << std::endl;
        tic();
        OT::Sample outputDesign = output(inputDesign, plugin[0], time_crb, online_tol, rbDim);
        toc("output design");

        OT::SaltelliSensitivityAlgorithm sensitivity(inputDesign, outputDesign, sampling_size);
        sensitivity.setUseAsymptoticDistribution( true );

        OT::Point firstOrder = sensitivity.getFirstOrderIndices();
        OT::Interval intervals = sensitivity.getFirstOrderIndicesInterval();
        OT::Point totalOrder = sensitivity.getTotalOrderIndices();
        OT::Interval totalIntervals = sensitivity.getTotalOrderIndicesInterval();

        Feel::cout << "Parameter names: " << tableRowHeader << std::endl;
        Feel::cout << "First order indices: " << firstOrder<< std::endl;
        Feel::cout << "Intervals : " << intervals << std::endl;
        Feel::cout << "Total order indices: " << totalOrder << std::endl;
        Feel::cout << "Intervals : " << totalIntervals << std::endl;
    }
    else    // if algo.poly
    {
        int nrun = ioption(_name="algo.nrun");
        std::vector<OT::Scalar> firstOrder(dim), firstOrderMin(dim), firstOrderMax(dim),
                                totalOrder(dim), totalOrderMin(dim), totalOrderMax(dim);
        std::fill(firstOrder.begin(), firstOrder.end(), 0.0);
        std::fill(totalOrder.begin(), totalOrder.end(), 0.0);
        std::fill(firstOrderMin.begin(), firstOrderMin.end(), 1.0);
        std::fill(totalOrderMin.begin(), totalOrderMin.end(), 1.0);
        std::fill(firstOrderMax.begin(), firstOrderMax.end(), 0.0);
        std::fill(totalOrderMax.begin(), totalOrderMax.end(), 0.0);
        OT::Scalar o1, ot;
        for (int r=0; r<nrun; ++r)
        {
            Feel::cout << tc::bold << tc::red << "Run " << r+1 << " over " << nrun << tc::reset << std::endl;
            OT::Sample input_sample = composed_distribution.getSample(sampling_size);
            OT::Sample output_sample = output(input_sample, plugin[0], time_crb, online_tol, rbDim);
            OT::FunctionalChaosAlgorithm polynomialChaosAlgorithm = OT::FunctionalChaosAlgorithm(input_sample, output_sample);

            polynomialChaosAlgorithm.run();
            OT::FunctionalChaosResult polynomialChaosResult = polynomialChaosAlgorithm.getResult();
            OT::FunctionalChaosSobolIndices sensitivityAnalysis = OT::FunctionalChaosSobolIndices(polynomialChaosResult);
            for (size_t i=0; i<dim; ++i)
            {
                o1 = sensitivityAnalysis.getSobolIndex(i);
                ot = sensitivityAnalysis.getSobolTotalIndex(i);
                firstOrder[i] += o1;
                totalOrder[i] += ot;
                if (o1 < firstOrderMin[i]) firstOrderMin[i] = o1;
                if (ot < totalOrderMin[i]) totalOrderMin[i] = ot;
                if (o1 > firstOrderMax[i]) firstOrderMax[i] = o1;
                if (ot > totalOrderMax[i]) totalOrderMax[i] = ot;
            }
        }

        for (size_t i=0; i<dim; ++i)
        {
            firstOrder[i] /= nrun;
            totalOrder[i] /= nrun;
        }

        Feel::cout << "Parameter names: " << tableRowHeader << std::endl;
        Feel::cout << "First order indices: " << firstOrder << std::endl;
        Feel::cout << "Total order indices: " << totalOrder << std::endl;
        Feel::cout << "FirstOrderIntervals" << std::endl;
        for (size_t i=0; i<dim; ++i)
            Feel::cout << "[" << firstOrderMin[i] << ", " << firstOrderMax[i] << "]" << std::endl;
        Feel::cout << "TotalOrderIntervals" << std::endl;
        for (size_t i=0; i<dim; ++i)
            Feel::cout << "[" << totalOrderMin[i] << ", " << totalOrderMax[i] << "]" << std::endl;

        exportValues("sensitivity.json", plugin[0]->parameterSpace()->parameterNames(),
            "polynomial-chaos", sampling_size,
            firstOrder, totalOrder, firstOrderMin, firstOrderMax, totalOrderMin, totalOrderMax);
    } // end if algo.poly
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

        ( "parameter", po::value<std::vector<std::string> >()->multitoken(), "database filename" )
        ( "sampling.size", po::value<int>()->default_value( 2000 ), "size of sampling" )
        ( "sampling.type", po::value<std::string>()->default_value( "random" ), "type of sampling" )
        ( "rb-dim", po::value<int>()->default_value( -1 ), "reduced basis dimension used (-1 use the max dim)" )
        ( "output_results.save.path", po::value<std::string>(), "output_results.save.path" )

        ( "algo.poly", po::value<bool>()->default_value(true), "use polynomial chaos" )
        ( "algo.nrun", po::value<int>()->default_value(5), "number to run algorithm" )

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

    OT::RandomGenerator::SetSeed( ::time(NULL) );
    plugin_ptr_t plugin = loadPlugin();
    // runCrbOnline( { plugin } );
    runSensitivityAnalysis( { plugin }, ioption(_name="sampling.size"), false );

    Feel::cout << tc::green << "Done ✓" << tc::reset << std::endl;
    return 0;
}
