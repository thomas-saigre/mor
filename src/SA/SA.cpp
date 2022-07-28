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

#if defined(FEELPP_HAS_MONGOCXX )
#include <bsoncxx/json.hpp>

#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>
#include <bsoncxx/builder/basic/document.hpp>
#include <bsoncxx/builder/basic/kvp.hpp>
#include <bsoncxx/builder/basic/array.hpp>
#endif

bool runCrbOnline( std::vector<std::shared_ptr<Feel::CRBPluginAPI>> plugin );

bool
runCrbOnline( std::vector<std::shared_ptr<Feel::CRBPluginAPI>> plugin )
{
    using namespace Feel;

    bool loadFiniteElementDatabase = boption(_name="crb.load-elements-database");

    Eigen::VectorXd/*typename crb_type::vectorN_type*/ time_crb;
    double online_tol = 1e-2;//Feel::doption(Feel::_name="crb.online-tolerance");
    bool print_rb_matrix = false;//boption(_name="crb.print-rb-matrix");
    auto muspace = plugin[0]->parameterSpace();

    std::ostringstream ostrmumin,ostrmumax;
    auto mumin=muspace->min();
    auto mumax=muspace->max();
    for ( uint16_type d=0;d<muspace->dimension();++d)
    {
        ostrmumin << mumin(d) << " ";
        ostrmumax << mumax(d) << " ";
    }
    std::cout << "dimension of parameter space : " << muspace->dimension() << "\n";
    std::cout << "min element in parameter space : "<< ostrmumin.str() << "\n";
    std::cout << "max element in parameter space : "<< ostrmumax.str() << "\n";


    auto mysampling = muspace->sampling();

    std::vector<double> inputParameter;
    if ( Environment::vm().count("parameter"))
    {
        auto inputParameterParsed = Environment::vm()["parameter"].as<std::vector<std::string> >();

        if ( inputParameterParsed.size() == 1 )
        {
            std::vector<std::string > stringParsedSplitted;
            boost::split( stringParsedSplitted, inputParameterParsed.front(), boost::is_any_of(" "), boost::token_compress_on );
            inputParameterParsed = stringParsedSplitted;
        }

        for ( std::string const& paramParsed : inputParameterParsed )
            inputParameter.push_back( std::stod(paramParsed) );
    }
    //inputParameter = Environment::vm()["parameter"].as<std::vector<double> >();
    if ( !inputParameter.empty() )
    {
        CHECK( inputParameter.size() == muspace->dimension() ) << "parameter has a wrong size : "<< inputParameter.size() << " but must be " << muspace->dimension() << ":"<<inputParameter;
        auto mu = muspace->element();
        for ( uint16_type d=0;d<muspace->dimension();++d)
            mu(d)=inputParameter[d];
        mysampling->push_back( mu );
    }
    else
    {
        int nSample = ioption(_name="sampling.size");
        std::string sampler = soption("sampling.type");
        mysampling->sample( nSample, sampler );
    }

    if ( loadFiniteElementDatabase )
        plugin[0]->initExporter();

    int rbDim = ioption(_name="rb-dim");
    int nSamples = mysampling->size();

    Feel::Table tableOutputResults;
    std::vector<std::string> tableRowHeader = muspace->parameterNames();
    tableRowHeader.push_back( "output");
    tableOutputResults.add_row( tableRowHeader );
    tableOutputResults.format().setFirstRowIsHeader( true );

    std::vector<double> tableRowValues(tableRowHeader.size());

    for ( int k=0;k<nSamples;++k )
    {
        auto const& mu = (*mysampling)[k];
        std::ostringstream ostrmu;
        for ( uint16_type d=0;d<muspace->dimension();++d)
            ostrmu << mu(d) << " ";
        // std::cout << "--------------------------------------\n";
        // std::cout << "mu["<<k<<"] : " << ostrmu.str() << "\n";
        //auto mu = crb->Dmu()->element();
        //std::cout << "input mu\n" << mu << "\n";
        for( auto const& p : plugin )
        {
            auto crbResult = p->run( mu, time_crb, online_tol, rbDim, print_rb_matrix);
            auto resOuptut = boost::get<0>( crbResult );
            auto resError = boost::get<0>( boost::get<6>( crbResult ) );
            //std::cout << "output " << resOuptut.back() << " " << resError.back() << "\n";

            int curRowValIndex = 0;
            for ( uint16_type d=0;d<muspace->dimension();++d)
                tableRowValues[curRowValIndex++] = mu(d);
            if ( !resOuptut.empty() )
                tableRowValues[curRowValIndex++] = resOuptut.back();
            tableOutputResults.add_row( tableRowValues );

            if ( loadFiniteElementDatabase )
            {
                p->exportField( (boost::format("sol-%1%")%k).str(), crbResult );
            }
        }
    }

    bool printResults = true;
    if ( printResults )
        std::cout << tableOutputResults << std::endl;
    bool saveResults = true;

    std::string outputResultPath = "output.csv";
    if ( Environment::vm().count("output_results.save.path") )
        outputResultPath = soption(_name="output_results.save.path");
    outputResultPath = Environment::expand( outputResultPath );
    if ( !fs::exists( fs::path(outputResultPath).parent_path() ) && !fs::path(outputResultPath).parent_path().empty() )
        fs::create_directories( fs::path(outputResultPath).parent_path() );

    if ( saveResults )
    {
        std::ofstream ofs( outputResultPath );
        tableOutputResults.exportCSV( ofs );
    }
    if ( loadFiniteElementDatabase )
        plugin[0]->saveExporter();

    return true;

}

std::shared_ptr<Feel::CRBPluginAPI>
loadPlugin()
{
    using namespace Feel;
    namespace dll=boost::dll;

    std::string crbmodelName = Environment::expand( soption(_name="crbmodel.name") );
    CRBModelDB crbmodelDB{ crbmodelName, uuids::nil_uuid() };

    std::string pluginname;
    std::string pluginlibname;
    std::string jsonPathStr;

    if ( countoption(_name="crbmodel.db.id") )
    {
        std::string crbmodel_dbid = Environment::expand( soption(_name="crbmodel.db.id") );
        crbmodelDB.updateIdFromId( crbmodel_dbid );
    }
    else if ( ioption(_name="crbmodel.db.last" ) )
    {
        crbmodelDB.updateIdFromDBLast( static_cast<crb::last>(ioption("crbmodel.db.last")) );
    }
    else
    {
        throw std::runtime_error( "no crbmodel selection, crbmodel.db.id or crbmodel.db.last should be defined" );
    }

    std::cout << "crbmodelDB.dbRepository()=" << crbmodelDB.dbRepository() << std::endl;

    fs::path jsonPath = fs::path(crbmodelDB.dbRepository())/fmt::format("{}.crb.json",crbmodelDB.name());
    jsonPathStr = jsonPath.string();
    if ( !fs::exists( jsonPath ) )
        throw std::runtime_error(fmt::format("crb db JSON file not found : {}",jsonPathStr));

    // TODO : move this part in CRBModelDB class
    std::ifstream ifs( jsonPathStr );
    nl::json jarg = nl::json::parse(ifs,nullptr,true,true);
    if ( jarg.contains( "crbmodel" ) )
    {
        auto const& j_crbmodel = jarg.at( "crbmodel" );
        std::string modelname = j_crbmodel.at( "name" ).get<std::string>();
        if ( j_crbmodel.contains( "plugin" ) )
        {
            auto const& j_plugin = j_crbmodel.at( "plugin" );
            pluginname = j_plugin.at( "name" ).get<std::string>();
            pluginlibname = j_plugin.at( "libname" ).get<std::string>();
        }
    }

    std::cout << "plugin name : " << pluginname << " libname : " << pluginlibname << std::endl;

    std::string pluginlibdir = Environment::expand( soption(_name="plugin.dir") );
    std::cout << "pluginlibdir="<<pluginlibdir<<std::endl;
    auto plugin = factoryCRBPlugin( pluginname, pluginlibname, pluginlibdir );
    std::cout << "Loaded the plugin " << plugin->name() << std::endl;

    bool loadFiniteElementDatabase = boption(_name="crb.load-elements-database");
    plugin->loadDB( jsonPathStr, (loadFiniteElementDatabase)? crb::load::all : crb::load::rb );

    return plugin;
}



inline Feel::AboutData makeAbout()
{
    Feel::AboutData about( "sensitibity_analysis",
                     "SA" ,
                     "0.1",
                     "Qensitibity analysis",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2022 Feel++ Consortium" );

    about.addAuthor( "Thomas Saigre", "developer", "saigre@math.unistra.fr", "" );
    return about;

}

int main( int argc, char** argv )
{
    using namespace Feel;
    po::options_description crbonlinerunoptions( "crb online run options" );
    crbonlinerunoptions.add_options()
        ( "plugin.dir", po::value<std::string>()->default_value(Info::libdir()) , "plugin directory" )
        ( "crbmodel.name", po::value<std::string>(), "CRB online code name" )
        ( "crbmodel.db.id", po::value<std::string>(), "CRB online code id" )
        ( "crbmodel.db.last", po::value<int>()->default_value( 2 ), "use last created(=1) or modified(=2) or not (=0)" )
        ( "crbmodel.db.root_directory", po::value<std::string>()->default_value( "${repository}/crbdb" ), "root directory of the CRB database " )

        ( "parameter", po::value<std::vector<std::string> >()->multitoken(), "database filename" )
        ( "sampling.size", po::value<int>()->default_value( 10 ), "size of sampling" )
        ( "sampling.type", po::value<std::string>()->default_value( "random" ), "type of sampling" )
        ( "rb-dim", po::value<int>()->default_value( -1 ), "reduced basis dimension used (-1 use the max dim)" )
        ( "output_results.save.path", po::value<std::string>(), "output_results.save.path" )

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

    runCrbOnline( { loadPlugin() } );

    OT::Point P(2);
    P.add(2.);
    Feel::cout << "OT::Point P :" << P << std::endl;

    OT::Normal N;
    OT::Sample S(N.getSample(10));

    Feel::cout << "OT::Sample S :" << S << std::endl;

    return 0;
}
