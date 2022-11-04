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

#include <omp.h>


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



int main( int argc, char** argv )
{
    using namespace Feel;
    po::options_description options( "options" );
    options.add_options()
        ( "app.test", po::value<int>()->default_value( 2 ), "test" )
        ;

    Environment env( _argc = argc, _argv = argv,
                     _desc = options,
                     _about = makeAbout() );

    assert( ioption("app.test") == 2 );

    return 0;
}
