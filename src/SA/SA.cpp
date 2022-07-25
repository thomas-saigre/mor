/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Thomas Saigre <saigre@math.unistra.fr>
       Date: 2022-07-25

  Copyright (C) 2022 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
    \file SA.cpp
    \author Thomas Saigre <@thomas-saigre>
    \date 2022-06-18
 */

#include <feel/options.hpp>
#include <openturns/OT.hxx>


#include <feel/feelmodels/heat/heat.hpp>


/** use Feel namespace */
using namespace Feel;

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline po::options_description makeOptions()
{
    po::options_description SAOptions( "Nirb options" );
    SAOptions.add_options()
    // meshes parameters
    ;

    return SAOptions.add( Feel::feel_options() ).add(Feel::toolboxes_options("heat"));
}

//-----------------------------------------
//-----------------------------------------
/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline AboutData makeAbout()
{
    AboutData about( "sensitibity_analysis",
                     "SA" ,
                     "0.1",
                     "Qensitibity analysis",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2022 Feel++ Consortium" );

    about.addAuthor( "Thomas Saigre", "developer", "saigre@math.unistra.fr", "" );
    return about;

}



//-----------------------------------------
//-----------------------------------------
/**
 * main function: entry point of the program
 */
int main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc = argc, _argv = argv,
                     _desc = makeOptions(),
                     _about = makeAbout() );
    

    Feel::cout << "Sensitivity analysis" << std::endl;


    OT::Point P(2);
    P.add(2.);
    Feel::cout << "OT::Point P :" << P << std::endl;

    OT::Normal N;
    OT::Sample S(N.getSample(10));

    Feel::cout << "OT::Sample S :" << S << std::endl;

}