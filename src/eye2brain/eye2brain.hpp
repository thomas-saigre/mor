/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_EYE2BRAIN_HPP
#define FEELPP_EYE2BRAIN_HPP 1

#include <feel/options.hpp>

#include <feel/feelmor/modelcrbbase.hpp>

namespace Feel
{

FEELPP_EXPORT po::options_description
makeEye2BrainOptions();
FEELPP_EXPORT AboutData
makeEye2BrainAbout( std::string const& str = "eye2brain" );

struct FEELPP_EXPORT Eye2BrainConfig
{
    typedef Mesh<Simplex<3>> mesh_type;
    typedef Pch_type<mesh_type,1> space_type;
};

class FEELPP_EXPORT Eye2Brain : public ModelCrbBase<ParameterSpace<>, Eye2BrainConfig::space_type >
{
    typedef ModelCrbBase<ParameterSpace<>, Eye2BrainConfig::space_type > super_type;

public:
    Eye2Brain();

    void initBetaQ();

    super_type::betaq_type computeBetaQ( parameter_type const& mu );


    void updateSpecificityModelPropertyTree( boost::property_tree::ptree & ptree ) const;

    //! initialisation of the model
    void initModel();

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);


    // parameter_type crbParametersFromUserParameters( feelpp4fastsim::UserParameters const& userParam ) const;
    // void updateFieldsExported( SourceCrbInterface * pvsource, element_type & feField, vtkSmartPointer<vtkUnstructuredGrid> vtkOutput );


};

}

#endif /* FEELPP_EYE2BRAIN_HPP */
