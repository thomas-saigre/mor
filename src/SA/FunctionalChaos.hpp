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
//! @file FunctionalChaos.hpp
//! @author Thomas Saigre <saigre@math.unistra.fr>
//! @date 04 Jan 2023
//! @copyright 2023 Feel++ Consortium
//! @brief Functional chaos sensitivity analysis adapted from https://openturns.github.io/openturns/latest/auto_meta_modeling/polynomial_chaos_metamodel/plot_chaos_sobol_confidence.html
//!


#include <openturns/OT.hxx>


OT::FunctionalChaosResult computeSparseLeastSquaresChaos( OT::Sample X, OT::Sample Y,
    OT::OrthogonalProductPolynomialFactory basis, OT::UnsignedInteger total_degree, OT::Distribution distribution)
{
    OT::LeastSquaresMetaModelSelectionFactory selectionAlgorithm;
    OT::LeastSquaresStrategy projectionStrategy(selectionAlgorithm);
    OT::EnumerateFunction enum_fonc = basis.getEnumerateFunction();
    OT::UnsignedInteger P = enum_fonc.getBasisSizeFromTotalDegree( total_degree );
    OT::FixedStrategy adaptiveStrategy( basis, P );

    OT::FunctionalChaosAlgorithm polynomialChaosAlgorithm(X, Y, distribution, adaptiveStrategy, projectionStrategy);
    polynomialChaosAlgorithm.run();

    OT::FunctionalChaosResult result = polynomialChaosAlgorithm.getResult();
    return result;
}

void checkMetaModel(OT::UnsignedInteger n_valid, OT::Sample X_test, OT::Sample Y_test, OT::Function metamodel)
{
    OT::MetaModelValidation validation(X_test, Y_test, metamodel);
    OT::Scalar Q2 = validation.computePredictivityFactor()[0];
    Feel::cout << Feel::tc::green << "Check of the metamodel : Q2 = " << Q2 << Feel::tc::reset << std::endl;
}


auto multiBootstrap(OT::Sample X, OT::Sample Y)
{
    OT::UnsignedInteger n = X.getSize();
    OT::Indices selection = OT::BootstrapExperiment::GenerateSelection(n, n);
    return std::make_tuple(X.select(selection), Y.select(selection));
}


auto computeChaosSensitivity( const OT::Sample X, const OT::Sample Y,
    OT::OrthogonalProductPolynomialFactory basis, OT::UnsignedInteger total_degree, OT::Distribution distribution)
{
    size_t dim_input = X.getDimension();
    OT::FunctionalChaosResult result = computeSparseLeastSquaresChaos(X, Y, basis, total_degree, distribution);
    OT::FunctionalChaosSobolIndices chaosSI(result);

    OT::Collection<OT::Scalar> first_order(dim_input);
    OT::Collection<OT::Scalar> total_order(dim_input);

    for (size_t i = 0; i < dim_input; ++i)
    {
        first_order[i] = chaosSI.getSobolIndex(i);
        total_order[i] = chaosSI.getSobolTotalIndex(i);
    }

    return std::make_tuple(first_order, total_order);
}

auto computeBootstrapChaosSobolIndices( const OT::Sample X, const OT::Sample Y,
    OT::OrthogonalProductPolynomialFactory basis, OT::UnsignedInteger total_degree, OT::Distribution distribution,
    size_t bootstrap_size, OT::Scalar eps = 1e-9)
{
    size_t dim_input = X.getDimension();
    OT::Sample fo_sample (0, dim_input);
    OT::Sample to_sample (0, dim_input);
    OT::Point low(dim_input), high(dim_input);
    for (size_t i = 0; i < dim_input; ++i)
    {
        low[i] = eps;
        high[i] = 1 - eps;
    }
    OT::Interval unit_eps(low, high);

    for (size_t i = 0; i < bootstrap_size; ++i)
    {
        auto [X_boot, Y_boot] = multiBootstrap(X, Y);

        auto [fo, to] = computeChaosSensitivity(X_boot, Y_boot, basis, total_degree, distribution);
        if (unit_eps.contains(fo) && unit_eps.contains(to))
        {
            fo_sample.add(fo);
            to_sample.add(to);
        }
    }
    return std::make_tuple(fo_sample, to_sample);
}


auto computeSobolIndicesConfidenceInterval(OT::Sample fo_sample, OT::Sample to_sample, OT::Scalar alpha=0.95)
{
    size_t dim_input = fo_sample.getDimension();
    OT::Point fo_lb(dim_input), fo_ub(dim_input), to_lb(dim_input), to_ub(dim_input);
    for (size_t i = 0; i < dim_input; ++i)
    {
        OT::Sample fo_i = fo_sample.getMarginal(i),
                   to_i = to_sample.getMarginal(i);
        OT::Scalar beta = (1. - alpha) / 2.;
        fo_lb[i] = fo_i.computeQuantile( beta )[0];
        fo_ub[i] = fo_i.computeQuantile( 1-beta )[0];
        to_lb[i] = to_i.computeQuantile( beta )[0];
        to_ub[i] = to_i.computeQuantile( 1-beta )[0];
    }

    OT::Interval fo_interval(fo_lb, fo_ub);
    OT::Interval to_interval(to_lb, to_ub);
    return std::make_tuple(fo_interval, to_interval);
}


OT::Graph computeAndDrawSobolIndices( Results &res, OT::Sample X, OT::Sample Y, OT::OrthogonalProductPolynomialFactory basis,
    OT::UnsignedInteger total_degree, OT::Distribution distribution, size_t bootstrap_size=500, OT::Scalar alpha = 0.95)
{
    size_t N = X.getSize();
    Feel::cout << "Compute bootstrap chaos Sobol indices" << std::endl;
    Feel::tic();
    auto  [fo_sample, to_sample] = computeBootstrapChaosSobolIndices(X, Y, basis, total_degree, distribution, bootstrap_size);
    Feel::toc("computeBootstrapChaosSobolIndices");
    Feel::cout << "Compute Sobol indices confidence interval" << std::endl;
    Feel::tic();
    auto [fo_interval, to_interval] = computeSobolIndicesConfidenceInterval(fo_sample, to_sample, alpha);
    Feel::toc("computeSobolIndicesConfidenceInterval");

    res.setIndices( fo_sample.computeMean(), 1);
    res.setIndices( to_sample.computeMean(), 2);
    res.setInterval( fo_interval, 1);
    res.setInterval( to_interval, 2);

    OT::Description input_names = X.getDescription();
    OT::Graph graph = OT::SobolIndicesAlgorithm::DrawSobolIndices( input_names,
        fo_sample.computeMean(), to_sample.computeMean(), fo_interval, to_interval);
    graph.setTitle(OT::String("Sobol indices for ") + std::to_string(N) + " samples");

    return graph;
}
