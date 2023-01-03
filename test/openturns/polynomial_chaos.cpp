#include <openturns/OT.hxx>
#include <vector>

/**
 * @brief Run the test adapted from https://openturns.github.io/openturns/latest/auto_reliability_sensitivity/sensitivity_analysis/plot_functional_chaos_sensitivity.html
 * 
 */
void run()
{
    // creation of the distribution
    size_t dimension = 8;
    OT::Description input_names = {"rw", "r", "Tu", "Hu", "Tl", "Hl", "L", "Kw"};
    OT::SymbolicFunction model = OT::SymbolicFunction(
        input_names, {"(2 * pi_ * Tu * (Hu-Hl)) / ( ln(r/rw) * (1+(2*L*Tu)/(ln(r/rw)*rw^2*Kw) + Tu/Tl) )"}
    );

     OT::Collection<OT::Distribution> coll = {
        OT::Normal(0.1, 0.0161812),
        OT::LogNormal(7.71, 1.0056),
        OT::Uniform(63070.0, 115600.0),
        OT::Uniform(990.0, 1110.0),
        OT::Uniform(63.1, 116.0),
        OT::Uniform(700.0, 820.0),
        OT::Uniform(1120.0, 1680.0),
        OT::Uniform(9855.0, 12045.0)
    };
    OT::Distribution distribution = OT::ComposedDistribution(coll);
    distribution.setDescription(input_names);

    // Freeze r, Tu, Tl from model to go faster
    OT::Indices selection = {1, 2, 4};
    OT::Indices complement = OT::Indices(selection).complement(dimension);
    distribution = distribution.getMarginal(complement);
    OT::ParametricFunction model2 = OT::ParametricFunction(model, selection, distribution.getMarginal(selection).getMean());
    input_names = distribution.getDescription();
    std::cout << "Input names = " << input_names << std::endl;
    dimension = distribution.getDimension();

    // design of experiments
    size_t size = 1000;
    OT::Sample X = distribution.getSample(size);
    OT::Sample Y = model2(X);

    // create a functional chaos model
    OT::FunctionalChaosAlgorithm algo(X, Y);
    algo.run();

    OT::FunctionalChaosResult result = algo.getResult();
    std::cout << "Residuals = " << result.getResiduals() << std::endl;
    std::cout << "Relative errors = " << result.getRelativeErrors() << std::endl;
    std::cout << std::endl;

    // sensitivity analysis
    OT::FunctionalChaosSobolIndices sensitivityAnalysis = OT::FunctionalChaosSobolIndices(result);
    // std::cout << sensitivityAnalysis << std::endl;

    // draw the Sobol' indices
    OT::Collection<OT::Scalar> first_order(dimension), total_order(dimension);
    for (size_t i = 0; i < dimension; ++i)
    {
        first_order[i] = sensitivityAnalysis.getSobolIndex(i);
        total_order[i] = sensitivityAnalysis.getSobolTotalIndex(i);
    }
    std::cout << "First order indices = " << first_order << std::endl;
    std::cout << "Total order indices = " << total_order << std::endl;
    
    OT::Graph graph = OT::SobolIndicesAlgorithm::DrawSobolIndices(input_names, first_order, total_order);
    graph.setTitle("Sobol' indices");
    graph.draw("sensitivity_analysis.png");

    for (size_t i = 0; i < dimension; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            std::cout <<  "Interaction indices (" << input_names[i] << ", " << input_names[j] << ") = " << sensitivityAnalysis.getSobolIndex({i, j}) << std::endl;
        }
    }

}


int main(int argc, char *argv[])
{
    run();
    return 0;
}