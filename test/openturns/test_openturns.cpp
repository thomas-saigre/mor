#include <openturns/OT.hxx>
#include <initializer_list>
#include <math.h>

#if 1
OT::ComposedDistribution my_composed_distribution()
{
    OT::Distribution rw = OT::Normal(0.1, 0.0161812);
    OT::Distribution r = OT::LogNormal(7.71, 1.0056);
    OT::Distribution Tu = OT::Uniform(63070.0, 115600.0);
    OT::Distribution Hu = OT::Uniform(990.0, 1110.0);
    OT::Distribution Tl = OT::Uniform(63.1, 116.0);
    OT::Distribution Hl = OT::Uniform(700.0, 820.0);
    OT::Distribution L = OT::Uniform(1120.0, 1680.0);
    OT::Distribution Kw = OT::Uniform(9855.0, 12045.0);

    std::initializer_list<OT::Distribution> distributions = {rw, r, Tu, Hu, Tl, Hl, L, Kw};

    return OT::ComposedDistribution( distributions );
}


OT::Sample output(OT::Sample input)
{
    size_t n = input.getSize();
    OT::Sample output(n, 1);
    for (int i=0; i<n; ++i)
    {
        OT::Point X = input[i];
        output[i] = OT::Point({ 2 * M_PI * X[2] * (X[3]-X[4]) / (log(X[1]/X[0]) * (1 + (2*X[6]*X[2]) / (log(X[1]/X[0]) * X[0]*X[0] * X[7]) + X[2]/X[6])) });
    }
    return output;
}
#endif




int main(int, char *[])
{
    OT::RandomGenerator::SetSeed(0);

    OT::ComposedDistribution C = my_composed_distribution();
    std::initializer_list<std::string> names = {"rw", "r", "Tu", "Hu", "Tl", "Hl", "L", "Kw"};
    OT::Description input_names =  names;
    C.setDescription(input_names);

    // OT::Sample T(C.getSample(10));
    // std::cout << "OT::Sample of C :\n" << T << std::endl;

    size_t size = 5000;
    bool computeSecondOrder = true;
    OT::SobolIndicesExperiment sobol(C, size, computeSecondOrder);
    OT::Sample inputDesign = sobol.generate();
    // std::cout << "inputDesign :\n" << inputDesign << std::endl;

    OT::Sample outputDesign = output(inputDesign);
    // std::cout << "outputDesign :\n" << outputDesign << std::endl;

    OT::SaltelliSensitivityAlgorithm sensitivity(inputDesign, outputDesign, size);
    sensitivity.setUseAsymptoticDistribution( true );

    auto firstOrder = sensitivity.getFirstOrderIndices();
    std::cout << firstOrder << std::endl;

    return 0;
}