import openturns as ot
from math import log, pi

ot.RandomGenerator.SetSeed(0)


dimension = 8
input_names = ['rw', 'r', 'Tu', 'Hu', 'Tl', 'Hl', 'L', 'Kw']

coll = [ot.Normal(0.1, 0.0161812),
        ot.LogNormal(7.71, 1.0056),
        ot.Uniform(63070.0, 115600.0),
        ot.Uniform(990.0, 1110.0),
        ot.Uniform(63.1, 116.0),
        ot.Uniform(700.0, 820.0),
        ot.Uniform(1120.0, 1680.0),
        ot.Uniform(9855.0, 12045.0)]
X = ot.ComposedDistribution(coll)
X.setDescription(input_names)

# T = X.getSample(10)
# print(T)

def computeOutput(X):
    s = 2 * pi * X[2] * (X[3]-X[4]) / (log(X[1]/X[0]) * (1 + (2*X[6]*X[2]) / (log(X[1]/X[0]) * X[0]**2 * X[7]) + X[2]/X[6]))
    return [s]

computeOutputOt = ot.PythonFunction(8, 1, computeOutput)
computeOutputOt.setDescription(input_names + ["Y"])

if False:
    params = X.getSample(3000)
    sample = computeOutputOt(params)


size = 5000
computeSecondOrder = True
sie = ot.SobolIndicesExperiment(X, size, computeSecondOrder)
inputDesign = sie.generate()
# print("inputDesign:\n", inputDesign)

outputDesign = computeOutputOt(inputDesign)
# print("outputDesign:\n", outputDesign)


sensitivityAnalysis = ot.SaltelliSensitivityAlgorithm(inputDesign, outputDesign, size)
sensitivityAnalysis.setUseAsymptoticDistribution(True)
# sensitivityAnalysis.draw()

print(sensitivityAnalysis.getFirstOrderIndices())
