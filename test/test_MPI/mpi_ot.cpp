#include <iostream>
#include <mpi.h>
#include <vector>
#include <openturns/OT.hxx>
#include <ctime>
#include <typeinfo>

using namespace std;

OT::Sample getSample(size_t size, size_t dim)
{
    OT::Collection<OT::Distribution> marginals( dim );
    for (size_t i=0; i<dim; ++i)
        marginals[i] = OT::Uniform(0.0, 1.0);
    OT::Distribution distribution = OT::ComposedDistribution( marginals );
    OT::Sample sample = distribution.getSample(size);
    return sample;
}


int gather(size_t sample_size, size_t dim)
{
    MPI_Init(nullptr, nullptr);
    int  world_size , world_rank;
    MPI_Comm_size(MPI_COMM_WORLD , &world_size );
    MPI_Comm_rank(MPI_COMM_WORLD , &world_rank );
    OT::RandomGenerator::SetSeed( ::time(NULL) * (world_rank + 1) );

    int master_rank = 0;
    if (world_rank == master_rank) cout << "Gather vectors" << endl;

    OT::Sample valsSend = getSample(sample_size, dim);
    std::cout << world_rank << " "  << valsSend << std::endl;

    OT::Sample valsRec_sample;
    double *valsRec;
    if (world_rank == master_rank)
    {
        valsRec = new double[world_size * sample_size * dim];
        valsRec_sample = OT::Sample(world_size * sample_size, dim);
    }
    
    MPI_Gather(valsSend.data(), sample_size * dim, MPI_DOUBLE, valsRec, sample_size * dim, MPI_DOUBLE, master_rank, MPI_COMM_WORLD);

    if (world_rank == master_rank){
        std::cout << "Gather: " << std::endl;
        for (size_t i=0; i<world_size * sample_size * dim; ++i)
        {
            cout << valsRec[i] << " ";
        }
        cout << endl;
        size_t size = world_size * sample_size;
        for (size_t i=0; i<size; ++i)
        {
            for (size_t j=0; j<dim; ++j)
            {
                valsRec_sample[i][j] = valsRec[i * dim + j];
            }
        }
        std::cout << valsRec_sample << std::endl;
        delete valsRec;
    }

    MPI_Finalize();
    return 0;
}


int main(int argc, char *argv[])
{
    size_t sample_size = (argc > 1) ? atoi(argv[1]) : 5;
    size_t dim = (argc > 2) ? atoi(argv[2]) : 2;
    gather(sample_size, dim);
}