#include <iostream>
#include <mpi.h>
#include <vector>
#include <openturns/OT.hxx>
#include <ctime>
#include <typeinfo>

using namespace std;

OT::Sample getSample(size_t size)
{
    OT::Distribution distribution = OT::Normal();
    OT::Sample sample = distribution.getSample(size);
    return sample;
}


int gather()
{
    MPI_Init(nullptr, nullptr);
    int  world_size , world_rank;
    MPI_Comm_size(MPI_COMM_WORLD , &world_size );
    MPI_Comm_rank(MPI_COMM_WORLD , &world_rank );
    OT::RandomGenerator::SetSeed( ::time(NULL) * (world_rank + 1) );

    int master_rank = 0;
    if (world_rank == master_rank) cout << "Gather vectors" << endl;

    OT::Sample valsSend = getSample(5);
    std::cout << world_rank << " "  << valsSend << std::endl;

    OT::Sample valsRec_sample;
    double *valsRec;
    if (world_rank == master_rank)
    {
        valsRec = new double[world_size * 5];
        valsRec_sample = OT::Sample(world_size * 5, 1);
    }
    
    MPI_Gather(valsSend.data(), 5, MPI_DOUBLE, valsRec, 5, MPI_DOUBLE, master_rank, MPI_COMM_WORLD);

    if (world_rank == master_rank){
        size_t size = world_size * 5;
        for (size_t i=0; i<size; ++i)
        {
            valsRec_sample[i] = OT::Point( {valsRec[i]} );
        }
        std::cout << valsRec_sample << std::endl;
        delete valsRec;
    }

    MPI_Finalize();
    return 0;
}


int main()
{
    gather();
}