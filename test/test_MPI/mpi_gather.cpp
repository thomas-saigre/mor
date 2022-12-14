#include <iostream>
#include <mpi.h>
#include <vector>

using namespace std;



int gather_array()
{
    MPI_Init(nullptr, nullptr);
    int  world_size , world_rank;
    MPI_Comm_size(MPI_COMM_WORLD , &world_size );
    MPI_Comm_rank(MPI_COMM_WORLD , &world_rank );

    int master_rank = 0;
    if (world_rank == master_rank) cout << "Gather vectors" << endl;
    
    int valsSend[5] = {1, 2, 3, 4, 5};

    int *valsRec;
    if (world_rank == master_rank)
        valsRec = new int[world_size * 5];

    MPI_Gather(valsSend, 5, MPI_INT, valsRec, 5, MPI_INT, master_rank, MPI_COMM_WORLD);

    if (world_rank == master_rank){
        size_t size = world_size * 5;
        for (size_t i = 0; i < size; i++)
            cout << valsRec[i] << " ";
        cout << endl;
        delete valsRec;
    }
    MPI_Finalize();
    return 0;
}

int gather_vectors()
{
    MPI_Init(nullptr, nullptr);
    int  world_size , world_rank;
    MPI_Comm_size(MPI_COMM_WORLD , &world_size );
    MPI_Comm_rank(MPI_COMM_WORLD , &world_rank );

    int master_rank = 0;
    if (world_rank == master_rank) cout << "Gather vectors" << endl;

    std::vector<int> valsSend = {1, 2, 3, 4, 5};

    std::vector<int> valsRec;
    if (world_rank == master_rank)
        valsRec.resize(world_size * 5);

    MPI_Gather(valsSend.data(), 5, MPI_INT, valsRec.data(), 5, MPI_INT, master_rank, MPI_COMM_WORLD);    

    if (world_rank == master_rank){
        size_t size = valsRec.size();
        for (size_t i = 0; i < size; i++)
            cout << valsRec[i] << " ";
        cout << endl;
    }
    MPI_Finalize();
    return 0;
}

int main()
{
#ifdef ARRAY
    gather_array();
#else
    gather_vectors();
#endif
}