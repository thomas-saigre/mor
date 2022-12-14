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
    
    int tab0[4] = {1 + 4*world_rank, 2 + 4*world_rank, 3 + 4*world_rank, 4 + 4*world_rank};
    int tab1[4] = {11 + 4*world_rank, 12 + 4*world_rank, 13 + 4*world_rank, 14 + 4*world_rank};

    cout << world_rank << "tab0: ";
    for (size_t i = 0; i < 4; i++)
        cout << tab0[i] << " ";
        cout << endl;
    cout << world_rank << "tab1: ";
    for (size_t i = 0; i < 4; i++)
        cout << tab1[i] << " ";
    cout << endl;

    int *tab0Gather, *tab1Gather;
    if (world_rank == master_rank)
    {
        tab0Gather = new int[world_size * 4];
        tab1Gather = new int[world_size * 4];
    }

    MPI_Gather(tab0, 4, MPI_INT, tab0Gather, 4, MPI_INT, master_rank, MPI_COMM_WORLD);
    MPI_Gather(tab1, 4, MPI_INT, tab1Gather, 4, MPI_INT, master_rank, MPI_COMM_WORLD);

    if (world_rank == master_rank){
        size_t size = world_size * 4;
        cout << "Gather0: ";
        for (size_t i = 0; i < size; i++)
            cout << tab0Gather[i] << " ";
        cout << endl;
        cout << "Gather1: ";
        for (size_t i = 0; i < size; i++)
            cout << tab1Gather[i] << " ";
        cout << endl;
        delete tab0Gather;
        delete tab1Gather;
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