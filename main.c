#include "hdf5.h"
#include "stdlib.h"

#define H5FILE_NAME     "h5p_file.h5"
#define DATASETNAME     "IntArray"
#define NX     256      /* dataset dimensions */
#define NY     5
#define RANK   2

/*
 * MPI variables
 */
int mpi_size, mpi_rank;
MPI_Comm comm;
MPI_Info info;

void write_data() {
    /*
     * HDF5 APIs definitions
     */
    hid_t file_id, dset_id;         /* file and dataset identifiers */
    hid_t filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t dimsf[2];               /* dataset dimensions */
    int *data;                      /* pointer to data buffer to write */
    hsize_t count[2];               /* hyperslab selection parameters */
    hsize_t offset[2];
    hid_t plist_id;                 /* property list identifier */
    int i;

    /*
     * Set up file access property list with parallel I/O access
     */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);

    /*
     * Create a new file collectively and release property list identifier.
     */
    file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    /*
     * Create the dataspace for the dataset.
     */
    dimsf[0] = NX;
    dimsf[1] = NY;
    filespace = H5Screate_simple(RANK, dimsf, NULL);

    /*
     * Create the dataset with default properties and close filespace.
     */
    dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_INT, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    count[0] = dimsf[0] / mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    memspace = H5Screate_simple(RANK, count, NULL);

    /*
     * Select hyperslab in the file.
     *
     * Specify the location in the dataset where we will save the data,
     * in this case it will be one or more rows depending on mpi_size and mpi_rank
     */
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    printf("Rank %d: Selecting offset=[%llu, %llu] count=[%llu, %llu]\n", mpi_rank, offset[0], offset[1], count[0], count[1]);

    /*
     * Initialize data buffer
     */
    data = (int *)malloc(sizeof(int) * count[0] * count[1]);
    for (i = 0; i < count[0] * count[1]; i++) {
        data[i] = mpi_rank * 7;
    }

    /*
     * Create property list for collective dataset write.
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                      plist_id, data);
    free(data);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Fclose(file_id);
}

int read_data() {
    /*
     * HDF5 APIs definitions
     */
    hid_t file_id, dset_id;         /* file and dataset identifiers */
    hid_t filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t dimsf[2] = {NX, NY};    /* dataset dimensions */
    int *data;                      /* pointer to data buffer to write */
    hsize_t count[2];               /* hyperslab selection parameters */
    hsize_t offset[2];
    int i;
    const int target_rank = (mpi_rank + 1) % mpi_size;
    int rc = 0;

    /*
     * Create a new file collectively and release property list identifier.
     */
    file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);

    /*
     * Create the dataset with default properties and close filespace.
     */
    dset_id = H5Dopen(file_id, DATASETNAME, H5P_DEFAULT);

    /*
     * Instead, we'll read the data written by the next rank
     */
    count[0] = dimsf[0] / mpi_size;
    count[1] = dimsf[1];
    offset[0] = target_rank * count[0];
    offset[1] = 0;
    memspace = H5Screate_simple(RANK, count, NULL);

    /*
     * Select hyperslab in the file.
     *
     * Specify the location in the dataset where we will save the data,
     * in this case it will be one or more rows depending on mpi_size and mpi_rank
     */
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    printf("Rank %d: Selecting offset=[%llu, %llu] count=[%llu, %llu]\n", mpi_rank, offset[0], offset[1], count[0], count[1]);

    /*
     * Initialize data buffer
     */
    data = (int *)malloc(sizeof(int) * count[0] * count[1]);

    /*
     * Read the selection into buffer
     */
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, data);

    for (i = 0; i < count[0] * count[1]; i++) {
        if (data[i] != target_rank * 7) {
            const hsize_t gx = offset[0] + i % count[0], gy = offset[1] + i % count[1];
            printf("Rank %d: Data at index %d was %d instead of %d! Global coords [%llu, %llu]\n",
                   mpi_rank, i, data[i], target_rank * 7, gx, gy);
            rc = 1;
        }
    }

    free(data);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Fclose(file_id);

    return rc;
}

int main(int argc, char **argv) {
    /*
     * Initialize MPI
     */
    comm = MPI_COMM_WORLD;
    info = MPI_INFO_NULL;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    if (mpi_rank == 0) printf("Writing\n");
    write_data();

    MPI_Barrier(comm);

    if (mpi_rank == 0) printf("\nReading\n");
    int rc = read_data();

    MPI_Finalize();

    return rc;
}
