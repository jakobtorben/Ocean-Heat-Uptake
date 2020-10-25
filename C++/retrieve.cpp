#include <iostream>
#include <netcdf>
#include <iterator>
#include <tuple>
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


// We are reading 2D data, a 6 x 12 grid.
static const int NX = 6;
static const int NY = 12;

// Return this in event of a problem.
static const int NC_ERR = 2;

int main()
{
    try
    {
        // This is the array we will read.
        int dataIn[NX][NY];

        // Open the file for read access
        NcFile dataFile("EN.4.2.1.f.analysis.g10.195001.nc", NcFile::read);

        // Retrieve the variable named "data"
        NcVar data=dataFile.getVar("temperature");
        if(data.isNull()) return NC_ERR;
        data.getVar(dataIn);

        // Check the values.
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
                if (dataIn[i][j] != i * NY + j)
                    return NC_ERR;

        // The netCDF file is automatically closed by the NcFile destructor
        //cout << "*** SUCCESS reading example file simple_xy.nc!" << endl;

        return 0;