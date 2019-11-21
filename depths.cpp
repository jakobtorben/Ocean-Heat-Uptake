#include <iostream>
#include <iterator>
#include <tuple>
using namespace std;

tuple<int, int> depth_ind(double depths[], int depths_size, double depth_from, double depth_to)
{
    int count = 0;
    int depth_from_i;
    int depth_to_i;
    for (int i = 0; i < depths_size; i++) {
        if ( (depths[i] >= depth_from) && (count == 0) )
        {
            depth_from_i = i;
            count += 1;
        }
        else if ( (depths[i] > depth_to) && (depths[i - 1] <= depth_to) ) {
            depth_to_i = i - 1; }
        }
    return make_tuple(depth_from_i, depth_to_i);
}

int main() {
    double depths[5] = {1.2, 2.2, 3.3, 4.4, 5.5};
    int depths_size = sizeof(depths)/ sizeof(depths[0]);
    double depth_from = 0;
    double depth_to = 4;

    auto [depth_from_i, depth_to_i] = depth_ind(depths, depths_size, depth_from, depth_to);

    cout << "Depth from " << depth_from_i << endl;
    cout << "Depth to " << depth_to_i << endl;

    return 0;
}