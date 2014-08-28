#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>

#include <cstdlib>

using namespace thrust;

int main(void)
{

	host_vector<int> h_vec(32<<20);
	generate(h_vec.begin(), h_vec.end(), rand);

	device_vector<int> d_vec = h_vec;

	sort(d_vec.begin(), d_vec.end());

	copy(d_vec.begin(), d_vec.end(), h_vec.begin());

	return 0;

}	
