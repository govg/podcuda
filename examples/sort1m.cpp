#include <cstdlib>
#include <vector>
#include <algorithm>

using namespace std;

int main()
{

	vector<int> a;
	a.resize(32<<20);

	sort(a.begin(), a.end());

	return 0;
}

