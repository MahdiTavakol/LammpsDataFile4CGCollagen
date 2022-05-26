#ifndef BOND_H
#define BOND_H
#include <fstream>
using namespace std;
class Bond
{
public:
	Bond(const int&, const int& );
	~Bond();
	void GetAtoms(int&, int&);
protected:
	int id1;
	int id2;
};
#endif
