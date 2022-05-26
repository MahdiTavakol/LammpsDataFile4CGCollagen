#include "Bonds.h"

Bond::Bond(const int& newAtom1, const int& newAtom2)
{
	id1 = newAtom1;
	id2 = newAtom2;
}


Bond::~Bond()
{
}

void Bond::GetAtoms(int& newAtom1, int& newAtom2)
{
	newAtom1 = id1;
	newAtom2 = id2;
}