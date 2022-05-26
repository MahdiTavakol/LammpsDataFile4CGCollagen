#include "Angles.h"



Angle::Angle(const int& atom1, const int& atom2, const int& atom3, const int& typeVal)
{
	id1 = atom1;
	id2 = atom2;
	id3 = atom3;
	type = typeVal;
}


Angle::~Angle()
{
}

void Angle::GetAtoms(int& atom1, int& atom2, int& atom3)
{
	atom1 = id1;
	atom2 = id2;
	atom3 = id3;
}

void Angle::GetType(int& typeVal)
{
	typeVal = this->type;
}