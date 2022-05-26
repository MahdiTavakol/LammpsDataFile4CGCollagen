#ifndef ANGLE_H
#define ANGLE_H
#include <fstream>
class Angle
{
public:
	Angle(const int&, const int&, const int&, const int& );
	~Angle();
	void GetAtoms(int&, int&, int&);
	void GetType(int& );
protected:
	int id1, id2, id3;
	int type;
};
#endif

