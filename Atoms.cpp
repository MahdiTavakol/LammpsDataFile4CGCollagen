#include "Atoms.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cctype>
using namespace std;

// class Atom definition
/*Atom::Atom()
{
	cout << "   " << endl;

}*/

Atom::Atom(int atomIDVal, string atomNameVal, string residueStr, char chainChar, 
           int resIDVal, float xVal, float yVal, float zVal, 
		   float betaVal, float occVal , char atomCharVal,
		   string atomTypeVal)
{
	atomID         = atomIDVal;
	atomName       = atomNameVal;
	residue        = residueStr;
	chain          = chainChar;
	residueID      = resIDVal;
	x              = xVal;
	y              = yVal;
	z              = zVal;
	beta           = betaVal;
	occupancy      = occVal;
	atomChar       = atomCharVal;
	atomType       = atomTypeVal;
	typeNumber     = 1;
}

Atom::Atom( const Atom& myAtom ):
	atomID( myAtom.atomID), atomName( myAtom.atomName),
	residue( myAtom.residue), chain( myAtom.chain),
	residueID(myAtom.residueID), 
	x( myAtom.x ), y( myAtom.y ), z( myAtom.z ),
	beta( myAtom.beta ), occupancy( myAtom.occupancy ),
	atomChar( myAtom.atomChar), typeNumber (myAtom.typeNumber)
{}

void Atom::SetAtom( int atomIDVal, string atomNameVal, string residueStr, char chainChar, 
           int resIDVal, float xVal, float yVal, float zVal, 
		   float betaVal, float occVal , char atomCharVal)
{
	atomID         = atomIDVal;
	atomName       = atomNameVal;
	residue        = residueStr;
	chain          = chainChar;
	residueID      = resIDVal;
	x              = xVal;
	y              = yVal;
	z              = zVal;
	beta           = betaVal;
	occupancy      = occVal;
	atomChar       = atomCharVal;
}

int Atom::ReadPDBLine(string line)
{
	if ( line.substr(0,3) == "END")  return 3;
	else if ( line.substr(0,4) != "ATOM" ) return 2;
	else if ( line.substr(0,4) == "ATOM")
	{
		stringstream SSR(" ");
		string atomIDStr = line.substr(4,7);
		SSR.clear();
		bool has_only_digits = (atomIDStr.find_first_not_of("0123456789") == std::string::npos);
		if (has_only_digits)
		{
			SSR << atomIDStr;
		}
		else
		{
			SSR << std::hex << atomIDStr;
		}
		SSR >> atomID;
		atomName = line.substr(13,2);
		residue  = line.substr(17,4);
		chain    = line[21];
		string residueIDStr = line.substr(22,4);
		SSR.clear();
		SSR.str(residueIDStr);
		SSR >> residueID;
		string xStr = line.substr(30,8);
		SSR.clear();
		SSR.str(xStr);
		SSR >> x;
		string yStr = line.substr(38,8);
		SSR.clear();
		SSR.str(yStr);
		SSR >> y;
		string zStr = line.substr(46,8);
		SSR.clear();
		SSR.str(zStr);
		SSR >> z;
		string betaStr = line.substr(55,5);
		SSR.clear();
		SSR.str(betaStr);
		SSR >> beta;
		/*string occStr = line.substr(61,5);
		SSR.clear();
		SSR.str(occStr);
		SSR >> occupancy;
		atomChar = line[72];*/
		return 1;
	}
	return 0;	
}

int Atom::ReadGROLine(const string& line)
{
	stringstream SSR(" ");
	string residueIDStr = line.substr(0, 5);
	SSR.clear();
	SSR.str(residueIDStr);
	SSR >> residueID;
	SSR.clear();
	residue = line.substr(5, 5);
	string atomNameStr = line.substr(12, 2);
	SSR.str(atomNameStr);
	SSR >> atomName;
	string atomIDStr = line.substr(15, 5);
	SSR.clear();
	SSR.str(atomIDStr);
	SSR >> atomID;
	string xStr = line.substr(20, 8);
	SSR.clear();
	SSR.str(xStr);
	SSR >> x;
	string yStr = line.substr(28, 8);
	SSR.clear();
	SSR.str(yStr);
	SSR >> y;
	string zStr = line.substr(36, 8);
	SSR.clear();
	SSR.str(zStr);
	SSR >> z;
	return 0;
}

void Atom::PrintAtomsXYZ(ofstream& file)
{
	
	file << " ";
	file.width(8);
	file << fixed << left << atomName ;
	file.width(10)                                           ;
	file << fixed << setprecision(6) << right << x           ;
	file << "      "                                         ;
	file.width(10)                                           ;
	file << fixed << setprecision(6) << right << y           ;
	file << "      "                                         ;
	file.width(10)                                           ;
	file << fixed << setprecision(6) << right << z           ;
	file << endl                                             ;
}

void Atom::PrintAtomsXYZ2(ofstream& file)
{
	this->SetAtomChar();
	file << " ";
	file.width(8);
	file << fixed << left << atomChar;
	file.width(10);
	file << fixed << setprecision(6) << right << x;
	file << "      ";
	file.width(10);
	file << fixed << setprecision(6) << right << y;
	file << "      ";
	file.width(10);
	file << fixed << setprecision(6) << right << z;
	file << endl;
}

void Atom::PrintPDBHeader(ofstream& file, const string& line)
{
	file << line << endl;
}

void Atom::PrintAtomsPDB(ofstream& file)
{
	file << "ATOM";
	file.width(7);
	file << fixed << right << this->atomID;
	file.width(4);
	file << fixed << left << this->atomName;
	file.width(4);
	file << fixed << left << this->residue;
	file.width(1);
	file << fixed << right << this->chain;
	file.width(4);
	file << fixed << right << this->residueID;
	file << "    ";
	file.width(8);
	file << fixed << right << setprecision(3) << this->x;
	file.width(8);
	file << fixed << right << setprecision(3) << this->y;
	file.width(8);
	file << fixed << right << setprecision(3) << this->z;
        file << " ";
        file.width(5);
	file << fixed << right << setprecision(2) << this->beta;
	file << " ";
    file.width(5);
	file << fixed << right << setprecision(2) << this->occupancy;
    file << "           ";
	file.width(1);
    file << fixed << right << this->atomChar;
	file << endl;
}

void Atom::PrintPDBFooter(ofstream& file)
{
	file << "END";
}

int Atom::GetAtomID()
{
	return atomID;
}

void  Atom::SetAtomID(const int& atomIDVal)
{
	this->atomID = atomIDVal;
}

string Atom::GetAtomName()
{
	return this->atomName;
}

string Atom::GetResidue()
{
	return this->residue;
}

int Atom::GetResidueID()
{
	return this->residueID;
}

void Atom::SetResidueID(const int& newID)
{
	this->residueID = newID;
}

void  Atom::GetCoordinates(float* coor )
{
	coor[0] = this->x;
	coor[1] = this->y;
	coor[2] = this->z;
}

void  Atom::SetCoordinates(float* coor)
{
	this->x = coor[0];
	this->y = coor[1];
	this->z = coor[2];
}

void Atom::GetFlag(bool &flagStatues)
{
	//flagStatues = flag;
}

void Atom::SetFlag(const bool &newflag)
{
	flag = newflag;
}

void Atom::GetBeta(float& betaVal)
{
	betaVal = this->beta;
}

void Atom::SetBeta(const float& betaVal)
{
	this->beta = betaVal;
}

void Atom::GetTypeNumber(int& typeNumberVal)
{
	typeNumberVal = this->typeNumber;
}

void Atom::SetTypeNumber(const int& typeNumberVal)
{
	this->typeNumber = typeNumberVal;
}

void Atom::SetType(const string &newtype)
{
	atomType = newtype;
}


void Atom::SetAtomChar()
{
    // Old Thing !!!!!!!!!!
}
