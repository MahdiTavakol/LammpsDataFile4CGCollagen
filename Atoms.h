#ifndef ATOM_H
#define ATOM_H
#include <fstream>
using namespace std;

// class Atom definition
class Atom
{
	public:
	    //Atom();
	    Atom(int = 0, string = "CA", string = "GLY", char = 'A', int = 0, float = 0.0f , 
		float = 0.0f , float = 0.0f, float = 0.0f, float = 0.0f, char = 'C', string = "CA"); // constructor
		Atom(const Atom&); // copy constructor
	    void   SetAtom( int , string , string , char , int , float , float , float , float , float , char);
		int    ReadPDBLine(string);
		int    ReadGROLine(const string&);
		void   PrintAtomsXYZ(ofstream& );
		void   PrintAtomsXYZ2(ofstream&); //just for debugging purposes
		void   PrintPDBHeader(ofstream&, const string& );
		void   PrintAtomsPDB(ofstream&);
		void   PrintPDBFooter(ofstream&);
        int    GetAtomID();
		void   SetAtomID(const int&);
		string GetAtomName();
		string GetResidue();
		int    GetResidueID();
		void   SetResidueID(const int&);
        void   GetCoordinates(float* );
		void   SetCoordinates(float*);
		void   GetFlag(bool& );
		void   SetFlag(const bool&);
		void   GetBeta(float&);
		void   SetBeta(const float&);
		void   SetType(const string&);
		void   SetAtomChar(); //just for debugging purposes
		void   GetTypeNumber(int&);
		void   SetTypeNumber(const int&);

	private:
		
		int atomID;
		string atomName;
		string residue;
		char chain;
		int residueID;
		float x;
		float y;
		float z;
		float beta;
		float occupancy;
		char  atomChar;
		int   flag = true; // For a later use;;;
		string atomType;
		int typeNumber;
}; // end class Atom

#endif
