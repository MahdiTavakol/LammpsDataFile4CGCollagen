#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include "Atoms.h"
#include "Bonds.h"
#include "Angles.h"



using namespace std;

float calculateAngle(float *a, float *b, float *c);
void writeDataFile( string fileName, string dataFileName, int type, bool mineralized, bool extraFibril, int fileType);
void isThereCrossLink(const string& dataFileName, vector<Atom>& atoms, int bonds);
int main(int argc, char* argv[])
{
	string dataFileName;
	string fileName;
	bool mineralized = false;
	bool extraFibril = false;
	int resp;
	for (int i = 1; i < argc-1; i++)
	{
		if (strcmp(argv[i],"-i") == 0)
			fileName = argv[i+1];	
		else if (strcmp(argv[i],"-d") == 0)
			dataFileName = argv[i+1];
		else if (strcmp(argv[i], "-m") == 0)
		{
			if (strcmp(argv[i + 1], "y") == 0)
				mineralized = true;
		}
	}
	if (argc == 1)
	{
		bool debug = true;
		string response;
		cout << "Enter the input file Name:" << endl;
		//fileName = "4-Microfibril-Circular-Modified.gro";
		cin >> fileName;
		dataFileName = fileName;
		if (fileName == "test")
		{
			fileName = "4-MT-Mineralized2-9.01";
			dataFileName = fileName;
			debug = false;
			resp = 2;
			fileName = fileName + ".pdb";
			mineralized = true;
			extraFibril = true;
		}
		if (debug)
		{
			cout << "Press 1 for gro file and 2 for pdb file: " << endl;
			cin >> resp;
			switch (resp)
			{
			case 1:
				fileName = fileName + ".gro";
				break;
			case 2:
				fileName = fileName + ".pdb";
				break;
			}
			cout << "Is it mineralized:" << endl;
			cin >> response;
			if (response == "y")
				mineralized = true;
			if (response == "Y")
				mineralized = true;
			if (mineralized)
			{
				cout << "Is there extra fibrillar minerals:" << endl;
				cin >> response;
				if (response == "y")
					extraFibril = true;
				if (response == "Y")
					extraFibril = true;
			}
		}
	}
	int type = 1;
	string dataFile2Name = dataFileName + "2.dat";
	string dataFile3Name = dataFileName + "3.dat";
	string dataFile1Name = dataFileName + "1.dat";
	dataFileName = dataFileName + ".dat";
	//writeDataFile(fileName, dataFile1Name, type, mineralized);
	//type = 2;
	//writeDataFile(fileName, dataFile2Name, type, mineralized);
	//type = 3;
	//writeDataFile(fileName, dataFile3Name, type, mineralized);
	writeDataFile(fileName, dataFileName, type, mineralized,extraFibril, resp);
	system("PAUSE");
	return 0;
}
void writeDataFile(string fileName, string dataFileName, int type, bool mineralized, bool extraFibril, int fileType)
{
	ofstream dataFile(dataFileName);
	ifstream coorFile(fileName);
	ifstream typeFile("bead-types.txt");
	vector<int> typeVector;
	//ofstream bondLengthFile("BondLength.txt");
	stringstream SSR(" ");
	int numberOfAtoms, numberOfBonds, numberOfAngles;
	int collagenBeads = 0;
	int HApBeads = 0;
	int HApBeadsIntra = 0;
	int HApBeadsExtra = 0;
	string line;
	vector<Atom> atoms;
	vector<Bond> bonds;
	vector<Angle> angles;
	vector<float> BLengths;
	float xlo; // = 126.0f; // -123.0f;
	float xhi; // = 373.0f; // 123.0f;
	float ylo; // = 134.0f; // -123.0f;
	float yhi; // = 372.0f; // 123.0f;
	float sideBox = 40.0f;
	float zlo = 0.0f;
	float zhi = 3442.0f;
	int angleTypes = 0;
	vector<float> angleValues;
	float mass = 1358.7f;
	float massMineral = 1548.0f;
	float epsilon = 6.87f;
	float epsilonMineral = 106.7f;
	float sigma = 14.72f;
	float sigmaMineral = 10.28f;
	float LJCut = 31.0f;
	float LJCutMineral = 13.85f;
	switch (fileType)
	{
	case 1:
		getline(coorFile, line); //comment line
		getline(coorFile, line);
		SSR.clear();
		SSR << line;
		SSR >> numberOfAtoms;
		break;
	case 2:
		cout << "Please enter the number of atoms: " << endl;
		cin >> numberOfAtoms;
		getline(coorFile, line); //comment line
		break;
	}


	int currResID = 0;
	int atomIDinRes = 0;
	Atom newAtom1;
	Atom newAtom2;
	Atom newAtom3;

	while (getline(typeFile, line))
	{
		stringstream iss(line);
		int typeVal;
		iss >> typeVal;
		typeVector.push_back(typeVal);
	}

	for (int i = 0; i < numberOfAtoms; i++)
	{
		getline(coorFile, line);
		Atom newAtom;
		switch (fileType)
		{
		case 1:
			newAtom.ReadGROLine(line);
			break;
		case 2:
			newAtom.ReadPDBLine(line);
			break;
		}
		atoms.push_back(newAtom);
		if (strcmp(newAtom.GetAtomName().c_str(), "HA") == 0 || strcmp(newAtom.GetAtomName().c_str(), "HA1") == 0 || strcmp(newAtom.GetAtomName().c_str(), "HA2") == 0)
		{
			HApBeadsIntra++;
			continue;
		}
		if (strcmp(newAtom.GetAtomName().c_str(), "EA") == 0 || strcmp(newAtom.GetAtomName().c_str(), "EA1") == 0 || strcmp(newAtom.GetAtomName().c_str(), "EA2") == 0)
		{
			HApBeadsExtra++;
			continue;
		}
		collagenBeads++;
		float beta;
		if (currResID != newAtom.GetResidueID())
		{
			atomIDinRes = 0;
			currResID = newAtom.GetResidueID();
			beta = 1.0f;
			atoms[i].SetBeta(beta);
			if (i)
				atoms[i - 1].SetBeta(beta);
		}
		else if (currResID == newAtom.GetResidueID())
		{
			atomIDinRes++;
			beta = 0.0f;
			atoms[i].SetBeta(beta);
		}
		if (atoms[i].GetAtomName().c_str() == "CL") 
			atoms[i].SetTypeNumber(typeVector[atomIDinRes]);
		


		if (atomIDinRes > 0)
		{
			Bond newBond(i-1,i);
			float coor1[3], coor2[3];
			atoms[i].GetCoordinates(coor1);
			atoms[i - 1].GetCoordinates(coor2);
			float x = coor1[0] - coor2[0];
			float y = coor1[1] - coor2[1];
			float z = coor1[2] - coor2[2];
			float length = sqrt(x*x + y*y + z*z);
			BLengths.push_back(length);
			bonds.push_back(newBond);
		}
		if (atomIDinRes > 1)
		{
			float coor1[3], coor2[3], coor3[3];
			int index;
			atoms[i - 2].GetCoordinates(coor1);
			atoms[i - 1].GetCoordinates(coor2);
			atoms[i].GetCoordinates(coor3);
			float newAngleValue;
			if ((coor3[2] > coor2[2]) && (coor2[2] > coor1[2]))
			{ 
				newAngleValue = calculateAngle(coor1, coor2, coor3);
			}
			else if (coor2[2] < coor1[2])
			{
				float coor2pbcz[3], coor3pbcz[3];
				coor2pbcz[0] = coor2[0];
				coor2pbcz[1] = coor2[1];
				coor2pbcz[2] = coor2[2] + zhi - zlo;
				coor3pbcz[0] = coor3[0];
				coor3pbcz[1] = coor3[1];
				coor3pbcz[2] = coor3[2] + zhi - zlo;
				newAngleValue = calculateAngle(coor1, coor2pbcz, coor3pbcz);
			}
			else if (coor3[2] < coor2[2])
			{

				float coor3pbcz[3];
				coor3pbcz[0] = coor3[0];
				coor3pbcz[1] = coor3[1];
				coor3pbcz[2] = coor3[2] + zhi - zlo;
				newAngleValue = calculateAngle(coor1, coor2, coor3pbcz);
			}
			vector<float>::iterator it = find(angleValues.begin(), angleValues.end(), newAngleValue);
			if (it == angleValues.end())
			{
				angleValues.push_back(newAngleValue);
				index = angleValues.size() - 1;
			}
			else
			{
				index = distance(angleValues.begin(), it);
			}
			Angle newAngle(i-2, i - 1, i, index+1);
			angles.push_back(newAngle);
		}
		int progress = (i + 1)*100 / numberOfAtoms;
		cout << "Calculating: " << progress << "%" << endl;
	}
	xlo = 999999.9f;
	ylo = 999999.9f;
	zlo = 999999.9f;
	xhi = -999999.9f;
	yhi = -999999.9f;
	zhi = -999999.9f;
	for (int i = 0; i < numberOfAtoms; i++)
	{
		float coor[3];
		atoms[i].GetCoordinates(coor);
		if (coor[0] < xlo)
			xlo = coor[0];
		if (coor[0] > xhi)
			xhi = coor[0];
		if (coor[1] < ylo)
			ylo = coor[1];
		if (coor[1] > yhi)
			yhi = coor[1];
		if (coor[2] < zlo)
			zlo = coor[2];
		if (coor[2] > zhi)
			zhi = coor[2];
	}
	xlo = 10.0f*xlo;
	xhi = 10.0f*xhi;
	ylo = 10.0f*ylo;
	yhi = 10.0f*yhi;
	zlo = 10.0f*zlo;
	zhi = 10.0f*zhi;

	char ans;
	cout << "Is the whole simulation box filled with extrafibrillar mineral?" << endl;
	cin >> ans;
	if (ans == 'y')
		sideBox = 16.322f / 4.0f;
	xlo -= sideBox;
	xhi += sideBox;
	ylo -= sideBox;
	yhi += sideBox;
	zhi += 14.0f/2.0f;


	HApBeads = numberOfAtoms - collagenBeads;
	float mineralizationPercent = 100.0f*((float)HApBeads)*massMineral / (((float)HApBeads)*massMineral + ((float)collagenBeads)*mass);
	float intraMineralizationPercent = 100.0f*((float)HApBeadsIntra)*massMineral / (((float)HApBeadsIntra)*massMineral + ((float)collagenBeads)*mass);
	float extraMineralizationPercent = 100.0f*((float)HApBeadsExtra)*massMineral / (((float)HApBeads)*massMineral + ((float)collagenBeads)*mass);

	numberOfBonds = bonds.size();
	numberOfAngles = angles.size();
	dataFile << "LAMMPS data file for CG model for mineralized collagen (" <<
		mineralizationPercent << "% Mineral , " <<
		intraMineralizationPercent << " % Intrafibrillar , "<<
		extraMineralizationPercent << " % Extrafibrillar ) written by Mahdi Tavakol (mahditavakol90@gmail.com)" << endl;
	cout << intraMineralizationPercent << "% Intrafibrillar and " <<
		extraMineralizationPercent << "% Extrafibrillar" << endl << "Are you happy?" << endl;
	char input;
	cin >> input;
	if (input == 'n' || input == 'N')
		return;
	cout << "Finishing up!!" << endl;
	dataFile << endl;
	dataFile << numberOfAtoms << " atoms" << endl;
	dataFile << numberOfBonds << " bonds" << endl;
	dataFile << numberOfAngles << " angles" << endl;
	dataFile << "0 dihedrals" << endl;
	dataFile << "0 impropers" << endl;
	dataFile << endl;
	int numTypes = 1;
	if (mineralized)
	{
		if (type == 1) numTypes = 2;
		else if (type == 2) numTypes = 5;
	}
	else
	{
		if (type == 1) numTypes = 1;
		else if (type == 2) numTypes = 4;
	}
	if (extraFibril) numTypes++;
	dataFile << numTypes << " atom types" << endl;
	dataFile << "1 bond types" << endl;
    dataFile << angleValues.size() << " angle types" << endl << endl;
	dataFile << xlo << " " << xhi << " xlo xhi" << endl
		<< ylo << " " << yhi << " ylo yhi" << endl
		<< zlo << " " << zhi << " zlo zhi" << endl << endl;
	dataFile << "Masses" << endl << endl;

	for (int i = 0; i < numTypes - 2; i++)
		dataFile << "  " << i + 1 << " " << mass << endl;
	if (extraFibril)
	{
		dataFile << "  " << numTypes - 1 << " " << massMineral << endl;
	}
	else
	{
		dataFile << "  " << numTypes - 1 << " " << mass << endl;
	}
	if (mineralized)
	{
		dataFile << "  " << numTypes << " " << massMineral << endl;
	}
	else
	{
		dataFile << "  " << numTypes << " " << mass << endl;
	}
	dataFile << endl;
	dataFile << "Bond Coeffs" << endl << endl;
	dataFile << "  1 8.565 14.00 62.780 18.20 00.000 21.00 21.00" << endl << endl;
	dataFile << "Angle Coeffs" << endl << endl;
	for (unsigned int i = 0; i < angleValues.size(); i++)
	{
			dataFile << "  " << i + 1 << " 7.49 " << angleValues[i] << endl;
	}
	dataFile << endl;
	dataFile << "Atoms" << endl << endl;

	for (int i = 0; i < numberOfAtoms; i++)
	{
		int molecule_tag = atoms[i].GetResidueID();
		int mineralType;
		int atom_type = 1;
		float q = 0.0f;
		float coor[3];
		atoms[i].GetCoordinates(coor);
		if (type == 1)
			mineralType = 2;
		else if (type == 2)
			mineralType = 5;

		if (strcmp(atoms[i].GetAtomName().c_str(), "HA") == 0 || strcmp(atoms[i].GetAtomName().c_str(), "HA1") == 0) atom_type = mineralType;
		if (strcmp(atoms[i].GetAtomName().c_str(), "EA") == 0 || strcmp(atoms[i].GetAtomName().c_str(), "EA1") == 0) atom_type = mineralType + 1;
		//else if (type == 2) atoms[i].GetTypeNumber(atom_type);
		//else if (type == 1) atom_type = 1;
		dataFile << "  " << i + 1 << " " << molecule_tag << " " << atom_type << " " << fixed << setprecision(2) << q;
		if (fileType == 1)
		dataFile << fixed << setprecision(6) << " " << 10.0f*coor[0] << " " << 10.0f*coor[1] << " " << 10.0f*coor[2] << endl;
		else if (fileType == 2)
		dataFile << fixed << setprecision(6) << " " << 10.0f*coor[0] << " " << 10.0f*coor[1] << " " << 10.0f*coor[2] << endl;
	}
	dataFile << endl;

	dataFile << "Bonds" << endl << endl;
	for (int i = 0; i < numberOfBonds; i++)
	{
		int atom1, atom2;
		bonds[i].GetAtoms(atom1, atom2);
		dataFile << "  " << i+1 << " 1 " << atom1+1 << " " << atom2+1  << endl;
	}
	dataFile << endl;
	dataFile << "Angles" << endl << endl;
	for (int i = 0; i < numberOfAngles; i++)
	{
		int atom1, atom2, atom3, type;
		type = 1;
		angles[i].GetAtoms(atom1, atom2,atom3);
		angles[i].GetType(type);
		dataFile << "  " << i + 1 << " "<< type << " " << atom1 + 1 << " " << atom2 + 1 << " " << atom3 + 1 << endl;
	}

	//for (unsigned int i = 0; i < BLengths.size(); i++)
		//bondLengthFile << BLengths[i] << endl;

	//isThereCrossLink(dataFileName, atoms, numberOfBonds);
}

float calculateAngle(float *a, float *b, float *c)
{
	float r1 = sqrt(pow((a[0] - b[0]), 2) + pow((a[1] - b[1]), 2) + pow((a[2] - b[2]), 2));
	float r2 = sqrt(pow((c[0] - b[0]), 2) + pow((c[1] - b[1]), 2) + pow((c[2] - b[2]), 2));
	float r3 = sqrt(pow((a[0] - c[0]), 2) + pow((a[1] - c[1]), 2) + pow((a[2] - c[2]), 2));
	//float cos = (r1*r1 + r2*r2 - r3*r3) / (2 * r1*r2);
	float r1Dotr2 = (a[0] - b[0])*(c[0] - b[0]) + (a[1] - b[1])*(c[1] - b[1]) + (a[2] - b[2])*(c[2] - b[2]);
	float r1r2 = r1*r2;
	float cos = r1Dotr2 / r1r2;
	float angle = 180*acos(cos)/3.141592f;
	float roundedAngle = roundf(angle * 1) / 1;
	return roundedAngle;
}

void isThereCrossLink(const string& dataFileName, vector<Atom>& atoms, int bonds)
{
	ifstream dataFile(dataFileName);
	string line;
	int atom1, atom2, id, type;
	int atom1Type, atom2Type;
	while (getline(dataFile, line))
	{
		if (line.size() > 4)
			if (strcmp(line.substr(0, 5).c_str(),"Bonds") == 0)
				break;
	}
	getline(dataFile, line);
	for (int i = 0; i < bonds; i++)
	{
		getline(dataFile, line);
		stringstream ss(line);
		ss >> id >> type >> atom1 >> atom2;
		atom1--;
		atom2--;
		atom1Type = atoms[atom1].GetResidueID();
		atom2Type = atoms[atom2].GetResidueID();
		if (i == 50)
		{
			cout << atom1Type << "==" << atom2Type << endl;
		}
		if (atom1Type != atom2Type)
		{
			cout << "Crossline " << endl << endl << endl << "FUCK!!" << endl;
			break;
		}
	}
}