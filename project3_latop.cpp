// to compile this program run  g++  -o RunFileName project_2.cpp -01 -larmadillo llapack -lblas  -Wall 
// the -Wall attribute can be used to show all the is used to show all the errors you may possible get in your score code while runnig
#include <armadillo> // the Armadillo linear algebra package of c++
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <string>
#include <sstream>
#include <complex> //for using complex numbers in the program
using namespace std;
int main()
{



//****************************************************************************************************************
// this part of the program construts a generalized Hamiltonian of an molecule specified in the file name seciton
//building the Huckel matrix using the Armadillo Package.
//****************************************************************************************************************
	const double pi = 3.1415926535897;
	string filename;
	int atoms;
	arma::mat connect;
	ifstream myfile("anthracene"); //the name of the molcule is specified here as the connectivity file name
	if (myfile.is_open())
	{
		string str1;
		getline(myfile,str1);
		istringstream ee(str1);
		ee >> atoms;
		connect=arma::zeros(atoms,3);
		string str;
		int row=0;
		while(getline(myfile,str))
		{	
			istringstream ss(str);
			int num;
			int col=0;
			while(ss >> num)
			{
				connect(row,col)=num;
				col++;
			}
			row++;
		}
	}
	myfile.close();
	arma::mat Huckel= arma::zeros(atoms,atoms);
	arma::mat Huckel2= arma::zeros(atoms,atoms); 
	for (int i=0; i<atoms; i++)
        {
		int n1=connect(i,0)-1;
		int n2=connect(i,1)-1;
		int n3=connect(i,2)-1;
		Huckel(n1,n2)=-1;
		Huckel(n2,n1)=-1;
		Huckel2(n1,n2)=2.94;
		Huckel2(n2,n1)=2.94; 
		if (n3 > 0)
		{
			Huckel(n1,n3)=-1;
			Huckel(n3,n1)=-1;
			Huckel2(n1,n3)=2.94;
			Huckel2(n1,n3)=2.94;
		}
	}	
	Huckel2.diag().fill(5.94);
	arma::vec eigenvalues;
	arma::mat eigenvectors;
	arma::eig_sym(eigenvalues, eigenvectors, Huckel);
	cout<<"Printing the eigenvalues of the Rounded-off Hamiltonian"<<endl;
	eigenvalues.print();
	arma::vec eigenvalues2;
	arma::mat eigenvectors2;
	arma::eig_sym(eigenvalues2, eigenvectors2, Huckel2);
	cout<<"Printing the orbital energies in eV"<<endl;
	eigenvalues2.print();
	complex <double> I(0.0,1.0);
	arma::cx_mat HB= arma::zeros<arma::cx_mat>(4,4);




//***************************************************************************************************************************************************
//(2-1) This part of the program is on the 1D Huckel band structure and applies it to the pi band of polyacetelyle and single walled carbon nanotube
//***************************************************************************************************************************************************

	arma::cx_mat dHuckel(2,2);
	arma::cx_vec eigenvalues22;
	arma::cx_mat eigenvectors22;
	fstream myfile1("bandgap.txt",fstream::out | fstream::trunc);
	for (double k=-pi;k<=pi;k=k+0.01)
	{
		dHuckel(0,0)=0;
		dHuckel(0,1)=-1.1-0.9*exp(-I*k);
		dHuckel(1,0)=-1.1-0.9*exp(I*k);
		dHuckel(1,1)=0;
		arma::eig_gen(eigenvalues22, eigenvectors22, dHuckel);
		arma::real(eigenvalues22).print();
		myfile1<<arma::real(eigenvalues22);
	}
		



//******************************************************************************************************************************************************
//(2-2) This part of the program is on the where we have metallic behaviour of transpolyacetelene due to equidistatn C-C bond
//******************************************************************************************************************************************************


	arma::cx_mat ddHuckel(2,2);
	arma::cx_vec eigenvalues222;
	arma::cx_mat eigenvectors222;
        fstream myfile22("bandgap1.txt", fstream::out|fstream::trunc);	
	for(double k=-pi; k<=pi; k+=0.01)
	{
		ddHuckel(0,0)=0;
		ddHuckel(0,1)=(-exp(-I*k))-1.0;
		ddHuckel(1,0)=-1.0-exp(I*k);
		ddHuckel(1,1)=0;
		arma::eig_gen(eigenvalues222,eigenvectors222,ddHuckel);
		myfile22<<arma::real(eigenvalues222);
	}


//*******************************************************************************************************










}	
