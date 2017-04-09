// to compile this program run  g++  -o RunFileName project_2.cpp -01 -larmadillo llapack -lblas  -Wall 
// the -Wall attribute can be used to show all the is used to show all the errors you may possible get in your score code while runnig
#include <armadillo> // the Armadillo linear algebra package of c++
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <string>
#include <complex> //for using complex numbers in the program
using namespace std;
int main()
{
//building the Huckel matrix using the Armadillo Package
	string filename;
	string line;
	int atoms=6;
	cout<<"Please enter the name of the molecule file"<<filename<<endl;
	cin>>filename;
	ifstream myfile(filename.c_str());
	if (myfile.is_open())
	{
		while (getline (myfile,line))
		{
			cout<<line<<endl;
		}
	myfile.close();
	}	
	arma::mat Huckel= arma::zeros(atoms,atoms);
	arma::vec subdiagonal= arma::zeros(atoms-1);
	subdiagonal.fill(-1);
	cout<<subdiagonal<<endl;
	Huckel.diag(-1)=subdiagonal;
	Huckel.diag(1)=subdiagonal;
	Huckel(5,0)=-1;
	Huckel(0,5)=-1;
	Huckel.print();
	arma::vec eigenvalues;
	arma::mat eigenvectors;
	arma::eig_sym(eigenvalues, eigenvectors, Huckel);
	eigenvalues.print();
	complex <double> I(0.0,1.0);
	arma::cx_mat HB= arma::zeros<arma::cx_mat>(4,4);

}	
