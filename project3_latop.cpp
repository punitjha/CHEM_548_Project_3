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
//building the Huckel matrix using the Armadillo Package
	string filename;
	int atoms;
	arma::mat connect;
//	cout<<"Please enter the name of the molecule file"<<filename<<endl;
//	cin>>filename;
	ifstream myfile("benzene");
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
				cout<<connect(row,col)<<endl;
				col++;
			}
			row++;
		}
	}
	myfile.close();
	connect.print();
	arma::mat Huckel= arma::zeros(atoms,atoms);
	for (int i=0; i<atoms; i++)
        {
		for (int j=0;j<3;j++)
		{
			if (connect(i,j) != 0)
			{
				Huckel()
		}
	}	

	Huckel.print();
	arma::vec eigenvalues;
	arma::mat eigenvectors;
	arma::eig_sym(eigenvalues, eigenvectors, Huckel);
	eigenvalues.print();
	complex <double> I(0.0,1.0);
	arma::cx_mat HB= arma::zeros<arma::cx_mat>(4,4);

}	
