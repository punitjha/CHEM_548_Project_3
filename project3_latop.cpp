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
	arma::mat connect(atoms,3);
//	cout<<"Please enter the name of the molecule file"<<filename<<endl;
//	cin>>filename;
	ifstream myfile("benzene");
	if (myfile.is_open())
	{
		string str1;
		getline(myfile,str1);
		istringstream ee(str1);
		ee >> atoms;
		string str;
		int row=0;
		while(getline(myfile,str))
		{	
			cout<<"this is the row"<<row<<endl;
			istringstream ss(str);
			int num;
			int col=0;
			while(ss >> num)
			{
				cout<<"this is the col"<<col<<endl;
				cout<<"this is num"<<num<<endl;
				connect(row,col)=num;
				cout<<connect(row,col)<<endl;
				col++;
			}
			row++;
		}
	}
//			for (int i=0; i<atoms; i++)
//			{	
//				for(int j=0; j<3;j++)//the loop variable for the connectivity file is =2
//				{
//					myfile>>connect(i,j);
//					cout<<"the elements of the matrix"<<i<<j<<connect(i,j)<<endl;
//				}
//			}
//		}	
	cout<<"this is the atom value"<<atoms<<endl;
	myfile.close();
//	connect.print();
	arma::mat Huckel= arma::zeros(atoms,atoms);
//	arma::vec subdiagonal= arma::zeros(atoms-1);
//	subdiagonal.fill(-1);
//	cout<<subdiagonal<<endl;
//	Huckel.diag(-1)=subdiagonal;
//	Huckel.diag(1)=subdiagonal;
//	Huckel(5,0)=-1;
//	Huckel(0,5)=-1;
	for (int i=0; i<atoms; i++)
        {
		for (int j=0;j<3;j++)
		{
			cout<<"this";
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
