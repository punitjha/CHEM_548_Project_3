//***************************************************************************************************************************************
//																        *
// to compile this program run  g++  -o RunFileName project_2.cpp -01 -larmadillo llapack -lblas  -Wall                                 *
// the -Wall attribute can be used to show all the is used to show all the errors you may possible get in your score code while runnig  *
//                															*
//***************************************************************************************************************************************

#include <armadillo> 			      			//the Armadillo linear algebra package of c++
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <string>
#include <sstream>
#include <complex> 						//for using complex numbers in the program
using namespace std;
int main()
{

//***************************************************************************************************************************************
//																	*
// this part of the program construts a generalized Hamiltonian of an molecule specified in the file name seciton building the Huckel   *
//matrix using the Armadillo Package.                                                                                                   *
//																	*
//***************************************************************************************************************************************
	
	const double pi = 3.1415926535897;			//defining the constant pi
	string filename;
	int atoms;						//storage for number of atoms
	arma::vec k1=arma::linspace(-pi/2, pi/2, 600);		//wavefactor k1 defined as a vector
	int k1_len=k1.n_elem;
	arma::vec k2=arma::linspace(-pi, pi, 100);		//wavefactor k2 defined as a vector
	complex <double> I(0.0,1.0);				//defining the complex no. I for use on the codes
	arma::mat connect;
	ifstream myfile("anthracene"); 				//the name of the molcule is specified here as the connectivity file name
	if (myfile.is_open())
	{
		string str1; 					//this gets the first line of the file that it read
		getline(myfile,str1);				//using the getline is to read the first line
		istringstream ee(str1);				//istring is used to extract all the words from the line
		ee >> atoms; 					// no. of atoms is assigned here
		connect=arma::zeros(atoms,3);			// the connectivity matrix-- has 3 rows as any carbon can at max connect to 3 
		string str;					// this str variable is for use in getline below
		int row=0;
		while(getline(myfile,str))
		{	
			istringstream ss(str);
			int num;
			int col=0;
			while(ss >> num)
			{
				connect(row,col)=num; 		//assigning the connectivity matrix elements form the file
				col++; 				// its done column wise
			}
			row++;
		}
	}
	myfile.close();
	arma::mat Huckel= arma::zeros(atoms,atoms); 		//the Huckel matrix for trans-polyacetylene as semiconductor
	arma::mat Huckel2= arma::zeros(atoms,atoms); 		//the H matrix for trans-polyacetylene alpha and beta assigned to it in electron-volts.
	for (int i=0; i<atoms; i++)
        {
		int n1=connect(i,0)-1;				//extracting the elements of the connectivity matrix
		int n2=connect(i,1)-1;
		int n3=connect(i,2)-1; 				// normally n3=0 for linear hydro-carbons
		Huckel(n1,n2)=-1;
		Huckel(n2,n1)=-1;
		Huckel2(n1,n2)=2.94;
		Huckel2(n2,n1)=2.94; 
		if (n3 > 0) 					//this means that the hydorcarbon has rings
		{
			Huckel(n1,n3)=-1;
			Huckel(n3,n1)=-1;
			Huckel2(n1,n3)=2.94;
			Huckel2(n1,n3)=2.94;
		}
	}	
	Huckel2.diag().fill(5.94);				//filling the diagonal matrix elements with the value of 5.94eV 
	arma::vec eigenvalues;
	arma::mat eigenvectors;
	arma::eig_sym(eigenvalues, eigenvectors, Huckel); 	//solving for the eigenvalues and the eigenvectors of the constructed Huckel matrix
	arma::vec eigenvalues2;
	arma::mat eigenvectors2;
	arma::eig_sym(eigenvalues2, eigenvectors2, Huckel2);
	arma::cx_mat HB= arma::zeros<arma::cx_mat>(4,4);




//************************************************************************************************************************************************
//																		 *
//(2-1) This part of the program is on the 1D Huckel band structure and applies it to the pi band of polyacetelyle showing its metallic behaviour*
//																		 *
//************************************************************************************************************************************************




	arma::cx_mat dHuckel(2,2);				//declaring the Huckel matrix/Hamiltonian
	arma::cx_vec eigenvalues22;				//the vector which will store the eigenvalues
	arma::cx_mat eigenvectors22;				// the matrix wich will store the eigenvectors
	fstream myfile1("bandgap.txt",fstream::out | fstream::trunc);
	for (int i=0; i<k1_len; i++)	 			//this for loop initializes the values to the Huckel Hamiltonian
	{
		dHuckel(0,0)=0;
		dHuckel(0,1)=-1.1-0.9*exp(-I*k1(i));
		dHuckel(1,0)=-1.1-0.9*exp(I*k1(i));
		dHuckel(1,1)=0;
		arma::eig_gen(eigenvalues22, eigenvectors22, dHuckel);
		myfile1<<arma::real(eigenvalues22);
	}
		



//***********************************************************************************************************************************************
//																		*
//(2-2) This part of the program is on the where we have metallic behaviour of transpolyacetelene due to equidistatn C-C bond                   *
//																		*
//***********************************************************************************************************************************************




	arma::cx_mat ddHuckel(2,2);						//declaring the Huckel matrix
	arma::cx_vec eigenvalues222;						//the vec which will store the eigenvalues
	arma::cx_mat eigenvectors222;						// the matrix which will store the eigenvectors
        fstream myfile22("bandgap1.txt", fstream::out|fstream::trunc);	
	for(int i=0; i<k1_len; i++)						//initializing the matrix elements of the Huckel matrix
	{
		ddHuckel(0,0)=0;
		ddHuckel(0,1)=(-exp(-I*k1(i)))-1.0;
		ddHuckel(1,0)=-1.0-exp(I*k1(i));
		ddHuckel(1,1)=0;
		arma::eig_gen(eigenvalues222,eigenvectors222,ddHuckel); 	//calculating the eigenvalues 
		myfile22<<arma::real(eigenvalues222);                   	//writing the eigenvalues of the matrix to a file 
	}





//************************************************************************************************************************************************
//																		 *
// this part of the program calculates the band of the 4 walled carbon nanotube                                                                  *
//																		 *
//************************************************************************************************************************************************

	int carbon=4;								//entering the number of carbon atoms
	fstream myfile33("carbon.txt", fstream::out|fstream::trunc);
	arma::cx_mat DHuckel(carbon,carbon);					//declaring the Hamiltonian matrix-- sum of 3 other matrices
	arma::cx_mat mat_1(carbon,carbon);					//this is the first matrix in the Hamiltonian
	arma::cx_mat mat_2(carbon,carbon);					//this is the second matrix in the Hamiltonian--has anti-diagonal terms
	arma::cx_mat mat_3(carbon,carbon);					//this is the third  matrix in the Hamiltonian
	arma::cx_vec eigencarbon;						//this strores the eigenvalues of the matrix obatined
	arma::cx_mat eigenvec_carbons;						//thsi stores the eigenvectors of the Hamiltonian Martix
	int diag_marker =3;							//this is used to create the anti-diagonal matrix elements
	mat_1(carbon-1,0)=-1;							//assigning the lowermost corner element as -1
	mat_3(0,carbon-1)=-1;							//assigning the topmost corner element as -1
	for(int i=0; i<4; i++)
	{
		for (int j=0; j<4; j++)						//initializing the matrix elements of the mat_1 and mat_2
		{
			if((j-i) == 1)
			{
				mat_1(i,j)=-1;
			}
			if ((i-j)==1)
			{
				mat_3(i,j)=-1;
			}
		}
		mat_2(i,diag_marker)=-1;
		diag_marker--;
	}
	for (int i=0; i<k1_len; i++)
	{
		DHuckel=mat_1*exp(-I*k1(i))+mat_2+mat_3*exp(I*k1(i));	      //creating the Hamiltonian matrix elements as the sum of mat_1, mat_2, mat_3
		arma::eig_gen(eigencarbon,eigenvec_carbons,DHuckel);	      //solving for the eigenvalues and the eigenvectors
		arma::vec eigencarbon1=arma::real(eigencarbon);		      //extracting the real parts of the eigenvalues
		arma::vec eigencarbon2=arma::sort(eigencarbon1);	      //sorting the real part of the eigenvalues
		myfile33<<real(eigencarbon2);				      //writing to file
	}


//***************************************************************************************************************************************************
//																		    *
// this part of the program calculates the energy bands of an wider single walled nano tube which has 10 carbons on the circumference		    *
//																		    *
//***************************************************************************************************************************************************



	int carbon10=10;							//entering the number of carbon atoms
	fstream myfile310("carbon1.txt", fstream::out|fstream::trunc);
//	fstream myfile333("checking.xlxs",fstream::out|fstream::trunc);
	arma::cx_mat Huckel_10(carbon10,carbon10);				//declaring the Hamiltonian matrix-- sum of 3 other matrices
	arma::cx_mat mat_10(carbon10,carbon10);					//this is the first matrix in the Hamiltonian
	arma::cx_mat mat_20(carbon10,carbon10);					//this is the second matrix in the Hamiltonian--has anti-diagonal terms
	arma::cx_mat mat_30(carbon10,carbon10);					//this is the third  matrix in the Hamiltonian
	arma::cx_vec eigencarbon10;						//this strores the eigenvalues of the matrix obatined
	arma::cx_mat eigenvec_carbons10;					//thsi stores the eigenvectors of the Hamiltonian Martix
	int diag_marker10 =9;							//this is used to create the anti-diagonal matrix elements(NOTE 9)
	mat_10(carbon10-1,0)=-1;						//assigning the lowermost corner element as -1
	mat_30(0,carbon10-1)=-1;						//assigning the topmost corner element as -1
	for(int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)					//initializing the matrix elements of the mat_1 and mat_2
		{
			if((j-i) == 1)
			{
				mat_10(i,j)=-1;
			}
			if ((i-j)==1)
			{
				mat_30(i,j)=-1;
			}
		}
		mat_20(i,diag_marker10)=-1;
		diag_marker10--;
	}
	mat_20.save("checking.xlxs",arma::arma_ascii);
	for (int i=0; i<k1_len; i++)
	{	
		Huckel_10=mat_10*exp(-I*k1(i))+mat_20+mat_30*exp(I*k1(i));       //creating the Hamiltonian matrix elements as the sum of mat_1, mat_2, mat_3
		arma::eig_gen(eigencarbon10,eigenvec_carbons10,Huckel_10);	 //solving for the eigenvalues and the eigenvectors
		arma::vec eigencarbon1=arma::real(eigencarbon10);		 //extracting the real parts of the eigenvalues
		arma::vec eigencarbon2=arma::sort(eigencarbon1);	         //sorting the real part of the eigenvalues
		for (int j=0; j<10 ; j++)
		{
			myfile310<<eigencarbon2(j)<<"	";			//writing the sorting eigenvalues to the file element wise
		}
		myfile310<<endl;						//line-beaking in the file
	}





//***************************************************************************************************************************************************
//																		    *
// 2D-Huckel band calculation of pi band of graphene												    *
//																		    *
//***************************************************************************************************************************************************

	arma::cx_mat H_k1(2,2); 						//this is the k1 dependent Hamiltonian matrix
	arma::cx_mat H_k2(2,2);							//this is the k2 dependent Hamiltonian matrix
	


















}	
