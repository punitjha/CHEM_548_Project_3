//**************************************************************************
//										               *
// To compile this program run:                                            *
//g++  -o RunFileName project_2.cpp -01 -larmadillo llapack -lblas  -Wall  * 
//									   			   *
//the -Wall attribute can be used to show all the is used to show all the  *
//errors you may possible get in your score code while runnig              *
//                							   		   *						
//**************************************************************************


#include <armadillo>   	     //the Armadillo linear algebra package of c++
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>          //for using complex numbers in the program 						
using namespace std;
int main()
{

//**************************************************************************
//												   *								
//this part of the program construts a generalized Hamiltonian of an 	   *
//molecule specified in the file name seciton building the Huckel          *
//matrix using the Armadillo Package.                                      *                                                             
//								           			   *					
//**************************************************************************


	const double pi = 3.1415926535897;			
	string filename;
	int atoms;							//storage for no. of atoms
	arma::vec k1=arma::linspace(-pi, pi, 100);	//wavefactor k1
	int k1_len=k1.n_elem;					//k1 and k2 len same
	arma::vec k2=arma::linspace(-pi, pi, 100);	//wavefactor k2
	complex <double> I(0.0,1.0);				//complex no. I  
	arma::mat connect;
	ifstream myfile("benzene"); 			//molcule connectivity file name
	if (myfile.is_open())
	{
		string str1; 					//gets the first lin
		getline(myfile,str1);				
		istringstream ee(str1);				//istring - extract words from the line
		ee >> atoms; 					//no. of atoms is assigned here
		connect=arma::zeros(atoms,3);			//connectivity matrix-- has 3 max
		string str;						//str for use in getline below
		int row=0;
		while(getline(myfile,str))
		{	
			istringstream ss(str);
			int num;
			int col=0;
			while(ss >> num)
			{
				connect(row,col)=num; 		//assigning the connectivity matrix 
				col++; 				// its done column wise
			}
			row++;
		}
	}
	myfile.close();
	arma::mat Huckel= arma::zeros(atoms,atoms); 	//Huckel matrix for TPA as semiconductor
	arma::mat Huckel2= arma::zeros(atoms,atoms); 	//H matrix  in electron-volts.
	for (int i=0; i<atoms; i++)
	{
		int n1=connect(i,0)-1;				//extracting the elements 
		int n2=connect(i,1)-1;
		int n3=connect(i,2)-1; 				//n3=0 for linear hydro-carbons
		Huckel(n1,n2)=-1;
		Huckel(n2,n1)=-1;
		Huckel2(n1,n2)=2.94;
		Huckel2(n2,n1)=2.94; 
		if (n3 > 0) 					//if hydorcarbon has rings
		{
			Huckel(n1,n3)=-1;
			Huckel(n3,n1)=-1;
			Huckel2(n1,n3)=2.94;
			Huckel2(n1,n3)=2.94;
		}
	}	
	Huckel2.diag().fill(5.94);				//diagonal elements = 5.94eV 
	arma::vec eigenvalues;
	arma::mat eigenvectors;
	arma::eig_sym(eigenvalues, eigenvectors, Huckel); 	//eigenvalues and the eigenvectors matrix
	eigenvalues.print();
//	eigenvectors.print();
	arma::vec eigenvalues2;
	arma::mat eigenvectors2;
	arma::eig_sym(eigenvalues2, eigenvectors2, Huckel2);
	eigenvalues2.print();
//	eigenvectors2.print();
//	arma::cx_mat HB= arma::zeros<arma::cx_mat>(4,4);
	


//****************************************************************************
//												     *								 
//(2-1) This part of the program is on the 1D Huckel band structure and	     *
// applies it to the pi band of polyacetelyle showing its metallic behaviour *
//												     *
//****************************************************************************



	arma::cx_mat dHuckel(2,2);				//Huckel matrix/Hamiltonian
	arma::cx_vec eigenvalues22;				//vector of eigenvalues
	arma::cx_mat eigenvectors22;				//matrix-eigenvectors
	fstream myfile1("bandgap.txt",fstream::out | fstream::trunc);
	for (int i=0; i<k1_len; i++)	 			//for loop initializes-Huckel Hamiltonian
	{
		dHuckel(0,0)=0;
		dHuckel(0,1)=-1.1-0.9*exp(-I*k1(i));
		dHuckel(1,0)=-1.1-0.9*exp(I*k1(i));
		dHuckel(1,1)=0;
		arma::eig_gen(eigenvalues22, eigenvectors22, dHuckel);
		myfile1<<arma::real(eigenvalues22);
	}




//*******************************************************************************
//											         	  *								
//(2-2) This part of the program is on the where we have metallic behaviour of  *
// transpolyacetelene due to equidistatn C-C bond                               *
//													  *					
//*******************************************************************************




	arma::cx_mat ddHuckel(2,2);					//Huckel matrix
	arma::cx_vec eigenvalues222;					//vec-eigenvalues
	arma::cx_mat eigenvectors222;					//matrix-eigenvectors
	fstream myfile22("bandgap1.txt", fstream::out|fstream::trunc);	
	for(int i=0; i<k1_len; i++)					//initializing -Huckel matrix
	{
		ddHuckel(0,0)=0;
		ddHuckel(0,1)=(-exp(-I*k1(i)))-1.0;
		ddHuckel(1,0)=-1.0-exp(I*k1(i));
		ddHuckel(1,1)=0;
		arma::eig_gen(eigenvalues222,eigenvectors222,ddHuckel); //calculating the eigenvalues 
		myfile22<<arma::real(eigenvalues222);                   //writing the eigenvalues to a file 
	}





//******************************************************************************
//									       			 *						 
// this part of the program calculates the band of the 4 walled carbon nanotube*                                                                  
//									       			 *					 
//****************************************************************************** 

	int carbon=4;							//number of carbon atoms
	fstream myfile33("carbon.txt", fstream::out|fstream::trunc);
	arma::cx_mat DHuckel(carbon,carbon);		 	//Hamiltonian matrix-- sum of 3 matrices
	arma::cx_mat mat_1(carbon,carbon);				//first matrix in the H
	arma::cx_mat mat_2(carbon,carbon);				//second matrix in the H--has anti-diagonal
	arma::cx_mat mat_3(carbon,carbon);				//third  matrix in the H
	arma::cx_vec eigencarbon;					//eigenvalues 
	arma::cx_mat eigenvec_carbons;				//eigenvectors
	int diag_marker =3;						//mark anti-diagonal matrix elements
	mat_1(carbon-1,0)=-1;						//lowermost corner element as -1
	mat_3(0,carbon-1)=-1;						//topmost corner element as -1
	for(int i=0; i<4; i++)
	{
		for (int j=0; j<4; j++)					//initializing mat_1 and mat_2
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
		DHuckel=mat_1*exp(-I*k1(i))+mat_2+mat_3*exp(I*k1(i));	      //H matrix as the sum of mat_1, mat_2, mat_3
		arma::eig_gen(eigencarbon,eigenvec_carbons,DHuckel);	      //eigenvalues and the eigenvectors
		arma::vec eigencarbon1=arma::real(eigencarbon);		      //extracting the real parts
		arma::vec eigencarbon2=arma::sort(eigencarbon1);	      //sorting the real parts
		for(int j=0; j<4;j++)
		{
			myfile33<<eigencarbon2(j)<<"	";		            //writing to file
		}
		myfile33<<endl;
	}


//******************************************************************************
//												       *								    
// this part of the program calculates the energy bands of an wider single     *
//walled nano tube which has 10 carbons on the circumference		       *
//									     			       *								    
//******************************************************************************



	int carbon10=10;						//number of carbon atoms
	fstream myfile310("carbon1.txt", fstream::out|fstream::trunc);
	arma::cx_mat Huckel_10(carbon10,carbon10);			//Hamiltonian matrix-- sum of 3
	arma::cx_mat mat_10(carbon10,carbon10);				//first matrix in the H
	arma::cx_mat mat_20(carbon10,carbon10);				//second matrix in the H--has a
	arma::cx_mat mat_30(carbon10,carbon10);				//third  matrix in the H
	arma::cx_vec eigencarbon10;						//eigenvalues 
	arma::cx_mat eigenvec_carbons10;					//eigenvectors
	int diag_marker10 =9;							//mark anti-diagonal matrix ele
	mat_10(carbon10-1,0)=-1;						//lowermost corner element as -
	mat_30(0,carbon10-1)=-1;						//topmost corner element as -1 	
	for(int i=0; i<10; i++)
	{
		for (int j=0; j<10; j++)				      //matrix elements-- mat_1 and mat_2
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
	for (int i=0; i<k1_len; i++)
	{	
		Huckel_10=mat_10*exp(-I*k1(i))+mat_20+mat_30*exp(I*k1(i));  //H matrix sum of mat_1, mat_2, mat_3
		arma::eig_gen(eigencarbon10,eigenvec_carbons10,Huckel_10);  //eigenvalues and the eigenvectors
		arma::vec eigencarbon1=arma::real(eigencarbon10);	      //extract real parts
		arma::vec eigencarbon2=arma::sort(eigencarbon1);	      //sort real part 
		for (int j=0; j<10 ; j++)
		{
			myfile310<<eigencarbon2(j)<<'\t';			      //eigenvalues to file 
		}
		myfile310<<endl;							      //line-beaking file
	}
	



//*****************************************************************************
//												      *
// 2D-Huckel band calculation of pi band of graphene				      *								
//									      			*							    
//*****************************************************************************

	arma::cx_mat H_total(2,2);
	arma::cx_mat H_1(2,2);					//H matrix which has -1 at (1,0)
	arma::cx_mat H_2(2,2);					//H matrix which has -1 at (0,1)
	arma::cx_mat H_antidiag(2,2);				//anti-diagonal matrix
	fstream myfile333("H2D.txt", fstream::out|fstream::trunc);
	H_1(0,0)=0;							//intializing H_1
	H_1(0,1)=0;
	H_1(1,0)=-1;
	H_1(1,1)=0;
	H_2(0,0)=0;
	H_2(0,1)=-1;
	H_2(1,0)=0;
	H_2(1,1)=0;
	H_antidiag(0,0)=0;
	H_antidiag(0,1)=-1;
	H_antidiag(1,0)=-1;
	H_antidiag(1,1)=0;
	for (int i=0; i<k1_len; i++)
	{
		for(int j=0; j<k1_len; j++)
		{
			H_total=H_1*exp(-I*k1(i))+H_2*exp(-I*k2(j))+H_antidiag+H_2*exp(I*k1(i))+H_1*exp(I*k2(j));
			arma::cx_vec eigenval33;
	    		arma::cx_mat eigenvec33;
			arma::eig_gen(eigenval33,eigenvec33,H_total);
			arma::vec rpart= arma::real(eigenval33);	
			myfile333<<k1(i)<<'\t'<<k2(j)<<'\t'<<rpart(0)<<'\t'<<rpart(1)<<endl;
		}
	}		

//******************************************************************************
// The Su-Schrieffer-Heeker model for benzene and benene cation and	       *
// Jahn-Teller effect							       *	
//****************************************************************************** 
	
	int atoms_b=0;
	arma::mat connect_b;
	ifstream myfile44("benzene"); 			//molecule connectivity file name
	if (myfile44.is_open())
	{
		string str11; 					//gets the first lin
		getline(myfile44,str11);				
		istringstream ee(str11);			//istring - extract words from the line
		ee >> atoms_b; 					//no. of atoms is assigned here
		connect_b=arma::zeros(atoms_b,3);			//connectivity matrix-- has 3 max
		string str;						//str for use in getline below
		int row=0;
		while(getline(myfile44,str))
		{	
			istringstream ss(str);
			int num;
			int col=0;
			while(ss >> num)
			{
				connect_b(row,col)=num; 		//assigning the connectivity matrix 
				col++; 				// its done column wise
			}
			row++;
		}
	}
	myfile44.close();
	arma::mat Huckel_b= arma::zeros(atoms_b,atoms_b); 	//Huckel matrix for TPA as semiconductor
	for (int i=0; i<atoms_b; i++)
	{
		int n1=connect_b(i,0)-1;				//extracting the elements 
		int n2=connect_b(i,1)-1;
		int n3=connect_b(i,2)-1; 				//n3=0 for linear hydro-carbons
		Huckel_b(n1,n2)=-1;
		Huckel_b(n2,n1)=-1;
		if (n3 > 0) 					//if hydorcarbon has rings
		{
			Huckel_b(n1,n3)=-1;
			Huckel_b(n3,n1)=-1;
		}
	}
	arma::vec delta=arma::linspace(-0.25,0.25,400);
	int delta_len=delta.n_elem;				//k1 and k2 len same
	arma::mat eigen_b=arma::zeros(delta_len,6);	//eigenvalue for diffrent 1
	for(int i=0; i<delta_len; i++)
	{	
		for(int row=0; row <6; row++)
		{
			for(int col =0; col<6; col++)
			{
				if (Huckel_b(row,col) !=0)
			      {
					if((row+col != 3) && (row+col != 9)) 			//look at H 
					{
						  Huckel_b(row,col)=-1.0+delta(i);
					}
					else
						  Huckel_b(row,col)=-1.0-2.0*delta(i);
				}
		      }
		}
		arma::vec eg_b;
		arma::mat ev_b;
		arma::eig_sym(eg_b, ev_b, Huckel_b); 	//eigenvalues
		//arma::vec eg_b1=arma::sort(eg_b);
		eigen_b.row(i)+=eg_b.t();
//	  	eigen_b.print();
//		cout<<endl;
	}
	//eigen_b.print();
	//delta.print();
	fstream myfile_11("ben_sym.txt", fstream::out|fstream::trunc);
	arma::vec ee_b=arma::zeros(delta_len);
	arma::vec ee_bc=arma::zeros(delta_len);
	for(int i=0; i<delta_len;i++)
	{
		ee_b(i)=2.0*eigen_b(i,0)+2.0*eigen_b(i,1)+2.0*eigen_b(i,2)+(0.5*5.0*(12*delta(i)*delta(i)));
		ee_bc(i)=2.0*eigen_b(i,0)+2.0*eigen_b(i,1)+eigen_b(i,2)+(0.5*5.0*(12*delta(i)*delta(i)));
		myfile_11<<delta(i)<<'\t'<<ee_b(i)<<'\t'<<ee_bc(i)<<endl;
	}


//*******************************************************************************
// The Su-Schrieffer-Heeker model with bond length variations for polyacetylene *
//and study of the Peierls theorem						*
//*******************************************************************************


	arma::cx_mat Hck_p(2,2);					//Huckel matrix
	arma::cx_vec ee_p;						//vec-eigenvalues
	arma::cx_mat ev_p;						//matrix-eigenvectors
	fstream myfile_p("poly.txt", fstream::out|fstream::trunc);
	arma::vec delta_1=arma::linspace(-0.25,0.25,100);
	int delta_1_len=delta_1.n_elem;
	arma::vec ee_poly=arma::zeros(delta_1_len);
	arma::vec del_k=arma::zeros(delta_1_len);
//	int marker=0;
	for(int j=0; j<delta_1_len;j++)
	{	
		double product=0.0;
		double sum=0.0;
		for(int i=0; i<k1_len; i++)					//initializing -Huckel matrix
		{
			Hck_p(0,0)=0.0;
			Hck_p(0,1)=((-1.0-delta_1(j))*exp(-I*k1(i)))+(-1.0+delta_1(j));
			Hck_p(1,0)=((-1.0-delta_1(j))*exp(I*k1(i)))+(-1.0+delta_1(j));
		//	cout<<Hck_p(1,0)<<endl;
			Hck_p(1,1)=0.0;
			arma::eig_gen(ee_p,ev_p,Hck_p); 			//calculating the eigenvalues 
			arma::vec eigen_ll=arma::real(ee_p);
			sum=sum+eigen_ll(1);
			product=product+0.5*k1(i)*2.0*delta_1(j)*delta_1(j);
			cout<<"this is the product individually"<<product<<endl;
//			ee_poly.row(marker)+=eigen_ll.t();
//			myfile_p<<eigen_ll(0)<<'\t'<<eigen_ll(1)<<endl;                    
//			marker++;
		}
	ee_poly(j)=sum;
	del_k(j)=product;
	cout<<ee_poly(j)<<'\t'<<del_k(j)<<endl;
	}
	//ee_poly.print();
	//cout<<endl;
	//del_k.print();
	arma::vec ee_hck=arma::zeros(delta_1_len);
	for(int i=0; i<delta_1_len; i++)
	{
		ee_hck(i)=(1/k1_len)*(ee_poly(i)+del_k(i));
		myfile_p<<delta_1(i)<<'\t'<<ee_hck(i)<<endl;
	}
	double a=-1+0.25;
	complex	<double> b=exp(I*3.14);
	double c=(-1.0-0.25);
	cout<<a<<b<<c<<endl;
	cout<<a*b+c;


















}	
