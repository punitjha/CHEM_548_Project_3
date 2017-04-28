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


	const double pi = 3.1415926535898;			
	string filename;
	int atoms;							//storage for no. of atoms
	arma::vec k1=arma::linspace(-pi, pi,100);	//wavefactor k1
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
	arma::vec eigenvalues2;
	arma::mat eigenvectors2;
	arma::eig_sym(eigenvalues2, eigenvectors2, Huckel2);
	eigenvalues2.print();
	


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
	mat_10(carbon10-1,0)=-1.0;						//lowermost corner element as -
	mat_30(0,carbon10-1)=-1.0;						//topmost corner element as -1 	
	mat_20(0,carbon10-1)=mat_20(carbon10-1,0)=-1.0;
	int diag=2;
	for(int i=0; i<10; i++)
	{
		if((diag-i)==1 && (diag!= 10))
		{
			mat_20(diag,i)=mat_20(i,diag)=-1.0;
			diag=diag+2;
		}
		for (int j=0; j<10; j++)				      //matrix elements-- mat_1 and mat_2
		{
			if((j-i) == 1)
			{
				mat_10(i,j)=-1.0;
			}
			if ((i-j)==1)
			{
				mat_30(i,j)=-1.0;
			}
		}
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
		eigen_b.row(i)+=eg_b.t();
	}
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
	fstream myfile_p("poly1.txt", fstream::out|fstream::trunc);
	arma::vec delta_1=arma::linspace(-0.25,0.25,1000);
	int delta_1_len=delta_1.n_elem;
	arma::vec ee_poly=arma::zeros(delta_1_len);
	arma::mat ee_store=arma::zeros(delta_1_len,2);
	for(int j=0; j<delta_1_len;j++)
	{	
		double sum=0.0;
		for(int i=0; i<k1_len; i++)					//initializing -Huckel matrix
		{
			Hck_p(0,0)=0.0;
			Hck_p(0,1)=((-1.0-delta_1(j))*exp(-I*k1(i)))+(-1.0+delta_1(j));
			Hck_p(1,0)=((-1.0-delta_1(j))*exp(I*k1(i)))+(-1.0+delta_1(j));
			Hck_p(1,1)=0.0;
			arma::eig_gen(ee_p,ev_p,Hck_p); 			//calculating the eigenvalues 
			arma::vec eigen_ll=arma::real(ee_p);
			sum=sum+eigen_ll(1);
			if(i==0)
			{
				ee_store.row(j)+=eigen_ll.t();
			}
		}
	ee_poly(j)=sum;
	}
	arma::vec ee_hck=arma::zeros(delta_1_len);
	for(int i=0; i<delta_1_len; i++)
	{
		ee_hck(i)=((1.0/k1_len)*ee_poly(i))+0.5*2.0*delta_1(i)*delta_1(i)*2.0;
		myfile_p<<delta_1(i)<<'\t'<<ee_hck(i)<<'\t'<<(ee_poly(i)/(k1_len))<<endl;
	}


//*******************************************************************************
//										*
// Using the Extended Huckel theory to demonstrate the origin of the Jahn-Teller*
// Based on Sohlber et all, J. Chem. Educ. 2013, 90 (4), 463–469.		*
//										*	
//*******************************************************************************

	const double R=(1.39/0.529);			//the bond lenght in Bohr
	const double Z=1.56;				//orbital exponent
	const double Z_eff=1.95;			//effective nuclear charge
	arma::vec w=arma::linspace(0.02,0.5,30);	//the distortion parameter in b
//	w.print();
	int w_len=w.n_elem;
	arma::mat a_val(w_len,2);
	arma::mat b_val(w_len,2);
	arma::mat x_val(w_len,2);
	for(int i=0; i<w_len; i++)
	{
		int col=0;
		double a_ang=(2.0/6.0)*pi;
		double x1_len=0.0;
		double x2_len=0.0;
		while(a_ang>(2.0/9.0)*pi)
		{
			double b_ang=(2.0/3.0)*pi-a_ang;
			x1_len=sqrt(((R-w(i))*(R-w(i)))/(2-2*cos(a_ang)));
			x2_len=sqrt(((R+w(i))*(R+w(i)))/(2-2*cos(b_ang)));
	//		cout<<b_ang<<" "<<a_ang<<" "<<x1_len<<" "<<x2_len<<endl;
			if(abs(x1_len-x2_len) < 0.0001 )
			{
				//cout<<"this is col"<<col<<endl;
				//cout<<a_ang<<b_ang<<x1_len<<endl;
				if (col == 0)
				{
					a_val(i,col)=a_ang;
					b_val(i,col)=b_ang;
					x_val(i,col)=x1_len;
					col++;

				}
				else
				{
					if (abs((x_val(i,col-1))-(x1_len))<0.00001 )
					{
					//	cout<<x_val(i,col-1)<<" "<<x1_len<<" "<<abs((x_val(i,col-1))-(x1_len))<<endl;
						a_val(i,col)=a_ang;
						b_val(i,col)=b_ang;
						x_val(i,col)=x1_len;
						col++;
					}
				}
			}
			a_ang=a_ang-0.00001;
		}
	}
//	a_val.print("this is a");
//	b_val.print("this is b"); 
//	x_val.print("this is x");
	arma::mat c_ang=(pi-a_val)/2;
	arma::mat d_ang=(pi-b_val)/2;
//	c_ang.print("this is c");
	arma::vec r12=R-w;
	arma::vec r23=R+w;
	arma::vec r13=arma::zeros(w_len);
	for(int i=0;i<w_len;i++)
	{
	       	r13(i)=sqrt(r12(i)*r12(i)+r23(i)*r23(i)-(2*r12(i)*r23(i))*cos(c_ang(i,0)+d_ang(i,0)));
	}
//	r13.print("this is r13");
	arma::vec r14=arma::zeros(w_len);
	for(int i=0;i<w_len;i++)
	{
		r14(i)=sqrt(2*x_val(i,0)*x_val(i,0)-2*x_val(i,0)*x_val(i,0)*cos(2*a_val(i,0)+b_val(i,0)));
	}
//	r14.print("this is r14");
	//setting the values of non-unique interaction distances according to symmetry
	arma::vec r34=r12;
	arma::vec r56=r12;
	arma::vec r45=r23;
	arma::vec r16=r23;
	arma::vec r24=r13;
	arma::vec r35=r13;
	arma::vec r46=r13;
	arma::vec r15=r13;
	arma::vec r26=r13;
	arma::vec r36=r14;
	arma::vec r25=r14;
//	r14.raw_print("this r14");
//	r25.raw_print("this is r25");
	arma::mat S=arma::zeros(6,6);
	arma::mat Ham=arma::zeros(6,6);
	Ham.diag().fill(-10.77/27.212);
	arma::vec Elec=arma::zeros(w_len);
	arma::vec Elec_cat=arma::zeros(w_len);
	arma::vec Enuc=arma::zeros(w_len);
	for(int i=0;i<w_len;i++)
	{
		Enuc(i)=Z_eff*Z_eff/r12(i)+Z_eff*Z_eff/r13(i)+Z_eff*Z_eff/r14(i)+Z_eff*Z_eff/r15(i)+
			Z_eff*Z_eff/r16(i)+Z_eff*Z_eff/r23(i)+Z_eff*Z_eff/r24(i)+Z_eff*Z_eff/r25(i)+
			Z_eff*Z_eff/r26(i)+Z_eff*Z_eff/r34(i)+Z_eff*Z_eff/r35(i)+Z_eff*Z_eff/r36(i)+
			Z_eff*Z_eff/r45(i)+Z_eff*Z_eff/r46(i)+Z_eff*Z_eff/r56(i);
	}
//	Ham.print();
	for(int i=0; i<w_len;i++)
	{
		S(0,2)=(1+Z*r13(i)+(2*Z*Z*r13(i)*r13(i))/5+(Z*Z*Z*r13(i)*r13(i)*r13(i))/15)*exp(-Z*r13(i));
		S(0,3)=(1+Z*r14(i)+(2*Z*Z*r14(i)*r14(i))/5+(Z*Z*Z*r14(i)*r14(i)*r14(i))/15)*exp(-Z*r14(i));
		S(0,4)=(1+Z*r15(i)+(2*Z*Z*r15(i)*r15(i))/5+(Z*Z*Z*r15(i)*r15(i)*r15(i))/15)*exp(-Z*r15(i));
		S(1,3)=(1+Z*r24(i)+(2*Z*Z*r24(i)*r24(i))/5+(Z*Z*Z*r24(i)*r24(i)*r24(i))/15)*exp(-Z*r24(i));
		S(1,4)=(1+Z*r25(i)+(2*Z*Z*r25(i)*r25(i))/5+(Z*Z*Z*r25(i)*r25(i)*r25(i))/15)*exp(-Z*r25(i));
		S(1,5)=(1+Z*r26(i)+(2*Z*Z*r26(i)*r26(i))/5+(Z*Z*Z*r26(i)*r26(i)*r26(i))/15)*exp(-Z*r26(i));
		S(2,4)=(1+Z*r35(i)+(2*Z*Z*r35(i)*r35(i))/5+(Z*Z*Z*r35(i)*r35(i)*r35(i))/15)*exp(-Z*r35(i));
		S(2,5)=(1+Z*r36(i)+(2*Z*Z*r36(i)*r36(i))/5+(Z*Z*Z*r36(i)*r36(i)*r36(i))/15)*exp(-Z*r36(i));
		S(3,5)=(1+Z*r46(i)+(2*Z*Z*r46(i)*r46(i))/5+(Z*Z*Z*r46(i)*r46(i)*r46(i))/15)*exp(-Z*r46(i));
	//setting the S matrix for the bonded atoms
		S(0,1)=(1+Z*r12(i)+(2*Z*Z*r12(i)*r12(i))/5+(Z*Z*Z*r12(i)*r12(i)*r12(i))/15)*exp(-Z*r12(i));
		S(2,3)=S(0,1);
		S(4,5)=S(0,1);
		S(1,2)=(1+Z*r23(i)+(2*Z*Z*r23(i)*r23(i))/5+(Z*Z*Z*r23(i)*r23(i)*r23(i))/15)*exp(-Z*r23(i));
		S(3,4)=S(1,2);
		S(0,5)=S(1,2);
		arma::mat A=S.t();
		arma::mat B=S;
		S=A+B;
	//	S.print("this is initial S");
		S.diag().fill(1);
	//	S.print(" this is S");
	//	cout<<endl;
	//loop to initialize the matrix elements of the Hamiltonian Marix
		for(int row=0;row<6;row++)
		{
			for(int col=0;col<6;col++)
			{
				if(row!=col)
				{
				//	cout<<row<<" "<<col<<"  "<<endl;
					Ham(row,col)=0.875*S(row,col)*(Ham(row,row)+Ham(col,col));
				}
			}
		}
		S.print("this is S");
		cout<<endl;
		Ham.print("this is H");
		cout<<endl;
		arma::mat temp_H=Ham-S;
		temp_H.print();
		arma::cx_vec ee;
		arma::cx_mat ev;
		arma::eig_pair(ee,ev,Ham,S);
		arma::vec e_real=arma::real(ee);
		e_real=arma::sort(e_real);
		e_real.print();
		cout<<endl;
		Elec(i)=2*e_real(0)+2*e_real(1)+2*e_real(2);
		Elec_cat(i)=2*e_real(0)+e_real(1)*2+e_real(2);
		S=arma::zeros(6,6);
	}
	Elec_cat.print("this is the Elec");
	Enuc.print("this is the Enuc");
	arma::vec Etot=Elec_cat+Enuc;
	Etot.raw_print();

//*******************************************************************************************************
//
//
//
//********************************************************************************************************
	
	
	
	arma::mat a_cote(w_len,2);
	arma::mat b_cote(w_len,2);
	arma::mat x_cote(w_len,2);
	for(int i=0; i<w_len; i++)
	{
		int col=0;
		double a_ang=(1/4.0)*pi;
		double x1_len=0.0;
		double x2_len=0.0;
		while(a_ang>(1.0/6.0)*pi)
		{
			double b_ang=(1.0/2.0)*pi-a_ang;
			x1_len=sqrt(((R-w(i))*(R-w(i)))/(2-2*cos(a_ang)));
			x2_len=sqrt(((R+w(i))*(R+w(i)))/(2-2*cos(b_ang)));
	//		cout<<b_ang<<" "<<a_ang<<" "<<x1_len<<" "<<x2_len<<endl;
			if(abs(x1_len-x2_len) < 0.0001 )
			{
				//cout<<"this is col"<<col<<endl;
				//cout<<a_ang<<b_ang<<x1_len<<endl;
				if (col == 0)
				{
					a_cote(i,col)=a_ang;
					b_cote(i,col)=b_ang;
					x_cote(i,col)=x1_len;
					col++;

				}
				else
				{
					if (abs((x_cote(i,col-1))-(x1_len))<0.00001 )
					{
					//	cout<<x_val(i,col-1)<<" "<<x1_len<<" "<<abs((x_val(i,col-1))-(x1_len))<<endl;
						a_cote(i,col)=a_ang;
						b_cote(i,col)=b_ang;
						x_cote(i,col)=x1_len;
						col++;
					}
				}
			}
			a_ang=a_ang-0.00001;
		}
	}
//	a_val.print("this is a");
//	b_val.print("this is b"); 
//	x_val.print("this is x");
//	arma::mat c_cote=(pi-a_cote)/2;
//	arma::mat d_cote=(pi-b_cote)/2;
//	c_ang.print("this is c");
	arma::vec c_r12=R-w;
	arma::vec c_r23=R+w;
	arma::vec c_r13=arma::zeros(w_len);
	for(int i=0;i<w_len;i++)
	{
	       	c_r13(i)=sqrt(2)*x_cote(i);
	}
//	r13.print("this is r13");
	arma::vec c_r14=arma::zeros(w_len);
	for(int i=0;i<w_len;i++)
	{
		c_r14(i)=sqrt(2*x_cote(i)*x_cote(i)-2*x_cote(i)*x_cote(i)*cos(2*a_cote(i)+b_cote(i)));
	}
	arma::vec c_r25=arma::zeros(w_len);
	for(int i=0;i<w_len;i++)
	{
		c_r25(i)=sqrt(2*x_cote(i)*x_cote(i)-2*x_cote(i)*x_cote(i)*cos(a_cote(i)+2*b_cote(i)));
	}
	arma::vec c_r15= arma::zeros(w_len);
	for(int i=0;i<w_len;i++)
	{
		c_r15(i)=2*x_cote(i);
	}

//	r14.print("this is r14");
	//setting the values of non-unique interaction distances according to symmetry
	arma::vec c_r34=c_r12;
	arma::vec c_r56=c_r12;
	arma::vec c_r78=c_r12;
	arma::vec c_r45=c_r23;
	arma::vec c_r67=c_r23;
	arma::vec c_r18=c_r23;
	arma::vec c_r24=c_r13;
	arma::vec c_r35=c_r13;
	arma::vec c_r46=c_r13;
	arma::vec c_r57=c_r13;
	arma::vec c_r68=c_r13;
	arma::vec c_r17=c_r13;
	arma::vec c_r28=c_r13;
	arma::vec c_r36=c_r14;
	arma::vec c_r58=c_r14;
	arma::vec c_r27=c_r14;
	arma::vec c_r47=c_r25;
	arma::vec c_r16=c_r25;
	arma::vec c_r38=c_r25;
	arma::vec c_r26=c_r15;
	arma::vec c_r37=c_r15;
	arma::vec c_r48=c_r15;

	arma::mat S=arma::zeros(8,8);
	arma::mat Ham=arma::zeros(8,8);
	Ham.diag().fill(-10.77/27.212);
	arma::vec Elec=arma::zeros(w_len);
	arma::vec Elec_cat=arma::zeros(w_len);
	arma::vec Enuc=arma::zeros(w_len);
	for(int i=0;i<w_len;i++)
	{
		Enuc(i)=Z_eff*Z_eff/r12(i)+Z_eff*Z_eff/r13(i)+Z_eff*Z_eff/r14(i)+Z_eff*Z_eff/r15(i)+
			Z_eff*Z_eff/r16(i)+Z_eff*Z_eff/r23(i)+Z_eff*Z_eff/r24(i)+Z_eff*Z_eff/r25(i)+
			Z_eff*Z_eff/r26(i)+Z_eff*Z_eff/r34(i)+Z_eff*Z_eff/r35(i)+Z_eff*Z_eff/r36(i)+
			Z_eff*Z_eff/r45(i)+Z_eff*Z_eff/r46(i)+Z_eff*Z_eff/r56(i);
	}
	Ham.print();
	for(int i=0; i<w_len;i++)
	{
		S(0,2)=(1+Z*r13(i)+(2*Z*Z*r13(i)*r13(i))/5+(Z*Z*Z*r13(i)*r13(i)*r13(i))/15)*exp(-Z*r13(i));
		S(0,3)=(1+Z*r14(i)+(2*Z*Z*r14(i)*r14(i))/5+(Z*Z*Z*r14(i)*r14(i)*r14(i))/15)*exp(-Z*r14(i));
		S(0,4)=(1+Z*r15(i)+(2*Z*Z*r15(i)*r15(i))/5+(Z*Z*Z*r15(i)*r15(i)*r15(i))/15)*exp(-Z*r15(i));
		S(1,3)=(1+Z*r24(i)+(2*Z*Z*r24(i)*r24(i))/5+(Z*Z*Z*r24(i)*r24(i)*r24(i))/15)*exp(-Z*r24(i));
		S(1,4)=(1+Z*r25(i)+(2*Z*Z*r25(i)*r25(i))/5+(Z*Z*Z*r25(i)*r25(i)*r25(i))/15)*exp(-Z*r25(i));
		S(1,5)=(1+Z*r26(i)+(2*Z*Z*r26(i)*r26(i))/5+(Z*Z*Z*r26(i)*r26(i)*r26(i))/15)*exp(-Z*r26(i));
		S(2,4)=(1+Z*r35(i)+(2*Z*Z*r35(i)*r35(i))/5+(Z*Z*Z*r35(i)*r35(i)*r35(i))/15)*exp(-Z*r35(i));
		S(2,5)=(1+Z*r36(i)+(2*Z*Z*r36(i)*r36(i))/5+(Z*Z*Z*r36(i)*r36(i)*r36(i))/15)*exp(-Z*r36(i));
		S(3,5)=(1+Z*r46(i)+(2*Z*Z*r46(i)*r46(i))/5+(Z*Z*Z*r46(i)*r46(i)*r46(i))/15)*exp(-Z*r46(i));
	//setting the S matrix for the bonded atoms
		S(0,1)=(1+Z*r12(i)+(2*Z*Z*r12(i)*r12(i))/5+(Z*Z*Z*r12(i)*r12(i)*r12(i))/15)*exp(-Z*r12(i));
		S(2,3)=S(0,1);
		S(4,5)=S(0,1);
		S(1,2)=(1+Z*r23(i)+(2*Z*Z*r23(i)*r23(i))/5+(Z*Z*Z*r23(i)*r23(i)*r23(i))/15)*exp(-Z*r23(i));
		S(3,4)=S(1,2);
		S(0,5)=S(1,2);
		arma::mat A=S.t();
		arma::mat B=S;
		S=A+B;
	//	S.print("this is initial S");
		S.diag().fill(1);
	//	S.print(" this is S");
	//	cout<<endl;
	//loop to initialize the matrix elements of the Hamiltonian Marix
		for(int row=0;row<6;row++)
		{
			for(int col=0;col<6;col++)
			{
				if(row!=col)
				{
				//	cout<<row<<" "<<col<<"  "<<endl;
					Ham(row,col)=0.875*S(row,col)*(Ham(row,row)+Ham(col,col));
				}
			}
		}
		S.print("this is S");
		cout<<endl;
		Ham.print("this is H");
		cout<<endl;
		arma::mat temp_H=Ham-S;
		temp_H.print();
		arma::cx_vec ee;
		arma::cx_mat ev;
		arma::eig_pair(ee,ev,Ham,S);
		arma::vec e_real=arma::real(ee);
		e_real=arma::sort(e_real);
		e_real.print();
		cout<<endl;
		Elec(i)=2*e_real(0)+2*e_real(1)+2*e_real(2);
		Elec_cat(i)=2*e_real(0)+e_real(1)*2+e_real(2);
		S=arma::zeros(6,6);
	}
	Elec_cat.print("this is the Elec");
	Enuc.print("this is the Enuc");
	arma::vec Etot=Elec_cat+Enuc;
	Etot.raw_print();
//


}	
