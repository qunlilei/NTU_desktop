#include<iostream>
#include<fstream>
#include<cstdlib>
#include <iomanip>
#include <math.h>

using namespace std;

int main()
{
	int n_seg=200;

	int fold=100;

	int nx_row=1;
	int ny_row=1;
	int nz_row=1;

	double chain_space_x=100;
	double chain_space_y=100;

	double box_lx=nx_row*chain_space_x;
	double box_ly=ny_row*chain_space_y;
	int n_chain=nx_row*ny_row*nz_row;

	int n_drug=n_chain*n_seg*0;


	double height=nz_row*fold;


	//##########################################################################
	ofstream outfile;
	outfile.open("configuration.txt",ios::out);
	if(!outfile)
	{
		cout<<"open error"<<endl;
		exit(1);
	}

	outfile<<setiosflags(ios::fixed)<<setprecision(4);
	outfile<<endl;
	outfile<<	n_seg*n_chain+n_drug                                                <<" atoms"   <<endl;
	outfile<<	(n_seg-1)*n_chain									                    <<" bonds"   <<endl;
	outfile<<	(n_seg-2)*n_chain		                                                <<" angles"   <<endl;
	outfile<<	0				     	                                                <<" dihedrals"<<endl;
	outfile<<endl;

	outfile<<	2	<<      " atom types"<<endl;
	outfile<<	1	<<		" bond types"<<endl;
	outfile<<	2	<< 		" angle types"<<endl;
	outfile<<	0	<< 		" dihedral types"<<endl;	

	outfile<<endl;
	outfile<<	-box_lx/2	<<" "<< 		box_lx/2	        <<   " xlo" << " xhi"<<endl;
	outfile<<	-box_ly/2	<<" "<<  		box_ly/2	        <<   " ylo" << " yhi"<<endl;
	outfile<<	-height/2	<<" "<< 		height/2			<<   " zlo" << " zhi"<<endl;
	outfile<<endl;

	outfile<<	"Masses"<<endl;
	outfile<<endl;
	outfile<<	1	<<" "<<    1.	<<endl;
	outfile<<	2	<<" "<<	   1.	<<endl;


//###################################################################################################
	outfile<<endl;
	outfile<<	"Atoms"<<endl;
	outfile<<endl;

	for(int i=1;i<=n_chain;i++)
	{
		double x_start=(i%nx_row+0.5)*chain_space_x-box_lx/2;
		double y_start=( (i-1)/nx_row+0.5)*chain_space_y-box_ly/2;
		double z_start= (i-1)/(nx_row*ny_row)*(n_seg)-height/2;


		//-----------------------------------------------------------------------------------------------------------------------
		int direct=1;
		int incre_z=0;
		int incre_x=0;
		for(int j=1;j<=n_seg;j++)
		{
			incre_z=incre_z+direct;
			if((j-1)%fold==0 && j!=1)  
			{
				direct=-1*direct;  
				incre_z=incre_z+direct;
				incre_x++;
			}
			outfile<< (i-1)*n_seg+j   <<" "<<  i  <<" "<<  1  <<" "<< -1 <<" "<<  x_start+incre_x      <<" "<<   y_start  <<" "<<    z_start+incre_z <<endl;
		}
		//-----------------------------------------------------------------------------------------------------------------------
	}

	for(int i=1;i<=n_drug;i++)
	{
		double x_start=(i%nx_row+0.9)*chain_space_x-(i-1)/(fold*nx_row*ny_row)-box_lx/2;
		double y_start=( (i-1)/nx_row+0.9)*chain_space_y-box_ly/2;
		double z_start=(((i-1)/(nx_row*ny_row))%fold)-height/2;


		//-----------------------------------------------------------------------------------------------------------------------
			outfile<< n_chain*n_seg+i   <<" "<<  n_chain+1  <<" "<<  2  <<" "<<  1 <<" "<<  x_start      <<" "<<   y_start  <<" "<<    z_start <<endl;
		//-----------------------------------------------------------------------------------------------------------------------
	}


	outfile<<endl;



//############################################################################################################
	outfile<<	"Bonds"<<endl;
	outfile<<endl;

	for(int i=1;i<=n_chain;i++)
	{
		for(int j=1;j<=n_seg-1;j++)
		{
			outfile<<	(i-1)*n_seg +j  <<" "<<  1 <<" "<< (i-1)*n_seg +j  <<" "<< (i-1)*n_seg +j+1 <<endl;
		}
	}

	outfile<<endl;
//############################################################################################################################
	outfile<<	"Angles"<<endl;
	outfile<<endl;

	for(int i=1;i<=n_chain;i++)
	{
		for(int j=1;j<=n_seg-2;j++)
		{
			outfile<<   (i-1)*n_seg +j  <<" "<<  1  <<" "<<  (i-1)*n_seg+j  <<" "<< (i-1)*n_seg+j+1 <<" "<< (i-1)*n_seg+j+2 <<endl;
		}
	}

	outfile<<endl;
//###############################################################################################################################################
//	outfile<<	"Dihedrals"<<endl;
//	outfile<<endl;

//	for(i=1;i<=n_chain;i++)
//	{
//		for(j=1;j<=n_seg-3;j++)
//		{
//			outfile<<   j  <<" "<<  1  <<" "<<  j  <<" "<< j+1 <<" "<< j+2 <<" "<< j+3 <<endl;
//		}
//	}

	return 0;
}
