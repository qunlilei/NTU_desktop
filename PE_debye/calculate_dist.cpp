#include<iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include<vector>
#include<math.h>
#define SQR2(x) ((x)*(x))

using namespace std;

const double PI  =3.141592653589793238463;
struct Helix
{ 
	double x,y,z;
	int   ix,iy,iz;
	int    type,id;
	double vect_x,vect_y,vect_z; 
};

ofstream outfile;
//////////////////////////////////////////////////////////////////////////
int main()
{
	int bin_number=500;
	double *main_chain, *side_chain1, *side_chain2;
	main_chain =new double[bin_number];
	side_chain1=new double[bin_number];
	side_chain2=new double[bin_number];
	for(int i=0;i<bin_number;i++)
	{
		 main_chain[i]  = 0.;
		side_chain1[i]  = 0.;
		side_chain2[i]  = 0.;
	}
//-----------------------------------------------------------------------------------------------------
	double temp;
	double box_xl,box_xh,box_yl,box_yh,box_zl,box_zh, Lx_box,Ly_box,Lz_box;
	int number0;

	ifstream infile("bottlebrush.5000000.lammpstrj",ios::in);
	if(!infile)
	{
	  cerr<<"open error!"<<endl;
	}



//-------------------------------------------------------------------------------------------------

	int frame=0;
	int frame_number=0;
	string line;
	getline(infile,line);
	while(infile)
	{
		infile>>frame;
		cout<<frame<<endl;
		for(int i=0;i<2;i++)  getline(infile,line);
		infile>>number0;
		cout<<"atom_number "<<number0<<endl;
		for(int i=0;i<2;i++)  getline(infile,line);
		infile>>box_xl>>box_xh;
		infile>>box_yl>>box_yh;
		infile>>box_zl>>box_zl;
		for(int i=0;i<2;i++)  getline(infile,line);
		Lx_box=box_xh-box_xl;
		Ly_box=box_yh-box_yl;
		Lz_box=box_zh-box_zl;
		
		//---------------------------------------------------------------------------------------
		for(int i=0;i<number0;i++)
		{
			double x,y,z,charge;
			int    ix,iy,iz,type,id,mol;
			infile>> id >> type >> x>> y>> z>> ix>> iy>>iz;
//			cout<< id<<endl;
			if(frame > 0*10000)
			{
				double r=sqrt( (x+0*ix*Lx_box)*(x+0*ix*Lx_box)+(y+0*iy*Ly_box)*(y+0*iy*Ly_box)+(z+0*iz*Lz_box)*(z+0*iz*Lz_box));
				int bin= int ( r/ (Lx_box/bin_number) );
				if(type==1)    main_chain[bin] +=1.0;
				if(type==2)   side_chain1[bin] +=1.0;
				if(type==3)   side_chain2[bin] +=1.0;
				frame_number++;
			}
		}
		for(int i=0;i<2;i++)  getline(infile,line);
		cout<<frame_number/number0<<" "<<Lx_box << endl;
	}

	//---------------------------------------------------------------------------------------
//		char *pPath = new char[50];
//		sprintf(pPath,"data/%04d.dat",frame);

//		ofstream component;
//		component.open(pPath,ios::out);
//		if(!component)
//		{
//			cout<<"open error"<<endl;
//			exit(1);
//		}

//		delete [] pPath;
	//---------------------------------------------------------------------------------------
	outfile.open("component_distribution.dat",ios::out);
	if(!outfile)
	{
		cout<<"open error"<<endl;
	}

	outfile<<setiosflags(ios::fixed)<<setprecision(4);

	for(int i=1;i<bin_number;i++)
	{
		double r=(i*Lx_box/bin_number);
		 main_chain[i]  /= r*r*4*3.1415*Lx_box/bin_number;
		side_chain1[i]  /= r*r*4*3.1415*Lx_box/bin_number;
		side_chain2[i]  /= r*r*4*3.1415*Lx_box/bin_number;
		 main_chain[i]  /= frame_number/number0;
		side_chain1[i]  /= frame_number/number0;
		side_chain2[i]  /= frame_number/number0;
		 main_chain[i]  *= sqrt(2.0)/2.0;
		side_chain1[i]  *= sqrt(2.0)/2.0;
		side_chain2[i]  *= sqrt(2.0)/2.0;
		outfile<<r/2.0<<" "<<main_chain[i]<<" "<<side_chain1[i]<<" "<<side_chain2[i]<<endl;
	}
	outfile.close();


return 0;
}





































//#####################################################################################################################
double perio_dist_2d(Helix *p_x1, Helix *p_x2, double Lx_box,double  Ly_box)
{
	double ds_x=p_x2->x-p_x1->x;
	if     (ds_x >  Lx_box/2.0 ) ds_x=ds_x-Lx_box;
	if     (ds_x < -Lx_box/2.0 ) ds_x=ds_x+Lx_box;
	double ds_y=p_x2->y-p_x1->y;
	if     (ds_y >  Ly_box/2.0 ) ds_y=ds_y-Ly_box;
	if     (ds_y < -Ly_box/2.0 ) ds_y=ds_y+Ly_box;
	double distance=sqrt(ds_x*ds_x+ds_y*ds_y);
	return distance;
}

//#####################################################################################################################
double perio_angle_2d(Helix *p_x1, Helix *p_x2,double  Lx_box,double  Ly_box)
{
	double ds_x=p_x2->x-p_x1->x;
	if     (ds_x >  Lx_box/2.0 ) ds_x=ds_x-Lx_box;
	if     (ds_x < -Lx_box/2.0 ) ds_x=ds_x+Lx_box;
	double ds_y=p_x2->y-p_x1->y;
	if     (ds_y >  Ly_box/2.0 ) ds_y=ds_y-Ly_box;
	if     (ds_y < -Ly_box/2.0 ) ds_y=ds_y+Ly_box;
	double angle=atan2(ds_y,ds_x);
	return angle;
}
//######################################################################################################################
void cal_correlation( Helix *x,double Lx_box, double Ly_box, int n)
{
	double MO1_cos=0.0;
	double MO1_sin=0.0;

	double BO2_cos=0.0;
	double BO2_sin=0.0;

	double BO3_cos_r=0.0;
	double BO3_sin_r=0.0;
	double BO3_cos_l=0.0;
	double BO3_sin_l=0.0;

	double BO4_cos=0.0;
	double BO4_sin=0.0;

	double BO6_cos=0.0;
	double BO6_sin=0.0;
	double BO6_inter_angle11=0.0;
	double BO6_inter_angle12=0.0;
	int BO6_n12=0;
	int BO6_n11=0;

	for(int i=1;i<n;i++)
	{
		MO1_cos=MO1_cos+x[i].vect_x;
		MO1_sin=MO1_sin+x[i].vect_y;
//---------------------------------------------------------------------------------------
		Helix *p_BO6_neigb[8];
		double neigb_dist[8];
		for(int k=0;k<7;k++) p_BO6_neigb[k]=x;
		neigb_dist[0]=0.0;
		for(int k=1;k<7;k++) neigb_dist[k]=Lx_box+Ly_box;
		for(int j=0;j<n;j++)
		{
			if(j!=i)
			{
				double perio_dist_2d(Helix *p_x1, Helix *p_x2, double Lx_box,double  Ly_box);
				double distance=perio_dist_2d( x+i, x+j, Lx_box, Ly_box);
				for(int k=1;k<=6;k++)
				{
					if( neigb_dist[k-1] <= distance  &&  distance < neigb_dist[k] )
					{
						for(int ii=6;ii>=k;ii--)
						{
							neigb_dist[ii+1]=neigb_dist[ii];
							p_BO6_neigb[ii+1]=p_BO6_neigb[ii];
						}
						neigb_dist[k]=distance;
						p_BO6_neigb[k]=x+j;
						break;
					}
				}
			}
		}
		for(int j=1;j<=2;j++)
		{
			double perio_angle_2d(Helix *p_x1, Helix *p_x2,double  Lx_box,double  Ly_box);
			double angle=perio_angle_2d( x+i, p_BO6_neigb[j], Lx_box, Ly_box);
			BO2_cos=BO2_cos+cos(2.0*angle);
			BO2_sin=BO2_sin+sin(2.0*angle);
		}
		for(int j=1;j<=3;j++)
		{
			double perio_angle_2d(Helix *p_x1, Helix *p_x2,double  Lx_box,double  Ly_box);
			double angle=perio_angle_2d( x+i, p_BO6_neigb[j], Lx_box, Ly_box);
			if(x[i].type==1) {BO3_cos_r=BO3_cos_r+cos(3.0*angle);   BO3_sin_r=BO3_sin_r+sin(3.0*angle);}
			if(x[i].type==2) {BO3_cos_l=BO3_cos_l+cos(3.0*angle);   BO3_sin_l=BO3_sin_l+sin(3.0*angle);}
		}
		for(int j=1;j<=4;j++)
		{
			double perio_angle_2d(Helix *p_x1, Helix *p_x2,double  Lx_box,double  Ly_box);
			double angle=perio_angle_2d( x+i, p_BO6_neigb[j], Lx_box, Ly_box);
			BO4_cos=BO4_cos+cos(4.0*angle);
			BO4_sin=BO4_sin+sin(4.0*angle);
		}
		for(int j=1;j<=6;j++)
		{
			double perio_angle_2d(Helix *p_x1, Helix *p_x2,double  Lx_box,double  Ly_box);
			double angle=perio_angle_2d( x+i, p_BO6_neigb[j], Lx_box, Ly_box);
			BO6_cos=BO6_cos+cos(6.0*angle);
			BO6_sin=BO6_sin+sin(6.0*angle);
			if(x[i].type==p_BO6_neigb[j]->type)
			{	
				BO6_n11++;
				double alpha=fabs(atan2(x[i].vect_y,x[i].vect_x)-atan2(p_BO6_neigb[j]->vect_y,p_BO6_neigb[j]->vect_x))/PI*180;
				if(alpha>180) alpha=360-alpha;
				BO6_inter_angle11=BO6_inter_angle11+alpha;
			}
			else
			{	
				BO6_n12++;
				double alpha=fabs(atan2(x[i].vect_y,x[i].vect_x)+atan2(p_BO6_neigb[j]->vect_y,p_BO6_neigb[j]->vect_x)-2.0*angle )/PI*180;
				if(alpha>180) alpha=360-alpha;
				BO6_inter_angle12=BO6_inter_angle12+ alpha;
			}
//			cout<<atan2(x[i].vect_y,x[i].vect_x)+atan2(p_BO6_neigb[j]->vect_y,p_BO6_neigb[j]->vect_x)<<" "<<2.0*angle<<endl;
		}
	}
//---------------------------------------------------------------------------------------------------
	double MO1_mod=sqrt(MO1_cos*MO1_cos+MO1_sin*MO1_sin)/n;
	double BO2_mod=sqrt(BO2_cos*BO2_cos+BO2_sin*BO2_sin)/n/2;
	double BO3_mod=( sqrt(SQR2(BO3_cos_r)+SQR2(BO3_sin_r))+sqrt(SQR2(BO3_cos_l)+SQR2(BO3_sin_l)) )/n/3;
	double BO4_mod=sqrt(BO4_cos*BO4_cos+BO4_sin*BO4_sin)/n/4;
	double BO6_mod=sqrt(BO6_cos*BO6_cos+BO6_sin*BO6_sin)/n/6;
	BO6_inter_angle11=BO6_inter_angle11/BO6_n11;
	if(BO6_n12!=0) BO6_inter_angle12=BO6_inter_angle12/BO6_n12;


	outfile<< MO1_mod<<" "<< BO2_mod<<" "<< BO3_mod<<" "<< BO4_mod<<" "<< BO6_mod <<" "<< BO6_inter_angle11 <<" "<< BO6_inter_angle12 <<endl;
	return;
}

