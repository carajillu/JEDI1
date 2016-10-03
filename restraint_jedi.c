/*
#
#
#                          JEDI
#        Just Evaluate Druggability at the Interface
#
#              Remi CUCHILLO & Julien MICHEL
#
#
#
*/


#include "metadyn.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

//#include <gsl/gsl_math.h>
//#include <gsl/gsl_eigen.h>


/* 
# This is the parser to get the atom numbers and cutoffs
#
# JEDI LIST <apolar> <polar> b_grid n_grid CUTOFF_close CUTOFF_far CUTOFF_enclosure CUTOFF_hull CUTOFF_surface CUTOFF_contact CUTOFF_hydro CUTOFF_con dump
#
*/

int PREFIX read_jedi(char **word, int count, t_plumed_input *input, FILE *fplog)
{

  int x, y, i, i_max, j, iw;
  double delta            = 0.0;
  double cutoff_close     = 0.0;
  double cutoff_far       = 0.0;
  double cutoff_hydro     = 0.0;
  double cutoff_enclosure = 0.0;
  double cutoff_hull      = 0.0;
  double cutoff_surface   = 0.0;
  double cutoff_contact   = 0.0;
  double cutoff_con       = 0.0;
  double dump_matrix      = 0.0;
  double n_pocket         = 0.0;
  double n_grid           = 0.0;
  double b_grid           = 0.0;
  double COM_X            = 0.0;
  double COM_Y            = 0.0;
  double COM_Z            = 0.0;

// THIS CV HAS PBC ON BY DEFAULT
//  colvar.cell_pbc[count]=1; 

// LOOK FOR '' LIST '' KEYWORD <apolar> <polar>
  iw = seek_word(word,"LIST");
  if(iw>=0) {
//------------ <apolar>
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
//------------ <polar>
             j=plumed_get_group(word[iw+2],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][1]=j;

  } else {plumed_error("NEEDED LIST KEYWORD FOR NEWCV");}

// LOOK FOR '' SIGMA '' KEYWORD  
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);  
             colvar.delta_r[count]  = (real) delta; }

  iw=seek_word(word,"COM_X");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &COM_X);  
             colvar.COM_X[count]  = (real) COM_X; }

  iw=seek_word(word,"COM_Y");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &COM_Y);  
             colvar.COM_Y[count]  = (real) COM_Y; }

  iw=seek_word(word,"COM_Z");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &COM_Z);  
             colvar.COM_Z[count]  = (real) COM_Z; }

  iw=seek_word(word,"b_grid");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &b_grid);  
             colvar.b_grid[count]  = (real) b_grid; }

  iw=seek_word(word,"n_grid");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &n_grid);  
             colvar.n_grid[count]  = (real) n_grid; }

  iw=seek_word(word,"n_pocket");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &n_pocket);  
             colvar.n_pocket[count]  = (real) n_pocket; }

  iw=seek_word(word,"CUTOFF_close");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &cutoff_close);  
             colvar.cutoff_close[count]  = (real) cutoff_close; }

  iw=seek_word(word,"CUTOFF_far");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &cutoff_far);  
             colvar.cutoff_far[count]  = (real) cutoff_far; }

  iw=seek_word(word,"CUTOFF_enclosure");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &cutoff_enclosure);  
             colvar.cutoff_enclosure[count]  = (real) cutoff_enclosure; }

  iw=seek_word(word,"CUTOFF_hull");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &cutoff_hull);  
             colvar.cutoff_hull[count]  = (real) cutoff_hull; }

  iw=seek_word(word,"CUTOFF_surface");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &cutoff_surface);  
             colvar.cutoff_surface[count]  = (real) cutoff_surface; }

  iw=seek_word(word,"CUTOFF_contact");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &cutoff_contact);  
             colvar.cutoff_contact[count]  = (real) cutoff_contact; }

  iw=seek_word(word,"CUTOFF_hydro");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &cutoff_hydro);  
             colvar.cutoff_hydro[count]  = (real) cutoff_hydro; }

  iw=seek_word(word,"CUTOFF_con");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &cutoff_con);  
             colvar.cutoff_con[count]  = (real) cutoff_con; }

  iw=seek_word(word,"dump");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &dump_matrix);  
             colvar.dump_matrix[count]  = (real) dump_matrix; }

//fill up neighbors list  from neighbors.txt

  FILE* m_File;
  int n[27] = {0};
  i = 0;
  m_File = fopen("neighbors.txt","rt");
  while(!feof(m_File) && i < colvar.n_grid[0])
  {
    j = 0;
    fscanf(m_File, "%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i", &n[0],&n[1],&n[2],&n[3],&n[4],&n[5],&n[6],&n[7],&n[8],&n[9],&n[10],&n[11],&n[12],&n[13],&n[14],&n[15],&n[16],&n[17],&n[18],&n[19],&n[20],&n[21],&n[22],&n[23],&n[24],&n[25],&n[26]);
    for(j=0;j<27;j++)
    {
      colvar.neighbors[i][j] = n[j];
    }
    i++;
  }

//fill up rays
  i = 0;
  m_File = fopen("rays.txt","rt");
  while(!feof(m_File) && i < colvar.n_grid[0]){
    for(j=0;j<44;j++){
      if (!fscanf(m_File, "%i", &colvar.rays[i][j]))
        break;
//      printf("test : %i\n",colvar.rays[i][j]);
    }
    i++;
  }


//fill up gridpoints coordinates from gridpoints_coord.txt

//  FILE* m_File;
  real n_d[3] = {0.0};
  i = 0;
  m_File = fopen("gridpoints_coord.txt","rt");
  while(!feof(m_File) && i < colvar.n_grid[0])
  {
    j = 0;
    fscanf(m_File, "%lf %lf %lf", &n_d[0],&n_d[1],&n_d[2]);
    for(j=0;j<3;j++)
    {
      colvar.grid_pos[i][j] = n_d[j];
    }
    i++;
  }

//  FILE* m_File;
  real m[1] = {0.0};
  i = 0;
  m_File = fopen("pocket.txt","rt");
  while(!feof(m_File) && i < colvar.n_grid[0])
  {
    fscanf(m_File, "%lf",&m[0]);
    colvar.pocket[i] = m[0];
//    printf("test : %f\n",colvar.pocket[i]);
    i++;
  }


// SET CV TYPE 0 is not used
  colvar.type_s[count]   = 0;

// ALLOCATE ARRAY FOR DERIVATIVES
  snew(colvar.myder[count], colvar.natoms[count]);

  return colvar.natoms[count];

}

//--------------------------------------------------------------------------------

// This is the routine where JEDI and its derivatives are calculated

void PREFIX jedi_restraint(int i_c, struct mtd_data_s *mtd_data)
{
// subroutines to calulate the switching functions and their derivatives

// equation 3
double s_on(double k,double v,double v_min,double delta)
{
  double s, m;
  m = (v-v_min)/(delta);
  if(m < 0.){
    s=0.;
  }else if(m > 1.){
    s=k;
  }else{
    s=k*(1.-(pow((1.-pow(m,2)),2))*((1.+2*(pow(m,2)))));
  }
  return s;
}
// equation 2
double s_off(double k,double v,double v_min,double delta)
{
  double s, m;
  m = (v-v_min)/(delta);
  if(m < 0.){
    s=k;
  }else if(m > 1.){
    s=0.;
  }else{
    s=k*(pow((1.-pow(m,2)),2)*((1.+2*(pow(m,2)))));
  }
  return s;
}
// equation 8 SI
double ds_off_dm(double k,double v,double v_min,double delta)
{
  double s, m;
  m = (v-v_min)/(delta);
  if(m < 0.){
    s = 0.;
  }else if(m > 1.){
    s = 0.;
  }else{
    s = 4*k*m*(pow((1.-pow(m,2)),2)) - 4*k*m*(1.-pow(m,2))*((1.+2*(pow(m,2))));
  }
  return s;
}
// equation 10 SI
double ds_off_dk(double k,double v,double v_min,double delta)
{
  double s, m;
  m = (v-v_min)/(delta);
  if(m < 0.){
    s = 1.;
  }else if(m > 1.){
    s = 0.;
  }else{
    s = (pow((1.-pow(m,2)),2)) * ((1.+2*(pow(m,2))));
  }
  return s;
}
// equation 13 SI
double ds_on_dm(double k,double v,double v_min,double delta)
{
  double s, m;
  m = (v-v_min)/(delta);
  if(m < 0.){
    s = 0.;
  }else if(m > 1.){
    s = 0.;
  }else{
    s = -4*k*m*(pow((1.-pow(m,2)),2)) + 4*k*m*(1.-pow(m,2))*((1.+2*(pow(m,2))));
  }
  return s;
}
// equation 14 SI
double ds_on_dk(double k,double v,double v_min,double delta)
{
  double s, m;
  m = (v-v_min)/(delta);
  if(m < 0.){
    s = 0.;
  }else if(m > 1.){
    s = 1.;
  }else{
    s = 1. - ( (pow((1.-pow(m,2)),2)) * ((1.+2*(pow(m,2)))) );
  }
  return s;
}


//declaration of all the variables
  int i, j, jj, k, l, N, gridpoint, protein, gridpoint1, gridpoint2;

  real COM_X, COM_Y, COM_Z;// coordinates of the center of mass intial structure
  real jedi, volume, a, enclosure, b, hydrophobicity, c, connectivity, d, constant, grd_x, grd_y, grd_z;
  real resolution, Vg, mod_rij, mod_rjk, mod_rik, pe, sum, min_dist, penalty_close, penalty_far;
  real ncoord, cutoff_hull, cutoff_close, cutoff_far, cutoff_enclosure, D_far, dump_matrix, n_grid, b_grid;
  real D_close, D_enclosure, cutoff_surface, cutoff_contact, cutoff_hydro, cutoff_con;
  real s1, s2, s3, s4, current, time1, time2, Va, Vmax, Vmin, D_Vmax, D_Vmin,dump_time, prev_snap;
  real connectivity_tot, ncontact, apolar, polar, D_hydro, hydrophobicity_tot;
  real  surface, hull, ncontact_hull, D_hull, ncontact_surface, D_surface, D_contact;

  rvec rij,rjk,rik;

  int size_grid;

  unsigned int size_protein;
  unsigned int size_apolar;
  unsigned int size_polar;
  unsigned int size_connect_map;

  n_grid           =  colvar.n_grid[i_c]; // total number of grid points (double)
  size_grid        =  n_grid; // total number of grid points !! for loops and array !!(integer)

  int active_grid[size_grid];// array with the index of active grid point (0 if inactive)

  double penalty_list[size_grid];// array with the activity of each grid point
  double enclosure_list[size_grid];// array with the penalty of distant contact
  double mind_list[size_grid];// array with the penalty of close contact
  double s3_m[size_grid];// array with the derivative of distant contact
  double s4_m[size_grid];// array with the derivative of close contact
  double sum_list[size_grid];// array with the sum of the distances between each grid point and all atoms in the CV
  double min_dist_list[size_grid];// array with the minimum distance between each grid point and all atoms in the CV
  double NP[size_grid];// array with the number of protein atoms within around each grid points according to s_on
  double new_x[size_grid];// array with update of x coordinates of each grid points according to translation/rotation
  double new_y[size_grid];// array with update of y coordinates of each grid points according to translation/rotation
  double new_z[size_grid];// array with update of z coordinates of each grid points according to translation/rotation
  double connectivity_list[size_grid];// array with the number of active grid points around each grid points according to s_on
  double hull_list[size_grid];// array with the hull score of each grid points
  double surface_list[size_grid];// array with the surface score of each grid points
  double ncontact_surface_list[size_grid];// array number of protein atoms in contact with each grid point (for derivatives)
  double hydrophobicity_list[size_grid];// array hydrophobicity score of each grid point
  double apolarity[size_grid];// array number of apolar protein atoms around each grid point
  double polarity[size_grid];// array number of polar protein atoms around each grid point

  double ray[size_grid];// 

  char *command[200];

//initialisation of (riables)

  volume      = 0.;// score volume : total number of active and partially active grid points
  N           = 0;// to modify
  D_close     = 0.05;// delta close contact
  D_far       = 0.05;// delta distant contact
  D_enclosure = 5.0;// delta NP contact ( for distant contact)
  resolution  = 0.15;// grid spacing
  s1          = 0.0;// minimum volume according our drug-like dataset
  s2          = 0.0;// maximum volume according our drug-like dataset
  Va          = 0.0;// penalty non drug-like pockets ( s1 * s2 )
  Vmin        = colvar.cutoff_hull[i_c];// minimum volume if volume < Vmin --> s1 = 0
  D_Vmin      = 10.0;// delta minimum volume volume > Vmin + D_Vmin --> s1 = 1
  Vmax        = 90.0;// minimum volume if volume < Vmax --> s2 = 1
  D_Vmax      = 50.0;// delta maximum volume volume > Vmax + D_Vmin --> s1 = 0
  enclosure   = 0.;// score enclosure : surface / hull
  hull        = 0.;// score hull
  surface     = 0.;// score surface
  D_hull      = 3;// delta hull
  D_surface   = 0.05;// delta surface
  D_contact   = 3;// delta contact for surface
  jedi        =  0.;// jedi score :-)
  a           =  0.059*Vmax;//just to normalise the volume according the trainning dataset Vmax
  b           =  24.29;//coefficient hydrophobicity derived from PLS
//  c           =  0.8032;//coefficient connectivity derived from PLS
  constant    =  -13.39;//constant derived from PLS
  Vg          = resolution*resolution*resolution; // grid point volume nm^3
  time1       = mtd_data->dt;// time step from gromacs
  time2       = mtd_data->istep;// step i from gromacs
  current     = time1*time2;// current snapshot (in ps)

/////////////////////////////////////// CHECK TIME

  dump_matrix      =  colvar.dump_matrix[i_c];//paramater to control the output of jedi (just to check, could be set to 0)
  COM_X            =  colvar.COM_X[i_c];//coordinate x COM first snapshot
  COM_Y            =  colvar.COM_Y[i_c];//coordinate y COM first snapshot
  COM_Z            =  colvar.COM_Z[i_c];//coordinate z COM first snapshot
  b_grid           =  colvar.b_grid[i_c];//atom number of the first grid point
  cutoff_close     =  colvar.cutoff_close[i_c];//cutoff close contact
  cutoff_far       =  colvar.cutoff_far[i_c];//cutoff distant contact
  cutoff_enclosure =  colvar.cutoff_enclosure[i_c];//cutoff NP for distant contact
//  cutoff_hull      =  colvar.cutoff_hull[i_c];//cutoff hull
//  cutoff_surface   =  colvar.cutoff_surface[i_c];//cutoff surface
//  cutoff_contact   =  colvar.cutoff_contact[i_c];//cutoff contact for surface
  cutoff_hydro     =  colvar.cutoff_hydro[i_c];//cutoff hydrophobicity
//  cutoff_con       =  colvar.cutoff_con[i_c];//cutoff connectivity
  size_protein     =  colvar.list[i_c][0] + colvar.list[i_c][1];//total number of protein atoms in the CV
  size_apolar      =  colvar.list[i_c][0];//total number of apolar atoms in the CV
  size_polar       =  colvar.list[i_c][1];//total number of polar atoms in the CV

  size_connect_map   = 27;//maximum number of grid point surrounding a grid point ( used to speed up derivatives calculations)
  connectivity_tot   = 0.;//sum of the connectivity scores of each grid point
  connectivity       = 0.;//score connectivity
  hydrophobicity     = 0.;//score hydrophobicity
  hydrophobicity_tot = 0.;//sum of the hydrophobicity scores of each grid point
  D_hydro            = 0.05;

  double it,beta;
  double connectivity_map[size_grid][27];
  int short_list_grid[size_grid][27];

  beta = 5.0;
//-----------------------------------------
//-----------------------------------------
//
//       STEP1  :  Protein Diffusion  
//
//-----------------------------------------
//-----------------------------------------

//get grid coordinates :-)

int start, max;
start = b_grid;
max   = n_grid;
int nbr[max];
for (i=0; i<max; i++)
{
  j = b_grid+i-1;
  nbr[i] = j;
}
// put as input
// Check the translation of the center of mass

// calculate COM from current snapshot
    double COM_x, sum_x;
    double COM_y, sum_y;
    double COM_z, sum_z;
    double mass_tot;

    mass_tot = 0.0;
    sum_x    = 0.0;
    sum_y    = 0.0;
    sum_z    = 0.0;
    for(j=0;j<colvar.list[i_c][0] + colvar.list[i_c][1];j++) {
      i = colvar.cvatoms[i_c][j];
      sum_x    += ( mtd_data->mass[i] ) * ( mtd_data->pos[i][0] );
      sum_y    += ( mtd_data->mass[i] ) * ( mtd_data->pos[i][1] );
      sum_z    += ( mtd_data->mass[i] ) * ( mtd_data->pos[i][2] );
      mass_tot += mtd_data->mass[i];
    }

    COM_x = sum_x/mass_tot;
    COM_y = sum_y/mass_tot;
    COM_z = sum_z/mass_tot;
   // printf("COM_x: %lf COM_y: %lf COM_z: %lf",COM_x,COM_y,COM_z);

/*
// This part of the code was a first approach to compute
// the rotation matrix. However, without the library libmatheval
//           www.gnu.org/software/libmatheval/
// it was not an efficient approach ... It is much easier to compile
// PLUMED 2.0 with this library. Therefore, the implementation of
// a methodolgy to modify the grid point coordinates (e.g. Kabsch algorithm) may 
// be of interest to speed up a bit the calculation. This current version of 
// the code is using gromacs to compute the rotation matrix of the heavy atoms of the protein.
 
 
 //matrix of inertia
    double xx, xy, xz, yx, yy, yz, zx, zy, zz;

    xx  = 0.0;
    xy  = 0.0;
    xz  = 0.0;

    yx  = 0.0;
    yy  = 0.0;
    yz  = 0.0;

    zx  = 0.0;
    zy  = 0.0;
    zz  = 0.0;


    for(j=0;j<colvar.list[i_c][0] + colvar.list[i_c][1];j++) {
      i  = colvar.cvatoms[i_c][j];

      xx  +=   (mtd_data->mass[i] * ((mtd_data->pos[i][1] * mtd_data->pos[i][1]) + (mtd_data->pos[i][2] * mtd_data->pos[i][2])));
      xy  += - (mtd_data->mass[i] * mtd_data->pos[i][0] * mtd_data->pos[i][1]);
      xz  += - (mtd_data->mass[i] * mtd_data->pos[i][0] * mtd_data->pos[i][2]);

      yx  += - (mtd_data->mass[i] * mtd_data->pos[i][0] * mtd_data->pos[i][1]);
      yy  +=   (mtd_data->mass[i] * ((mtd_data->pos[i][0] * mtd_data->pos[i][0]) + (mtd_data->pos[i][2] * mtd_data->pos[i][2])));
      yz  += - (mtd_data->mass[i] * mtd_data->pos[i][0] * mtd_data->pos[i][2]);

      zx  += - (mtd_data->mass[i] * mtd_data->pos[i][0] * mtd_data->pos[i][2]);
      zy  += - (mtd_data->mass[i] * mtd_data->pos[i][1] * mtd_data->pos[i][2]);
      zz  +=   (mtd_data->mass[i] * ((mtd_data->pos[i][0] * mtd_data->pos[i][0]) + (mtd_data->pos[i][1] * mtd_data->pos[i][1]))); 

    }

    xx = xx - ( COM_y*COM_y + COM_z*COM_z ) * mass_tot;
    xy = xy + COM_x*COM_y*mass_tot;
    xz = xz + COM_x*COM_z*mass_tot;
    yx = yx + COM_x*COM_y*mass_tot;
    yy = yy - ( COM_x*COM_x + COM_z*COM_z ) * mass_tot;
    yz = yz + COM_y*COM_z*mass_tot;
    zx = zx + COM_x*COM_z*mass_tot;
    zy = zy + COM_y*COM_z*mass_tot;
    zz = zz - ( COM_x*COM_x + COM_y*COM_y ) * mass_tot;

printf("xx : %f\n",xx);
printf("xy : %f\n",xy);
printf("xz : %f\n",xz);

printf("yx : %f\n",yx);
printf("yy : %f\n",yy);
printf("yz : %f\n",yz);

printf("zx : %f\n",zx);
printf("zy : %f\n",zy);
printf("zz : %f\n",zz);

printf("number of atoms : %i\n",colvar.list[i_c][0] + colvar.list[i_c][1]);
printf("%f %f %f\n",xx,xy,xz);
printf("%f %f %f\n",yx,yy,yz);
printf("%f %f %f\n",zx,zy,zz);


  double lambda = 0.0;
//  double det    = 0.0;

  double const_det = (xx*yy*zz)-(xz*zx*yy)+(xz*yx*zy)+(xy*yz*zx)-(xy*yx*zz)-(yz*zy*xx);
//  det = -(lambda*lambda*lambda) + (lambda*lambda)*(yy+xx+zz) + lambda*(-(xx*yy)-(yy*zz)-(xx*zz)+(xz*zx)+(xy*yx)+(yz*zy)) + const_det;

  double a_eq,b_eq,c_eq,d_eq;
  double f, g, h, root1, root2, root3;
  double num_i, num_j, num_k, num_l, num_m, num_n, num_p;

  a_eq = -1.0;
  b_eq = yy+xx+zz;
  c_eq = -(xx*yy)-(yy*zz)-(xx*zz)+(xz*zx)+(xy*yx)+(yz*zy);
  d_eq = const_det;

  f=( (3.0*c_eq/a_eq)-((b_eq*b_eq)/(a_eq*a_eq)) ) / (3.0);
  g=( ( (2.0*b_eq*b_eq*b_eq)/(a_eq*a_eq*a_eq) ) - ((9.0*b_eq*c_eq)/(a_eq*a_eq)) + (27.0*d_eq/a_eq) )/(27.0);
  h=( (g*g)/4.0 )+ ( (f*f*f)/27.0 );
  if( h > 0.0){
    root1 = 0.0;
    root2 = 0.0;
    root3 = 0.0;
  }else if (h == 0.0 && f == 0.0 && g == 0.0){
    root1 = 0.0;
    root2 = 0.0;
    root3 = 0.0;    
  }else{
    num_i = sqrt(( (g*g)/4.0 ) - h);
    num_j = pow(num_i,(1/3.0));
    num_k = acos(-(g/(2.0*num_i)));
    num_l = -num_j;
    num_m = cos(num_k/3.0);
    num_n = sqrt(3.0) * sin(num_k/3.0);
    num_p = b_eq/(-3.0*a_eq);
    root1 = (2.0*num_j) * cos(num_k/3.0) - (b_eq/(3.0*a_eq));
    root2 = num_l*(num_m+num_n)+num_p;
    root3 = num_l*(num_m-num_n)+num_p;
  }

  printf("\nroot1 : %f\n",root1);
  printf("\nroot2 : %f\n",root2);
  printf("\nroot3 : %f\n",root3);

*/

// generate a pdb file of the last snapshot of the trajectory (prev_snap)
// to compute the rotation matrix.
    if (current > 0.0)
    {
    prev_snap = current-dump_matrix;
    sprintf(command,"echo 1 | trjconv-jedi -f md.xtc -o last.gro -b %f -e %f -s md.tpr > /dev/null 2>&1",prev_snap,prev_snap);
    system(command);
    }

// compute the rotation matrix using gromacs between the first and the last snapshots of the trajectory
// using the heavy atom of the entire protein (group 3 in the default index file of gromacs)
  system("echo 3 | g_rotmat-jedi -f last.gro -s md.tpr > /dev/null 2>&1");

  remove("rotmat.xvg");
  remove("coord.xvg");
  remove("#last.gro.1#");

  float score[9] = {10.}; 
  FILE* file = NULL;
  file = fopen("rotmat-jedi", "r");
 
  if (file != NULL)
  {
    fscanf(file, "%f %f %f %f %f %f %f %f %f\n", &score[0], &score[1], &score[2], &score[3], &score[4], &score[5], &score[6], &score[7], &score[8]);
    fclose(file);
  }

//consider grid point coordinates do not change during the simulation, we may optimise the code by storing information
  //printf("new_x new_y new_z\n");
  for(i=0;i<max;i++) {
    gridpoint = i;
    new_x[i]  = ( score[0]*(colvar.grid_pos[gridpoint][0]) + score[3]*(colvar.grid_pos[gridpoint][1]) + score[6]*(colvar.grid_pos[gridpoint][2]) ) + (COM_x-COM_X) ;
    new_y[i]  = ( score[1]*(colvar.grid_pos[gridpoint][0]) + score[4]*(colvar.grid_pos[gridpoint][1]) + score[7]*(colvar.grid_pos[gridpoint][2]) ) + (COM_y-COM_Y) ;
    new_z[i]  = ( score[2]*(colvar.grid_pos[gridpoint][0]) + score[5]*(colvar.grid_pos[gridpoint][1]) + score[8]*(colvar.grid_pos[gridpoint][2]) ) + (COM_z-COM_Z) ;
/*  printf("%5f %5f %5f", new_x[i], new_y[i], new_z[i]);*/
/*  printf("\n");*/
  }
  
/*
//consider grid point coordinates do not change during the simulation, we may optimise the code by storing information
  float score[9] = {10.}; 
  for(i=0;i<max;i++) {
    gridpoint = i;
    new_x[i]  = colvar.grid_pos[gridpoint][0];
    new_y[i]  = colvar.grid_pos[gridpoint][1];
    new_z[i]  = colvar.grid_pos[gridpoint][2];
  }
*/

//-----------------------------------------
//-----------------------------------------
//
//       STEP2  :  JEDI score  
//
//-----------------------------------------
//-----------------------------------------

  double pe_list[26] = {0.0};
  double ray_score = .0;

//----------> Compute activity of grid points and also VOLUME
 // printf("%s %s %s\n","mind_dist","cutoff_close","D_close");
  for(i=0;i<max;i++) {
    N += 1;
    ncoord = 0.;
    sum = 0.;
    penalty_close = 0.;
    penalty_far = 0.;
    min_dist = 0.;
    grd_x = new_x[i];
    grd_y = new_y[i];
    grd_z = new_z[i];
    k = 0;
    // Store array coordinates that do no change in the inner loop calculation in a variable that is in the stack
    // e.g. xgrid_i = new_x[i] ; etc ...
//-----> apolar check
    for(j=0;j<colvar.list[i_c][0];j++) {
      protein   = colvar.cvatoms[i_c][j];
      rij[0]    = grd_x-mtd_data->pos[protein][0];
      rij[1]    = grd_y-mtd_data->pos[protein][1];
      rij[2]    = grd_z-mtd_data->pos[protein][2];
      // No PBC check !!  Code may need changes for PBC 
      mod_rij   = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
      if (mod_rij < 0.03){
        mod_rij = 0.03;
      }
/*      printf("grd x %f y %f z %f prot x %f y %f z %f \n",grd_x,grd_y,grd_z,mtd_data->pos[protein][0],mtd_data->pos[protein][1],mtd_data->pos[protein][2]);*/
/*      printf("i %d j %d mod_rij %f\n",i,j,mod_rij);*/
      // calculate mindist between grid points and protein atoms //
      // also get number of neighbors in the same pass ...
      ncoord += s_off(1.0,mod_rij,cutoff_far,D_far);
      pe = exp(beta/mod_rij);
      sum += pe;
    }
//-----> polar check

    for(j=colvar.list[i_c][0];j<colvar.list[i_c][0] + colvar.list[i_c][1];j++) {
      protein   = colvar.cvatoms[i_c][j];
      rij[0]    = grd_x-mtd_data->pos[protein][0];
      rij[1]    = grd_y-mtd_data->pos[protein][1];
      rij[2]    = grd_z-mtd_data->pos[protein][2];
      mod_rij   = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
      if (mod_rij < 0.03){
        mod_rij = 0.03;
      }
/*      printf("grd x %f y %f z %f prot x %f y %f z %f \n",grd_x,grd_y,grd_z,mtd_data->pos[protein][0],mtd_data->pos[protein][1],mtd_data->pos[protein][2]);*/
/*      printf("i %d j %d mod_rij %f\n",i,j,mod_rij);*/
      // calculate mindist between grid points and protein atoms //
      ncoord += s_off(1.0,mod_rij,cutoff_far,D_far);
      pe = exp(beta/mod_rij);
      sum += pe;
    }
    sum_list[i] = sum;
    min_dist = beta/log(sum);
    min_dist_list[i] = min_dist;
   //  printf("sum %f\n",sum);
    NP[i] = ncoord;
    mind_list[i]=s_on(1.0,min_dist,cutoff_close,D_close);
/*    printf("%lf %lf %lf\n",min_dist,cutoff_close,D_close);*/
/*    exit(-1);*/
    ncoord = 0.;
  }
  N = 0;
  //printf("%s %s %s\n","mind_list","enclosure_list","colvar.pocket");
  for(i=0;i<max;i++) {
    sum = .0;
    ray_score = 0.0;
    grd_x = new_x[i];
    grd_y = new_y[i];
    grd_z = new_z[i];
    for(j=0;j<44;j++){
      k = colvar.rays[i][j] - b_grid;
      if (colvar.rays[i][j] > 0){
        rij[0]    = grd_x-new_x[k];
        rij[1]    = grd_y-new_y[k];
        rij[2]    = grd_z-new_z[k];
        // No PBC check !!  Code may need changes for PBC 
        mod_rij   = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
        if (mod_rij < 0.03){
          mod_rij = 0.03;
        }
        ray_score += s_off(1.0,min_dist_list[k],0.15,cutoff_far) * s_on(1.0,mod_rij,0.25,0.05) * s_off(1.0,mod_rij,0.3,0.05);
      }else{
        break;
      }
    }
    ray[i]=ray_score;
    penalty_far = (s_on(1.0,ray_score,10.0,cutoff_enclosure));
    enclosure_list[i]  = (s_on(1.0,ray_score,10.0,cutoff_enclosure));
    penalty_list[i] = mind_list[i]*enclosure_list[i]*colvar.pocket[i];
    //printf("%lf %lf %lf\n",mind_list[i],enclosure_list[i],colvar.pocket[i]);
    
    
    

    if (0.0 < penalty_list[i]){
      active_grid[N] = i;
      N += 1;
    }
    volume += penalty_list[i];
/*    printf("The volume is %f\n",volume);*/
/*    printf("i %d penalty_list[i] %f mind_list[i] %f enclosure_list[i] %f pocket[i] %f volume %f\n",i,penalty_list[i],mind_list[i],enclosure_list[i],colvar.pocket[i],volume);*/
  }

// JCN Sep 2016. Printing activities of grid points
    FILE *actPoint;
    actPoint=fopen("grid_analysis.txt","w");
   
    fprintf(actPoint,"Size of penalty_list is: %i \n",sizeof(penalty_list));
    fprintf(actPoint,"Size of active_grid is: %i \n",sizeof(active_grid));
    fprintf(actPoint,"Value of size_grid: %i \n",size_grid);
     
    int isactive_grid=0;
    int isactive_penalty=0;
    for (int k=0; k<sizeof(active_grid); k++)
        {
         if(active_grid[k]>0)
            {
             isactive_grid++;
            } 
        }
    fprintf(actPoint,"%i values of active_grid, out of %i, are higher than 0\n",isactive_grid,sizeof(active_grid));

    for (int k=0; k<sizeof(penalty_list); k++)
        {
         if(active_grid[k]>0)
            {
             isactive_penalty++;
            }
        }
    fprintf(actPoint,"%i values of active_grid, out of %i, are higher than 0\n",isactive_penalty,sizeof(penalty_list));
    fclose(actPoint);
    

//----------> "Drug-like" volume
  s1 = (s_off(1.0,volume,Vmax,D_Vmax));
  s2 = (s_on(1.0,volume,Vmin,D_Vmin));
  Va = s1 * s2;

//  printf("%i\n",N);

//----------> Connectivity
/*
// Think how to make this work if the spacing was not 1 Angstrom
// Could loop only on active points 
  for(i=0;i<N;i++){
    ncontact = 0.;
    // No need to test distances as only looking at neighbors
    for(j=0;j<27;j++){
      k = colvar.neighbors[active_grid[i]][j] - b_grid;
      if(colvar.neighbors[active_grid[i]][j] > 0){
        rij[0]    = new_x[active_grid[i]]-new_x[k];
        rij[1]    = new_y[active_grid[i]]-new_y[k];
        rij[2]    = new_z[active_grid[i]]-new_z[k];
        mod_rij   = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
        connectivity_map[active_grid[i]][j] = mod_rij;
        short_list_grid[active_grid[i]][j] = k;
        ncontact += (s_off(1.0,mod_rij,0.15,0.15))*penalty_list[k];
      }else{
        break;
      }
    }
    connectivity_list[active_grid[i]] = ncontact;
    connectivity_tot += (s_on(1.0,ncontact,1.0,cutoff_con))*penalty_list[active_grid[i]];
  }
  if(volume > 0.){
    connectivity = connectivity_tot/volume;
  }else{
    connectivity = 0.;
  }
*///uncomment for connectivity calculation

//----------> Hydrophobicity

  for(i=0;i<N;i++) {
    apolar = 0.;
    polar  = 0.;
    for(j=0;j<colvar.list[i_c][0];j++) {
      protein   = colvar.cvatoms[i_c][j];
      rij[0]    = new_x[active_grid[i]]-mtd_data->pos[protein][0];
      rij[1]    = new_y[active_grid[i]]-mtd_data->pos[protein][1];
      rij[2]    = new_z[active_grid[i]]-mtd_data->pos[protein][2];
      mod_rij   = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
      apolar += s_off(1.0,mod_rij,cutoff_hydro,D_hydro);
    }
    apolarity[active_grid[i]] = apolar;
    for(j=colvar.list[i_c][0];j<colvar.list[i_c][0]+colvar.list[i_c][1];j++) {
      protein   = colvar.cvatoms[i_c][j];
      rij[0]    = new_x[active_grid[i]]-mtd_data->pos[protein][0];
      rij[1]    = new_y[active_grid[i]]-mtd_data->pos[protein][1];
      rij[2]    = new_z[active_grid[i]]-mtd_data->pos[protein][2];
      mod_rij   = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
      polar += s_off(1.0,mod_rij,cutoff_hydro,D_hydro);
    }
    polarity[active_grid[i]] = polar;
    if(polar + apolar > 0.){
      hydrophobicity_list[active_grid[i]]=apolar/(apolar+polar);   
    }else{
      hydrophobicity_list[active_grid[i]]=0.;   
    }
    hydrophobicity_tot+=hydrophobicity_list[active_grid[i]]*penalty_list[active_grid[i]];
  }
  if(volume > 0.){
    hydrophobicity = hydrophobicity_tot/volume;   
  }else{
    hydrophobicity = 0.0;
  }
/*printf("Va a volume Vmax b hydrophobicity constant:\n %f %f %f %f %f %f %f\n",Va,a,volume,Vmax,b,hydrophobicity,constant);*/

//----------> JEDI SCORE
  connectivity = 0.0;
  
//  jedi=Va*(a*volume/Vmax+b*hydrophobicity+c*connectivity+constant);//JEDI score with connectivity
  jedi=Va*(a*volume/Vmax+b*hydrophobicity+constant);//JEDI score without connectivity
  colvar.ss0[i_c] = jedi;

//-------> output descriptors and grid update every dumb_matrix (in ps)

  if( fmod(current,dump_matrix) == 0.0 )
  {
      FILE* output_derivatives = NULL;
      FILE* grid_dump = NULL;
//      FILE* prot_dump = NULL; 

// write output file output.txt
      output_derivatives = fopen("output.txt", "a");
      if (output_derivatives != NULL)
      {
          fprintf(output_derivatives,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",current, jedi, Va, hydrophobicity, volume/Vmax, COM_x, COM_y, COM_z, score[0], score[1], score[2], score[3], score[4], score[5], score[6], score[7], score[8]);
//          fprintf(output_derivatives,"%f %f %i\n",volume,hydrophobicity,N);
          fclose(output_derivatives);
      }
  }
/*
      char *name[20];
      i = time2 ;
      sprintf(name,"grid_JEDI.pdb",i);
      grid_dump = fopen(name, "w");
      for(i=0;i<max;i++) {
        fprintf(grid_dump,"%f %f %f %f\n",mind_list[i],penalty_list[i],connectivity_list[i],ray[i]);
      }
      fclose(grid_dump);

      char *name[20];
      i = time2 ;
      sprintf(name,"grid%d.gro",i);
      grid_dump = fopen(name, "w");
      fprintf(grid_dump,"Grid generated during the JEDI run\n %d\n",max);
      for(i=0;i<max;i++) {
        fprintf(grid_dump,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",i,"GRD","PR",i,new_x[i],new_y[i],new_z[i]);
      }
      fprintf(grid_dump,"   0.00000   0.00000   0.00000\n");

      char *name2[20];
      i = time2 ;
      sprintf(name2,"prot%d.gro",i);
      prot_dump = fopen(name2, "w");
      fprintf(prot_dump,"Protein generated during the JEDI run\n %d\n",colvar.list[i_c][0]+colvar.list[i_c][1]);
      for(i=0;i<colvar.list[i_c][0]+colvar.list[i_c][1];i++) {
        protein   = colvar.cvatoms[i_c][i];
        fprintf(prot_dump,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",i,"GRD","PR",i,mtd_data->pos[protein][0],mtd_data->pos[protein][1],mtd_data->pos[protein][2]);
      }
      fprintf(prot_dump,"   0.00000   0.00000   0.00000\n");
*/
/*exit(0);*/
//-----------------------------------------
//-----------------------------------------
//
//       STEP 3  :  Derivatives  
//
//-----------------------------------------
//-----------------------------------------



// Derivatives should be computed at the same time as potential 

  double dV_dx;
  double ds1_dx, ds1_dy, ds1_dz;//derivatives Vmin
  double ds2_dx, ds2_dy, ds2_dz;//derivatives Vmax
  double da_dx, da_dy, da_dz;//derivatives activity
  double ds2_dm, dAp_dm, dAp_dx, dP_dx;

  double dVa_dx_apolar[colvar.list[i_c][0]];
  double dVa_dx_polar[colvar.list[i_c][1]];

  double dV_dx_apolar[colvar.list[i_c][0]];
  double dV_dx_polar[colvar.list[i_c][1]];

  double dVa_dy_apolar[colvar.list[i_c][0]];
  double dVa_dy_polar[colvar.list[i_c][1]];

  double dV_dy_apolar[colvar.list[i_c][0]];
  double dV_dy_polar[colvar.list[i_c][1]];

  double dVa_dz_apolar[colvar.list[i_c][0]];
  double dVa_dz_polar[colvar.list[i_c][1]];

  double dV_dz_apolar[colvar.list[i_c][0]];
  double dV_dz_polar[colvar.list[i_c][1]];

  double dH_dx_apolar[colvar.list[i_c][0]];
  double dH_dx_polar[colvar.list[i_c][1]];

  double dH_dy_apolar[colvar.list[i_c][0]];
  double dH_dy_polar[colvar.list[i_c][1]];

  double dH_dz_apolar[colvar.list[i_c][0]];
  double dH_dz_polar[colvar.list[i_c][1]];

  double dHa_dx, dHa_dy, dHa_dz, dP_dm, ds1_dm;

  int index;
/*
  for(j=0;j<colvar.list[i_c][0];j++) {
    dVa_dx_apolar[j] = 0.0;
    dV_dx_apolar[j]  = 0.0;
    dVa_dy_apolar[j] = 0.0;
    dV_dy_apolar[j]  = 0.0;
    dVa_dz_apolar[j] = 0.0;
    dV_dz_apolar[j]  = 0.0;
    dH_dx_apolar[j]  = 0.0;
    dH_dy_apolar[j]  = 0.0;
    dH_dz_apolar[j]  = 0.0;
  }
  for(j=colvar.list[i_c][0];j<colvar.list[i_c][0]+colvar.list[i_c][1];j++) {
    dVa_dx_polar[j] = 0.0;
    dV_dx_polar[j]  = 0.0;
    dVa_dy_polar[j] = 0.0;
    dV_dy_polar[j]  = 0.0;
    dVa_dz_polar[j] = 0.0;
    dV_dz_polar[j]  = 0.0;
    dH_dx_polar[j]  = 0.0;
    dH_dy_polar[j]  = 0.0;
    dH_dz_polar[j]  = 0.0;
  }
*/
  double da_dx_apolar[1200][size_apolar];// !!!! max 2000 active grid points and 500 apolar atom in the CV
  double da_dy_apolar[1200][size_apolar];// !!!! max 2000 active grid points and 500 apolar atom in the CV
  double da_dz_apolar[1200][size_apolar];// !!!! max 2000 active grid points and 500 apolar atom in the CV

  double da_dx_polar[1200][size_polar];// !!!! max 2000 active grid points and 500 polar atom in the CV
  double da_dy_polar[1200][size_polar];// !!!! max 2000 active grid points and 500 polar atom in the CV
  double da_dz_polar[1200][size_polar];// !!!! max 2000 active grid points and 500 polar atom in the CV
//  printf("start derivatives\n");

/*
  for (int i = 0; i < 1500; i++) {
    for (int j = 0; j < size_apolar; j++) {
      da_dx_apolar[i][j] = 10.0;
      da_dy_apolar[i][j] = 10.0;
      da_dz_apolar[i][j] = 10.0;
    }
    for (int j = 0; j < size_polar; j++) {
      da_dx_polar[i][j] = 10.0;
      da_dy_polar[i][j] = 10.0;
      da_dz_polar[i][j] = 10.0;
//      printf("da_dx_apolar : %f\n",da_dx_apolar[i][j]);
//      printf("da_dy_apolar : %f\n",da_dy_apolar[i][j]);
//      printf("da_dz_apolar : %f\n",da_dz_apolar[i][j]);
//      printf("da_dx_polar : %f\n",da_dx_polar[i][j]);
//      printf("da_dy_polar : %f\n",da_dy_polar[i][j]);
//      printf("da_dz_polar : %f\n",da_dz_polar[i][j]);
    }
  }
*/
  double dSenclosure_dx, dSmind_dx, dSenclosure_dy, dSmind_dy, dSenclosure_dz, dSmind_dz;//dSmin: just for CC pemalty
  double dsmind_dx, dsmind_dy, dsmind_dz;//dsmin: just for solvent exposed pemalty
  double dmind_dx, dmind_dy, dmind_dz, dmind_dr, dr_dx, dr_dy, dr_dz,dSmind_dm,dSenclosure_dm,s1_ray,s2_ray;

  double dsum_enclosure_x, dsum_enclosure_y, dsum_enclosure_z, ds_off_dx, ds_on_dx, ds_off_dy, ds_on_dy, ds_off_dz, ds_on_dz;
//if more than 2000 active grid points --> dynamic allocation ?
/*
  double **da_dx_apolar;
  double **da_dy_apolar;
  double **da_dz_apolar;

  double **da_dx_polar;
  double **da_dy_polar;
  double **da_dz_polar;

  da_dx_apolar = (double**)malloc(size_grid * sizeof(double*));
  if (da_dx_apolar == NULL)
  {
    exit(0);
  }
  for (int i = 0; i < size_grid; i++) {
    da_dx_apolar[i] = (double*)malloc(size_apolar * sizeof(double));
    if (da_dx_apolar[i] == NULL)
    {
      exit(0);
    }
  }
  for (int i = 0; i < size_grid; i++) {
    for (int j = 0; j < size_apolar; j++) {
    	da_dx_apolar[i][j] = 0.;
    } 
  }


  da_dy_apolar = (double**)malloc(size_grid * sizeof(double*));
  if (da_dy_apolar == NULL)
  {
    exit(0);
  }
  for (int i = 0; i < size_grid; i++) {
    da_dy_apolar[i] = (double*)malloc(size_apolar * sizeof(double));
  }
  for (int i = 0; i < size_grid; i++) {
    for (int j = 0; j < size_apolar; j++) {
    	da_dy_apolar[i][j] = 0.;
    } 
  }


  da_dz_apolar = (double**)malloc(size_grid * sizeof(double*));
  if (da_dz_apolar == NULL)
  {
    exit(0);
  }
  for (int i = 0; i < size_grid; i++) {
    da_dz_apolar[i] = (double*)malloc(size_apolar * sizeof(double));
  }
  for (int i = 0; i < size_grid; i++) {
    for (int j = 0; j < size_apolar; j++) {
    	da_dz_apolar[i][j] = 0.;
    } 
  }


  da_dx_polar = (double**)malloc(size_grid * sizeof(double*));
  if (da_dx_polar == NULL)
  {
    exit(0);
  }
  for (int i = 0; i < size_grid; i++) {
    da_dx_polar[i] = (double*)malloc(size_polar * sizeof(double));
  }
  for (int i = 0; i < size_grid; i++) {
    for (int j = 0; j < size_polar; j++) {
    	da_dx_polar[i][j] = 0.;
    } 
  }


  da_dy_polar = (double**)malloc(size_grid * sizeof(double*));
  if (da_dy_polar == NULL)
  {
    exit(0);
  }
  for (int i = 0; i < size_grid; i++) {
    da_dy_polar[i] = (double*)malloc(size_polar * sizeof(double));
  }
  for (int i = 0; i < size_grid; i++) {
    for (int j = 0; j < size_polar; j++) {
    	da_dy_polar[i][j] = 0.;
    } 
  }

  da_dz_polar = (double**)malloc(size_grid * sizeof(double*));
  if (da_dz_polar == NULL)
  {
    exit(0);
  }
  for (int i = 0; i < size_grid; i++) {
    da_dz_polar[i] = (double*)malloc(size_polar * sizeof(double));
  }
  for (int i = 0; i < size_grid; i++) {
    for (int j = 0; j < size_polar; j++) {
    	da_dz_polar[i][j] = 0.;
    } 
  }
*/



//  output_derivatives = fopen("derivatives.txt", "a");

//----------> da_dx, da_dy, da_dz ( activity )

  for(i=0;i<N;i++) {
    grd_x=new_x[active_grid[i]];
    grd_y=new_y[active_grid[i]];
    grd_z=new_z[active_grid[i]];
// da_dx_apolar
      for(j=0;j<colvar.list[i_c][0];j++) {
//        printf("%i %i\n\n\n",i,j);
        protein   = colvar.cvatoms[i_c][j];
        rij[0]    = new_x[active_grid[i]]-mtd_data->pos[protein][0];
        rij[1]    = new_y[active_grid[i]]-mtd_data->pos[protein][1];
        rij[2]    = new_z[active_grid[i]]-mtd_data->pos[protein][2];
        mod_rij   = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
        if (mod_rij < 0.03){
          mod_rij = 0.03;
        }
//        fprintf(output_derivatives,"rij[0]    : %f\n",rij[0]);
//        fprintf(output_derivatives,"rij[1]    : %f\n",rij[1]);
//        fprintf(output_derivatives,"rij[2]    : %f\n",rij[2]);
//        fprintf(output_derivatives,"mod_rij   : %f\n",mod_rij);
        dsum_enclosure_x = 0.0;

        dmind_dr = ((beta*beta)*(exp(beta/mod_rij))) / ((mod_rij*mod_rij) * (sum_list[active_grid[i]]) * pow((log(sum_list[active_grid[i]])),2));
        dSmind_dm = ds_on_dm(1.0,min_dist_list[active_grid[i]],cutoff_close,D_close) * 1.0/D_close;

        dr_dx  = - (new_x[active_grid[i]]-mtd_data->pos[protein][0]) / (mod_rij);
        dmind_dx  = dmind_dr*(dr_dx);
        dSmind_dx = dSmind_dm * dmind_dx;
        dsum_enclosure_y = 0.0;

        dr_dy  = - (new_y[active_grid[i]]-mtd_data->pos[protein][1]) / (mod_rij);
        dmind_dy  = dmind_dr*(dr_dy);
        dSmind_dy = dSmind_dm * dmind_dy;
        dsum_enclosure_z = 0.0;

        dr_dz  = - (new_z[active_grid[i]]-mtd_data->pos[protein][2]) / (mod_rij);
        dmind_dz  = dmind_dr*(dr_dz);
        dSmind_dz = dSmind_dm * dmind_dz;

//        fprintf(output_derivatives,"dr_dx      : %f\n",dr_dx);
//        fprintf(output_derivatives,"dmind_dx   : %f\n",dmind_dx);
//        fprintf(output_derivatives,"ddSmind_dx : %f\n",dSmind_dx);

//        fprintf(output_derivatives,"dr_dy      : %f\n",dr_dy);
//        fprintf(output_derivatives,"dmind_dy   : %f\n",dmind_dy);
//        fprintf(output_derivatives,"ddSmind_dy : %f\n",dSmind_dy);

//        fprintf(output_derivatives,"dr_dz      : %f\n",dr_dz);
//        fprintf(output_derivatives,"dmind_dz   : %f\n",dmind_dz);
//        fprintf(output_derivatives,"ddSmind_dz : %f\n",dSmind_dz);
        for(l=0;l<44;l++){
          if (colvar.rays[active_grid[i]][l] > 0){
            k = colvar.rays[active_grid[i]][l] - b_grid;
            rjk[0]    = new_x[k]-mtd_data->pos[protein][0];//protein j and grid point k
            rjk[1]    = new_y[k]-mtd_data->pos[protein][1];//protein j and grid point k
            rjk[2]    = new_z[k]-mtd_data->pos[protein][2];//protein j and grid point k
            rik[0]    = new_x[k]-grd_x;
            rik[1]    = new_y[k]-grd_y;
            rik[2]    = new_z[k]-grd_z;
           // No PBC check !!  Code may need changes for PBC 
            mod_rjk   = sqrt(rjk[0]*rjk[0]+rjk[1]*rjk[1]+rjk[2]*rjk[2]);//protein j and grid point k
            mod_rik   = sqrt(rik[0]*rik[0]+rik[1]*rik[1]+rik[2]*rik[2]);//grid point i and grid point k
            if (mod_rjk < 0.03){
              mod_rjk = 0.03;
            }
            if (mod_rik < 0.03){
              mod_rik = 0.03;
            }
            s1_ray    = s_on(1.0,mod_rik,0.25,0.05);
            s2_ray    = s_off(1.0,mod_rik,0.3,0.05);
//            s1_ray    = 1.0;
//            s2_ray    = 1.0;
            dmind_dr  = ((beta*beta)*(exp(beta/mod_rjk))) / ((mod_rjk*mod_rjk) * (sum_list[k]) * pow((log(sum_list[k])),2));
            dSmind_dm = ds_off_dm(1.0,min_dist_list[k],0.15,cutoff_far) * 1.0/cutoff_far ;
//with respect to x
            dr_dx     = - (new_x[k]-mtd_data->pos[protein][0]) / (mod_rjk);
            dmind_dx  = dmind_dr*(dr_dx);
            dsmind_dx = dSmind_dm * dmind_dx;
            dsum_enclosure_x += (dsmind_dx * ( s1_ray * s2_ray ));

//            fprintf(output_derivatives,"mod_rjk          : %f\n",mod_rjk);
//            fprintf(output_derivatives,"dr_dx            : %f\n",dr_dx);
//            fprintf(output_derivatives,"dmind_dx         : %f\n",dmind_dx);
//            fprintf(output_derivatives,"sum_list[k]      : %f\n",(sum_list[k]));
//            fprintf(output_derivatives,"min_dist_list    : %f\n",min_dist_list[k]);
//            fprintf(output_derivatives,"cutoff_far       : %f\n",cutoff_far);
//            fprintf(output_derivatives,"dmind_dx         : %f\n",dmind_dx);
//            fprintf(output_derivatives,"dsmind_dx        : %f\n",dsmind_dx);
//            fprintf(output_derivatives,"dsum_enclosure_x : %f\n",dsum_enclosure_x);
//            fprintf(output_derivatives,"ds_off_dm        : %f\n",ds_off_dm(1.0,min_dist_list[k],0.15,cutoff_far));
//            fprintf(output_derivatives,"mod_rik          : %f\n",mod_rik);
//            fprintf(output_derivatives,"s_on             : %f\n",s_on(1.0,mod_rik,0.25,0.05));
//            fprintf(output_derivatives,"s_off            : %f\n\n",s_off(1.0,mod_rik,0.3,0.05));
//with respect to y
            dr_dy     = - (new_y[k]-mtd_data->pos[protein][1]) / (mod_rjk);
            dmind_dy  = dmind_dr*(dr_dy);
            dsmind_dy = dSmind_dm * dmind_dy;
            dsum_enclosure_y += (dsmind_dy * ( s1_ray * s2_ray ));
//            fprintf(output_derivatives,"dr_dy            : %f\n",dr_dy);
//            fprintf(output_derivatives,"min_dist_list    : %f\n",min_dist_list[k]);
//            fprintf(output_derivatives,"cutoff_far       : %f\n",cutoff_far);
//            fprintf(output_derivatives,"mod_rjk          : %f\n",mod_rjk);
//            fprintf(output_derivatives,"sum_list[k]      : %f\n",(sum_list[k]));
//            fprintf(output_derivatives,"dmind_dx         : %f\n",dmind_dx);
//            fprintf(output_derivatives,"dsmind_dx        : %f\n",dsmind_dx);
//            fprintf(output_derivatives,"dsum_enclosure_x : %f\n",dsum_enclosure_x);
//            fprintf(output_derivatives,"ds_off_dm        : %f\n",ds_off_dm(1.0,min_dist_list[k],0.15,cutoff_far));
//            fprintf(output_derivatives,"mod_rik          : %f\n",mod_rik);
//            fprintf(output_derivatives,"s_on             : %f\n",s_on(1.0,mod_rik,0.25,0.05));
//            fprintf(output_derivatives,"s_off            : %f\n\n",s_off(1.0,mod_rik,0.3,0.05));

//with respect to z
            dr_dz     = - (new_z[k]-mtd_data->pos[protein][2]) / (mod_rjk);
            dmind_dz  = dmind_dr*(dr_dz);
            dsmind_dz = dSmind_dm * dmind_dz;
            dsum_enclosure_z += (dsmind_dz * ( s1_ray * s2_ray ));
          }else{
            break;
          }
        }
        dSenclosure_dm     = ds_on_dm(1.0,ray[active_grid[i]],10.0,cutoff_enclosure);
        dSenclosure_dx     = dSenclosure_dm * 1.0/cutoff_enclosure * dsum_enclosure_x;
        da_dx_apolar[i][j] = ( enclosure_list[active_grid[i]] * dSmind_dx + mind_list[active_grid[i]] * dSenclosure_dx) * colvar.pocket[active_grid[i]];
//        fprintf(output_derivatives,"dsum_enclosure_x    : %f\n",dsum_enclosure_x);
//        fprintf(output_derivatives,"ds_on_dm    : %f\n",ds_on_dm(1.0,ray[active_grid[i]],10.0,cutoff_enclosure));
//        fprintf(output_derivatives,"enclosure_list    : %f\n",enclosure_list[active_grid[i]]);
//        fprintf(output_derivatives,"dSmind_dx         : %f\n",dSmind_dx);
//        fprintf(output_derivatives,"mind_list         : %f\n",mind_list[active_grid[i]]);
//        fprintf(output_derivatives,"dSenclosure_dx    : %f\n",dSenclosure_dx);
//        fprintf(output_derivatives,"colvar.pocket     : %f\n\n",colvar.pocket[active_grid[i]]);
//        fprintf(output_derivatives,"da_dx_apolar : %f\n",da_dx_apolar[i][j]);
        dSenclosure_dy     = dSenclosure_dm * 1.0/cutoff_enclosure * dsum_enclosure_y;
        da_dy_apolar[i][j] = ( enclosure_list[active_grid[i]] * dSmind_dy + mind_list[active_grid[i]] * dSenclosure_dy) *colvar.pocket[active_grid[i]];
//        fprintf(output_derivatives,"enclosure_list    : %f\n",enclosure_list[active_grid[i]]);
//        fprintf(output_derivatives,"dSmind_dy         : %f\n",dSmind_dy);
//        fprintf(output_derivatives,"mind_list         : %f\n",mind_list[active_grid[i]]);
//        fprintf(output_derivatives,"dSenclosure_dy    : %f\n",dSenclosure_dy);
//        fprintf(output_derivatives,"colvar.pocket     : %f\n\n",colvar.pocket[active_grid[i]]);
//        fprintf(output_derivatives,"da_dy_apolar : %f\n",da_dy_apolar[i][j]);
        dSenclosure_dz     = dSenclosure_dm * 1.0/cutoff_enclosure * dsum_enclosure_z;
        da_dz_apolar[i][j] = ( enclosure_list[active_grid[i]] * dSmind_dz + mind_list[active_grid[i]] * dSenclosure_dz) *colvar.pocket[active_grid[i]];
        dsum_enclosure_x   = 0.;
        dsum_enclosure_y   = 0.;
        dsum_enclosure_z   = 0.;
//        fprintf(output_derivatives,"da_dx_apolar : %f\n",( enclosure_list[active_grid[i]] * dSmind_dx + mind_list[active_grid[i]] * dSenclosure_dx) * colvar.pocket[active_grid[i]]);
//        fprintf(output_derivatives,"da_dy_apolar : %f\n",( enclosure_list[active_grid[i]] * dSmind_dy + mind_list[active_grid[i]] * dSenclosure_dy) * colvar.pocket[active_grid[i]]);
//        fprintf(output_derivatives,"da_dz_apolar : %f\n",( enclosure_list[active_grid[i]] * dSmind_dz + mind_list[active_grid[i]] * dSenclosure_dz) * colvar.pocket[active_grid[i]]);
//        fprintf(output_derivatives,"da_dx_apolar : %f\n",da_dx_apolar[i][j]);
//        fprintf(output_derivatives,"da_dy_apolar : %f\n",da_dy_apolar[i][j]);
//        fprintf(output_derivatives,"da_dz_apolar : %f\n",da_dz_apolar[i][j]);
      }
// da_dx_polar
      for(j=colvar.list[i_c][0];j<colvar.list[i_c][0]+colvar.list[i_c][1];j++) {
        jj = j - colvar.list[i_c][0];
        protein   = colvar.cvatoms[i_c][j];
        rij[0]    = new_x[active_grid[i]]-mtd_data->pos[protein][0];
        rij[1]    = new_y[active_grid[i]]-mtd_data->pos[protein][1];
        rij[2]    = new_z[active_grid[i]]-mtd_data->pos[protein][2];
        mod_rij   = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
        if (mod_rij < 0.03){
          mod_rij = 0.03;
        }
        dmind_dr = ((beta*beta)*(exp(beta/mod_rij))) / ((mod_rij*mod_rij) * (sum_list[active_grid[i]]) * pow((log(sum_list[active_grid[i]])),2));
        dSmind_dm = ds_on_dm(1.0,min_dist_list[active_grid[i]],cutoff_close,D_close) * 1.0/D_close;
        dsum_enclosure_x = 0.0;
        dr_dx  = - (new_x[active_grid[i]]-mtd_data->pos[protein][0]) / (mod_rij);
        dmind_dx  = dmind_dr*(dr_dx);
        dSmind_dx = dSmind_dm * dmind_dx;
        dsum_enclosure_y = 0.0;
        dr_dy  = - (new_y[active_grid[i]]-mtd_data->pos[protein][1]) / (mod_rij);
        dmind_dy  = dmind_dr*(dr_dy);
        dSmind_dy = dSmind_dm * dmind_dy;
        dsum_enclosure_z = 0.0;
        dr_dz     = - (new_z[active_grid[i]]-mtd_data->pos[protein][2]) / (mod_rij);
        dmind_dz  = dmind_dr*(dr_dz);
        dSmind_dz = dSmind_dm * dmind_dz;
        for(l=0;l<44;l++){
          k = colvar.rays[active_grid[i]][l] - b_grid;
          if (colvar.rays[active_grid[i]][l] > 0){
            rjk[0]    = new_x[k]-mtd_data->pos[protein][0];//protein j and grid point k
            rjk[1]    = new_y[k]-mtd_data->pos[protein][1];//protein j and grid point k
            rjk[2]    = new_z[k]-mtd_data->pos[protein][2];//protein j and grid point k
            rik[0]    = new_x[k]-grd_x;
            rik[1]    = new_y[k]-grd_y;
            rik[2]    = new_z[k]-grd_z;
           // No PBC check !!  Code may need changes for PBC 
            mod_rjk   = sqrt(rjk[0]*rjk[0]+rjk[1]*rjk[1]+rjk[2]*rjk[2]);//protein j and grid point k
            mod_rik   = sqrt(rik[0]*rik[0]+rik[1]*rik[1]+rik[2]*rik[2]);//grid point i and grid point k
            if (mod_rjk < 0.03){
              mod_rjk = 0.03;
            }
            s1_ray    = s_on(1.0,mod_rik,0.25,0.05);
            s2_ray    = s_off(1.0,mod_rik,0.3,0.05);
//            s1_ray    = 1.0;
//            s2_ray    = 1.0;
            dmind_dr  = ((beta*beta)*(exp(beta/mod_rjk))) / ((mod_rjk*mod_rjk) * (sum_list[k]) * pow((log(sum_list[k])),2));
            dSmind_dm = ds_off_dm(1.0,min_dist_list[k],0.15,cutoff_far)* 1.0/cutoff_far;
//with respect to x
            dr_dx     = - (new_x[k]-mtd_data->pos[protein][0]) / (mod_rjk);
            dmind_dx  = dmind_dr*(dr_dx);
            dsmind_dx = dSmind_dm * dmind_dx;
            dsum_enclosure_x += (dsmind_dx * ( s1_ray * s2_ray ));
//with respect to y
            dr_dy     = - (new_y[k]-mtd_data->pos[protein][1]) / (mod_rjk);
            dmind_dy  = dmind_dr*(dr_dy);
            dsmind_dy = dSmind_dm * dmind_dy;
            dsum_enclosure_y += (dsmind_dy * ( s1_ray * s2_ray ));
//with respect to z
            dr_dz     = - (new_z[k]-mtd_data->pos[protein][2]) / (mod_rjk);
            dmind_dz  = dmind_dr*(dr_dz);
            dsmind_dz = dSmind_dm * dmind_dz;
            dsum_enclosure_z += (dsmind_dz * ( s1_ray * s2_ray ));
          }else{
            break;
          }
        }
        dSenclosure_dm     = ds_on_dm(1.0,ray[active_grid[i]],10.0,cutoff_enclosure) * 1.0/cutoff_enclosure;
        dSenclosure_dx     = dSenclosure_dm * dsum_enclosure_x;
        da_dx_polar[i][jj] = ( enclosure_list[active_grid[i]] * dSmind_dx + mind_list[active_grid[i]] * dSenclosure_dx) * colvar.pocket[active_grid[i]];
        dSenclosure_dy     = dSenclosure_dm * dsum_enclosure_y;
        da_dy_polar[i][jj] = ( enclosure_list[active_grid[i]] * dSmind_dy + mind_list[active_grid[i]] * dSenclosure_dy) *colvar.pocket[active_grid[i]];
        dSenclosure_dz     = dSenclosure_dm * dsum_enclosure_z;
        da_dz_polar[i][jj] = ( enclosure_list[active_grid[i]] * dSmind_dz + mind_list[active_grid[i]] * dSenclosure_dz) *colvar.pocket[active_grid[i]];
      //  printf("i %d jj %d da_dz_polar %f enclosure_list[active_grid[i]] %f dSmind_dz %f mind_list[active_grid[i]] %f dSenclosure_dz %f colvar.pocket[active_grid[i]] %f \n",i,jj,da_dz_polar[i][jj],enclosure_list[active_grid[i]],dSmind_dz,mind_list[active_grid[i]],dSenclosure_dz,colvar.pocket[active_grid[i]]);
        dsum_enclosure_x   = 0.;
        dsum_enclosure_y   = 0.;
        dsum_enclosure_z   = 0.;
//        printf("da_dx_polar : %f\n",( enclosure_list[active_grid[i]] * dSmind_dx + mind_list[active_grid[i]] * dSenclosure_dx) * colvar.pocket[active_grid[i]]);
//        printf("da_dy_polar : %f\n",( enclosure_list[active_grid[i]] * dSmind_dy + mind_list[active_grid[i]] * dSenclosure_dy) * colvar.pocket[active_grid[i]]);
//        printf("da_dz_polar : %f\n",( enclosure_list[active_grid[i]] * dSmind_dz + mind_list[active_grid[i]] * dSenclosure_dz) * colvar.pocket[active_grid[i]]);
//        printf("da_dx_polar : %f\n",da_dx_polar[i][jj]);
//        printf("da_dy_polar : %f\n",da_dy_polar[i][jj]);
//        printf("da_dz_polar : %f\n",da_dz_polar[i][jj]);
      }
    }

//----------> dC_dx, dC_dy, dC_dz ( connectivity )


//  real dCi_dx, dsCON_dm, dCON_dx, dCON_dy, dCON_dz;
//  real ds5_dak;
//  real ds5_da;
//  real dCa_dx, dCa_dy, dCa_dz;

//  double dC_dx_apolar[colvar.list[i_c][0]];
//  double dC_dx_polar[colvar.list[i_c][1]];

//  double dC_dy_apolar[colvar.list[i_c][0]];
//  double dC_dy_polar[colvar.list[i_c][1]];

//  double dC_dz_apolar[colvar.list[i_c][0]];
//  double dC_dz_polar[colvar.list[i_c][1]];

//  double dCON_dx_apolar[colvar.list[i_c][0]];
//  double dCON_dx_polar[colvar.list[i_c][1]];

//  double dCON_dy_apolar[colvar.list[i_c][0]];
//  double dCON_dy_polar[colvar.list[i_c][1]];

//  double dCON_dz_apolar[colvar.list[i_c][0]];
//  double dCON_dz_polar[colvar.list[i_c][1]];


//  real dsCON_dx, dsCON_dy, dsCON_dz;

//----------> dH_dx, dH_dy, dH_dz, dE_dx, dE_dy, dE_dz ( hydrophobicity / enclosure )
//  FILE* output_derivatives = NULL;
//  output_derivatives = fopen("deivatives.txt", "a");
//  fprintf(output_derivatives,"APOLAR ATOMS\n");

//apolar
  for(j=0;j<colvar.list[i_c][0];j++) {
    jj      = j ;
    dHa_dx  = 0.;
    dHa_dy  = 0.;
    dHa_dz  = 0.;
    da_dx   = 0.;
    da_dy   = 0.;
    da_dz   = 0.;
//      printf("dHa_dx : %f\n",dHa_dx);
//      printf("dHa_dy : %f\n",dHa_dy);
//      printf("dHa_dz : %f\n",dHa_dz);
    protein = colvar.cvatoms[i_c][j];
    for(i=0;i<N;i++) {
        gridpoint = active_grid[i];
        da_dx += da_dx_apolar[i][jj];
        da_dy += da_dy_apolar[i][jj];
        da_dz += da_dz_apolar[i][jj];
//        fprintf(output_derivatives,"da_dx_apolar : %f\n",da_dx_apolar[i][jj]);
//        fprintf(output_derivatives,"da_dy_apolar : %f\n",da_dy_apolar[i][jj]);
//        fprintf(output_derivatives,"da_dz_apolar : %f\n",da_dz_apolar[i][jj]);

//        dsCON_dm = ds_on_dm(penalty_list[gridpoint],connectivity_list[gridpoint],1.0,cutoff_con);
//        dCON_dx = 0.0;
//        dCON_dy = 0.0;
//        dCON_dz = 0.0;
//        for(k=0;k<27;k++) { // !!!!!!! Check the indexes, it is not the same when static allocation ...
//          if(colvar.neighbors[gridpoint][k] > 0){
//            index = colvar.neighbors[gridpoint][k] - b_grid;
//            ds5_da = ds_off_dk(penalty_list[index],connectivity_map[gridpoint][k],0.15,0.15);
//            dCON_dx += ds5_da * da_dx_apolar[index][jj];
//            dCON_dy += ds5_da * da_dy_apolar[index][jj];
//            dCON_dz += ds5_da * da_dz_apolar[index][jj];
//          }else{
//            break;
//          }
//        }
//        dCi_dx = dsCON_dm * (1.0/cutoff_con) * dCON_dx + ds_off_dk(penalty_list[gridpoint],connectivity_list[gridpoint],1.0,cutoff_con) * da_dx_apolar[gridpoint][jj];
//        dCa_dx += dCi_dx;
//        dCi_dx = dsCON_dm * (1.0/cutoff_con) * dCON_dy + ds_off_dk(penalty_list[gridpoint],connectivity_list[gridpoint],1.0,cutoff_con) * da_dy_apolar[gridpoint][jj];
//        dCa_dy += dCi_dx;
//        dCi_dx = dsCON_dm * (1.0/cutoff_con) * dCON_dz + ds_off_dk(penalty_list[gridpoint],connectivity_list[gridpoint],1.0,cutoff_con) * da_dz_apolar[gridpoint][jj];
//        dCa_dz += dCi_dx;
      // Does this ever happen??
        if( (polarity[gridpoint]+apolarity[gridpoint]) > 0.){
          rij[0]    = new_x[gridpoint]-mtd_data->pos[protein][0];
          rij[1]    = new_y[gridpoint]-mtd_data->pos[protein][1];
          rij[2]    = new_z[gridpoint]-mtd_data->pos[protein][2];
          mod_rij   = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
//          fprintf(output_derivatives,"da_dx     : %f\n",da_dx);
//          fprintf(output_derivatives,"rij[0]    : %f\n",rij[0]);
//          fprintf(output_derivatives,"rij[1]    : %f\n",rij[1]);
//          fprintf(output_derivatives,"rij[2]    : %f\n",rij[2]);
//          fprintf(output_derivatives,"mod_rij   : %f\n",mod_rij);
//    fprintf(output_derivatives,"mod_rij 2    : %f\n",mod_rij);
          dAp_dm  = ds_off_dm(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro);//same value for x,y,z
          dAp_dx  = ((dAp_dm) * (1.0/D_hydro) * (-(rij[0])/(mod_rij)) ) + ( ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro) * da_dx_apolar[i][jj]);//!!! dynamic allocation
          dHa_dx += (1.0/volume) * ( (dAp_dx*(polarity[gridpoint] + apolarity[gridpoint]) - (apolarity[gridpoint])*dAp_dx ) / ((polarity[gridpoint] + apolarity[gridpoint]) * (polarity[gridpoint] + apolarity[gridpoint])));
//          fprintf(output_derivatives,"dAp_dm    : %f\n",dAp_dm);
//          fprintf(output_derivatives,"dAp_dx    : %f\n",dAp_dx);
//          fprintf(output_derivatives,"dHa_dx    : %f\n",dHa_dx);
//          fprintf(output_derivatives,"dha_dx : %f\n",(1.0/volume) * ( (dAp_dx*(polarity[gridpoint] + apolarity[gridpoint]) - (apolarity[gridpoint])*dAp_dx ) / ((polarity[gridpoint] + apolarity[gridpoint]) * (polarity[gridpoint] + apolarity[gridpoint]))));
//          printf("apolarity : %f\n",apolarity[gridpoint]);
//          printf("polarity : %f\n",polarity[gridpoint]);
          dAp_dx  = ((dAp_dm) * (1.0/D_hydro) * (-(rij[1])/(mod_rij)) ) + ( ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro) * da_dy_apolar[i][jj]);//!!! dynamic allocation
//          printf("dAp_dm       : %f\n",dAp_dm);
//          printf("1.0/D_hydro  : %f\n",1.0/D_hydro);
//          printf("rij[1]       : %f\n",rij[1]);
//          printf("mod_rij      : %f\n",mod_rij);
//          printf("ds_off_dk    : %f\n",ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro));
//          printf("da_dy_apolar : %f\n",da_dy_apolar[i][jj]);
          dHa_dy += (1.0/volume) * ( (dAp_dx*(polarity[gridpoint] + apolarity[gridpoint]) - (apolarity[gridpoint])*dAp_dx ) / ((polarity[gridpoint] + apolarity[gridpoint]) * (polarity[gridpoint] + apolarity[gridpoint])));
//          printf("dHa_dx    : %f\n",dHa_dx);

          dAp_dx  = ((dAp_dm) * (1.0/D_hydro) * (-(rij[2])/(mod_rij)) ) + ( ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro) * da_dz_apolar[i][jj]);//!!! dynamic allocation
          dHa_dz += (1.0/volume) * ( (dAp_dx*(polarity[gridpoint] + apolarity[gridpoint]) - (apolarity[gridpoint])*dAp_dx ) / ((polarity[gridpoint] + apolarity[gridpoint]) * (polarity[gridpoint] + apolarity[gridpoint])));

        }
    }

//    dC_dx_apolar[jj] = (1.0/volume) * dCa_dx;
//    dC_dy_apolar[jj] = (1.0/volume) * dCa_dy;
//    dC_dz_apolar[jj] = (1.0/volume) * dCa_dz;
//    dCa_dx = 0.0;
//    dCa_dy = 0.0;
//    dCa_dz = 0.0;
//      printf("jj     : %i\n",jj);
//      fprintf(output_derivatives,"dHa_dx : %f\n",dHa_dx);
//      fprintf(output_derivatives,"dHa_dy : %f\n",dHa_dy);
//      fprintf(output_derivatives,"dHa_dz : %f\n",dHa_dz);

    dH_dx_apolar[jj] = dHa_dx;
    dH_dy_apolar[jj] = dHa_dy;
    dH_dz_apolar[jj] = dHa_dz;

    ds1_dm =ds_off_dm(1.0,volume,Vmax,D_Vmax) * (1.0/D_Vmax) * (Vg) ;

    dV_dx_apolar[jj] = (Vg) * (da_dx);
    dV_dy_apolar[jj] = (Vg) * (da_dy);
    dV_dz_apolar[jj] = (Vg) * (da_dz);

    ds1_dx = (ds1_dm) * (da_dx);
    ds1_dy = (ds1_dm) * (da_dy);
    ds1_dz = (ds1_dm) * (da_dz);

    ds2_dm = ds_on_dm(1.0,volume,Vmin,D_Vmin) * (1.0/D_Vmin) * (Vg) ;

    ds2_dx = (ds2_dm) * (da_dx);
    ds2_dy = (ds2_dm) * (da_dy);
    ds2_dz = (ds2_dm) * (da_dz);

    dVa_dx_apolar[jj] = s1 * ds2_dx + s2 * ds1_dx;
    dVa_dy_apolar[jj] = s1 * ds2_dy + s2 * ds1_dy;
    dVa_dz_apolar[jj] = s1 * ds2_dz + s2 * ds1_dz;
//-----  derivatives including connectivity
//    colvar.myder[i_c][j][0] = jedi * ((1.0/Va)*(dVa_dx_apolar[jj]) + (1.0/(a*volume+b*hydrophobicity+c*connectivity+constant)) * ( a*(dV_dx_apolar[jj]) + b*(dH_dx_apolar[jj]) + c*(dC_dx_apolar[jj]) ));
//    colvar.myder[i_c][j][1] = jedi * ((1.0/Va)*(dVa_dy_apolar[jj]) + (1.0/(a*volume+b*hydrophobicity+c*connectivity+constant)) * ( a*(dV_dy_apolar[jj]) + b*(dH_dy_apolar[jj]) + c*(dC_dy_apolar[jj]) ));
//    colvar.myder[i_c][j][2] = jedi * ((1.0/Va)*(dVa_dz_apolar[jj]) + (1.0/(a*volume+b*hydrophobicity+c*connectivity+constant)) * ( a*(dV_dz_apolar[jj]) + b*(dH_dz_apolar[jj]) + c*(dC_dz_apolar[jj]) ));
//-----  derivatives without connectivity

    colvar.myder[i_c][j][0] = jedi * ((1.0/Va)*(dVa_dx_apolar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dx_apolar[jj]) + b*(dH_dx_apolar[jj]) ));
/*     printf("jedi Va dVa_dx_apolar[jj] a volume b hydrophobicity constant dV_dx_apolar[jj] dH_dx_apolar[jj]: %f %f %f %f %f %f %f %f %f %f\n",jedi,Va,dVa_dx_apolar[jj],a,volume,b,hydrophobicity,constant,dV_dx_apolar[jj],dH_dx_apolar[jj]);*/
/*     printf("%f %f %f %f %f %f %f %f %f %f\n",Jedi,Va,dVa_dx_apolar[jj],a,volume,b,hydrophobicity,constant,dV_dx_apolar[jj],dH_dx_apolar[jj])*/
    colvar.myder[i_c][j][1] = jedi * ((1.0/Va)*(dVa_dy_apolar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dy_apolar[jj]) + b*(dH_dy_apolar[jj]) ));
    colvar.myder[i_c][j][2] = jedi * ((1.0/Va)*(dVa_dz_apolar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dz_apolar[jj]) + b*(dH_dz_apolar[jj]) ));
/*    printf("Apolar j %i is colvar.myder[i_c][j][0]: %f\n",j,colvar.myder[i_c][j][0]);*/
/*    printf("Apolar j %i is colvar.myder[i_c][j][1]: %f\n",j,colvar.myder[i_c][j][1]);*/
/*    printf("Apolar j %i is colvar.myder[i_c][j][2]: %f\n",j,colvar.myder[i_c][j][2]);*/
/*    exit(0);*/
  }

// polar
  for(j=colvar.list[i_c][0];j<colvar.list[i_c][0]+colvar.list[i_c][1];j++) {
    jj = j - colvar.list[i_c][0];
    protein   = colvar.cvatoms[i_c][j];
    dHa_dx = 0.;
    dHa_dy = 0.;
    dHa_dz = 0.;
    da_dx = 0.;
    da_dy = 0.;
    da_dz = 0.;
    for(i=0;i<N;i++) {
        gridpoint = active_grid[i];
        da_dx += da_dx_polar[i][jj];
        da_dy += da_dy_polar[i][jj];
        da_dz += da_dz_polar[i][jj];
//        dsCON_dm = ds_on_dm(penalty_list[gridpoint],connectivity_list[gridpoint],1.0,cutoff_con);
//        dCON_dx = 0.0;
//        dCON_dy = 0.0;
//        dCON_dz = 0.0;
//        for(k=0;k<27;k++) {
//          if(colvar.neighbors[gridpoint][k] > 0){
//            index = colvar.neighbors[gridpoint][k] - b_grid;
//            ds5_da = ds_off_dk(penalty_list[index],connectivity_map[gridpoint][k],0.15,0.15);
//            dCON_dx += ds5_da * da_dx_polar[index][jj];
//            dCON_dy += ds5_da * da_dy_polar[index][jj];
//            dCON_dz += ds5_da * da_dz_polar[index][jj];
//          }
//        }
//        dCi_dx = dsCON_dm * (1.0/cutoff_con) * dCON_dx + ds_off_dk(penalty_list[gridpoint],connectivity_list[gridpoint],1.0,cutoff_con) * da_dx_polar[gridpoint][jj];
//        dCa_dx += dCi_dx;
//        dCi_dx = dsCON_dm * (1.0/cutoff_con) * dCON_dy + ds_off_dk(penalty_list[gridpoint],connectivity_list[gridpoint],1.0,cutoff_con) * da_dy_polar[gridpoint][jj];
//        dCa_dy += dCi_dx;
//        dCi_dx = dsCON_dm * (1.0/cutoff_con) * dCON_dz + ds_off_dk(penalty_list[gridpoint],connectivity_list[gridpoint],1.0,cutoff_con) * da_dz_polar[gridpoint][jj];
//        dCa_dz += dCi_dx;
        if((polarity[gridpoint]+apolarity[gridpoint]) > 0.){
          rij[0]    = new_x[gridpoint]-mtd_data->pos[protein][0];
          rij[1]    = new_y[gridpoint]-mtd_data->pos[protein][1];
          rij[2]    = new_z[gridpoint]-mtd_data->pos[protein][2];
          mod_rij   = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
          max(mod_rij, 0.05);
          dP_dm   = ds_off_dm(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro);
          dP_dx   = ( (dP_dm) * (1.0/D_hydro) * (-(rij[0])/(mod_rij)) ) + ( ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro) * da_dx_polar[i][jj] );
          dHa_dx += (1.0/volume) * (-( ( apolarity[gridpoint] * dP_dx )/( (polarity[gridpoint]+apolarity[gridpoint])*(polarity[gridpoint]+apolarity[gridpoint]) ) ));
          dP_dx   = ( (dP_dm) * (1.0/D_hydro) * (-(rij[1])/(mod_rij)) ) + ( ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro) * da_dx_polar[i][jj] );
          dHa_dy += (1.0/volume) * (-( ( apolarity[gridpoint] * dP_dx )/( (polarity[gridpoint]+apolarity[gridpoint])*(polarity[gridpoint]+apolarity[gridpoint]) ) ));
          dP_dx   = ( (dP_dm) * (1.0/D_hydro) * (-(rij[2])/(mod_rij)) ) + ( ds_off_dk(penalty_list[gridpoint],mod_rij,cutoff_hydro,D_hydro) * da_dx_polar[i][jj] );
          dHa_dz += (1.0/volume) * (-( ( apolarity[gridpoint] * dP_dx )/( (polarity[gridpoint]+apolarity[gridpoint])*(polarity[gridpoint]+apolarity[gridpoint]) ) ));
        }
    }
//    dC_dx_polar[jj] = (1.0/volume) * dCa_dx;
//    dC_dy_polar[jj] = (1.0/volume) * dCa_dy;
//    dC_dz_polar[jj] = (1.0/volume) * dCa_dz;
//    dCa_dx = 0.0;
//    dCa_dy = 0.0;
//    dCa_dz = 0.0;

    dH_dx_polar[jj] = dHa_dx;
    dH_dy_polar[jj] = dHa_dy;
    dH_dz_polar[jj] = dHa_dz;

    jj = j - colvar.list[i_c][0];
    ds1_dm =ds_off_dm(1.0,volume,Vmax,D_Vmax) * (1.0/D_Vmax) * (Vg);

    dV_dx_polar[jj] = (Vg) * (da_dx);
    dV_dy_polar[jj] = (Vg) * (da_dy);
    dV_dz_polar[jj] = (Vg) * (da_dz);

    ds1_dx = (ds1_dm) * (da_dx);
    ds1_dy = (ds1_dm) * (da_dy);
    ds1_dz = (ds1_dm) * (da_dz);

    ds2_dm = ds_on_dm(1.0,volume,Vmin,D_Vmin) * (1.0/D_Vmin) * (Vg);

    ds2_dx = (ds2_dm) * (da_dx);
    ds2_dy = (ds2_dm) * (da_dy);
    ds2_dz = (ds2_dm) * (da_dz);

    dVa_dx_polar[jj] = s1 * ds2_dx + s2 * ds1_dx;
    dVa_dy_polar[jj] = s1 * ds2_dy + s2 * ds1_dy;
    dVa_dz_polar[jj] = s1 * ds2_dz + s2 * ds1_dz;
//-----  derivatives including connectivity  
//    colvar.myder[i_c][j][0] = jedi * ((1.0/Va)*(dVa_dx_polar[jj]) + (1.0/(a*volume+b*hydrophobicity+c*connectivity+constant)) * ( a*(dV_dx_polar[jj]) + b*(dH_dx_polar[jj]) + c*(dC_dx_polar[jj]) ));
//    colvar.myder[i_c][j][1] = jedi * ((1.0/Va)*(dVa_dy_polar[jj]) + (1.0/(a*volume+b*hydrophobicity+c*connectivity+constant)) * ( a*(dV_dy_polar[jj]) + b*(dH_dy_polar[jj]) + c*(dC_dy_polar[jj]) ));
//    colvar.myder[i_c][j][2] = jedi * ((1.0/Va)*(dVa_dz_polar[jj]) + (1.0/(a*volume+b*hydrophobicity+c*connectivity+constant)) * ( a*(dV_dz_polar[jj]) + b*(dH_dz_polar[jj]) + c*(dC_dz_polar[jj]) ));
//-----  derivatives without connectivity
    colvar.myder[i_c][j][0] = jedi * ((1.0/Va)*(dVa_dx_polar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dx_polar[jj]) + b*(dH_dx_polar[jj]) ));
    colvar.myder[i_c][j][1] = jedi * ((1.0/Va)*(dVa_dy_polar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dy_polar[jj]) + b*(dH_dy_polar[jj]) ));
    colvar.myder[i_c][j][2] = jedi * ((1.0/Va)*(dVa_dz_polar[jj]) + (1.0/(a*volume+b*hydrophobicity+constant)) * ( a*(dV_dz_polar[jj]) + b*(dH_dz_polar[jj]) ));
/*    printf("Polar j is colvar.myder[i_c][j][0]: %i %f\n",j,colvar.myder[i_c][j][0]);*/
/*    printf("Polar j is colvar.myder[i_c][j][1]: %i %f\n",j,colvar.myder[i_c][j][1]);*/
/*    printf("Polar j is colvar.myder[i_c][j][2]: %i %f\n",j,colvar.myder[i_c][j][2]);*/
//    colvar.myder[i_c][j][0] = jedi ;
//    colvar.myder[i_c][j][1] = jedi ;
//    colvar.myder[i_c][j][2] = jedi ;
  }

/*
//  output_derivatives = fopen("derivatives.txt", "a");
//  if (output_derivatives != NULL)
//  {

  fprintf(output_derivatives,"                  FINAL   DERIVATIVES                \n");
  fprintf(output_derivatives,"POLAR ATOMS\n");

  for(j=colvar.list[i_c][0];j<colvar.list[i_c][0]+colvar.list[i_c][1];j++) {
          jj = j - colvar.list[i_c][0];
          fprintf(output_derivatives,"%d\n",j);
          fprintf(output_derivatives,"dJEDI_dx     : %f\n",colvar.myder[i_c][j][0]);
          fprintf(output_derivatives,"dJEDI_dy     : %f\n",colvar.myder[i_c][j][1]);
          fprintf(output_derivatives,"dJEDI_dz     : %f\n",colvar.myder[i_c][j][2]);
          fprintf(output_derivatives,"dVa_dz_polar : %f\n",dVa_dz_polar[jj]);
          fprintf(output_derivatives,"dV_dz_polar  : %f\n",dV_dz_polar[jj]);
//          fprintf(output_derivatives,"dE_dz_polar  : %f\n",dE_dz_polar[jj]);
          fprintf(output_derivatives,"dH_dz_polar  : %f\n",dH_dz_polar[jj]);
//          fprintf(output_derivatives,"dC_dz_polar  : %f\n",dC_dz_polar[jj]);
          fprintf(output_derivatives,"\n");
  }

  fprintf(output_derivatives,"APOLAR ATOMS\n");
  for(j=0;j<colvar.list[i_c][0];j++) {
          fprintf(output_derivatives,"%d\n",j);
          fprintf(output_derivatives,"dJEDI_dx     : %f\n",colvar.myder[i_c][j][0]);
          fprintf(output_derivatives,"dJEDI_dy     : %f\n",colvar.myder[i_c][j][1]);
          fprintf(output_derivatives,"dJEDI_dz     : %f\n",colvar.myder[i_c][j][2]);
          fprintf(output_derivatives,"dVa_dx_apolar : %f\n",dVa_dx_apolar[j]);
          fprintf(output_derivatives,"dV_dx_apolar  : %f\n",dV_dx_apolar[j]);
//          fprintf(output_derivatives,"dE_dx_apolar  : %f\n",dE_dx_apolar[j]);
          fprintf(output_derivatives,"dH_dx_apolar  : %f\n",dH_dx_apolar[j]);
//          fprintf(output_derivatives,"dC_dx_apolar  : %f\n",dC_dx_apolar[j]);
//          fprintf(output_derivatives,"dVa_dy_apolar : %f\n",dVa_dy_apolar[j]);
//          fprintf(output_derivatives,"dV_dy_apolar  : %f\n",dV_dy_apolar[j]);
//          fprintf(output_derivatives,"dE_dy_apolar  : %f\n",dE_dy_apolar[j]);
//          fprintf(output_derivatives,"dH_dy_apolar  : %f\n",dH_dy_apolar[j]);
//          fprintf(output_derivatives,"dC_dy_apolar  : %f\n",dC_dy_apolar[j]);
//          fprintf(output_derivatives,"dVa_dz_apolar : %f\n",dVa_dy_apolar[j]);
//          fprintf(output_derivatives,"dV_dz_apolar  : %f\n",dV_dy_apolar[j]);
//          fprintf(output_derivatives,"dE_dz_apolar  : %f\n",dE_dy_apolar[j]);
//          fprintf(output_derivatives,"dH_dz_apolar  : %f\n",dH_dy_apolar[j]);
//          fprintf(output_derivatives,"dC_dz_apolar  : %f\n",dC_dy_apolar[j]);
          fprintf(output_derivatives,"\n");
  }
//  fclose(output_derivatives);
//  }


//--------- free the memory when dynamic allocation

  for (int i = 0; i < size_grid; i++) {
    free(da_dx_apolar[i]);
    free(da_dy_apolar[i]);
    free(da_dz_apolar[i]);
    free(da_dx_polar[i]);
    free(da_dy_polar[i]);
    free(da_dz_polar[i]);
  }
  free(da_dx_apolar);
  free(da_dy_apolar);
  free(da_dz_apolar);
  free(da_dx_polar);
  free(da_dy_polar);
  free(da_dz_polar);
*/
//  fclose(output_derivatives);
}

