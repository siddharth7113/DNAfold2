/*
 * A_state.c
 * -----------------------------------------------------------------------------
 * Legacy scorer utility used to extract selected conformations from conf_0.dat
 * into min_conf.dat based on line indices listed in min.dat.
 *
 * Data contract:
 *   - ch.dat      : used only to determine number of bead lines per conformation
 *   - min.dat     : selected conformation indices from Umin.py
 *   - conf_0.dat  : concatenated conformations
 *   - min_conf.dat: extracted conformations for downstream PDB conversion
 *
 * Note: this behavior is now also implemented in Python
 * (dnafold2.stage_tools.extract_min_conformations).
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define e 0.26
#define Beta 0.60          //AU=Beta*GC
#define BetaGC 0.93
#define gamma 0.1          //The energy threshold of base-pairing for statistics (碱基配对的能量Ubp0<gamma*Ubp算作配对形成)
// hydrogen bond
#define dNN1 8.55
#define dNN2 9.39
#define A -3.5             //base-pairing strength of GC
#define N_thread 15
FILE *input,*fc,*output,*inputmin;
int N,rf,bp,s[1000][1000],c[10000],NF=0,NS1=0,NS2=0,NU=0,N_other=0;
char type[10000];
float x[5000],y[5000],z[5000],Q[5000],f[5000],R[5000];
float uN,temperature[20]={25.0,31.0,37.0,43.0,49.0,55.0,63.0,71.0,78.0,86.0,94.0,102.0,110.0,120.0,130.0};
int Nstep,N0;
int main()
{
   	int id1,t1,i1,j1;
   	float HB();
   	void BasePairing(),Stem1_2();
   	fc=fopen("ch.dat","r+");
   	inputmin=fopen("min.dat","r+");
   	output=fopen("min_conf.dat","w+");
   	int i=0;
   	while(!feof(fc))
   	{
          	fscanf(fc,"%d %d %s %f %f %f %f %f %f\n",&t1,&id1,&type[i],&x[i],&y[i],&z[i],&R[i],&Q[i],&f[i]);
           	i++;
   	}
   	fclose(fc);
   	N=i;
   	int num;
   	i=0;
   	int Umin[1000];
   	while(!feof(inputmin))
   	{
   		fscanf(inputmin,"%d\n",&Umin[i]);
   		i++;
   	}
   	num=i;
	printf("----- %d\n",num);
	input=fopen("conf_0.dat","r+");

   	for(rf=0;rf<N_thread;rf++)
   	{
         	i1=0;j1=0;NF=0;NS1=0;NS2=0;NU=0;N_other=0;
         	while(!feof(input))
         	{
                	for(i1=j1*N+1;i1<=j1*N+N;i1++)
                	{
                
                       		fscanf(input,"%d %d %s %f %f %f %f %f %f\n",&t1,&id1,&type[i1-j1*N],&x[i1-j1*N],&y[i1-j1*N],&z[i1-j1*N],&R[i1-j1*N],&Q[i1-j1*N],&f[i1-j1*N]);
                       		for(i=0;i<num;i++)
                       		{
                       			if(j1==Umin[i])
                       			{
                       				fprintf(output,"%d %d %c %f %f %f %f %f %f\n",t1,id1,type[i1-j1*N],x[i1-j1*N],y[i1-j1*N],z[i1-j1*N],R[i1-j1*N],Q[i1-j1*N],f[i1-j1*N]);
                       			}
                       		}
                 	} // 读入构象；i1: No. of nt; j1: No. of conf.;           
                 	j1++;     
                 	
        	}
    	}
    	fclose(input);fclose(output);fclose(inputmin);
      	return 0;
}
