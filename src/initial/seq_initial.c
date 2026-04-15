/*
 * seq_initial.c
 * -----------------------------------------------------------------------------
 * Build the initial coarse-grained conformation file (ch.dat) from:
 *   1) seq.dat      : one DNA sequence line
 *   2) initial.dat  : template bead coordinates
 *
 * For every nucleotide there are three beads in order P-S-Base. This program
 * relabels each third bead with the sequence base (A/T/C/G) while preserving
 * template coordinates and auxiliary columns.
 *
 * Input columns in initial.dat:
 *   bead_id residue_id bead_type x y z R Q f
 *
 * Output:
 *   ch.dat with the same columns and sequence-correct base labels.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include<ctype.h>
#define MAX_LENGTH 10000
char rna_sequence[MAX_LENGTH];
//int n1 =0, n2 = 0;
using namespace std;
int main()
{
	FILE *fp1, *fp2, *fp3;
	fp1=fopen("seq.dat","r+");
	fp2=fopen("initial.dat","r+");
	fp3=fopen("ch.dat","w+");
	 int N;
	if (fgets(rna_sequence, MAX_LENGTH, fp1) == NULL) 
    	{
        	perror("Error reading RNA sequence");
        	return EXIT_FAILURE;
    	}
    	for (int i = 0; rna_sequence[i] != '\0'; i++) 
    	{
        	rna_sequence[i] = toupper(rna_sequence[i]);
    	}
    	//printf("%s",rna_sequence);
	N=strlen(rna_sequence)-1;
	//printf("%d\n",N);
	int a;
	int n1[10000], n2[10000];
	char type[10000];
	float x[10000], y[10000], z[10000], R[10000], Q[10000], f[10000];
	
	int i=0;
	while(!feof(fp2))
	{
				
		fscanf(fp2,"%d %d %s %f %f %f %f %f %f\n",&n1[i],&n2[i],&type[i],&x[i],&y[i],&z[i],&R[i],&Q[i],&f[i]);
		i++;
	}
	for(i=0;i<=3*N;i++)
	{
			if(fmod(i+1,3)==0)
			{
				fprintf(fp3,"%d %d %c %f %f %f %f %f %f\n",n1[i],n2[i],rna_sequence[i/3],x[i],y[i],z[i],R[i],Q[i],f[i]);
			}
			else
			{
				fprintf(fp3,"%d %d %c %f %f %f %f %f %f\n",n1[i],n2[i],type[i],x[i],y[i],z[i],R[i],Q[i],f[i]);
			}
	}
	
	printf("=================================================\n");
        printf("Successfully generated initial conformation ! ! !\n");
        printf("=================================================\n");
	fclose(fp1);fclose(fp3);fclose(fp2);	
	return 0;
}
