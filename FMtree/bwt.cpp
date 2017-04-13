
#include<stdio.h>
#include<stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "bwt.h"
#include<stdint.h>
#include<ctype.h>
///#include <nmmintrin.h>










void SACA_K(unsigned char *s, unsigned int *SA, unsigned int n, 
           unsigned int K, unsigned int m, int level);



void indenpendent_get_sa_fromFILE(unsigned int **sa, unsigned int cc, FILE *_ih_fp, char *refer)
{
	unsigned int i;
	unsigned int n;
	n = cc-1;

	fprintf(stdout, "r_n=%u\n", n);

	n++; // append the virtual sentinel
	fprintf(stdout, "Allocating input and output space: %u bytes = %.2lf MB", 5 * n, (double)5 * n / 1024 / 1024);
	unsigned char *s_ch = new unsigned char[n];
	//unsigned char *s_ch;
	unsigned int *SA = new unsigned int[n];
	//unsigned int *SA;
	SA = (unsigned int*)malloc(sizeof(unsigned int)*n+10);
	if (s_ch == NULL || SA == NULL) {
		delete[] s_ch; delete[] SA;
		fprintf(stdout, "\nInsufficient memory, exit!");
		return;
	}


	// read the string into buffer.
	fprintf(stdout, "\nReading input string...");
	fseek(stdin, 0, SEEK_SET);
	for(i=0;i<n;i++) s_ch[i]=refer[i];
	// set the virtual sentinel
	s_ch[n - 1] = 0;








	clock_t start, finish;
	double  duration;
	start = clock();

	fprintf(stdout, "\nConstructing the suffix array...");
	SACA_K(s_ch, SA, n, 256, n, 0);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;

	fprintf(stdout, "\nSize: %u bytes, Time: %5.3f seconds\n", n - 1, duration);

	SA[0] = n - 1;
	printf("sa_n=%d\n",n);
	*sa = SA;
}


unsigned int indenpendent_creadte_index(unsigned int text_length, char** input_refer, unsigned int compress_sa, char* filename)
{


	unsigned int compress_occ = 448, compress_SA_flag = 224;
	typedef uint64_t bwt_string_type;
	typedef uint32_t SA_flag_string_type;
	unsigned int bwt_warp_number = sizeof(bwt_string_type)* 8 / 2;
	unsigned int SA_flag_warp_number = sizeof(SA_flag_string_type)* 8;
	uint64_t* long_SA_flag = NULL;
	unsigned int single_occ_bwt = sizeof(bwt_string_type) / sizeof(unsigned);
	unsigned int occ_words = 4 / single_occ_bwt;
	unsigned int acctuall_bwt_gap = compress_occ + sizeof(unsigned int)* 4 * 8 / 2;
	unsigned int acctuall_SA_flag_gap = compress_SA_flag + sizeof(SA_flag_string_type)* 8;
	bwt_string_type mode_4[4];
	bwt_string_type mode = (bwt_string_type)-1;
	bwt_string_type mode_high_1 = (bwt_string_type)1 << (SA_flag_warp_number - 1);
	bwt_string_type mode_16 = (bwt_string_type)65535;
	bwt_string_type mode_32 = ((bwt_string_type)-1) >> 32;
	bwt_string_type mode_high;
	bwt_string_type mode_low;
	SA_flag_string_type mode_SA_flag = (SA_flag_string_type)(((SA_flag_string_type)-1) << (sizeof(SA_flag_string_type)* 8 - 1));

	bwt_string_type pop_count_mode[4];
	unsigned int bwt_count_hash_table_bit = 16;
	unsigned int SA_length;
	unsigned int na, nc, ng, nt;
	unsigned int nacgt[4][258];
	unsigned int bwt_step = 1;
	unsigned int ctoi[256];

	unsigned int shapline;
	char filenames[100], filenameo[100], filenameb[100];
	char *refer = (*input_refer);
	bwt_string_type *bwt;
	unsigned int **occ;
	unsigned int *sa;
	char ch;
	unsigned int i, j;
	FILE *f1, *f2, *fs, *fb;
	ctoi['A'] = 0;
	ctoi['C'] = 1;
	ctoi['G'] = 2;
	ctoi['T'] = 3;
	ctoi['a'] = 0;
	ctoi['c'] = 1;
	ctoi['g'] = 2;
	ctoi['t'] = 3;




	unsigned int occ_line_number = 0;
	occ_line_number = (text_length + 1) / compress_occ + 1;
	bwt_warp_number = sizeof(bwt_string_type)* 8 / 2;
	unsigned int bwt_length = (text_length + 1) / bwt_warp_number + 1;
	unsigned int occ_byte_length = (occ_line_number * 4 * sizeof(unsigned int))
		/ (sizeof(bwt_string_type)) + 1;
	bwt =
		(bwt_string_type *)malloc(sizeof(bwt_string_type)*(bwt_length + occ_byte_length + 1));

	unsigned int SA_flag_length = (text_length + 1) / SA_flag_warp_number + 1;
	unsigned int SA_flag_occ_length = (text_length + 1) / compress_SA_flag + 1;
	unsigned int SA_occ_byte_length = (SA_flag_occ_length * 1 * sizeof(unsigned int))
		/ (sizeof(SA_flag_string_type)) + 1;



	SA_flag_string_type* SA_flag =
		(SA_flag_string_type *)malloc(sizeof(SA_flag_string_type)*(SA_flag_length + SA_occ_byte_length + 1));



	for (i = 0; i < SA_flag_length + SA_occ_byte_length + 1; i++)
	{
		SA_flag[i] = (SA_flag_string_type)0;
	}




	for (i = 0; i < bwt_length + occ_byte_length + 1; i++)
	{
		bwt[i] = (bwt_string_type)0;
	}


	fprintf(stdout, "bwt_length=%llu,occ_byte_length=%llu, bwt_length+occ_byte_length=%llu\n",
		bwt_length, occ_byte_length, bwt_length + occ_byte_length);





	SA_length = text_length + 1;
	for (i = text_length; i >0; i--)
	{
		ch = refer[i - 1];
		refer[i] = ctoi[ch] + 1;
	}
	refer[0] = 0;



	strcpy(filename + strlen(filename), ".index");
	strcpy(filenames, filename);
	strcpy(filenameo, filename);
	strcpy(filenameb, filename);
	strcpy(filenames + strlen(filename), ".sa");
	strcpy(filenameo + strlen(filename), ".occ");
	strcpy(filenameb + strlen(filename), ".bwt");
	f2 = fopen(filename, "w");
	fs = fopen(filenames, "w");
	fb = fopen(filenameb, "w");




	indenpendent_get_sa_fromFILE(&sa, SA_length, f1, &refer[1]);


	printf("SA has been generated!\n");
	na = nc = ng = nt = 0;
	for (i = 0; i<bwt_step; i++)
	for (j = 0; j <= (1 << (i + i + 2)); j++)
		nacgt[i][j] = 0;


	unsigned int bwt_iterater = 0;
	i = 0;
	unsigned int shift_length = 0;
	bwt_string_type tmp_bwt = (bwt_string_type)0;


	unsigned int tmp_occ = 0;
	unsigned int occ_bwt_number = (sizeof(unsigned int)* 4) / sizeof(bwt_string_type);
	for (bwt_iterater = 0; bwt_iterater < occ_bwt_number; bwt_iterater++)
	{
		bwt[bwt_iterater] = (bwt_string_type)0;
	}



	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 258; j++)
		{
			nacgt[i][j] = 0;
		}
	}


	long long A_number[4] = { 0 };


	int flag = 0;

	i = 0;




	while (1)
	{

		tmp_bwt = (bwt_string_type)0;

		if (i >= SA_length)
		{
			break;
		}

		ch = refer[sa[i]] - 1;
		shift_length = (bwt_warp_number - i % bwt_warp_number - 1) * 2;
		tmp_bwt = (bwt_string_type)abs(ch);
		tmp_bwt = tmp_bwt << shift_length;





		bwt[bwt_iterater] = bwt[bwt_iterater] | tmp_bwt;





		tmp_bwt = (bwt_string_type)0;



		nacgt[0][ch + 1]++;



		A_number[abs(ch)]++;






		if (ch == -1 && i<SA_length)
		{
			shapline = i;
		}





		i++;



		if (i%bwt_warp_number == 0)
		{
			bwt_iterater++;
		}


		if (i % compress_occ == 0)
		{



			unsigned int single_bwt_occ = sizeof(bwt_string_type) / sizeof(unsigned int);




			for (j = 0; j < 4; j++)
			{

				tmp_bwt = (bwt_string_type)0;

				shift_length = (single_bwt_occ - j % single_bwt_occ - 1) * sizeof(unsigned int)* 8;



				tmp_bwt = (bwt_string_type)nacgt[0][j + 1];



				tmp_bwt = tmp_bwt << shift_length;





				bwt[bwt_iterater + j / single_bwt_occ] =
					bwt[bwt_iterater + j / single_bwt_occ] | tmp_bwt;


			}

			bwt_iterater = bwt_iterater + single_bwt_occ;




		}


		tmp_bwt = (bwt_string_type)0;


	}


	


	if (i%bwt_warp_number != 0 &&
		i % compress_occ != 0)
	{
		bwt_iterater++;
	}


	printf("Occ and BWT have been built!\n");

	

	fprintf(stdout, "write SA_length=%u, shapline=%u\n", SA_length, shapline);





	


	SA_flag_string_type tmp_SA_flag = (SA_flag_string_type)0;

	unsigned int SA_flag_iterater = 0;

	unsigned int SA_number = 0;

	i = 0;

	SA_flag[SA_flag_iterater] = (SA_flag_string_type)0;

	SA_flag_iterater++;

	while (1)
	{

		tmp_SA_flag = (SA_flag_string_type)0;

		if (i >= SA_length)
		{
			break;
		}



		if (sa[i] % compress_sa == 0)
		{
			tmp_SA_flag = (SA_flag_string_type)1;
			SA_number++;
		}
		else
		{
			tmp_SA_flag = (SA_flag_string_type)0;
		}


		shift_length = SA_flag_warp_number - i % SA_flag_warp_number - 1;

		tmp_SA_flag = tmp_SA_flag << shift_length;



		SA_flag[SA_flag_iterater] = SA_flag[SA_flag_iterater] | tmp_SA_flag;


		tmp_SA_flag = (SA_flag_string_type)0;



		i++;



		if (i%SA_flag_warp_number == 0)
		{
			SA_flag_iterater++;
		}


		if (i % compress_SA_flag == 0)
		{

			SA_flag[SA_flag_iterater] = SA_number;

			SA_flag_iterater++;

		}


		tmp_SA_flag = (SA_flag_string_type)0;


	}



	if (i%SA_flag_warp_number != 0 &&
		i % compress_SA_flag != 0)
	{
		SA_flag_iterater++;
	}












	fwrite(&bwt_iterater, sizeof(unsigned int), 1, fb);
	fwrite(bwt, sizeof(bwt_string_type), bwt_iterater, fb);
	printf("BWT has been writed!\n");

	fwrite(&SA_length, sizeof(SA_length), 1, f2);
	fwrite(&shapline, sizeof(shapline), 1, f2);

	printf("cc&sharp_line has been writed!\n");


	for (i = 0; i<bwt_step; i++)
	{

		nacgt[i][0] = 1;


		///fprintf(f2, "1");
		fwrite(&nacgt[i][0], sizeof(unsigned int), 1, f2);





		for (j = 1; j <= (1 << (i + i + 2)); j++)
		{
			nacgt[i][j] = nacgt[i][j] + nacgt[i][j - 1];

			fwrite(&nacgt[i][j], sizeof(unsigned int), 1, f2);
			fprintf(stdout, "nacgt[%u][%u]=%u\n", i, j, nacgt[i][j]);
			fflush(stdout);
		}
		///fprintf(f2, "\n");
	}


	fwrite(&SA_number, sizeof(SA_number), 1, fs);





	unsigned int ijkijkijkijk = 0;

	i = 0;

	while (1)
	{

		if (i >= SA_length)
		{
			break;
		}

		if (sa[i] % compress_sa == 0)
		{
			fwrite(&sa[i], sizeof(sa[i]), 1, fs);
			ijkijkijkijk++;
		}

		i++;

	}

	fprintf(stdout, "SA_number=%llu\n", SA_number);
	fprintf(stdout, "i=%llu, ijkijkijkijk=%llu\n", i, ijkijkijkijk);




	fprintf(stdout, "SA_flag_iterater=%llu\n", SA_flag_iterater);


	fwrite(&SA_flag_iterater, sizeof(SA_flag_iterater), 1, fs);
	fwrite(SA_flag, sizeof(SA_flag_string_type), SA_flag_iterater, fs);





	fwrite(&compress_sa, sizeof(unsigned int), 1, f2);
	fwrite(&compress_occ, sizeof(unsigned int), 1, f2);





	fclose(f2);

	fprintf(stdout, "Sucess!\n");
	fflush(stdout);

	fclose(fs);
	fclose(fb);


	return 0;


}



