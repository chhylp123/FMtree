// This is a demo program for showing how to call SACA_K.

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

bwt_index bitmapper_index_params;






void SACA_K(unsigned char *s, unsigned int *SA, unsigned int n, 
           unsigned int K, unsigned int m, int level);



unsigned int locate(unsigned int sp, unsigned int ep, unsigned int* locates)
{


	unsigned int ijkijk = 0;
	unsigned int t;
	for (t = sp; t<ep; t++, ijkijk++)
	{
		locates[ijkijk] = get_sa(bitmapper_index_params.SA_flag, t, bitmapper_index_params.bwt, bitmapper_index_params.sa,
			bitmapper_index_params.high_occ_table);
	}

}




void indenpendent_get_sa_fromFILE(unsigned int **sa, unsigned int cc, char *refer)
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


	///这里SA要改
	#define SA_counter_length 32

	///这里要改
	///unsigned int compress_occ = 448, compress_SA_flag = 224;
	unsigned int compress_occ = 256, high_compress_occ=65536,compress_SA_flag = 224;
	typedef uint64_t bwt_string_type;

	///这里SA要改
	///typedef uint32_t SA_flag_string_type;
	typedef uint64_t SA_flag_string_type;


	typedef uint32_t high_occ_table_type;
	unsigned int bwt_warp_number = sizeof(bwt_string_type)* 8 / 2;
	unsigned int SA_flag_warp_number = sizeof(SA_flag_string_type)* 8;
	uint64_t* long_SA_flag = NULL;
	///这里要改
	///unsigned int single_occ_bwt = sizeof(bwt_string_type) / sizeof(unsigned);
	unsigned int single_occ_bwt = sizeof(bwt_string_type) / sizeof(unsigned short);
	unsigned int occ_words = 4 / single_occ_bwt;
	unsigned int acctuall_bwt_gap = compress_occ + sizeof(unsigned int)* 4 * 8 / 2;

	///这里SA要改
	///unsigned int acctuall_SA_flag_gap = compress_SA_flag + sizeof(SA_flag_string_type)* 8;
	unsigned int acctuall_SA_flag_gap = compress_SA_flag + SA_counter_length;
	///这里SA要改
	unsigned int SA_counter_shift_length = sizeof(SA_flag_string_type)* 8 - SA_counter_length;


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
	///这里要改
	FILE *f1, *f2, *fs, *fb, *fo;
	ctoi['A'] = 0;
	ctoi['C'] = 1;
	ctoi['G'] = 2;
	ctoi['T'] = 3;
	ctoi['a'] = 0;
	ctoi['c'] = 1;
	ctoi['g'] = 2;
	ctoi['t'] = 3;




	unsigned int occ_line_number = 0;
	///this is the number of occ line
	occ_line_number = (text_length + 1) / compress_occ + 1;
	///this is the number of character which could be saved in one bwt_string_type
	bwt_warp_number = sizeof(bwt_string_type)* 8 / 2;
	///this is the number of bwt_string_type representing bwt string (does not include occ)
	unsigned int bwt_length = (text_length + 1) / bwt_warp_number + 1;
	///this is the number of bwt_string_type representing occ line
	///这里要改
	/**
	unsigned int occ_byte_length = (occ_line_number * 4 * sizeof(unsigned int))
		/ (sizeof(bwt_string_type)) + 1;
	**/
	unsigned int occ_byte_length = (occ_line_number * 4 * sizeof(unsigned short))
		/ (sizeof(bwt_string_type)) + 1;
	


	bwt =
		(bwt_string_type *)malloc(sizeof(bwt_string_type)*(bwt_length + occ_byte_length + 1));





	///这里SA要改
	/**
	unsigned int SA_flag_length = (text_length + 1) / SA_flag_warp_number + 1;
	unsigned int SA_flag_occ_length = (text_length + 1) / compress_SA_flag + 1;
	unsigned int SA_occ_byte_length = (SA_flag_occ_length * 1 * sizeof(unsigned int))
		/ (sizeof(SA_flag_string_type)) + 1;
	SA_flag_string_type* SA_flag =
		(SA_flag_string_type *)malloc(sizeof(SA_flag_string_type)*(SA_flag_length + SA_occ_byte_length + 1));
	**/
	unsigned int SA_flag_length = text_length + 1;
	unsigned int SA_flag_occ_length = (text_length + 1) / compress_SA_flag + 1;
	unsigned int SA_flag_byte_length = (SA_flag_length + SA_flag_occ_length * SA_counter_length + 1)
		/ SA_flag_warp_number + 1;
	SA_flag_string_type* SA_flag = (SA_flag_string_type *)malloc
		(sizeof(SA_flag_string_type)*SA_flag_byte_length);








	///这里要改
	high_occ_table_type* high_occ_table;
	unsigned int high_occ_table_length = ((text_length + 1) / high_compress_occ + 1)*4+1;
	high_occ_table = (high_occ_table_type*)malloc(sizeof(high_occ_table_type)*high_occ_table_length);

	///这里SA要改
	///for (i = 0; i < SA_flag_length + SA_occ_byte_length + 1; i++)
	for (i = 0; i < SA_flag_byte_length; i++)
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
	fo = fopen(filenameo, "w");




	indenpendent_get_sa_fromFILE(&sa, SA_length, &refer[1]);


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
	///这里要改
	/**
	unsigned int occ_bwt_number = (sizeof(unsigned int)* 4) / sizeof(bwt_string_type);
	**/
	unsigned int occ_bwt_number = (sizeof(unsigned short)* 4) / sizeof(bwt_string_type);
	for (bwt_iterater = 0; bwt_iterater < occ_bwt_number; bwt_iterater++)
	{
		bwt[bwt_iterater] = (bwt_string_type)0;
	}

	///这里要改
	unsigned int high_occ_table_iterater = 0;
	for (high_occ_table_iterater = 0; high_occ_table_iterater < 4; high_occ_table_iterater++)
	{
		high_occ_table[high_occ_table_iterater] = 0;
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

		
		if (i % high_compress_occ == 0)
		{
			for (j = 0; j < 4; j++)
			{
				high_occ_table[high_occ_table_iterater] = nacgt[0][j + 1];
				high_occ_table_iterater++;
			}
		}


		if (i % compress_occ == 0)
		{


			///这里要改,single_bwt_occ=4
			/**
			unsigned int single_bwt_occ = sizeof(bwt_string_type) / sizeof(unsigned int);
			**/
			unsigned int single_bwt_occ = sizeof(bwt_string_type) / sizeof(unsigned short);



			for (j = 0; j < 4; j++)
			{

				tmp_bwt = (bwt_string_type)0;
				
				///这里要改
				///shift_length = (single_bwt_occ - j % single_bwt_occ - 1) * sizeof(unsigned int)* 8;
				///j=0, shift_length = 48
				///j=1, shift_length = 32
				///j=2, shift_length = 16
				///j=3, shift_length = 0
				shift_length = (single_bwt_occ - j % single_bwt_occ - 1) * sizeof(unsigned short)* 8;
				/**
				fprintf(stderr, "j=%llu\n", j);
				fprintf(stderr, "shift_length=%llu\n", shift_length);
				**/

				tmp_bwt = (bwt_string_type)nacgt[0][j + 1];

				/**
				fprintf(stderr, "tmp_bwt=%llu\n", tmp_bwt);
				**/

				///这里要改
				tmp_bwt = tmp_bwt - high_occ_table[(i / high_compress_occ) * 4 + j];

				/**
				fprintf(stderr, "(i / high_compress_occ) * 4 + j=%llu\n", (i / high_compress_occ) * 4 + j);

				fprintf(stderr, "tmp_bwt=%llu\n", tmp_bwt);
				**/



				tmp_bwt = tmp_bwt << shift_length;





				bwt[bwt_iterater + j / single_bwt_occ] =
					bwt[bwt_iterater + j / single_bwt_occ] | tmp_bwt;


			}


			///这里要改
			///bwt_iterater = bwt_iterater + single_bwt_occ;
			bwt_iterater = bwt_iterater + occ_words;
			/**
			fprintf(stderr, "single_bwt_occ=%llu\n", single_bwt_occ);
			fprintf(stderr, "occ_words=%llu\n", occ_words);
			**/




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

	///这里SA要改
	unsigned int number_of_SA_flag_bits = 0;



	i = 0;

	SA_flag[SA_flag_iterater] = (SA_flag_string_type)0;


	///这里SA要改
	///注意这里SA_flag_iterater并不能++
	///SA_flag_iterater++;
	number_of_SA_flag_bits = number_of_SA_flag_bits + SA_counter_length;

	/**
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
	**/


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



		///这里SA要改
		///shift_length = SA_flag_warp_number - i % SA_flag_warp_number - 1;
		shift_length = SA_flag_warp_number - number_of_SA_flag_bits % SA_flag_warp_number - 1;

		tmp_SA_flag = tmp_SA_flag << shift_length;



		SA_flag[SA_flag_iterater] = SA_flag[SA_flag_iterater] | tmp_SA_flag;


		tmp_SA_flag = (SA_flag_string_type)0;



		///这里SA要改
		number_of_SA_flag_bits++;

		i++;

		///这里SA要改
		/**
		if (i%SA_flag_warp_number == 0)
		{
		SA_flag_iterater++;
		}


		if (i % compress_SA_flag == 0)
		{

		SA_flag[SA_flag_iterater] = SA_number;

		SA_flag_iterater++;

		}
		**/
		if (number_of_SA_flag_bits%SA_flag_warp_number == 0)
		{
			SA_flag_iterater++;
		}
		if (i % compress_SA_flag == 0)
		{

			if (number_of_SA_flag_bits % SA_flag_warp_number != 0)
			{
				fprintf(stderr, "ERROR! \n");
			}



			SA_flag[SA_flag_iterater] = SA_number;

			SA_flag[SA_flag_iterater] = (SA_flag_string_type)(SA_flag[SA_flag_iterater] << SA_counter_shift_length);

			number_of_SA_flag_bits = number_of_SA_flag_bits + SA_counter_length;

		}


		tmp_SA_flag = (SA_flag_string_type)0;


	}


	///这里SA要改
	/**
	if (i%SA_flag_warp_number != 0 &&
	i % compress_SA_flag != 0)
	{
	SA_flag_iterater++;
	}
	**/
	if (number_of_SA_flag_bits%SA_flag_warp_number != 0)
	{
		SA_flag_iterater++;
	}

	///这里SA要改,为了防止后面计算时溢出
	SA_flag_iterater++;



	///这里要改
	fwrite(&high_occ_table_iterater, sizeof(unsigned int), 1, fo);
	fwrite(high_occ_table, sizeof(high_occ_table_type), high_occ_table_iterater, fo);
	printf("Occ has been writed!\n");



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

	unsigned int tmp_site = 0;






	while (1)
	{

		if (i >= SA_length)
		{
			break;
		}

		if (sa[i] % compress_sa == 0)
		{

			///if (compress_sa >= 4)
			if (compress_sa >= 0)
			{
				///如果前一位是$，那么在sa必然是0，这个绝壁是不符合要求
				//所以不用记录，只是要判断是不是为0，为0就丢弃
				//不过这东西应该是被保存成010000000000000000000这个数
				//tmp_site = (unsigned int)0;
				ch = refer[sa[i]] - 1;
				tmp_site = (unsigned int)abs(ch);

				tmp_site = tmp_site << 30;

				tmp_site = tmp_site | (sa[i] / compress_sa);


				fwrite(&tmp_site, sizeof(tmp_site), 1, fs);
			}
			else
			{
				fwrite(&sa[i], sizeof(sa[i]), 1, fs);
			}

			


			
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
	///这里要改
	fwrite(&high_compress_occ, sizeof(unsigned int), 1, f2);
	






	fprintf(stdout, "Sucess!\n");
	fflush(stdout);



	fflush(fo);
	fflush(f2);
	fflush(fs);
	fflush(fb);

	fclose(fo);
	fclose(f2);
	fclose(fs);
	fclose(fb);

	return 0;


}



unsigned int init_bitmapper_index_params()
{
	unsigned int i;

	for (i = 0; i < 256; i++)
	{
		bitmapper_index_params.ctoi[i] = 4;
	}
	bitmapper_index_params.ctoi['A'] = 0;
	bitmapper_index_params.ctoi['C'] = 1;
	bitmapper_index_params.ctoi['G'] = 2;
	bitmapper_index_params.ctoi['T'] = 3;
	bitmapper_index_params.ctoi['a'] = 0;
	bitmapper_index_params.ctoi['c'] = 1;
	bitmapper_index_params.ctoi['g'] = 2;
	bitmapper_index_params.ctoi['t'] = 3;
	bitmapper_index_params.cut_thr = 1;

}







unsigned int load_index(char* filename_prefix)
{
	char filename[100], filename1[100], filenames[100], filenameo[100], filenameb[100];
	long long i, j, t;

	FILE *f1, *f2, *fs, *fb, *fout, *fo;



	bitmapper_index_params.pop_count_mode[0] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.pop_count_mode[0]
			= (bwt_string_type)(bitmapper_index_params.pop_count_mode[0] << 2) | (bwt_string_type)3;
	}


	bitmapper_index_params.pop_count_mode[1] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.pop_count_mode[1] =
			(bwt_string_type)(bitmapper_index_params.pop_count_mode[1] << 2) | (bwt_string_type)2;
	}


	bitmapper_index_params.pop_count_mode[2] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.pop_count_mode[2] =
			(bwt_string_type)(bitmapper_index_params.pop_count_mode[2] << 2) | (bwt_string_type)1;
	}

	bitmapper_index_params.pop_count_mode[3] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.pop_count_mode[3] =
			(bwt_string_type)(bitmapper_index_params.pop_count_mode[3] << 2) | (bwt_string_type)0;
	}

	bitmapper_index_params.mode_high = (bwt_string_type)0;
	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.mode_high =
			(bwt_string_type)(bitmapper_index_params.mode_high << 2) | (bwt_string_type)2;
	}


	bitmapper_index_params.mode_low = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.mode_low =
			(bwt_string_type)(bitmapper_index_params.mode_low << 2) | (bwt_string_type)1;
	}




	strcpy(filename, filename_prefix);
	strcpy(filename + strlen(filename), ".index");
	strcpy(filenames, filename);
	strcpy(filenameo, filename);
	strcpy(filenameb, filename);
	strcpy(filenames + strlen(filename), ".sa");
	strcpy(filenameo + strlen(filename), ".occ");
	strcpy(filenameb + strlen(filename), ".bwt");


	///这里要改
	fo = fopen(filenameo, "r");
	if (fo == NULL)
	{
		fprintf(stdout, "Failed to open %s! FMtree will exit ...\n", filenameo);
		return 1;
	}

	fs = fopen(filenames, "r");
	if (fs == NULL)
	{
		fprintf(stdout, "Failed to open %s! FMtree will exit ...\n", filenames);
		return 1;
	}
	fb = fopen(filenameb, "r");
	if (fb == NULL)
	{
		fprintf(stdout, "Failed to open %s! FMtree will exit ...\n", filenameb);
		return 1;
	}
	strcpy(filename1, filename);
	f2 = fopen(filename1, "r");

	if (f2 == NULL)
	{
		fprintf(stdout, "Failed to open %s! FMtree will exit ...\n", filename1);
		return 1;
	}


	fread(&bitmapper_index_params.SA_length, sizeof(bitmapper_index_params.SA_length), 1, f2);
	fread(&bitmapper_index_params.shapline, sizeof(bitmapper_index_params.shapline), 1, f2);

	fprintf(stdout, "shapline=%llu\n", bitmapper_index_params.shapline);



	for (j = 0; j <= 4; j++)
		fread(&bitmapper_index_params.nacgt[j], sizeof(bitmapper_index_params.nacgt[j]), 1, f2);




	fread(&bitmapper_index_params.compress_sa, sizeof(unsigned int), 1, f2);
	fread(&bitmapper_index_params.compress_occ, sizeof(unsigned int), 1, f2);
	fread(&bitmapper_index_params.high_compress_occ, sizeof(unsigned int), 1, f2);



	bitmapper_index_params.mode_4[0] = (bwt_string_type)-1;
	bitmapper_index_params.mode_4[1] = (bwt_string_type)-1;
	bitmapper_index_params.mode_4[2] = (bwt_string_type)-1;
	bitmapper_index_params.mode_4[3] = (bwt_string_type)0;


	fread(&bitmapper_index_params.bwt_length, sizeof(bitmapper_index_params.bwt_length), 1, fb);
	bitmapper_index_params.bwt = (bwt_string_type *)malloc(sizeof(bwt_string_type)*bitmapper_index_params.bwt_length);
	fread(bitmapper_index_params.bwt, sizeof(bwt_string_type), bitmapper_index_params.bwt_length, fb);
	printf("BWT has been loaded!\n");



	printf("SA_length=%u\n", bitmapper_index_params.SA_length);
	fread(&bitmapper_index_params.sparse_suffix_array_length,
		sizeof(bitmapper_index_params.sparse_suffix_array_length), 1, fs);
	printf("sparse_suffix_array_length=%u\n", bitmapper_index_params.sparse_suffix_array_length);
	bitmapper_index_params.sa
		= (unsigned int *)malloc(sizeof(unsigned int)*(bitmapper_index_params.sparse_suffix_array_length));
	fread(bitmapper_index_params.sa, sizeof(unsigned int), bitmapper_index_params.sparse_suffix_array_length, fs);








	fread(&bitmapper_index_params.SA_flag_iterater,
		sizeof(bitmapper_index_params.SA_flag_iterater), 1, fs);

	fprintf(stdout, "SA_flag_iterater=%llu\n", bitmapper_index_params.SA_flag_iterater);

	bitmapper_index_params.SA_flag =
		(SA_flag_string_type*)malloc(sizeof(SA_flag_string_type)*bitmapper_index_params.SA_flag_iterater);

	fread(bitmapper_index_params.SA_flag, sizeof(SA_flag_string_type),
		bitmapper_index_params.SA_flag_iterater, fs);


	fread(&bitmapper_index_params.high_occ_table_length, sizeof(bitmapper_index_params.high_occ_table_length), 1, fo);
	fprintf(stdout, "high_occ_table_length=%llu\n", bitmapper_index_params.high_occ_table_length);
	bitmapper_index_params.high_occ_table = (high_occ_table_type*)malloc(sizeof(high_occ_table_type)
		*bitmapper_index_params.high_occ_table_length);
	fread(bitmapper_index_params.high_occ_table, sizeof(high_occ_table_type),
		bitmapper_index_params.high_occ_table_length, fo);





	fclose(f2);
	fclose(fs);
	fclose(fb);
	fclose(fo);

	init_bitmapper_index_params();

	return 1;



}







