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
#include<stdint.h>
#include<ctype.h>
#include "bwt.h"
///#include <nmmintrin.h>
///这里要改
///unsigned int compress_sa = 8, compress_occ = 448, compress_SA_flag=224;
unsigned int compress_sa = 8, compress_occ = 256, high_compress_occ = 65536, compress_SA_flag = 224;

///这里SA要改
#define SA_counter_length 32

typedef uint64_t bwt_string_type;
///这里SA要改
///typedef uint32_t SA_flag_string_type;
typedef uint64_t SA_flag_string_type;

typedef uint32_t high_occ_table_type;
///typedef unsigned int bwt_string_type;
unsigned int bwt_warp_number = sizeof(bwt_string_type)* 8 / 2;
unsigned int SA_flag_warp_number = sizeof(SA_flag_string_type)* 8;
uint64_t tmp_SA_flag = (uint64_t)0;
uint64_t* long_SA_flag = NULL;
///这里要改
///unsigned int single_occ_bwt = sizeof(bwt_string_type)/sizeof(unsigned);
unsigned int single_occ_bwt = sizeof(bwt_string_type) / sizeof(unsigned short);
unsigned int occ_words = 4 / single_occ_bwt;
///这里要改
///unsigned int acctuall_bwt_gap = compress_occ + sizeof(unsigned int)* 4 * 8 / 2;
unsigned int acctuall_bwt_gap = compress_occ + sizeof(unsigned short)* 4 * 8 / 2;

///这里SA要改
///unsigned int acctuall_SA_flag_gap = compress_SA_flag + sizeof(SA_flag_string_type)*8;
unsigned int acctuall_SA_flag_gap = compress_SA_flag + SA_counter_length;

///这里SA要改
unsigned int SA_counter_shift_length = sizeof(SA_flag_string_type)* 8 - SA_counter_length;

bwt_string_type mode_4[4];
bwt_string_type mode = (bwt_string_type)-1;
bwt_string_type mode_high_1 = (bwt_string_type)1 << (SA_flag_warp_number-1);
bwt_string_type mode_16 = (bwt_string_type)65535;
bwt_string_type mode_32 = ((bwt_string_type)-1)>>32;
bwt_string_type mode_high;
bwt_string_type mode_low;
SA_flag_string_type mode_SA_flag = (SA_flag_string_type)(((SA_flag_string_type)-1) << (sizeof(SA_flag_string_type)* 8 - 1));
///这里SA要改
SA_flag_string_type SA_pop_count_mode = (SA_flag_string_type)(((SA_flag_string_type)-1) >> SA_counter_length);




bwt_string_type pop_count_mode[4];
///unsigned int cc;
unsigned int bwt_count_hash_table_bit = 16;
unsigned int text_length;
unsigned int SA_length;
unsigned int na, nc, ng, nt;
unsigned int nacgt[4][258];
unsigned int bwt_step = 1;
unsigned int ctoi[256];
char itoc[4] = { 'A', 'C', 'G', 'T' };

unsigned int total_gap;
unsigned int num_r = 10000, len_r = 12;
unsigned int SA_header_mode = (unsigned int)(((unsigned int)-1) >> 2);
unsigned int SA_header_mode_reverse = (unsigned int)(((unsigned int)-1) << 30);
FILE* _rg_fp;

unsigned int cnt_table[65536][4];
u_int64_t cnt_aux[4]={0xffffffffffffffffull,0xaaaaaaaaaaaaaaaaull,0x5555555555555555,0x0ull};

unsigned int shapline;


struct timeval occ_tv_begin, occ_tv_end;
int64_t timecc = 0;

const int Max = 200001;
//int num[Max];
//int r[Max * 3], sa[Max * 3];

unsigned int *wa, *wb, *wv, *wd;
#define F(x) ((x) / 3 + ((x) % 3 == 1 ? 0 : tb))
#define G(x) ((x) < tb ? (x) * 3 + 1 : ((x) - tb) * 3 + 2)
#define SEQ_MAX_LENGTH		1000			// Seq Max Length


typedef struct _rg_name_l
{
	char _rg_chrome_name[1000];
	unsigned int _rg_chrome_length;
	unsigned int start_location;
	unsigned int end_location;
	//struct _rg_name_l *next;
}_rg_name_l;



#ifndef uchar
#define uchar unsigned char
#endif
#ifndef ulong
#define ulong unsigned int
#endif


void SACA_K(unsigned char *s, unsigned int *SA, unsigned int n, 
           unsigned int K, unsigned int m, int level);


double Get_T(void)
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec + t.tv_usec / 1000000.0;
}


void dec_bit(bwt_string_type input)
{
	unsigned int bit_2[64];
	int i = 0;
	for (i = 0; i < 64; i++)
	{

		///fprintf(stdout, "input=%llu\n", input);

		bwt_string_type tmp = (bwt_string_type)((bwt_string_type)input&(bwt_string_type)1);

		input = input >> 1;


		///fprintf(stdout, "tmp=%llu\n", tmp);

		if (tmp == (bwt_string_type)0)
		{
			bit_2[i] = 0;
		}
		else
		{
			bit_2[i] = 1;
		}
	}

	for (i = 63; i >= 0; i--)
	{

		///fprintf(stdout, "%u\n", i);


		fprintf(stderr, "%u", bit_2[i]);
	}

	fprintf(stderr, "\n");

}



// uncomment the below line to verify the result SA
//#define _verify_sa

// output values:
//  1: s1<s2
//  0: s1=s2
// -1: s1>s2
int sless(unsigned char *s1, unsigned char *s2, unsigned int n) {
	for(unsigned int i=0; i<n; i++) {
		if (s1[i] < s2[i]) return 1;
		if (s1[i] > s2[i]) return -1;
	}
	return 0;
} 

// test if SA is sorted for the input string s
bool isSorted(unsigned int *SA, unsigned char *s, unsigned int n) {
	for(unsigned int i = 0;  i < n-1;  i++) {

		if (i%1000==0)
		{
			fprintf(stdout, "i=%u\n", i);

		}

	  unsigned int d=SA[i]<SA[i+1]?n-SA[i+1]:n-SA[i];
		int rel=sless(s+SA[i], s+SA[i+1], d);
		if(rel==-1 || rel==0 && SA[i+1]>SA[i])
			return false;
	}
	return true;  
}


void output_bit_data(uint64_t input_value)
{
	int bit_values[64];

	uint64_t mode_1 = (uint64_t)1;
	
	for (int i = 0; i < 64; i++)
	{
		bit_values[i] = input_value&mode_1;
		input_value = input_value >> 1;
	}

	for (int i = 63; i >= 0; i--)
	{
		fprintf(stdout, "%u", bit_values[i]);
	}

	fprintf(stdout, "\n");

}


int Vertified_Hash_Index(uchar* ref, ulong ref_length, uchar *pattern, ulong length)
{

	unsigned int *locs = NULL;

	/**

	fprintf(stdout, "reference length = %u!\n", ref_length);
	fprintf(stdout, "seed = !%s!\n", pattern);

	fprintf(stdout, "length of seed = %d\n", length);

	**/	

	unsigned int i = 0;
	int j = 0;
	int avali = 0;
	unsigned int list_length = 0;
	int flag = 0;
	for (i = 0; i <= ref_length - length; ++i)
	{

		for (j = 0; j<length; j++)
		{
			if (ref[i + j] != pattern[j])
			{
				flag = 0;
				break;
			}
			flag = 1;
			//list_length++;
		}
		if (flag)
		{
			list_length++;
			flag = 0;
		}


	}
	///fprintf(stdout, "i=%u!\n", i);
	fprintf(stdout, "times(scan)=%u!\n", list_length);






}

void convert_to_pizzachili_format()
{

	FILE* _ih_fp = fopen("patterns.txt", "r");

	unsigned int numocc;
	unsigned int query_length = 0;

	fread(&(query_length), sizeof(unsigned int), 1, _ih_fp);
	fread(&(numocc), sizeof(unsigned int), 1, _ih_fp);

	fprintf(stdout, "query_length=%u\n", query_length);
	fprintf(stdout, "numocc=%u\n", numocc);


	char* patterns = (char *)malloc((query_length*numocc)*(sizeof(char)));

	fread(patterns, sizeof(char), (query_length*numocc), _ih_fp);

	FILE* _out_fp = fopen("pizzachili_format_patterns.txt", "w");

	char haha[3];

	haha[0] = 'S';
	haha[1] = 'B';
	haha[2] = '\0';

	fprintf(_out_fp, "# number=%i length=%i file=%s forbidden=%s\n",
		numocc, query_length, haha, "");

	fwrite(patterns, sizeof(char), (query_length*numocc), _out_fp);




		/**
	if (fprintf(ofile, "# number=%i length=%i file=%s forbidden=%s\n",
	J, m, argv[1],
	forbid == NULL ? "" : (char *)forbid) <= 0)
	**/
}





void Generate_patterns() {

	char filename[100];



	fprintf(stdout, "Please input the text file name (patterns will be randomly generated from file):\n");
	scanf("%s", filename);


	fprintf(stdout, "\n\n\n*********************Warning*********************\n\n");



















	fprintf(stdout,
		"We strongly recommend that input text (genome) should only consist of the characters in {a, c, g, t, A, C, G, T}!\n");
	fprintf(stdout,
		"If not, the generated patterns may include the characters which do not belong to {a, c, g, t, A, C, G, T}. \n");


	fprintf(stdout,
		"You could preprocess input text (genome) using the program named 'preprocess'. \n");

	fprintf(stdout, "\n*********************Warning*********************\n\n\n");







	FILE *_ih_fp = fopen(filename, "r");

	fseek(_ih_fp, 0, SEEK_END);
	unsigned int n = ftell(_ih_fp);

	///fprintf(stdout, "r_n=%u\n", n);

	fseek(_ih_fp, 0, SEEK_SET);


	unsigned char* text = (unsigned char *)malloc((n)*(sizeof(unsigned char)));

	fread(text, sizeof(unsigned char), n, _ih_fp);


	fclose(_ih_fp);

	_ih_fp = fopen("patterns.txt", "w");

	unsigned int query_number=0;
	unsigned int query_length = 0;


	fprintf(stdout, "Please input the required pattern length:\n");
	scanf("%u", &query_length);

	fwrite(&query_length, sizeof(unsigned int), 1, _ih_fp);

	fprintf(stdout, "Please input the required pattern nunmber\n");
	scanf("%u", &query_number);

	fwrite(&query_number, sizeof(unsigned int), 1, _ih_fp);
















	
	unsigned int occur = 0;

	long long tmp_loc = 0;




	occur = 0;





	srandom((unsigned int)time(NULL));

	while (occur<query_number)
	{
		


		tmp_loc=random();

		if (tmp_loc + query_length<=n)
		{
			occur++;

			fwrite(text + tmp_loc, sizeof(unsigned char), query_length, _ih_fp);

			///fprintf(_ih_fp, "\n");

		}
		
	}





}










int c0(unsigned int *r, unsigned int a, unsigned int b)
{
	return r[a] == r[b] && r[a + 1] == r[b + 1] && r[a + 2] == r[b + 2];
}
int c12(int k, unsigned int *r, unsigned int a, unsigned int b)
{
	if (k == 2) return r[a] < r[b] || r[a] == r[b] && c12(1, r, a + 1, b + 1);
	else return r[a]<r[b] || r[a] == r[b] && wv[a + 1]<wv[b + 1];
}
void count_sort(unsigned int *r, unsigned int *a, unsigned int *b, unsigned int n, unsigned int m)
{
	int i;
	for (i = 0; i<n; i++) wv[i] = r[a[i]];
	for (i = 0; i<m; i++) wd[i] = 0;
	//printf("%d %d\n",n,m);
	for (i = 0; i<n; i++)
	{
		//if(m==262148)
		//printf("%d %d %d\n",a[0],i,wv[i]);

		wd[wv[i]]++;
	}
	//printf("%d %d\n",n,m);
	for (i = 1; i<m; i++) wd[i] += wd[i - 1];
	//printf("%d %d\n",n,m);
	for (i = n - 1; i >= 0; i--) b[--wd[wv[i]]] = a[i];

	return;
}



void qquicksort(int l, int r, unsigned int *sa, char *refer)
{
	unsigned int x, y, i, j, t;
	/*
	printf("%s\n",refer);
	for(t=0;t<cc;t++)
	printf("%d ",sa[t]);
	printf("\n");
	*/
	i = l;
	j = r;
	x = sa[(l + r) / 2];
	do{
		//printf("l:%d,r:%d,i:%d,j:%d\n",l,r,i,j);
		while (strcmp(&refer[sa[i]], &refer[x])<0) i++;
		while (strcmp(&refer[sa[j]], &refer[x])>0) j--;
		if (i + 10 <= j + 10)
		{
			//printf("%d %d\n",i,j);
			y = sa[i];
			sa[i] = sa[j];
			sa[j] = y;
			//printf("%d %d\n",i,j);
			i++;
			j--;
		}
		//printf("l:%d,r:%d,i:%d,j:%d\n",l,r,i,j);
	} while (i + 10 <= j + 10);

	if (l + 10<j + 10) qquicksort(l, j, sa, refer);
	if (i + 10<r + 10) qquicksort(i, r, sa, refer);
}

void put_sa_toFILE(unsigned int *sa, unsigned int cc)
{
	int64_t i;
	FILE *fp;
	fp = fopen("woyaosa", "w");
	for (i = 0; i<cc; i++)
	{
		fprintf(fp, "%d\n", sa[i]);
	}

	fclose(fp);
}

void get_sa_fromFILE(unsigned int **sa, unsigned int cc, FILE *_ih_fp, char *refer)
{
	unsigned int i;
	/**
	fseek(_ih_fp, 0, SEEK_SET);

	fseek(_ih_fp, 0, SEEK_END);
	
	unsigned int n = ftell(_ih_fp);
	**/
	unsigned int n;
	n = cc-1;

	///fread(&(suffix_array.n), sizeof(unsigned int), 1, _ih_fp);

	fprintf(stdout, "r_n=%u\n", n);

	///fseek(_ih_fp, 0, SEEK_SET);



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
	//fread((unsigned char *)s_ch, 1, n - 1, _ih_fp);
	//s_ch=(unsigned char *)refer;
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
	//for(i=0;i<10;i++)
	//	printf("SA[%d]=%d\n",i,SA[i]);

	*sa = SA;

	/**
	fprintf(stdout, "write bwa_sa\n", n);

	fwrite(&n, sizeof(n), 1, _ih_fp);
	fwrite(SA, sizeof(unsigned int), n, _ih_fp);


	fclose(_ih_fp);

	delete[] SA;
	**/

	
}



void loadRefGenome(char **refGen, _rg_name_l **refGenName, int* refChromeCont, unsigned int *refGenOff)
{
	char ch;
	unsigned int _rg_contGen = 0;
	unsigned int _rg_contChromeLength = 0;
	int _rg_contChrome = 0;

	(*refGen) = (char*)malloc(sizeof(char)*(*refGenOff));
	*refGenName = (_rg_name_l*)malloc(sizeof(_rg_name_l)* 1);


	while (fscanf(_rg_fp, "%c", &ch) > 0)
	{


		if (ch == '>')
		{
			fprintf(stdout, "%c\n", ch);
			fflush(stdout);
			_rg_name_l*  tmp_rg_name = (struct _rg_name_l*)malloc(sizeof(struct _rg_name_l)* (++_rg_contChrome));

			int i = 0;

			for (i = 0; i < (_rg_contChrome - 1); i++)
			{
				strcpy(tmp_rg_name[i]._rg_chrome_name, (*refGenName)[i]._rg_chrome_name);
				tmp_rg_name[i]._rg_chrome_length = (*refGenName)[i]._rg_chrome_length;
			}

			char *tmp;
			tmp = fgets(tmp_rg_name[_rg_contChrome - 1]._rg_chrome_name, SEQ_MAX_LENGTH, _rg_fp);
			if ((_rg_contChrome - 2) >= 0)
			{
				tmp_rg_name[_rg_contChrome - 2]._rg_chrome_length = _rg_contChromeLength;
			}
			if (tmp == NULL)
				fprintf(stdout, "Error reading the reference.\n");
			free((*refGenName));
			(*refGenName) = tmp_rg_name;
			int k;
			for (k = 0; k<strlen((*refGenName)[_rg_contChrome - 1]._rg_chrome_name); k++)
			{
				if ((*refGenName)[_rg_contChrome - 1]._rg_chrome_name[k] == ' ' || (*refGenName)[_rg_contChrome - 1]._rg_chrome_name[k] == '\n')
				{
					(*refGenName)[_rg_contChrome - 1]._rg_chrome_name[k] = '\0';
					break;
				}
			}
			fprintf(stdout, "Chrome Name =%s ********\n", (*refGenName)[_rg_contChrome - 1]._rg_chrome_name);
			_rg_contChromeLength = 0;
		}
		else if (!isspace(ch) && (ch != '\n'))
		{

			ch = toupper(ch);
			(*refGen)[_rg_contGen++] = ch;
			_rg_contChromeLength++;
		}
	}




	(*refGenName)[_rg_contChrome - 1]._rg_chrome_length = _rg_contChromeLength;

	(*refGen)[_rg_contGen] = '\0';
	*refChromeCont = _rg_contChrome;
	*refGenOff = _rg_contGen;
	fclose(_rg_fp);
}



unsigned int initLoadingRefGenome(char *fileName)
{
	_rg_fp = fopen(fileName, "r");

	if (_rg_fp == NULL)
	{
		return 0;
	}
	unsigned int file_length = 0;
	fseek(_rg_fp, 0, SEEK_END);
	file_length = ftell(_rg_fp);
	fseek(_rg_fp, 0, SEEK_SET);
	return file_length;
}



unsigned int creadte_index(unsigned int text_length, char** input_refer, unsigned int compress_sa, char* filename)
{
	///char filename[100], filenames[100], filenameo[100], filenameb[100];
	char filenames[100], filenameo[100], filenameb[100];
	char *refer = (*input_refer);
	bwt_string_type *bwt;
	unsigned int **occ;
	unsigned int *sa;
	char ch;
	unsigned int i, j;
	FILE *f1, *f2, *fs, *fb;
	unsigned int ctoi[256];
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
	///fo = fopen(filenameo, "w");
	fb = fopen(filenameb, "w");




	get_sa_fromFILE(&sa, SA_length, f1, &refer[1]);


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

		///ch = tmp_refer[sa[i]] - 1;
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


		fwrite(&nacgt[i][0], sizeof(unsigned int), 1, f2);





		for (j = 1; j <= (1 << (i + i + 2)); j++)
		{
			nacgt[i][j] = nacgt[i][j] + nacgt[i][j - 1];

			fwrite(&nacgt[i][j], sizeof(unsigned int), 1, f2);
			///fprintf(f2, " %d", nacgt[i][j]);
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

	fprintf(stdout, "**********\n");
	fflush(stdout);

	fclose(fs);
	fclose(fb);


	return 0;


}


int build_index()
{
	_rg_name_l *refChromeName = NULL;
	int refChromeCont = 0;
	char filename[100], filenames[100], filenameo[100], filenameb[100];
	char *refer;
	bwt_string_type *bwt;
	unsigned int **occ;
	unsigned int *sa;
	char ch;
	unsigned int i, j;
	int64_t occ_roll[4];

	///FILE *f1, *f2, *fs, *fo, *fb;
	FILE *f1, *f2, *fs, *fb;




	fprintf(stdout, "\n\n\n*********************Warning*********************\n\n");


	fprintf(stdout,
		"The input text must only consist of the characters in {a, c, g, t, A, C, G, T}!\n");
	fprintf(stdout,
		"If not, you could preprocess it using the program named 'preprocess'. \n");

	fprintf(stdout, "\n*********************Warning*********************\n\n\n");



	fprintf(stdout,
		"Input file name:\n");
	scanf("%s", filename);

	
	fprintf(stdout, "Please input the sampling distance D of SA (1 < D < 9):\n");
	scanf("%u", &compress_sa);

	if (compress_sa>=9 || compress_sa<=1)
	{
		fprintf(stdout, "Do not input allowed sampling distance D! FMtree will exit ...\n");
		return 1;
	}

	
	ctoi['A'] = 0;
	ctoi['C'] = 1;
	ctoi['G'] = 2;
	ctoi['T'] = 3;
	ctoi['a'] = 0;
	ctoi['c'] = 1;
	ctoi['g'] = 2;
	ctoi['t'] = 3;
	text_length = 0;
	f1 = fopen(filename, "r");


	if (NULL == f1)
	{
		fprintf(stdout, "Failed to open %s! FMtree will exit ...\n", filename);
		return 1;
	}


	unsigned int number_of_total_characters = 0;

	while (!feof(f1))
	{
		fscanf(f1, "%c", &ch);
		if ((ch == 'A') || (ch == 'C') || (ch == 'G') || (ch == 'T')||
			(ch == 'a') || (ch == 'c') || (ch == 'g') || (ch == 't'))
			text_length++;
		else
		{
			
			fprintf(stdout, "%u-th character does not belong to {a, c, g, t, A, C, G, T}. Its  ASCII Code is %u, and it is %c.\n",
				number_of_total_characters, ch, ch);
		}

		number_of_total_characters++;

	}

	printf("text_length=%llu\n", text_length);

	fclose(f1);


	refer = (char *)malloc(sizeof(char)*(text_length + 1));




	f1 = fopen(filename, "r");
	for (i = 0; i <text_length; i++)
	{
		fscanf(f1, "%c", &ch);

		refer[i] = ch;
	}





















	/**

	unsigned int ctoi[256];

	ctoi['A'] = 0;
	ctoi['C'] = 1;
	ctoi['G'] = 2;
	ctoi['T'] = 3;
	ctoi['a'] = 0;
	ctoi['c'] = 1;
	ctoi['g'] = 2;
	ctoi['t'] = 3;

	i = 0;

	unsigned int four_bit_i = 0;

	unsigned char tmp_char = 0;

	while (i <text_length)
	{
		tmp_char = (unsigned char)0;

		int inner_i_output = 3;



		for (inner_i_output = 3; inner_i_output >=0; inner_i_output--)
		{
			tmp_char = tmp_char|
				(unsigned char)((unsigned char)ctoi[refer[i]] << (inner_i_output * 2));
			
			///fprintf(stderr, "i=%u,refer[i]=%c\n", i, refer[i]);

			///dec_bit(tmp_char);
			

			i++;

			if (i >= text_length)
			{
				break;

			}

		}

		refer[i / 4-1] = tmp_char;

		
		///fprintf(stderr, "refer[i / 4]= \n");
		///dec_bit(refer[i / 4]);
		

		four_bit_i = i / 4;
	}


	
	for (i = 0; i < four_bit_i; i++)
	{
		fprintf(stderr, "%c", refer[i]);
		///tmp_char = refer[i];
		///fprintf(stderr, "i=%u\n",i);
		///dec_bit(tmp_char);
	}
	

	return 0;

	**/







	indenpendent_creadte_index(text_length, &refer, compress_sa, filename);


	return 0;

}

static inline int __occ_auxx(u_int64_t y, int c)
{
	// reduce nucleotide counting to bits counting
	y = ((c & 2) ? y : ~y) >> 1 & ((c & 1) ? y : ~y) & 0x5555555555555555ull;
	// count the number of 1s in y
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

static inline int __occ_aux(u_int64_t y, int c)
{
	// reduce nucleotide counting to bits counting
	y = y ^ cnt_aux[c];
	y = (y >> 1) & y & 0x5555555555555555ull;
	// count the number of 1s in y
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

unsigned int jhp_find_occ(unsigned int line, unsigned int *occ, int delta, bwt_string_type *bwt)
{

	unsigned int i;
	unsigned int b;

	unsigned int ans = nacgt[0][delta];


    if(line%compress_occ==0) 
		return ans + occ[(line / compress_occ) * 4 + delta];


	ans += occ[(line / compress_occ) * 4 + delta];



	//printf("line%d %d\n",ans,delta);
	for (i = (line / compress_occ)*compress_occ; i + 16<line; i = i + 16)
	{
		///if(((bwt[i/4]>>(((~i)&3)<<1))&3)==delta) ans++;
		b = bwt[i / 16];
		//printf("%d %d %d %d\n",b&0xffff,b>>16&0xffff,cnt_table[b>>16&0xffff][delta],cnt_table[b&0xffff][delta]);
		//printf("b%d\n",b);
		ans += (cnt_table[b & 0xffff][delta] + cnt_table[b >> 16 & 0xffff][delta]) & 0xffff;
		//printf("line%d\n",i);
	}



	//printf("line%d %d\n",ans,delta);
	b = bwt[i / 16] & ~((1ull << ((~(line - 1) & 15) << 1)) - 1);
	//printf("%d %d %d %d\n",b&0xffff,b>>16&0xffff,cnt_table[b>>16&0xffff][delta],cnt_table[b&0xffff][delta]);
	//printf("%d %d %d\n",b,b&0xffff,cnt_table[0][0]);
	ans += (cnt_table[b & 0xffff][delta] + cnt_table[b >> 16 & 0xffff][delta]) & 0xffff;
	//printf("line%d %d\n",ans,delta);
	if (delta == 0) ans -= ~(line - 1) & 15;
	//printf("line%d %d\n",ans,delta);
	if ((delta == 1) && (shapline >= (line / compress_occ)*compress_occ) && (shapline<line)) ans--;


	//printf("return%d\n",ans);
	//gettimeofday(&occ_tv_end, NULL);
	//timecc+=((occ_tv_end.tv_sec-occ_tv_begin.tv_sec)*1000000+occ_tv_end.tv_usec-occ_tv_begin.tv_usec);
	return ans;
}


unsigned int find_occ(unsigned int line, int delta, bwt_string_type *bwt, high_occ_table_type* high_occ_table)
{

	///fprintf(stdout, "delta=%u\n", delta);
	///fprintf(stdout, "line=%u\n", line);

	//gettimeofday(&occ_tv_begin, NULL);
	unsigned int i, j;
	unsigned int b;

	bwt_string_type ans = nacgt[0][delta];

	///这里是增加的
	ans = ans + high_occ_table[(line / high_compress_occ) * 4 + delta];

	/**
	fprintf(stdout, "delta=%llu\n",
		delta);
	fflush(stdout);


	fprintf(stdout, "compress_occ=%llu, acctuall_bwt_gap=%llu, bwt_warp_number=%llu\n", 
		compress_occ, acctuall_bwt_gap, bwt_warp_number);
	fflush(stdout);
	**/

	unsigned int actually_line = ((line / compress_occ)*acctuall_bwt_gap) / bwt_warp_number;

	/**
	fprintf(stdout, "actually_line=%llu, line=%llu\n", actually_line, line);
	fflush(stdout);


	fprintf(stdout, "ans=%llu\n", ans);
	fflush(stdout);
	**/
	
	///这里要改
	///ans += (bwt[actually_line + delta / single_occ_bwt] >> ((single_occ_bwt - 1 - delta%single_occ_bwt) * 32))&mode_32;
	///single_occ_bwt=4
	ans += (bwt[actually_line] >> ((single_occ_bwt - 1 - delta) * 16))&mode_16;

	/**
	fprintf(stdout, "bwt[actually_line + delta / single_occ_bwt]=%llu\n", 
		bwt[actually_line + delta / single_occ_bwt]);


	dec_bit(bwt[actually_line + delta / single_occ_bwt]);

	fprintf(stdout, "ans=%llu\n", ans);
	fflush(stdout);
	**/

	///ans += occ[(line / compress_occ) * 4 + delta];


	if (line%compress_occ == 0)
		return ans;

	unsigned int need_line = line % compress_occ;

	actually_line = actually_line + 4 / single_occ_bwt;


	bwt_string_type tmp_bwt;

	bwt_string_type mode_low_2 = (bwt_string_type)3;

	
	bwt_string_type P_A, P_B;


	i = 0;

	while (i + bwt_warp_number <= need_line)
	{
		P_A = bwt[actually_line ++];
		P_B = P_A^pop_count_mode[delta];
		P_A = P_B >> 1;
		P_A = P_A & P_B;
		P_A = P_A & mode_low;

		ans = ans + __builtin_popcountll(P_A);

		i = i + bwt_warp_number;

	}


	need_line = need_line - i;


	if (need_line!=0)
	{



		P_A = bwt[actually_line++];

		P_B = mode << ((bwt_warp_number - need_line) * 2);

		P_A = P_A & P_B;


		if (delta==0)
		{
			P_B = ~P_B;

			P_A = P_A | P_B;
		}


		P_B = P_A^pop_count_mode[delta];
		P_A = P_B >> 1;
		P_A = P_A & P_B;
		P_A = P_A & mode_low;


		ans = ans + __builtin_popcountll(P_A);


	}



	


	/**
	for (i = 0; i < need_line; i++)
	{
		tmp_bwt = (bwt[actually_line + i / bwt_warp_number] >>
			((bwt_warp_number - (i % bwt_warp_number) - 1) * 2))&mode_low_2;

		if (tmp_bwt == delta)
		{
			ans++;
		}
	}
	**/









	if ((delta == 1) && (shapline >= (line / compress_occ)*compress_occ) && (shapline<line)) 
		ans--;

	return ans;
}


void find_occ_all(unsigned int line, bwt_string_type *bwt,
	unsigned int* ans_A, unsigned int* ans_C, unsigned int* ans_G, unsigned int* ans_T)
{


	unsigned int i;

	(*ans_A) = 0;
	(*ans_C) = 0;
	(*ans_G) = 0;
	(*ans_T) = 0;



	unsigned int actually_line = ((line / compress_occ)*acctuall_bwt_gap) / bwt_warp_number;




	(*ans_A) += (bwt[actually_line]>> 32);

	(*ans_C) += bwt[actually_line] & mode_32;

	actually_line++;

	(*ans_G) += (bwt[actually_line] >> 32);

	(*ans_T) += bwt[actually_line] & mode_32;

	actually_line++;


	if (line%compress_occ == 0)
	{
		(*ans_A) += nacgt[0][0];
		(*ans_C) += nacgt[0][1];
		(*ans_G) += nacgt[0][2];
		(*ans_T) += nacgt[0][3];
		return;
	}
		

	unsigned int need_line = line % compress_occ;



	bwt_string_type P_tmp, P_A, P_B;


	i = 0;

	while (i + bwt_warp_number <= need_line)
	{

		P_tmp = bwt[actually_line++];


		P_A = P_tmp;
		P_B = P_A^pop_count_mode[1];
		P_A = P_B >> 1;
		P_A = P_A & P_B;
		P_A = P_A & mode_low;
		(*ans_C) = (*ans_C) + __builtin_popcountll(P_A);


		P_A = P_tmp;
		P_B = P_A^pop_count_mode[2];
		P_A = P_B >> 1;
		P_A = P_A & P_B;
		P_A = P_A & mode_low;
		(*ans_G) = (*ans_G) + __builtin_popcountll(P_A);


		P_A = P_tmp;
		P_B = P_A^pop_count_mode[3];
		P_A = P_B >> 1;
		P_A = P_A & P_B;
		P_A = P_A & mode_low;
		(*ans_T) = (*ans_T) + __builtin_popcountll(P_A);








		i = i + bwt_warp_number;

	}


	need_line = need_line - i;

	if (need_line != 0)
	{

		P_tmp = bwt[actually_line++];
		P_B = mode << ((bwt_warp_number - need_line) * 2);
		P_tmp = P_tmp & P_B;



		P_A = P_tmp;
		P_B = P_A^pop_count_mode[1];
		P_A = P_B >> 1;
		P_A = P_A & P_B;
		P_A = P_A & mode_low;
		(*ans_C) = (*ans_C) + __builtin_popcountll(P_A);


		P_A = P_tmp;
		P_B = P_A^pop_count_mode[2];
		P_A = P_B >> 1;
		P_A = P_A & P_B;
		P_A = P_A & mode_low;
		(*ans_G) = (*ans_G) + __builtin_popcountll(P_A);


		P_A = P_tmp;
		P_B = P_A^pop_count_mode[3];
		P_A = P_B >> 1;
		P_A = P_A & P_B;
		P_A = P_A & mode_low;
		(*ans_T) = (*ans_T) + __builtin_popcountll(P_A);







	}




	///(*ans_A) = line - (*ans_C) - (*ans_G) - (*ans_T);


	

	if ((shapline >= (line / compress_occ)*compress_occ) && (shapline<line))
		(*ans_C)--;

	(*ans_A) = line - (*ans_C) - (*ans_G) - (*ans_T);


	if (line>shapline)
	{
		(*ans_A)--;
	}






	(*ans_A) += nacgt[0][0];
	(*ans_C) += nacgt[0][1];
	(*ans_G) += nacgt[0][2];
	(*ans_T) += nacgt[0][3];




	return;
}










void find_occ_all_sp_ep(unsigned int sp, unsigned int ep, bwt_string_type *bwt,
	unsigned int* ans_A_sp, unsigned int* ans_C_sp, unsigned int* ans_G_sp, unsigned int* ans_T_sp,
	unsigned int* ans_A_ep, unsigned int* ans_C_ep, unsigned int* ans_G_ep, unsigned int* ans_T_ep)
{


	unsigned int i_sp;
	unsigned int i_ep;

	bwt_string_type P_tmp, P_A, P_B;

	(*ans_A_sp) = 0;
	(*ans_C_sp) = 0;
	(*ans_G_sp) = 0;
	(*ans_T_sp) = 0;


	(*ans_A_ep) = 0;
	(*ans_C_ep) = 0;
	(*ans_G_ep) = 0;
	(*ans_T_ep) = 0;



	unsigned int actually_line_sp = ((sp / compress_occ)*acctuall_bwt_gap) / bwt_warp_number;
	unsigned int actually_line_ep = ((ep / compress_occ)*acctuall_bwt_gap) / bwt_warp_number;



	(*ans_A_sp) += (bwt[actually_line_sp] >> 32);

	(*ans_C_sp) += bwt[actually_line_sp] & mode_32;

	actually_line_sp++;

	(*ans_G_sp) += (bwt[actually_line_sp] >> 32);

	(*ans_T_sp) += bwt[actually_line_sp] & mode_32;

	actually_line_sp++;



	(*ans_A_ep) += (bwt[actually_line_ep] >> 32);

	(*ans_C_ep) += bwt[actually_line_ep] & mode_32;

	actually_line_ep++;

	(*ans_G_ep) += (bwt[actually_line_ep] >> 32);

	(*ans_T_ep) += bwt[actually_line_ep] & mode_32;

	actually_line_ep++;





	unsigned int need_line_sp = sp % compress_occ;
	unsigned int need_line_ep = ep % compress_occ;



	if (need_line_sp != 0)
	{




		i_sp = 0;

		while (i_sp + bwt_warp_number <= need_line_sp)
		{

			P_tmp = bwt[actually_line_sp++];


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_C_sp) = (*ans_C_sp) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[2];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_G_sp) = (*ans_G_sp) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[3];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_T_sp) = (*ans_T_sp) + __builtin_popcountll(P_A);








			i_sp = i_sp + bwt_warp_number;

		}


		need_line_sp = need_line_sp - i_sp;

		if (need_line_sp != 0)
		{

			P_tmp = bwt[actually_line_sp++];
			P_B = mode << ((bwt_warp_number - need_line_sp) * 2);
			P_tmp = P_tmp & P_B;



			P_A = P_tmp;
			P_B = P_A^pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_C_sp) = (*ans_C_sp) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[2];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_G_sp) = (*ans_G_sp) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[3];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_T_sp) = (*ans_T_sp) + __builtin_popcountll(P_A);







		}




		///(*ans_A) = line - (*ans_C) - (*ans_G) - (*ans_T);




		if ((shapline >= (sp / compress_occ)*compress_occ) && (shapline<sp))
			(*ans_C_sp)--;

		(*ans_A_sp) = sp - (*ans_C_sp) - (*ans_G_sp) - (*ans_T_sp);


		if (sp>shapline)
		{
			(*ans_A_sp)--;
		}


	}



	(*ans_A_sp) += nacgt[0][0];
	(*ans_C_sp) += nacgt[0][1];
	(*ans_G_sp) += nacgt[0][2];
	(*ans_T_sp) += nacgt[0][3];


















	if (need_line_ep != 0)
	{




		i_ep = 0;

		while (i_ep + bwt_warp_number <= need_line_ep)
		{

			P_tmp = bwt[actually_line_ep++];


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[2];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[3];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);








			i_ep = i_ep + bwt_warp_number;

		}


		need_line_ep = need_line_ep - i_ep;

		if (need_line_ep != 0)
		{

			P_tmp = bwt[actually_line_ep++];
			P_B = mode << ((bwt_warp_number - need_line_ep) * 2);
			P_tmp = P_tmp & P_B;



			P_A = P_tmp;
			P_B = P_A^pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[2];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[3];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);







		}




		///(*ans_A) = line - (*ans_C) - (*ans_G) - (*ans_T);




		if ((shapline >= (ep / compress_occ)*compress_occ) && (shapline<ep))
			(*ans_C_ep)--;

		(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);


		if (ep>shapline)
		{
			(*ans_A_ep)--;
		}


	}



	(*ans_A_ep) += nacgt[0][0];
	(*ans_C_ep) += nacgt[0][1];
	(*ans_G_ep) += nacgt[0][2];
	(*ans_T_ep) += nacgt[0][3];








	return;
}







void find_occ_all_sp_ep_optimal(unsigned int sp, unsigned int ep, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* ans_A_sp, unsigned int* ans_C_sp, unsigned int* ans_G_sp, unsigned int* ans_T_sp,
	unsigned int* ans_A_ep, unsigned int* ans_C_ep, unsigned int* ans_G_ep, unsigned int* ans_T_ep)
{


	unsigned int i_sp;
	unsigned int i_ep;

	unsigned int flag = 0;

	bwt_string_type P_tmp, P_A, P_B;


	///这里要改
	/**
	(*ans_A_sp) = 0;
	(*ans_C_sp) = 0;
	(*ans_G_sp) = 0;
	(*ans_T_sp) = 0;


	(*ans_A_ep) = 0;
	(*ans_C_ep) = 0;
	(*ans_G_ep) = 0;
	(*ans_T_ep) = 0;
	**/

	unsigned int high_occ_table_line = (sp / high_compress_occ) * 4;

	(*ans_A_sp) = high_occ_table[high_occ_table_line];
	(*ans_C_sp) = high_occ_table[high_occ_table_line + 1];
	(*ans_G_sp) = high_occ_table[high_occ_table_line + 2];
	(*ans_T_sp) = high_occ_table[high_occ_table_line + 3];

	high_occ_table_line = (ep / high_compress_occ) * 4;
	(*ans_A_ep) = high_occ_table[high_occ_table_line];
	(*ans_C_ep) = high_occ_table[high_occ_table_line + 1];
	(*ans_G_ep) = high_occ_table[high_occ_table_line + 2];
	(*ans_T_ep) = high_occ_table[high_occ_table_line + 3];



	unsigned int actually_line_sp = ((sp / compress_occ)*acctuall_bwt_gap) / bwt_warp_number;
	unsigned int actually_line_ep = ((ep / compress_occ)*acctuall_bwt_gap) / bwt_warp_number;

	if (actually_line_sp == actually_line_ep)
	{
		flag = 1;
	}




	

	///这里要改
	/**
	(*ans_A_sp) += (bwt[actually_line_sp] >> 32);

	(*ans_C_sp) += bwt[actually_line_sp] & mode_32;

	actually_line_sp++;

	(*ans_G_sp) += (bwt[actually_line_sp] >> 32);

	(*ans_T_sp) += bwt[actually_line_sp] & mode_32;

	actually_line_sp++;



	(*ans_A_ep) += (bwt[actually_line_ep] >> 32);

	(*ans_C_ep) += bwt[actually_line_ep] & mode_32;

	actually_line_ep++;

	(*ans_G_ep) += (bwt[actually_line_ep] >> 32);

	(*ans_T_ep) += bwt[actually_line_ep] & mode_32;

	actually_line_ep++;
	**/


	(*ans_A_sp) += (bwt[actually_line_sp] >> 48);

	(*ans_C_sp) += (bwt[actually_line_sp] >> 32) & mode_16;

	(*ans_G_sp) += (bwt[actually_line_sp] >> 16) & mode_16;

	(*ans_T_sp) += bwt[actually_line_sp] & mode_16;

	actually_line_sp++;



	(*ans_A_ep) += (bwt[actually_line_ep] >> 48);

	(*ans_C_ep) += (bwt[actually_line_ep] >> 32) & mode_16;

	(*ans_G_ep) += (bwt[actually_line_ep] >> 16) & mode_16;

	(*ans_T_ep) += bwt[actually_line_ep] & mode_16;

	actually_line_ep++;



	unsigned int need_line_sp = sp % compress_occ;
	unsigned int need_line_ep = ep % compress_occ;



	if (need_line_sp != 0)
	{




		i_sp = 0;

		while (i_sp + bwt_warp_number <= need_line_sp)
		{

			P_tmp = bwt[actually_line_sp++];


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_C_sp) = (*ans_C_sp) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[2];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_G_sp) = (*ans_G_sp) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[3];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_T_sp) = (*ans_T_sp) + __builtin_popcountll(P_A);








			i_sp = i_sp + bwt_warp_number;

		}


		if (flag == 1)
		{
			i_ep = i_sp;

			(*ans_C_ep) = (*ans_C_sp);

			(*ans_G_ep) = (*ans_G_sp);

			(*ans_T_ep) = (*ans_T_sp);

			actually_line_ep = actually_line_sp;

			flag = 2;
		}


		need_line_sp = need_line_sp - i_sp;

		if (need_line_sp != 0)
		{

			P_tmp = bwt[actually_line_sp++];
			P_B = mode << ((bwt_warp_number - need_line_sp) * 2);
			P_tmp = P_tmp & P_B;



			P_A = P_tmp;
			P_B = P_A^pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_C_sp) = (*ans_C_sp) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[2];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_G_sp) = (*ans_G_sp) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^pop_count_mode[3];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & mode_low;
			(*ans_T_sp) = (*ans_T_sp) + __builtin_popcountll(P_A);







		}




		///(*ans_A) = line - (*ans_C) - (*ans_G) - (*ans_T);




		if ((shapline >= (sp / compress_occ)*compress_occ) && (shapline<sp))
			(*ans_C_sp)--;

		(*ans_A_sp) = sp - (*ans_C_sp) - (*ans_G_sp) - (*ans_T_sp);


		if (sp>shapline)
		{
			(*ans_A_sp)--;
		}


	}



	(*ans_A_sp) += nacgt[0][0];
	(*ans_C_sp) += nacgt[0][1];
	(*ans_G_sp) += nacgt[0][2];
	(*ans_T_sp) += nacgt[0][3];


















	if (need_line_ep != 0)
	{


		if (flag == 2)
		{



			while (i_ep + bwt_warp_number <= need_line_ep)
			{

				P_tmp = bwt[actually_line_ep++];


				P_A = P_tmp;
				P_B = P_A^pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & mode_low;
				(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^pop_count_mode[2];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & mode_low;
				(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^pop_count_mode[3];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);








				i_ep = i_ep + bwt_warp_number;

			}


			need_line_ep = need_line_ep - i_ep;

			if (need_line_ep != 0)
			{

				P_tmp = bwt[actually_line_ep++];
				P_B = mode << ((bwt_warp_number - need_line_ep) * 2);
				P_tmp = P_tmp & P_B;



				P_A = P_tmp;
				P_B = P_A^pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & mode_low;
				(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^pop_count_mode[2];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & mode_low;
				(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^pop_count_mode[3];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);







			}




			///(*ans_A) = line - (*ans_C) - (*ans_G) - (*ans_T);




			if ((shapline >= (ep / compress_occ)*compress_occ) && (shapline<ep))
				(*ans_C_ep)--;

			(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);


			if (ep>shapline)
			{
				(*ans_A_ep)--;
			}












		}
		else
		{



			i_ep = 0;

			while (i_ep + bwt_warp_number <= need_line_ep)
			{

				P_tmp = bwt[actually_line_ep++];


				P_A = P_tmp;
				P_B = P_A^pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & mode_low;
				(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^pop_count_mode[2];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & mode_low;
				(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^pop_count_mode[3];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);








				i_ep = i_ep + bwt_warp_number;

			}


			need_line_ep = need_line_ep - i_ep;

			if (need_line_ep != 0)
			{

				P_tmp = bwt[actually_line_ep++];
				P_B = mode << ((bwt_warp_number - need_line_ep) * 2);
				P_tmp = P_tmp & P_B;



				P_A = P_tmp;
				P_B = P_A^pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & mode_low;
				(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^pop_count_mode[2];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & mode_low;
				(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^pop_count_mode[3];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);







			}




			///(*ans_A) = line - (*ans_C) - (*ans_G) - (*ans_T);




			if ((shapline >= (ep / compress_occ)*compress_occ) && (shapline<ep))
				(*ans_C_ep)--;

			(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);


			if (ep>shapline)
			{
				(*ans_A_ep)--;
			}

		}




	}



	(*ans_A_ep) += nacgt[0][0];
	(*ans_C_ep) += nacgt[0][1];
	(*ans_G_ep) += nacgt[0][2];
	(*ans_T_ep) += nacgt[0][3];








	return;
}




unsigned int get_sa(SA_flag_string_type* SA_flag, unsigned int line, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int *sa)
{
	unsigned int l = line;
	bwt_string_type delta;
	unsigned int i = 0;


	if (line == shapline) return 0;

	/**
	fprintf(stdout, "l=%llu, compress_SA_flag=%llu, SA_flag_warp_number=%llu, acctuall_SA_flag_gap=%llu\n", 
		l, compress_SA_flag, SA_flag_warp_number, acctuall_SA_flag_gap);
	**/

	unsigned int actually_line 
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;

	unsigned int last = l % compress_SA_flag;

	/**
	fprintf(stdout, "actually_line=%llu, last=%llu\n",
		actually_line, last);
	**/

	actually_line = actually_line + last / SA_flag_warp_number;

	/**
	fprintf(stdout, "actually_line=%llu\n",
		actually_line);
	**/
	


	///while ((l%compress_sa) != 0)
	while (((SA_flag[actually_line] << (last % SA_flag_warp_number))&mode_SA_flag)==0)
	{

		actually_line = ((l / compress_occ)*acctuall_bwt_gap) + l % compress_occ + 64;

		delta = (bwt[actually_line / bwt_warp_number] 
			>> ((bwt_warp_number - actually_line%bwt_warp_number - 1) * 2))
			&(bwt_string_type)3;

		///fprintf(stdout, "delta=%u\n", delta);

		l = find_occ(l, delta, bwt, high_occ_table);

		///fprintf(stdout, "l=%u\n", l);

		//printf("%d!\n",l);
		i++;

		if (l == shapline) return i;



		actually_line
			= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;

		last = l % compress_SA_flag;

		actually_line = actually_line + last / SA_flag_warp_number;

	}


	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	unsigned int ans = SA_flag[actually_line];
	

	if (last!=0)
	{
		actually_line++;

		unsigned int j = 0;


		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcount(SA_flag[actually_line++]);

			j = j + SA_flag_warp_number;

		}


		last = last - j;


		if (last != 0)
		{
			mode << (SA_flag_warp_number - last);
			ans = ans + 
				__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line++] 
				& (mode << (SA_flag_warp_number - last))));
		}


	}
	/**
	fprintf(stdout, "ans=%llu\n",
		ans);
	**/

	return sa[ans] + i;
}











unsigned int get_sa_restrict_steps_more_than_3
(SA_flag_string_type* SA_flag, unsigned int line, 
bwt_string_type *bwt, high_occ_table_type* high_occ_table,
unsigned int *sa, unsigned int need_step, unsigned int* accessed_sa)
{
	unsigned int l = line;
	bwt_string_type delta;
	unsigned int i = 0;


	if (line == shapline) return 0;


	///这里SA要改
	///unsigned int actually_line= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;
	unsigned int actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	///这里SA要改
	///unsigned int last = l % compress_SA_flag;
	unsigned int last = l % compress_SA_flag + SA_counter_length;


	actually_line = actually_line + last / SA_flag_warp_number;




	while (((SA_flag[actually_line] << (last % SA_flag_warp_number))&mode_SA_flag) == 0)
	{

		if (i>=need_step)
		{
			return 0;
		}


		///似乎只有这里要改，因为这里涉及到取BWT了
		///这里做了修改
		///actually_line = ((l / compress_occ)*acctuall_bwt_gap) + l % compress_occ + 64;
		actually_line = ((l / compress_occ)*acctuall_bwt_gap) + l % compress_occ + 32;
		delta = (bwt[actually_line / bwt_warp_number]
			>> ((bwt_warp_number - actually_line%bwt_warp_number - 1) * 2))
			&(bwt_string_type)3;


		l = find_occ(l, delta, bwt, high_occ_table);


		i++;

		if (l == shapline)
		{
			(*accessed_sa) = i;

			return i;
		}
			
			
			


		///这里SA要改
		///actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;
		actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

		///这里SA要改
		///last = l % compress_SA_flag;
		last = l % compress_SA_flag + SA_counter_length;

		actually_line = actually_line + last / SA_flag_warp_number;

	}


	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	///这里SA要改
	///unsigned int ans = SA_flag[actually_line];
	unsigned int ans = SA_flag[actually_line] >> SA_counter_shift_length;


	///这里SA要改
	/**
	if (last != 0)
	{
		actually_line++;



		unsigned int j = 0;




		///while (j + SA_flag_warp_number <= last)
		while (j + 64 <= last)
		{
			tmp_SA_flag = (uint64_t)0;
			tmp_SA_flag = SA_flag[actually_line++];
			tmp_SA_flag = tmp_SA_flag << SA_flag_warp_number;
			tmp_SA_flag = tmp_SA_flag | SA_flag[actually_line++];
			ans = ans + __builtin_popcountll(tmp_SA_flag);

			j = j + 64;

		}


		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcount(SA_flag[actually_line++]);

			j = j + SA_flag_warp_number;

		}





		last = last - j;


		if (last != 0)
		{
			///mode << (SA_flag_warp_number - last);
			ans = ans +
				__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line++]
				& (mode << (SA_flag_warp_number - last))));
		}


	}
	**/
	if (last != 0)
	{
		///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
		unsigned int j = 0;
		SA_flag_string_type tmp_SA_pop_count;
		tmp_SA_pop_count = SA_flag[actually_line] & SA_pop_count_mode;




		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (mode << (SA_flag_warp_number - last))));
		}


	}

	(*accessed_sa) = (sa[ans] & SA_header_mode)*compress_sa + i;

	return 1;
}





unsigned int get_sa_restrict_steps
(SA_flag_string_type* SA_flag, unsigned int line,
bwt_string_type *bwt, high_occ_table_type* high_occ_table,
unsigned int *sa, unsigned int need_step, unsigned int* accessed_sa)
{
	unsigned int l = line;
	bwt_string_type delta;
	unsigned int i = 0;


	if (line == shapline) return 0;


	///这里SA要改
	///unsigned int actually_line= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;
	unsigned int actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	///这里SA要改
	///unsigned int last = l % compress_SA_flag;
	unsigned int last = l % compress_SA_flag + SA_counter_length;


	actually_line = actually_line + last / SA_flag_warp_number;




	while (((SA_flag[actually_line] << (last % SA_flag_warp_number))&mode_SA_flag) == 0)
	{

		if (i >= need_step)
		{
			return 0;
		}


		///似乎只有这里要改，因为这里涉及到取BWT了
		///这里做了修改
		///actually_line = ((l / compress_occ)*acctuall_bwt_gap) + l % compress_occ + 64;
		actually_line = ((l / compress_occ)*acctuall_bwt_gap) + l % compress_occ + 32;
		delta = (bwt[actually_line / bwt_warp_number]
			>> ((bwt_warp_number - actually_line%bwt_warp_number - 1) * 2))
			&(bwt_string_type)3;


		l = find_occ(l, delta, bwt, high_occ_table);


		i++;

		if (l == shapline)
		{
			(*accessed_sa) = i;

			return i;
		}





		///这里SA要改
		///actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;
		actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

		///这里SA要改
		///last = l % compress_SA_flag;
		last = l % compress_SA_flag + SA_counter_length;

		actually_line = actually_line + last / SA_flag_warp_number;

	}


	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	///这里SA要改
	///unsigned int ans = SA_flag[actually_line];
	unsigned int ans = SA_flag[actually_line] >> SA_counter_shift_length;


	///这里SA要改
	if (last != 0)
	{
		///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
		unsigned int j = 0;
		SA_flag_string_type tmp_SA_pop_count;
		tmp_SA_pop_count = SA_flag[actually_line] & SA_pop_count_mode;




		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (mode << (SA_flag_warp_number - last))));
		}


	}


	(*accessed_sa) = sa[ans] + i;

	return 1;

}


/**

unsigned int get_sa_restrict_steps
(SA_flag_string_type* SA_flag, unsigned int line,
bwt_string_type *bwt, high_occ_table_type* high_occ_table, unsigned int *sa, unsigned int need_step, unsigned int* accessed_sa)
{
	unsigned int l = line;
	bwt_string_type delta;
	unsigned int i = 0;


	if (line == shapline) return 0;



	unsigned int actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;

	unsigned int last = l % compress_SA_flag;



	actually_line = actually_line + last / SA_flag_warp_number;




	while (((SA_flag[actually_line] << (last % SA_flag_warp_number))&mode_SA_flag) == 0)
	{

		if (i >= need_step)
		{
			return 0;
		}





		///似乎只有这里要改，因为这里涉及到取BWT了
		///这里做了修改



		actually_line = ((l / compress_occ)*acctuall_bwt_gap) + l % compress_occ + 32;
		delta = (bwt[actually_line / bwt_warp_number]
			>> ((bwt_warp_number - actually_line%bwt_warp_number - 1) * 2))
			&(bwt_string_type)3;


		l = find_occ(l, delta, bwt, high_occ_table);


		i++;

		if (l == shapline)
		{
			(*accessed_sa) = i;

			return i;
		}






		actually_line
			= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;

		last = l % compress_SA_flag;

		actually_line = actually_line + last / SA_flag_warp_number;

	}


	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	unsigned int ans = SA_flag[actually_line];


	if (last != 0)
	{
		actually_line++;



		unsigned int j = 0;




		///while (j + SA_flag_warp_number <= last)
		while (j + 64 <= last)
		{
			tmp_SA_flag = (uint64_t)0;
			tmp_SA_flag = SA_flag[actually_line++];
			tmp_SA_flag = tmp_SA_flag << SA_flag_warp_number;
			tmp_SA_flag = tmp_SA_flag | SA_flag[actually_line++];
			ans = ans + __builtin_popcountll(tmp_SA_flag);

			j = j + 64;

		}


		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcount(SA_flag[actually_line++]);

			j = j + SA_flag_warp_number;

		}





		last = last - j;


		if (last != 0)
		{
			///mode << (SA_flag_warp_number - last);
			ans = ans +
				__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line++]
				& (mode << (SA_flag_warp_number - last))));
		}


	}


	(*accessed_sa) = sa[ans] + i;

	return 1;
}
**/





unsigned int get_sa_restrict_zero_steps_more_than_3
(SA_flag_string_type* SA_flag, unsigned int *sa, unsigned int sp, unsigned int* start)
{
	unsigned int l;
	unsigned int i = 0;
	///uint64_t tmp_SA_flag;
	unsigned int last;

	unsigned int actually_line;

	unsigned int ans;

	unsigned int last_wods;

	SA_flag_string_type tmp_SA_pop_count;


	/**
	l = sp;


	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;

	last = l % compress_SA_flag;

	actually_line = actually_line + last / SA_flag_warp_number;

	if (((SA_flag[actually_line] << (last % SA_flag_warp_number))&mode_SA_flag) == 0)
		return 0;



	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;



	ans = SA_flag[actually_line];




	if (last != 0)
	{
		actually_line++;

		unsigned int j = 0;


		while (j + 64 <= last)
		{
			tmp_SA_flag = (uint64_t)0;
			tmp_SA_flag = SA_flag[actually_line++];
			tmp_SA_flag = tmp_SA_flag << SA_flag_warp_number;
			tmp_SA_flag = tmp_SA_flag | SA_flag[actually_line++];
			ans = ans + __builtin_popcountll(tmp_SA_flag);

			j = j + 64;

		}


		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcount(SA_flag[actually_line++]);

			j = j + SA_flag_warp_number;

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line++]
				& (mode << (SA_flag_warp_number - last))));
		}


	}
	**/



	l = sp;


	actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	///这里SA要改
	///last = l % compress_SA_flag;
	last = l % compress_SA_flag + SA_counter_length;

	///这里SA要改
	///ans = SA_flag[actually_line];
	ans = SA_flag[actually_line] >> SA_counter_shift_length;


	///这里SA要改
	if (last != 0)
	{
		///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
		unsigned int j = 0;

		tmp_SA_pop_count = SA_flag[actually_line] & SA_pop_count_mode;




		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (mode << (SA_flag_warp_number - last))));
		}


	}


	///(*start) = sa[ans];
	(*start) = (sa[ans] & SA_header_mode)*compress_sa;


	return 1;

}




unsigned int get_sa_restrict_zero_steps
(SA_flag_string_type* SA_flag, unsigned int *sa, unsigned int sp, unsigned int* start)
{
	unsigned int l;
	unsigned int i = 0;
	///uint64_t tmp_SA_flag;
	unsigned int last;

	unsigned int actually_line;

	unsigned int ans;

	unsigned int last_wods;

	SA_flag_string_type tmp_SA_pop_count;


	/**
	l = sp;


	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;

	last = l % compress_SA_flag;

	actually_line = actually_line + last / SA_flag_warp_number;

	if (((SA_flag[actually_line] << (last % SA_flag_warp_number))&mode_SA_flag) == 0)
		return 0;



	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;



	ans = SA_flag[actually_line];




	if (last != 0)
	{
		actually_line++;

		unsigned int j = 0;


		while (j + 64 <= last)
		{
			tmp_SA_flag = (uint64_t)0;
			tmp_SA_flag = SA_flag[actually_line++];
			tmp_SA_flag = tmp_SA_flag << SA_flag_warp_number;
			tmp_SA_flag = tmp_SA_flag | SA_flag[actually_line++];
			ans = ans + __builtin_popcountll(tmp_SA_flag);

			j = j + 64;

		}


		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcount(SA_flag[actually_line++]);

			j = j + SA_flag_warp_number;

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line++]
				& (mode << (SA_flag_warp_number - last))));
		}


	}
	**/



	l = sp;


	actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	///这里SA要改
	///last = l % compress_SA_flag;
	last = l % compress_SA_flag + SA_counter_length;

	///这里SA要改
	///ans = SA_flag[actually_line];
	ans = SA_flag[actually_line] >> SA_counter_shift_length;


	///这里SA要改
	if (last != 0)
	{
		///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
		unsigned int j = 0;

		tmp_SA_pop_count = SA_flag[actually_line] & SA_pop_count_mode;




		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (mode << (SA_flag_warp_number - last))));
		}


	}


	(*start) = sa[ans];

	return 1;

}









void direct_get_sa_interval
(SA_flag_string_type* SA_flag, unsigned int sp, unsigned int ep, unsigned int* start, unsigned int* length)
{
	unsigned int l;
	unsigned int i = 0;

	unsigned int last;

	unsigned int actually_line;

	unsigned int ans;


	l = sp;


	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	last = l % compress_SA_flag;

	ans = SA_flag[actually_line];


	if (last != 0)
	{
		actually_line++;

		unsigned int j = 0;


		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcount(SA_flag[actually_line++]);

			j = j + SA_flag_warp_number;

		}


		last = last - j;


		if (last != 0)
		{
			mode << (SA_flag_warp_number - last);
			ans = ans +
				__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line++]
				& (mode << (SA_flag_warp_number - last))));
		}


	}


	(*start) = ans;



	///return sa[ans] + i;














	l = ep;


	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	last = l % compress_SA_flag;

	ans = SA_flag[actually_line];


	if (last != 0)
	{
		actually_line++;

		unsigned int j = 0;


		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcount(SA_flag[actually_line++]);

			j = j + SA_flag_warp_number;

		}


		last = last - j;


		if (last != 0)
		{
			mode << (SA_flag_warp_number - last);
			ans = ans +
				__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line++]
				& (mode << (SA_flag_warp_number - last))));
		}


	}




	(*length) = ans - (*start);

}



void direct_get_sa_interval_long
(SA_flag_string_type* SA_flag, unsigned int sp, unsigned int ep, unsigned int* start, unsigned int* length)
{
	unsigned int l;
	unsigned int i = 0;
	///uint64_t tmp_SA_flag;
	unsigned int last;

	unsigned int actually_line;

	unsigned int ans;

	SA_flag_string_type tmp_SA_pop_count;


	l = sp;


	actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	///这里SA要改
	///last = l % compress_SA_flag;
	last = l % compress_SA_flag + SA_counter_length;

	///这里SA要改
	///ans = SA_flag[actually_line];
	ans = SA_flag[actually_line] >> SA_counter_shift_length;


	///这里SA要改
	/**
	if (last != SA_counter_length)
	{
		actually_line++;

		unsigned int j = 0;

		


		///while (j + SA_flag_warp_number <= last)
		while (j + 64 <= last)
		{
			tmp_SA_flag = (uint64_t)0;
			tmp_SA_flag = SA_flag[actually_line++];
			tmp_SA_flag = tmp_SA_flag << SA_flag_warp_number;
			tmp_SA_flag = tmp_SA_flag | SA_flag[actually_line++];
			ans = ans + __builtin_popcountll(tmp_SA_flag);

			j = j + 64;

		}


		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcount(SA_flag[actually_line++]);

			j = j + SA_flag_warp_number;

		}


		last = last - j;


		if (last != 0)
		{
			///mode << (SA_flag_warp_number - last);
			ans = ans +
				__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line++]
				& (mode << (SA_flag_warp_number - last))));
		}


	}
	**/
	if (last != 0)
	{
		///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
		unsigned int j = 0;
		
		tmp_SA_pop_count = SA_flag[actually_line] & SA_pop_count_mode;




		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (mode << (SA_flag_warp_number - last))));
		}


	}




	(*start) = ans;



	///return sa[ans] + i;














	l = ep;

	/**
	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	last = l % compress_SA_flag;

	ans = SA_flag[actually_line];


	if (last != 0)
	{
		actually_line++;

		unsigned int j = 0;




		while (j + 64 <= last)
		{
			tmp_SA_flag = (uint64_t)0;
			tmp_SA_flag = SA_flag[actually_line++];
			tmp_SA_flag = tmp_SA_flag << SA_flag_warp_number;
			tmp_SA_flag = tmp_SA_flag | SA_flag[actually_line++];
			ans = ans + __builtin_popcountll(tmp_SA_flag);

			j = j + 64;

		}


		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcount(SA_flag[actually_line++]);

			j = j + SA_flag_warp_number;

		}


		last = last - j;


		if (last != 0)
		{
			////mode << (SA_flag_warp_number - last);
			ans = ans +
				__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line++]
				& (mode << (SA_flag_warp_number - last))));
		}


	}
	**/



	actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	///这里SA要改
	///last = l % compress_SA_flag;
	last = l % compress_SA_flag + SA_counter_length;

	///这里SA要改
	///ans = SA_flag[actually_line];
	ans = SA_flag[actually_line] >> SA_counter_shift_length;


	///这里SA要改
	/**
	if (last != SA_counter_length)
	{
	actually_line++;

	unsigned int j = 0;




	///while (j + SA_flag_warp_number <= last)
	while (j + 64 <= last)
	{
	tmp_SA_flag = (uint64_t)0;
	tmp_SA_flag = SA_flag[actually_line++];
	tmp_SA_flag = tmp_SA_flag << SA_flag_warp_number;
	tmp_SA_flag = tmp_SA_flag | SA_flag[actually_line++];
	ans = ans + __builtin_popcountll(tmp_SA_flag);

	j = j + 64;

	}


	while (j + SA_flag_warp_number <= last)
	{

	ans = ans + __builtin_popcount(SA_flag[actually_line++]);

	j = j + SA_flag_warp_number;

	}


	last = last - j;


	if (last != 0)
	{
	///mode << (SA_flag_warp_number - last);
	ans = ans +
	__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line++]
	& (mode << (SA_flag_warp_number - last))));
	}


	}
	**/
	if (last != 0)
	{
		///主要last本身已经加过SA_counter_length了，所以这里j=0；否则就应该把j = SA_counter_length
		unsigned int j = 0;

		tmp_SA_pop_count = SA_flag[actually_line] & SA_pop_count_mode;




		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (mode << (SA_flag_warp_number - last))));
		}


	}




	(*length) = ans - (*start);

}

void direct_get_sa_interval_long_optimal
(SA_flag_string_type* SA_flag, unsigned int sp, unsigned int ep, unsigned int* start, unsigned int* length)
{
	
	unsigned int i = 0;
	uint64_t tmp_SA_flag;
	unsigned int last_sp;

	unsigned int last_ep;

	unsigned int actually_line_sp;

	unsigned int actually_line_ep;

	unsigned int ans_sp;

	unsigned int ans_ep;

	unsigned int j_sp = 0;

	unsigned int j_ep = 0;

	unsigned int flag = 0;




	actually_line_sp
		= ((sp / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	last_sp = sp % compress_SA_flag;

	ans_sp = SA_flag[actually_line_sp];




	actually_line_ep
		= ((ep / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	last_ep = ep % compress_SA_flag;

	ans_ep = SA_flag[actually_line_ep];

	if (actually_line_sp == actually_line_ep)
	{
		flag = 1;
	}




	if (last_sp != 0)
	{
		actually_line_sp++;

		


		while (j_sp + 64 <= last_sp)
		{
			tmp_SA_flag = (uint64_t)0;
			tmp_SA_flag = SA_flag[actually_line_sp++];
			tmp_SA_flag = tmp_SA_flag << SA_flag_warp_number;
			tmp_SA_flag = tmp_SA_flag | SA_flag[actually_line_sp++];
			ans_sp = ans_sp + __builtin_popcountll(tmp_SA_flag);

			j_sp = j_sp + 64;

		}

		
		if (flag==1)
		{
			actually_line_ep = actually_line_sp;

			j_ep = j_sp;

			ans_ep = ans_sp;

			flag = 2;


		}
		


		while (j_sp + SA_flag_warp_number <= last_sp)
		{

			ans_sp = ans_sp + __builtin_popcount(SA_flag[actually_line_sp++]);

			j_sp = j_sp + SA_flag_warp_number;

		}


		last_sp = last_sp - j_sp;


		if (last_sp != 0)
		{
			ans_sp = ans_sp +
				__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line_sp++]
				& (mode << (SA_flag_warp_number - last_sp))));
		}


	}


	(*start) = ans_sp;



















	


	if (last_ep != 0)
	{


		if (flag == 2)
		{

			while (j_ep + 64 <= last_ep)
			{
				tmp_SA_flag = (uint64_t)0;
				tmp_SA_flag = SA_flag[actually_line_ep++];
				tmp_SA_flag = tmp_SA_flag << SA_flag_warp_number;
				tmp_SA_flag = tmp_SA_flag | SA_flag[actually_line_ep++];
				ans_ep = ans_ep + __builtin_popcountll(tmp_SA_flag);

				j_ep = j_ep + 64;

			}


			while (j_ep + SA_flag_warp_number <= last_ep)
			{

				ans_ep = ans_ep + __builtin_popcount(SA_flag[actually_line_ep++]);

				j_ep = j_ep + SA_flag_warp_number;

			}


			last_ep = last_ep - j_ep;


			if (last_ep != 0)
			{
				ans_ep = ans_ep +
					__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line_ep++]
					& (mode << (SA_flag_warp_number - last_ep))));
			}


		}
		else
		{


			actually_line_ep++;

			while (j_ep + 64 <= last_ep)
			{
				tmp_SA_flag = (uint64_t)0;
				tmp_SA_flag = SA_flag[actually_line_ep++];
				tmp_SA_flag = tmp_SA_flag << SA_flag_warp_number;
				tmp_SA_flag = tmp_SA_flag | SA_flag[actually_line_ep++];
				ans_ep = ans_ep + __builtin_popcountll(tmp_SA_flag);

				j_ep = j_ep + 64;

			}


			while (j_ep + SA_flag_warp_number <= last_ep)
			{

				ans_ep = ans_ep + __builtin_popcount(SA_flag[actually_line_ep++]);

				j_ep = j_ep + SA_flag_warp_number;

			}


			last_ep = last_ep - j_ep;


			if (last_ep != 0)
			{
				ans_ep = ans_ep +
					__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line_ep++]
					& (mode << (SA_flag_warp_number - last_ep))));
			}
		}


	}




	(*length) = ans_ep - (*start);
	

}



void direct_get_sa_interval_direct_long
(SA_flag_string_type* SA_flag, unsigned int sp, unsigned int ep, unsigned int* start, unsigned int* length)
{
	unsigned int l;
	unsigned int i = 0;

	unsigned int j = 0;
	unsigned int long_j;

	unsigned int last;

	unsigned int actually_line;

	unsigned int ans;


	l = sp;


	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	last = l % compress_SA_flag;

	ans = SA_flag[actually_line];


	if (last != 0)
	{
		actually_line++;

		
		j = 0;
		long_j = 0;
		
		long_SA_flag = (uint64_t*)(SA_flag+actually_line);

		///while (j + SA_flag_warp_number <= last)
		while (j + 64 <= last)
		{

			ans = ans + __builtin_popcountll(long_SA_flag[long_j++]);

			j = j + 64;

		}

		actually_line = long_j * 2;

		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcount(SA_flag[actually_line++]);

			j = j + SA_flag_warp_number;

		}


		last = last - j;


		if (last != 0)
		{
			mode << (SA_flag_warp_number - last);
			ans = ans +
				__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line++]
				& (mode << (SA_flag_warp_number - last))));
		}


	}


	(*start) = ans;



	///return sa[ans] + i;














	l = ep;


	actually_line
		= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number;

	last = l % compress_SA_flag;

	ans = SA_flag[actually_line];


	if (last != 0)
	{
		actually_line++;

		j = 0;
		long_j = 0;

		long_SA_flag = (uint64_t*)(SA_flag + actually_line);

		while (j + 64 <= last)
		{

			ans = ans + __builtin_popcountll(long_SA_flag[long_j++]);

			j = j + 64;

		}

		actually_line = long_j * 2;

		while (j + SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcount(SA_flag[actually_line++]);

			j = j + SA_flag_warp_number;

		}


		last = last - j;


		if (last != 0)
		{
			mode << (SA_flag_warp_number - last);
			ans = ans +
				__builtin_popcount((SA_flag_string_type)(SA_flag[actually_line++]
				& (mode << (SA_flag_warp_number - last))));
		}


	}




	(*length) = ans - (*start);

}





void accesss_SA(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length)
{
	unsigned int SA_start, SA_length, t;

	if (ep - sp>1)
	{
		///这里SA要改
		direct_get_sa_interval_long(SA_flag, sp, ep, &SA_start, &SA_length);




		for (t = 0; t<SA_length; t++)
		{
			locates[(*tmp_SA_length)++] = sa[t + SA_start] + occ;
		}
	}
	else if (ep - sp == 1)
	{
		///这里SA要改
		if (get_sa_restrict_zero_steps(SA_flag,
			sa, sp, (&(locates[(*tmp_SA_length)]))) == 1)
		{
			locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] + occ;
			(*tmp_SA_length)++;

		}
	}

}


void accesss_pre_SA_debug(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep,
	unsigned int* tmp_SA_length, unsigned int delta, bwt_string_type *bwt)
{
	unsigned int SA_start, SA_length, t, tmp_SA, tmp_SA_tail, aviable_BWT_single;

	unsigned int* aviable_BWT_index = NULL;

	unsigned int actually_line, l, last, aviable_BWT_index_i = 0;
	bwt_string_type tmp_delta;

	if (ep - sp>1)
	{

		direct_get_sa_interval_long(SA_flag, sp, ep, &SA_start, &SA_length);


		///aviable_BWT_index = (unsigned int*)malloc(sizeof(unsigned int)*SA_length);

		for (t = sp; t<ep; t++)
		{
			l = t;


			actually_line
				= ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;

			last = l % compress_SA_flag;

			actually_line = actually_line + last / SA_flag_warp_number;


			if (((SA_flag[actually_line] << (last % SA_flag_warp_number))&mode_SA_flag) == 1)
			{
				actually_line = ((l / compress_occ)*acctuall_bwt_gap) + l % compress_occ + 64;

				tmp_delta = (bwt[actually_line / bwt_warp_number]
					>> ((bwt_warp_number - actually_line%bwt_warp_number - 1) * 2))
					&(bwt_string_type)3;


				if (tmp_delta == delta)
				{
					tmp_SA = sa[aviable_BWT_index_i + SA_start];

					tmp_SA_tail = (tmp_SA & SA_header_mode)*compress_sa;

					if (tmp_SA_tail != 0)
					{
						locates[(*tmp_SA_length)++] = tmp_SA_tail - 1;
					}
				}

				aviable_BWT_index_i++;
				
			}
		}



		/**
		for (t = 0; t<SA_length; t++)
		{
			tmp_SA = sa[t + SA_start];
			tmp_SA_tail = (tmp_SA & SA_header_mode)*compress_sa;



			if (tmp_SA_tail != 0 &&
				((SA_header_mode_reverse&tmp_SA) >> 30) == delta)
			{
				locates[(*tmp_SA_length)++] = tmp_SA_tail - 1;
			}


			///locates[(*tmp_SA_length)++] = sa[t + SA_start] & SA_header_mode - 1;

		}


		free(aviable_BWT_index);
		**/
	}
	else if (ep - sp == 1)
	{


		


		if (get_sa_restrict_zero_steps(SA_flag,
			sa, sp, (&(locates[(*tmp_SA_length)]))) == 1)
		{



			l = sp;

			actually_line = ((l / compress_occ)*acctuall_bwt_gap) + l % compress_occ + 64;

			tmp_delta = (bwt[actually_line / bwt_warp_number]
				>> ((bwt_warp_number - actually_line%bwt_warp_number - 1) * 2))
				&(bwt_string_type)3;


			if (tmp_delta == delta)
			{

				/**
				tmp_SA = sa[aviable_BWT_index_i + SA_start];

				tmp_SA_tail = (tmp_SA & SA_header_mode)*compress_sa;

				**/

				tmp_SA = locates[(*tmp_SA_length)];
				tmp_SA_tail = (tmp_SA & SA_header_mode)*compress_sa;


				if (tmp_SA_tail != 0)
				{
					locates[(*tmp_SA_length)] = tmp_SA_tail - 1;
					(*tmp_SA_length)++;
				}
			}



			/**
			tmp_SA = locates[(*tmp_SA_length)];
			tmp_SA_tail = (tmp_SA & SA_header_mode)*compress_sa;



			if (tmp_SA_tail != 0 &&
				((SA_header_mode_reverse&tmp_SA) >> 30) == delta)
			{
				locates[(*tmp_SA_length)] = tmp_SA_tail - 1;
				(*tmp_SA_length)++;
			}
			**/


			///locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] & SA_header_mode - 1;
			///(*tmp_SA_length)++;

		}
	}

}




void accesss_pre_SA(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep,
	unsigned int* tmp_SA_length, unsigned int delta)
{
	unsigned int SA_start, SA_length, t, tmp_SA, tmp_SA_tail;

	if (ep - sp>1)
	{
		///这里SA要改
		direct_get_sa_interval_long(SA_flag, sp, ep, &SA_start, &SA_length);




		for (t = 0; t<SA_length; t++)
		{
			tmp_SA = sa[t + SA_start];
			tmp_SA_tail = (tmp_SA & SA_header_mode)*compress_sa;



			if (tmp_SA_tail != 0&&
				((SA_header_mode_reverse&tmp_SA) >> 30) == delta)
			{
				locates[(*tmp_SA_length)++] = tmp_SA_tail - 1;
			}

				
			///locates[(*tmp_SA_length)++] = sa[t + SA_start] & SA_header_mode - 1;
			
		}
	}
	else if (ep - sp == 1)
	{




		///这里SA要改
		if (get_sa_restrict_zero_steps(SA_flag,
			sa, sp, (&(locates[(*tmp_SA_length)]))) == 1)
		{


			tmp_SA = locates[(*tmp_SA_length)];
			tmp_SA_tail = (tmp_SA & SA_header_mode)*compress_sa;



			if (tmp_SA_tail != 0 &&
				((SA_header_mode_reverse&tmp_SA) >> 30) == delta)
			{
				locates[(*tmp_SA_length)] = tmp_SA_tail - 1;
				(*tmp_SA_length)++;
			}



			///locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] & SA_header_mode - 1;
			///(*tmp_SA_length)++;

		}
	}

}



void accesss_SA_more_than_3(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length)
{
	unsigned int SA_start, SA_length, t;

	if (ep - sp>1)
	{


		///这里SA要改
		direct_get_sa_interval_long(SA_flag, sp, ep, &SA_start, &SA_length);




		for (t = 0; t<SA_length; t++)
		{
			locates[(*tmp_SA_length)++] = (sa[t + SA_start] & SA_header_mode)*compress_sa + occ;
		}
	}
	else if (ep - sp == 1)
	{
		///这里SA要改
		if (get_sa_restrict_zero_steps_more_than_3(SA_flag,
			sa, sp, (&(locates[(*tmp_SA_length)]))) == 1)
		{
			locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] + occ;
			(*tmp_SA_length)++;

		}

	}

}






void accesss_SA_cur_more_than_3(unsigned int *sa, SA_flag_string_type* SA_flag, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length, unsigned int need_step)
{

	unsigned int t;
	for (t = sp; t<ep; t++)
	{
		///这里SA要改
		if (get_sa_restrict_steps_more_than_3(SA_flag, t, bwt, high_occ_table, sa,
			need_step, (&(locates[(*tmp_SA_length)]))) == 1)
		{
			locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] + occ;
			(*tmp_SA_length)++;
		}

	}

}



void accesss_SA_cur(unsigned int *sa, SA_flag_string_type* SA_flag, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length, unsigned int need_step)
{

	unsigned int t;
	for (t = sp; t<ep; t++)
	{

		///这里SA要改
		if (get_sa_restrict_steps(SA_flag, t, bwt, high_occ_table, sa,
			need_step, (&(locates[(*tmp_SA_length)]))) == 1)
		{

			locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] + occ;
			(*tmp_SA_length)++;

		}

	}

}


int unsigned_int_compareEntrySize(const void *a, const void *b) {
	unsigned int a_list = *(unsigned int *)a;
	unsigned int b_list = *(unsigned int *)b;
	if (a_list>b_list)
		return 1;
	else if (a_list<b_list)
		return -1;
	else
		return 0;

}




void search_from_bwt_more_than_3(unsigned int *sa, SA_flag_string_type* SA_flag, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	int na, int nc, int ng, int nt, int num_reads, FILE *f1)
{
	long long i, j;
	FILE *fout;
	int length_read;
	char* reads;
	int delta;
	unsigned int top, bot, t;
	unsigned int pre_top, pre_bot;
	unsigned int top_A, top_C, top_G, top_T;
	unsigned int bot_A, bot_C, bot_G, bot_T;

	ctoi['A'] = 0;
	ctoi['C'] = 1;
	ctoi['G'] = 2;
	ctoi['T'] = 3;



	FILE* _ih_fp = fopen("patterns.txt", "r");

	unsigned int numocc;
	unsigned int query_length = 0;

	fread(&(query_length), sizeof(unsigned int), 1, _ih_fp);
	fread(&(numocc), sizeof(unsigned int), 1, _ih_fp);

	fprintf(stdout, "query_length=%u\n", query_length);
	fprintf(stdout, "numocc=%u\n", numocc);


	char* patterns = (char *)malloc((query_length*numocc)*(sizeof(char)));

	fread(patterns, sizeof(char), (query_length*numocc), _ih_fp);


	reads = patterns;

	reads = patterns;

	length_read = query_length;

	num_reads = numocc;

	unsigned int* locates;

	double start = clock();

	double interval_start = clock();
	double interval_time = 0;

	double locate_start = clock();
	double locate_time = 0;

	long long number_of_locations = 0;

	unsigned int SA_start;
	unsigned int SA_length;
	unsigned int total_top;
	unsigned int total_bot;

	unsigned int tmp_SA_length = 0;
	unsigned int ijkijk = 0;

	unsigned int tree_nodes_number = (pow(4, compress_sa) - 1) / 3;

	fprintf(stdout, "tree_nodes_number=%llu\n", tree_nodes_number);

	unsigned int* sp_tree = (unsigned int*)malloc(sizeof(unsigned int)* tree_nodes_number);
	unsigned int* ep_tree = (unsigned int*)malloc(sizeof(unsigned int)* tree_nodes_number);
	unsigned int tree_index = 0;
	unsigned int tree_layer_length;
	unsigned int tree_layers;
	unsigned int tree_layers_i;

	unsigned int rare_interval = 0;

	unsigned int father_sp;
	unsigned int father_ep;

	unsigned int cut_thr = 1;
	unsigned int need_step = 0;



	double SA_flag_interval_length = 0;
	double SA_flag_interval_length_optimal = 0;

	long long number_of_hits = 0;

	unsigned int occ;

	unsigned int perform_interval = 0;


	struct  timeval  start_timeval;
	struct  timeval  end_timeval;
	unsigned long timer;
	gettimeofday(&start_timeval, NULL);

	long long total_debug_number = 0;



	for (i = 1; i <= num_reads; i++)
	{

		if (i % 10000 == 0)
		{
			fprintf(stdout, "i=%llu\n", i);
		}





		j = length_read - 1;
		delta = ctoi[reads[j]];






		top = nacgt[0][delta];
		bot = nacgt[0][delta + 1];


		///is_exist = 1;


		for (j=j-1; j >= 1; j--)
		{
			
			delta = ctoi[reads[j]];
			top = find_occ(top, delta, bwt,high_occ_table);
			bot = find_occ(bot, delta, bwt, high_occ_table);


			if (bot <= top)
			{
				break;
			}

		}



		pre_top = top;
		pre_bot = bot;


		for (; j >= 0; j--)
		{
			if (bot <= top)
			{
				break;
			}


			delta = ctoi[reads[j]];
			top = find_occ(top, delta, bwt, high_occ_table);
			bot = find_occ(bot, delta, bwt, high_occ_table);




		}


		if (bot <= top)
		{
			reads = reads + length_read;
			continue;
		}






		number_of_hits = bot - top;


		///这个地方要删了
		///total_debug_number = total_debug_number + number_of_hits;



		


		locates = (unsigned int *)malloc((bot - top)*sizeof(unsigned int));


		total_top = top;
		total_bot = bot;
		tmp_SA_length = 0;
		ijkijk = 0;


		tree_index = 0;

		sp_tree[tree_index] = top;
		ep_tree[tree_index] = bot;
		///SA_interval_tree_occ[tree_index] = 0;

		tree_layers = 0;

		need_step = compress_sa - tree_layers - 1;

		occ = 0;

		if (ep_tree[tree_index] - sp_tree[tree_index] <= cut_thr)
		{
			///这个函数内部要改
			///这里SA要改
			accesss_SA_cur_more_than_3(sa, SA_flag, bwt, high_occ_table, 
				locates,
				sp_tree[tree_index], ep_tree[tree_index], occ,
				&tmp_SA_length, need_step);

			sp_tree[tree_index] = ep_tree[tree_index] = (unsigned int)-1;
		}
		else
		{
		    ///这里不要改
			///这里SA要改
			accesss_SA_more_than_3(sa, SA_flag, locates,
				sp_tree[tree_index], ep_tree[tree_index], occ,
				&tmp_SA_length);
		}



		
		if (length_read>=2&&
			tmp_SA_length != number_of_hits)
		{

			delta = ctoi[reads[0]];

			///这里不要改
			///这里SA要改
			accesss_pre_SA(sa, SA_flag, locates,
				pre_top, pre_bot,
				&tmp_SA_length,delta);
			

		}
		





		tree_layer_length = 1;







		tree_layer_length = 1;

		tree_index = 1;

		//for (tree_layers = 1; tree_layers < compress_sa; tree_layers++)
		for (tree_layers = 1; tree_layers < compress_sa - 1; tree_layers++)
		{


			if (tmp_SA_length == number_of_hits)
			{
				break;
			}


			//need_step = compress_sa - tree_layers - 1;
			need_step = compress_sa - tree_layers - 2;

			occ = tree_layers;


			for (tree_layers_i = 0; tree_layers_i < tree_layer_length; tree_layers_i++)
			{

				father_sp = sp_tree[(tree_index - 1) / 4];
				father_ep = ep_tree[(tree_index - 1) / 4];





				if (father_sp == (unsigned int)-1)
				{
					sp_tree[tree_index] = sp_tree[tree_index + 1]
						= sp_tree[tree_index + 2] = sp_tree[tree_index + 3]
						= (unsigned int)-1;
				}
				else
				{
					///这里要改
					find_occ_all_sp_ep_optimal(father_sp, father_ep, bwt, high_occ_table,
						&sp_tree[tree_index], &sp_tree[tree_index + 1],
						&sp_tree[tree_index + 2], &sp_tree[tree_index + 3],
						&ep_tree[tree_index], &ep_tree[tree_index + 1],
						&ep_tree[tree_index + 2], &ep_tree[tree_index + 3]);





					if (ep_tree[tree_index] - sp_tree[tree_index] <= cut_thr)
					{
						///这里要改
						accesss_SA_cur_more_than_3(sa, SA_flag, bwt, high_occ_table,
							locates,
							sp_tree[tree_index], ep_tree[tree_index], occ,
							&tmp_SA_length, need_step);

						sp_tree[tree_index] = ep_tree[tree_index] = (unsigned int)-1;
					}
					else
					{
						///这里不要改
						accesss_SA_more_than_3(sa, SA_flag, locates,
							sp_tree[tree_index], ep_tree[tree_index], occ,
							&tmp_SA_length);
					}


					if (tmp_SA_length == number_of_hits)
					{

						break;
					}


					if (ep_tree[tree_index + 1] - sp_tree[tree_index + 1] <= cut_thr)
					{

						///这里要改
						accesss_SA_cur_more_than_3(sa, SA_flag, bwt, high_occ_table,
							locates,
							sp_tree[tree_index + 1], ep_tree[tree_index + 1], occ,
							&tmp_SA_length, need_step);

						sp_tree[tree_index + 1] = ep_tree[tree_index + 1] = (unsigned int)-1;
					}
					else
					{
						///这里不要改
						accesss_SA_more_than_3(sa, SA_flag, locates,
							sp_tree[tree_index + 1], ep_tree[tree_index + 1], occ,
							&tmp_SA_length);
					}


					if (tmp_SA_length == number_of_hits)
					{
						break;
					}


					if (ep_tree[tree_index + 2] - sp_tree[tree_index + 2] <= cut_thr)
					{

						///这里要改
						accesss_SA_cur_more_than_3(sa, SA_flag, bwt, high_occ_table,
							locates,
							sp_tree[tree_index + 2], ep_tree[tree_index + 2], occ,
							&tmp_SA_length, need_step);

						sp_tree[tree_index + 2] = ep_tree[tree_index + 2] = (unsigned int)-1;
					}
					else
					{
						///这里不要改
						accesss_SA_more_than_3(sa, SA_flag, locates,
							sp_tree[tree_index + 2], ep_tree[tree_index + 2], occ,
							&tmp_SA_length);
					}


					if (tmp_SA_length == number_of_hits)
					{
						break;
					}

					if (ep_tree[tree_index + 3] - sp_tree[tree_index + 3] <= cut_thr)
					{

						///这里要改
						accesss_SA_cur_more_than_3(sa, SA_flag, bwt, high_occ_table,
							locates,
							sp_tree[tree_index + 3], ep_tree[tree_index + 3], occ,
							&tmp_SA_length, need_step);

						sp_tree[tree_index + 3] = ep_tree[tree_index + 3] = (unsigned int)-1;
					}
					else
					{
						///这里不要改
						accesss_SA_more_than_3(sa, SA_flag, locates,
							sp_tree[tree_index + 3], ep_tree[tree_index + 3], occ,
							&tmp_SA_length);
					}


					if (tmp_SA_length == number_of_hits)
					{
						break;
					}

				}

				tree_index = tree_index + 4;
			}

			tree_layer_length = tree_layer_length * 4;
		}


		
	/**	
		qsort(locates, tmp_SA_length, sizeof(unsigned int), unsigned_int_compareEntrySize);

		for (t = 0; t<tmp_SA_length; t++)
		{

			fprintf(stderr, "i=%llu, site=%u\n", i, locates[t]);
			fflush(stderr);
		}
		
	**/	



		free(locates);



		number_of_locations = number_of_locations + tmp_SA_length;

		


		reads = reads + length_read;

	}

	gettimeofday(&end_timeval, NULL);
	timer = 1000000 * (end_timeval.tv_sec - start_timeval.tv_sec) + end_timeval.tv_usec - start_timeval.tv_usec;


	double finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;


	fprintf(stdout, "\n\n\n*********************Result*********************\n");


	fprintf(stdout, "searching Time: %ld microsecond\n", timer);

	fprintf(stdout, "searching Time: %5.3f seconds\n", duration);

	fprintf(stdout, "number of matched locations=%ld\n", number_of_locations);


	fprintf(stdout, "************************************************\n");

	///fprintf(stdout, "total_debug_number=%ld\n", total_debug_number);
	


}






void search_from_bwt(unsigned int *sa, SA_flag_string_type* SA_flag, bwt_string_type *bwt, high_occ_table_type* high_occ_table, int na, int nc, int ng, int nt, int num_reads, FILE *f1)
{
	long long i, j;
	FILE *fout;
	int length_read;
	char* reads;
	int delta;
	unsigned int top, bot, t;
	unsigned int top_A, top_C, top_G, top_T;
	unsigned int bot_A, bot_C, bot_G, bot_T;

	ctoi['A'] = 0;
	ctoi['C'] = 1;
	ctoi['G'] = 2;
	ctoi['T'] = 3;



	FILE* _ih_fp = fopen("patterns.txt", "r");

	unsigned int numocc;
	unsigned int query_length = 0;

	fread(&(query_length), sizeof(unsigned int), 1, _ih_fp);
	fread(&(numocc), sizeof(unsigned int), 1, _ih_fp);

	fprintf(stdout, "query_length=%u\n", query_length);
	fprintf(stdout, "numocc=%u\n", numocc);


	char* patterns = (char *)malloc((query_length*numocc)*(sizeof(char)));

	fread(patterns, sizeof(char), (query_length*numocc), _ih_fp);


	reads = patterns;

	reads = patterns;

	length_read = query_length;

	num_reads = numocc;

	unsigned int* locates;

	double start = clock();

	double interval_start = clock();
	double interval_time = 0;

	double locate_start = clock();
	double locate_time = 0;

	long long number_of_locations = 0;

	unsigned int SA_start;
	unsigned int SA_length;
	unsigned int total_top;
	unsigned int total_bot;

	unsigned int tmp_SA_length = 0;
	unsigned int ijkijk = 0;

	unsigned int tree_nodes_number = (pow(4, compress_sa)-1)/3;

	fprintf(stdout, "tree_nodes_number=%llu\n", tree_nodes_number);

	unsigned int* sp_tree = (unsigned int*)malloc(sizeof(unsigned int)* tree_nodes_number);
	unsigned int* ep_tree = (unsigned int*)malloc(sizeof(unsigned int)* tree_nodes_number);
	unsigned int tree_index = 0;
	unsigned int tree_layer_length;
	unsigned int tree_layers;
	unsigned int tree_layers_i;

	unsigned int rare_interval = 0;

	unsigned int father_sp;
	unsigned int father_ep;

	unsigned int cut_thr = 1;
	unsigned int need_step = 0;


	
	double SA_flag_interval_length = 0;
	double SA_flag_interval_length_optimal = 0;

	long long number_of_hits = 0;

	unsigned int occ;
	
	unsigned int perform_interval = 0;


	struct  timeval  start_timeval;
	struct  timeval  end_timeval;
	unsigned long timer;
	gettimeofday(&start_timeval, NULL);
	



	for (i = 1; i <= num_reads; i++)
	{
		
		if (i%10000==0)
		{
			fprintf(stdout, "i=%llu\n", i);
		}
		




		j = length_read - 1;
		delta = ctoi[reads[j]];


		

		top = nacgt[0][delta];
		bot = nacgt[0][delta + 1];


		j--;


		for (; j >= 0; j--)
		{
			delta = ctoi[reads[j]];
			top = find_occ(top, delta, bwt, high_occ_table);
			bot = find_occ(bot, delta, bwt, high_occ_table);
		}



		number_of_hits = bot - top;
		
		locates = (unsigned int *)malloc((bot - top)*sizeof(unsigned int));

		
		total_top = top;
		total_bot = bot;
		tmp_SA_length = 0;
		ijkijk = 0;


		tree_index = 0;

		sp_tree[tree_index] = top;
		ep_tree[tree_index] = bot;
		///SA_interval_tree_occ[tree_index] = 0;

		tree_layers = 0;

		need_step = compress_sa - tree_layers-1;

		occ = 0;

		if (ep_tree[tree_index] - sp_tree[tree_index] <= cut_thr)
		{
			///这里SA要改
			accesss_SA_cur(sa, SA_flag, bwt, high_occ_table,
				locates,
				sp_tree[tree_index], ep_tree[tree_index], occ,
				&tmp_SA_length, need_step);

			sp_tree[tree_index] = ep_tree[tree_index] = (unsigned int)-1;
		}
		else
		{
			///这里SA要改
			accesss_SA(sa, SA_flag, locates,
				sp_tree[tree_index], ep_tree[tree_index], occ,
				&tmp_SA_length);
		}


		


		tree_layer_length = 1;



			



		tree_layer_length = 1;

		tree_index = 1;

		for (tree_layers = 1; tree_layers < compress_sa; tree_layers++)
		{


			if (tmp_SA_length == number_of_hits)
			{
				break;
			}


			need_step = compress_sa - tree_layers - 1;

			occ = tree_layers;


			for (tree_layers_i = 0; tree_layers_i < tree_layer_length; tree_layers_i++)
			{

				father_sp = sp_tree[(tree_index - 1) / 4];
				father_ep = ep_tree[(tree_index - 1) / 4];





				if (father_sp==(unsigned int)-1)
				{
					sp_tree[tree_index] = sp_tree[tree_index + 1]
						= sp_tree[tree_index + 2] = sp_tree[tree_index + 3]
						= (unsigned int)-1;
				}
				else
				{
					find_occ_all_sp_ep_optimal(father_sp, father_ep, bwt, high_occ_table,
						&sp_tree[tree_index], &sp_tree[tree_index + 1],
						&sp_tree[tree_index + 2], &sp_tree[tree_index + 3],
						&ep_tree[tree_index], &ep_tree[tree_index + 1],
						&ep_tree[tree_index + 2], &ep_tree[tree_index + 3]);





					if (ep_tree[tree_index] - sp_tree[tree_index] <= cut_thr)
					{
						accesss_SA_cur(sa, SA_flag, bwt, high_occ_table,
							locates,
							sp_tree[tree_index], ep_tree[tree_index], occ,
							&tmp_SA_length, need_step);

						sp_tree[tree_index] = ep_tree[tree_index] = (unsigned int)-1;
					}
					else
					{
						accesss_SA(sa, SA_flag, locates,
							sp_tree[tree_index], ep_tree[tree_index], occ,
							&tmp_SA_length);
					}


					if (tmp_SA_length == number_of_hits)
					{

						break;
					}


					if (ep_tree[tree_index+1] - sp_tree[tree_index+1] <= cut_thr)
					{
						accesss_SA_cur(sa, SA_flag, bwt, high_occ_table,
							locates,
							sp_tree[tree_index+1], ep_tree[tree_index+1], occ,
							&tmp_SA_length, need_step);

						sp_tree[tree_index+1] = ep_tree[tree_index+1] = (unsigned int)-1;
					}
					else
					{
						accesss_SA(sa, SA_flag, locates,
							sp_tree[tree_index+1], ep_tree[tree_index+1], occ,
							&tmp_SA_length);
					}


					if (tmp_SA_length == number_of_hits)
					{
						break;
					}


					if (ep_tree[tree_index+2] - sp_tree[tree_index+2] <= cut_thr)
					{
						accesss_SA_cur(sa, SA_flag, bwt, high_occ_table,
							locates,
							sp_tree[tree_index+2], ep_tree[tree_index+2], occ,
							&tmp_SA_length, need_step);

						sp_tree[tree_index+2] = ep_tree[tree_index+2] = (unsigned int)-1;
					}
					else
					{
						accesss_SA(sa, SA_flag, locates,
							sp_tree[tree_index+2], ep_tree[tree_index+2], occ,
							&tmp_SA_length);
					}


					if (tmp_SA_length == number_of_hits)
					{
						break;
					}

					if (ep_tree[tree_index+3] - sp_tree[tree_index+3] <= cut_thr)
					{
						accesss_SA_cur(sa, SA_flag, bwt, high_occ_table,
							locates,
							sp_tree[tree_index+3], ep_tree[tree_index+3], occ,
							&tmp_SA_length, need_step);

						sp_tree[tree_index+3] = ep_tree[tree_index+3] = (unsigned int)-1;
					}
					else
					{
						accesss_SA(sa, SA_flag, locates,
							sp_tree[tree_index+3], ep_tree[tree_index+3], occ,
							&tmp_SA_length);
					}


					if (tmp_SA_length == number_of_hits)
					{
						break;
					}

				}

				tree_index = tree_index + 4;
			}

			tree_layer_length = tree_layer_length * 4;
		}

		/**
		qsort(locates, tmp_SA_length, sizeof(unsigned int), unsigned_int_compareEntrySize);

		for (t = 0; t<tmp_SA_length; t++)
		{

			fprintf(stderr, "i=%llu, site=%u\n", i, locates[t]);
			fflush(stderr);
		}
		**/
		
		
		

		
		free(locates);

		

		number_of_locations = number_of_locations + tmp_SA_length;


		

		reads = reads + length_read;

	}

	gettimeofday(&end_timeval, NULL);
	timer = 1000000 * (end_timeval.tv_sec - start_timeval.tv_sec) + end_timeval.tv_usec - start_timeval.tv_usec;


	double finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;


	fprintf(stdout, "\n\n\n*********************Result*********************\n");
	

	fprintf(stdout, "searching Time: %ld microsecond\n", timer);

	fprintf(stdout, "searching Time: %5.3f seconds\n", duration);

	fprintf(stdout, "number of matched locations=%ld\n", number_of_locations);


	fprintf(stdout, "************************************************\n");
	



}






int search_bwt()
{
	char filename[100], filename1[100], filenames[100], filenameo[100], filenameb[100];
	bwt_string_type *bwt;
	unsigned int bwt_length;
	unsigned int *sa;
	unsigned int *occ_whole;
	unsigned int tt;
	int num_reads;
	char ch;
	long long i, j, t;
	u_int32_t x;
	u_int8_t xx;
	struct timeval tv_begin, tv_end;

	FILE *f1, *f2, *fs, *fb, *fout, *fo;


	pop_count_mode[0] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)*8 / 2; i++)
	{
		pop_count_mode[0] = (bwt_string_type)(pop_count_mode[0] << 2) | (bwt_string_type)3;
	}


	pop_count_mode[1] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		pop_count_mode[1] = (bwt_string_type)(pop_count_mode[1] << 2) | (bwt_string_type)2;
	}


	pop_count_mode[2] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		pop_count_mode[2] = (bwt_string_type)(pop_count_mode[2] << 2) | (bwt_string_type)1;
	}

	pop_count_mode[3] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		pop_count_mode[3] = (bwt_string_type)(pop_count_mode[3] << 2) | (bwt_string_type)0;
	}
	/**
	dec_bit(pop_count_mode[0]);
	dec_bit(pop_count_mode[1]);
	dec_bit(pop_count_mode[2]);
	dec_bit(pop_count_mode[3]);
	**/

	mode_high=(bwt_string_type)0;
	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		mode_high = (bwt_string_type)(mode_high << 2) | (bwt_string_type)2;
	}


	mode_low = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		mode_low = (bwt_string_type)(mode_low << 2) | (bwt_string_type)1;
	}

	/**
	dec_bit(mode_high);
	dec_bit(mode_low);
	**/


	fprintf(stdout, "Please input the prefix of index name:\n");
	scanf("%s", filename);
	strcpy(filename + strlen(filename), ".index");
	strcpy(filenames, filename);
	strcpy(filenameo, filename);
	strcpy(filenameb, filename);
	strcpy(filenames + strlen(filename), ".sa");
	strcpy(filenameo + strlen(filename), ".occ");
	strcpy(filenameb + strlen(filename), ".bwt");





	fprintf(stdout, "\n\n\n*********************Warning*********************\n\n");


	fprintf(stdout, "All patterns in 'patterns.txt' cannot include the characters which do not belong {a,c, g, t, A, C, G, T}. \n");
	fprintf(stdout, "If some patterns consist of such characters, the results of FMtree will be incorrect. \n");


	fprintf(stdout, "\n*********************Warning*********************\n\n\n");


	///这里要改
	fo = fopen(filenameo, "r");
	if (fo == NULL)
	{
		fprintf(stdout, "Failed to open %s! FMtree will exit ...\n", filenameo);
		return 1;
	}

	fs = fopen(filenames, "r");
	if (fs==NULL)
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


	fread(&SA_length, sizeof(SA_length), 1, f2);
	fread(&shapline, sizeof(shapline), 1, f2);

	fprintf(stdout, "shapline=%llu\n", shapline);




	for (i = 0; i<bwt_step; i++)
	{
		for (j = 0; j <= (1 << (i + i + 2)); j++)
			fread(&nacgt[i][j], sizeof(nacgt[i][j]), 1, f2);
	}

	fread(&compress_sa, sizeof(unsigned int), 1, f2);
	fread(&compress_occ, sizeof(unsigned int), 1, f2);
	///这里要改
	fread(&high_compress_occ, sizeof(unsigned int), 1, f2);

	

	fprintf(stdout, "compress_sa=%llu\n", compress_sa);
	fprintf(stdout, "compress_occ=%llu\n", compress_occ);
	fprintf(stdout, "high_compress_occ=%llu\n", high_compress_occ);



	na = nacgt[0][0];
	nc = nacgt[0][1];
	ng = nacgt[0][2];
	nt = nacgt[0][3];




	mode_4[0] = (bwt_string_type)-1;
	mode_4[1] = (bwt_string_type)-1;
	mode_4[2] = (bwt_string_type)-1;
	mode_4[3] = (bwt_string_type)0;


	fread(&bwt_length, sizeof(bwt_length), 1, fb);
	bwt = (bwt_string_type *)malloc(sizeof(bwt_string_type)*bwt_length);
	fread(bwt, sizeof(bwt_string_type), bwt_length, fb);
	printf("BWT has been loaded!\n");


	printf("SA_length=%u\n", SA_length);
	unsigned int sparse_suffix_array_length = 0;
	fread(&sparse_suffix_array_length, sizeof(sparse_suffix_array_length), 1, fs);
	printf("sparse_suffix_array_length=%u\n", sparse_suffix_array_length);
	sa = (unsigned int *)malloc(sizeof(unsigned int)*(sparse_suffix_array_length));
	fread(sa, sizeof(unsigned int), sparse_suffix_array_length, fs);


	unsigned int SA_flag_iterater = 0;

	fread(&SA_flag_iterater, sizeof(SA_flag_iterater), 1, fs);

	fprintf(stdout, "SA_flag_iterater=%llu\n", SA_flag_iterater);
	
	SA_flag_string_type*  SA_flag=NULL;

	SA_flag = (SA_flag_string_type*)malloc(sizeof(SA_flag_string_type)*SA_flag_iterater);

	fread(SA_flag, sizeof(SA_flag_string_type), SA_flag_iterater, fs);



	unsigned int high_occ_table_length = 0;
	fread(&high_occ_table_length, sizeof(high_occ_table_length), 1, fo);
	fprintf(stdout, "high_occ_table_length=%llu\n", high_occ_table_length);
	high_occ_table_type* high_occ_table;
	high_occ_table = (high_occ_table_type*)malloc(sizeof(high_occ_table_type)*high_occ_table_length);
	fread(high_occ_table, sizeof(high_occ_table_type), high_occ_table_length, fo);

	


	fclose(f2);
	fclose(fs);
	fclose(fb);
	fclose(fo);


	if (compress_sa >= 4)
	{
		search_from_bwt_more_than_3(sa, SA_flag, bwt, high_occ_table, na, nc, ng, nt, num_reads, f1);
	}
	else
	{
		search_from_bwt(sa, SA_flag, bwt, high_occ_table, na, nc, ng, nt, num_reads, f1);
	}

	

	return 0;
}









int main()
{
	int mode;
	fprintf(stdout, "input mode:\n1 for index\n2 for search\n3 make patterns:");
	scanf("%d", &mode);


	if (mode == 1) build_index();
	else if (mode == 2) search_bwt();
	else if (mode == 3) Generate_patterns();
	else
	{
		fprintf(stdout, "Do not input required options! FMtree will exit ...\n");
	}

	return 0;
}
