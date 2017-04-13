
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







typedef uint32_t SA_flag_string_type;



unsigned int compress_sa = 4, compress_occ = 448, compress_SA_flag=224;;
typedef uint64_t bwt_string_type;
///typedef unsigned int bwt_string_type;
unsigned int bwt_warp_number = sizeof(bwt_string_type)* 8 / 2;
unsigned int single_occ_bwt = sizeof(bwt_string_type)/sizeof(unsigned);
unsigned int occ_words = 4 / single_occ_bwt;
unsigned int acctuall_bwt_gap = compress_occ + sizeof(unsigned int)* 4 * 8 / 2;
bwt_string_type mode_4[4];
bwt_string_type mode = (bwt_string_type)-1;
bwt_string_type mode_16 = (bwt_string_type)65535;
bwt_string_type mode_32 = ((bwt_string_type)-1)>>32;
bwt_string_type mode_high;
bwt_string_type mode_low;

bwt_string_type pop_count_mode[4];
unsigned int SA_flag_warp_number = sizeof(SA_flag_string_type)* 8;
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


		fprintf(stdout, "%u", bit_2[i]);
	}

	fprintf(stdout, "\n");

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

	unsigned int query_number = 0;
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



		tmp_loc = random();

		if (tmp_loc + query_length <= n)
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

		wd[wv[i]]++;
	}
	for (i = 1; i<m; i++) wd[i] += wd[i - 1];
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
	fseek(_ih_fp, 0, SEEK_SET);

	fseek(_ih_fp, 0, SEEK_END);
	unsigned int n = ftell(_ih_fp);

	n = cc - 1;

	///fread(&(suffix_array.n), sizeof(unsigned int), 1, _ih_fp);

	fprintf(stdout, "r_n=%u\n", n);

	fseek(_ih_fp, 0, SEEK_SET);



	n++; // append the virtual sentinel
	fprintf(stdout, "Allocating input and output space: %u bytes = %.2lf MB", 5 * n, (double)5 * n / 1024 / 1024);
	unsigned char *s_ch = new unsigned char[n];
	//unsigned char *s_ch;
	unsigned int *SA = new unsigned int[n];
	//unsigned int *SA;
	SA = (unsigned int*)malloc(sizeof(unsigned int)*n + 10);
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
	for (i = 0; i<n; i++) s_ch[i] = refer[i];
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
	printf("sa_n=%d\n", n);
	*sa = SA;

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
	scanf("%llu", &compress_sa);



	if (compress_sa >= 9 || compress_sa <= 1)
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
		if ((ch == 'A') || (ch == 'C') || (ch == 'G') || (ch == 'T') ||
			(ch == 'a') || (ch == 'c') || (ch == 'g') || (ch == 't'))
			text_length++;
		else
		{

			fprintf(stdout, "%u-th character does not belong to {a, c, g, t, A, C, G, T}. Its  ASCII Code is %u, and it is %c.\n",
				number_of_total_characters, ch, ch);
		}

		number_of_total_characters++;

	}


	fclose(f1);
	

	f1 = fopen(filename, "r");

	unsigned int occ_line_number = 0;
	occ_line_number = (text_length+1) / compress_occ + 1;



	
	refer = (char *)malloc(sizeof(char)*(text_length + 1));
	bwt_warp_number = sizeof(bwt_string_type) * 8 / 2;






	unsigned int bwt_length = (text_length + 1) / bwt_warp_number + 1;





	unsigned int occ_byte_length = (occ_line_number * 4 * sizeof(unsigned int)) 
		/ (sizeof(bwt_string_type))+1;






	bwt = (bwt_string_type *)malloc(sizeof(bwt_string_type)*(bwt_length + occ_byte_length+1));






	for (i = 0; i < bwt_length + occ_byte_length+1; i++)
	{
		bwt[i] = (bwt_string_type)0;
	}


	fprintf(stdout, "bwt_length=%llu,occ_byte_length=%llu, bwt_length+occ_byte_length=%llu\n",
		bwt_length, occ_byte_length, bwt_length + occ_byte_length);

	SA_length = 1;
	while (fscanf(f1, "%c", &ch)!=EOF)
	{
		
		if ((ch == 'A') || (ch == 'C') || (ch == 'G') || (ch == 'T') ||
			(ch == 'a') || (ch == 'c') || (ch == 'g') || (ch == 't'))
		{

			if (SA_length == 38627 || SA_length == 38628 || SA_length == 38629)
			{
				fprintf(stdout, "SA_length=%llu, ch=%c\n", SA_length, ch);
			}

			refer[SA_length++] = ctoi[ch] + 1;

			if (SA_length == 38627 || SA_length == 38628 || SA_length == 38629)
			{
				fprintf(stdout, "refer[38628]=%llu, refer[38627]=%llu,refer[38629]=%llu\n", 
					refer[38628], refer[38627], refer[38629]);
			}

		}
	}


	printf("SA_length=%llu\n", SA_length);



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

	printf("SA_length=%llu\n", SA_length);




	FILE* suffix = fopen("suffix.txt", "w");

	get_sa_fromFILE(&sa, SA_length, f1, &refer[1]);


	fwrite(&SA_length, sizeof(unsigned int), 1, suffix);
	fwrite(sa, sizeof(unsigned int), SA_length, suffix);

	fclose(suffix);






	



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
	unsigned int occ_bwt_number = (sizeof(unsigned int)*4)/sizeof(bwt_string_type);
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


	long long A_number[4] = {0};



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

				shift_length = (single_bwt_occ - j % single_bwt_occ - 1) * sizeof(unsigned int) * 8;



				tmp_bwt = (bwt_string_type)nacgt[0][j + 1];



				tmp_bwt = tmp_bwt << shift_length;




				bwt[bwt_iterater + j / single_bwt_occ] =
					bwt[bwt_iterater + j / single_bwt_occ] | tmp_bwt;


			}
			
			bwt_iterater = bwt_iterater + single_bwt_occ;



		}


		tmp_bwt = (bwt_string_type)0;
		

	}


	printf("Occ has been writed!\n");

	fprintf(stdout, "bwt_iterater=%u\n", bwt_iterater);


	fprintf(stdout, "A_number[0]=%u\n", A_number[0]);
	fprintf(stdout, "A_number[1]=%u\n", A_number[1]);
	fprintf(stdout, "A_number[2]=%u\n", A_number[2]);
	fprintf(stdout, "A_number[3]=%u\n", A_number[3]);

	fprintf(stdout, "nacgt[0][0]=%u\n", nacgt[0][0]);
	fprintf(stdout, "nacgt[0][1]=%u\n", nacgt[0][1]);
	fprintf(stdout, "nacgt[0][2]=%u\n", nacgt[0][2]);
	fprintf(stdout, "nacgt[0][3]=%u\n", nacgt[0][3]);
	fprintf(stdout, "nacgt[0][4]=%u\n", nacgt[0][4]);




	if (i%bwt_warp_number != 0 &&
		i % compress_occ != 0)
	{
		bwt_iterater++;
	}

	fwrite(&bwt_iterater, sizeof(unsigned int), 1, fb);
	fwrite(bwt, sizeof(bwt_string_type), bwt_iterater, fb);
	printf("BWT has been writed!\n");

	fwrite(&SA_length, sizeof(SA_length), 1, f2);
	fwrite(&shapline, sizeof(shapline), 1, f2);

	printf("cc&sharp_line has been writed!\n");


	fprintf(stdout, "write SA_length=%u, shapline=%u\n", SA_length, shapline);
	

	for (i = 0; i<bwt_step; i++)
	{
		nacgt[i][0] = 1;
		fwrite(&nacgt[i][0], sizeof(unsigned int), 1, f2);


		


		for (j = 1; j <= (1 << (i + i + 2)); j++)
		{
			nacgt[i][j] = nacgt[i][j] + nacgt[i][j - 1];
			fwrite(&nacgt[i][j], sizeof(unsigned int), 1, f2);
			fprintf(stdout, "nacgt[%u][%u]=%u\n", i,j, nacgt[i][j]);
			fflush(stdout);
		}
	}



	unsigned int SA_number = SA_length / compress_sa + ((SA_length % compress_sa) != 0);
	fwrite(&SA_number, sizeof(unsigned int), 1, fs);
	fprintf(stdout, "SA_number=%u\n", SA_number);
	fflush(stdout);
	for (i = 0; i<SA_length; i = i + compress_sa)
	{
		fwrite(&sa[i], sizeof(sa[i]), 1, fs);
	}



	



	fwrite(&compress_sa, sizeof(unsigned int), 1, f2);
	fwrite(&compress_occ, sizeof(unsigned int), 1, f2);




	






	fclose(f1);

	


	fclose(f2);

	fprintf(stdout, "Sucess!\n");
	fflush(stdout);

	fclose(fs);
	fclose(fb);

	

	return 0;
}

static inline int __occ_auxx(u_int64_t y, int c)
{
	y = ((c & 2) ? y : ~y) >> 1 & ((c & 1) ? y : ~y) & 0x5555555555555555ull;
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

static inline int __occ_aux(u_int64_t y, int c)
{
	y = y ^ cnt_aux[c];
	y = (y >> 1) & y & 0x5555555555555555ull;
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



	for (i = (line / compress_occ)*compress_occ; i + 16<line; i = i + 16)
	{
		b = bwt[i / 16];
		ans += (cnt_table[b & 0xffff][delta] + cnt_table[b >> 16 & 0xffff][delta]) & 0xffff;
	}



	b = bwt[i / 16] & ~((1ull << ((~(line - 1) & 15) << 1)) - 1);
	ans += (cnt_table[b & 0xffff][delta] + cnt_table[b >> 16 & 0xffff][delta]) & 0xffff;
	if (delta == 0) ans -= ~(line - 1) & 15;
	if ((delta == 1) && (shapline >= (line / compress_occ)*compress_occ) && (shapline<line)) ans--;


	return ans;
}


unsigned int find_occ(unsigned int line, int delta, bwt_string_type *bwt)
{


	unsigned int i, j;
	unsigned int b;

	bwt_string_type ans = nacgt[0][delta];

	unsigned int actually_line = ((line / compress_occ)*acctuall_bwt_gap) / bwt_warp_number;



	ans += 
		(bwt[actually_line + delta / single_occ_bwt] 
		>> ((single_occ_bwt - 1 - delta%single_occ_bwt) * 32))&mode_32;



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
		P_A = P_A & mode_low;
		P_B = P_B & mode_low;
		P_A = P_A & P_B;

		ans = ans + __builtin_popcountll(P_A);


	}



	








	if ((delta == 1) && (shapline >= (line / compress_occ)*compress_occ) && (shapline<line)) 
		ans--;

	return ans;
}


unsigned int get_sa(unsigned int line, bwt_string_type *bwt, unsigned int *sa)
{
	unsigned int l = line;
	bwt_string_type delta;
	unsigned int i = 0;

	unsigned int actually_line;

	if (line == shapline) return 0;
	while ((l%compress_sa) != 0)
	{

		actually_line = ((l / compress_occ)*acctuall_bwt_gap) + l % compress_occ + 64;

		delta = (bwt[actually_line / bwt_warp_number] 
			>> ((bwt_warp_number - actually_line%bwt_warp_number - 1) * 2))
			&(bwt_string_type)3;


		l = find_occ(l, delta, bwt);


		i++;
		if (l == shapline) return i;
	}

	return sa[l / compress_sa] + i;
}

void search_from_bwt(unsigned int *sa, bwt_string_type *bwt, int na, int nc, int ng, int nt, int num_reads, FILE *f1)
{
	long long i, j;
	FILE *fout;
	int length_read;
	char* reads;
	int delta;
	unsigned int top, bot, t;

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
			top = find_occ(top, delta, bwt);
			bot = find_occ(bot, delta, bwt);			
		}




		locates = (unsigned int *)malloc((bot - top)*sizeof(unsigned int));
		
		unsigned int ijkijk = 0;
		for (t = top; t<bot; t++, ijkijk++)
		{
			locates[ijkijk] = get_sa(t, bwt, sa);
				
			
		}


		/**
		qsort(locates, bot - top, sizeof(unsigned int), unsigned_int_compareEntrySize);

		for (t = 0; t<bot - top; t++)
		{

			fprintf(stderr, "i=%llu, site=%u\n", i, locates[t]);
			fflush(stderr);
		}
		**/

		

		
		free(locates);


		number_of_locations = number_of_locations + (bot - top);

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


	FILE *f1, *f2, *fs, *fb, *fout;
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





	///printf("input fa file name:\n");
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

	fprintf(stdout, "compress_sa=%llu\n", compress_sa);
	fprintf(stdout, "compress_occ=%llu\n", compress_occ);




	na = nacgt[0][0];
	nc = nacgt[0][1];
	ng = nacgt[0][2];
	nt = nacgt[0][3];


	unsigned int bwt_count_hash_length = pow(2, bwt_count_hash_table_bit);
	for (i = 0; i < bwt_count_hash_length; ++i)
	{
		for (j = 0; j < 4; ++j)
		{
			x = ((((i >> 8) & 3) == j) 
				+ (((i >> 10) & 3) == j) 
				+ (((i >> 12) & 3) == j) 
				+ (((i >> 14) & 3) == j) 
				+ ((i & 3) == j) 
				+ (((i >> 2) & 3) == j) 
				+ (((i >> 4) & 3) == j) 
				+ (((i >> 6) & 3) == j));
			cnt_table[i][j] = x;
		}
	}

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




















	
	fclose(f2);
	fclose(fs);
	fclose(fb);


	search_from_bwt(sa, bwt, na, nc, ng, nt, num_reads, f1);

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
