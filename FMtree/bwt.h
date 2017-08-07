
#ifndef __BWT__
#define __BWT__




#define SA_counter_length 32
typedef uint64_t bwt_string_type;
typedef uint64_t SA_flag_string_type;
typedef uint32_t high_occ_table_type;

typedef struct bwt_index
{
	unsigned int *sa;
	SA_flag_string_type* SA_flag;
	bwt_string_type *bwt;
	high_occ_table_type* high_occ_table;




	unsigned int SA_length;
	unsigned int bwt_length;
	unsigned int sparse_suffix_array_length;
	unsigned int SA_flag_iterater;
	unsigned int high_occ_table_length;


	unsigned int shapline;



	unsigned int compress_sa = 8;
	unsigned int compress_occ = 256;
	unsigned int high_compress_occ = 65536;
	unsigned int compress_SA_flag = 224;
	unsigned int bwt_warp_number = sizeof(bwt_string_type)* 8 / 2;
	unsigned int SA_flag_warp_number = sizeof(SA_flag_string_type)* 8;
	uint64_t tmp_SA_flag = (uint64_t)0;
	uint64_t* long_SA_flag = NULL;
	unsigned int single_occ_bwt = sizeof(bwt_string_type) / sizeof(unsigned short);
	unsigned int occ_words = 4 / single_occ_bwt;
	unsigned int acctuall_bwt_gap = compress_occ + sizeof(unsigned short)* 4 * 8 / 2;
	unsigned int acctuall_SA_flag_gap = compress_SA_flag + SA_counter_length;
	unsigned int SA_counter_shift_length = sizeof(SA_flag_string_type)* 8 - SA_counter_length;
	bwt_string_type mode_4[4];
	bwt_string_type mode = (bwt_string_type)-1;
	bwt_string_type mode_high_1 = (bwt_string_type)1 << (SA_flag_warp_number - 1);
	bwt_string_type mode_16 = (bwt_string_type)65535;
	bwt_string_type mode_32 = ((bwt_string_type)-1) >> 32;
	bwt_string_type mode_high;
	bwt_string_type mode_low;
	SA_flag_string_type mode_SA_flag = (SA_flag_string_type)(((SA_flag_string_type)-1) << (sizeof(SA_flag_string_type)* 8 - 1));
	SA_flag_string_type SA_pop_count_mode = (SA_flag_string_type)(((SA_flag_string_type)-1) >> SA_counter_length);
	bwt_string_type pop_count_mode[4];
	unsigned int bwt_count_hash_table_bit = 16;
	unsigned int text_length;
	unsigned int na, nc, ng, nt;
	unsigned int nacgt[5];
	unsigned int bwt_step = 1;
	unsigned int ctoi[256];
	char itoc[4] = { 'A', 'C', 'G', 'T' };

	unsigned int total_gap;
	unsigned int num_r = 10000, len_r = 12;
	unsigned int SA_header_mode = (unsigned int)(((unsigned int)-1) >> 2);
	unsigned int SA_header_mode_reverse = (unsigned int)(((unsigned int)-1) << 30);
	FILE* _rg_fp;


	int delta;



	unsigned int tree_nodes_number;
	unsigned int* sp_tree;
	unsigned int* ep_tree;
	unsigned int tree_index;
	unsigned int tree_layers;
	unsigned int need_step;
	unsigned int cut_thr=1;
	unsigned int tree_layer_length = 1;
	unsigned int* FMtree_queue;
	///char* FMtree_queue_int;

	unsigned int FMtree_queue_start_point = 0;
	unsigned int FMtree_queue_end_point = 0;




} bwt_index;


extern bwt_index bitmapper_index_params;




/**
unsigned int locate(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences);
**/


unsigned int locate(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences);

unsigned int locate_debug(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences);

unsigned int locate_less_than_4(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences);

unsigned int init_bitmapper_index_params();

unsigned int indenpendent_creadte_index(unsigned int text_length, char** input_refer, unsigned int compress_sa, char* filename);

/**
unsigned int count(char* pattern, unsigned int length, 
	unsigned int* sp, unsigned int* ep, unsigned int* sp1, unsigned int* ep1);
**/



///unsigned int find_occ_fm_index(unsigned int line, int delta, bwt_string_type *bwt, high_occ_table_type* high_occ_table);
/**
void bwt_accesss_SA_cur_more_than_3(unsigned int *sa, SA_flag_string_type* SA_flag, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length, unsigned int need_step);

**/

/**
unsigned int bwt_get_sa_restrict_steps_more_than_3
(SA_flag_string_type* SA_flag, unsigned int line,
bwt_string_type *bwt, high_occ_table_type* high_occ_table,
unsigned int *sa, unsigned int need_step, unsigned int* accessed_sa);
**/

/**
void bwt_direct_get_sa_interval_long
(SA_flag_string_type* SA_flag, unsigned int sp, unsigned int ep, unsigned int* start, unsigned int* length);
**/
/**
void bwt_accesss_SA_more_than_3(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length);
**/

/**
unsigned int bwt_get_sa_restrict_zero_steps_more_than_3
(SA_flag_string_type* SA_flag, unsigned int *sa, unsigned int sp, unsigned int* start);
**/


/**
void bwt_accesss_pre_SA(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep,
	unsigned int* tmp_SA_length, unsigned int delta);
**/

/**
unsigned int bwt_get_sa_restrict_zero_steps
(SA_flag_string_type* SA_flag, unsigned int *sa, unsigned int sp, unsigned int* start);
**/

/**
void bwt_find_occ_all_sp_ep_optimal(unsigned int sp, unsigned int ep, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* ans_A_sp, unsigned int* ans_C_sp, unsigned int* ans_G_sp, unsigned int* ans_T_sp,
	unsigned int* ans_A_ep, unsigned int* ans_C_ep, unsigned int* ans_G_ep, unsigned int* ans_T_ep);
**/

unsigned int load_index(char* filename_prefix);






inline unsigned int find_occ_fm_index(unsigned int line, int delta, bwt_string_type *bwt, high_occ_table_type* high_occ_table)
{

	unsigned int i, j;
	unsigned int b;

	bwt_string_type ans = bitmapper_index_params.nacgt[delta];

	///这里是增加的
	///ans = ans + high_occ_table[(line / bitmapper_index_params.high_compress_occ) * 4 + delta];
	ans = ans + high_occ_table[((line >> 16) << 2) + delta];

	/**
	unsigned int actually_line = ((line / bitmapper_index_params.compress_occ)*bitmapper_index_params.acctuall_bwt_gap)
		/ bitmapper_index_params.bwt_warp_number;

	ans += (bwt[actually_line] >> ((bitmapper_index_params.single_occ_bwt - 1 - delta) * 16))
		&bitmapper_index_params.mode_16;
	**/

	unsigned int actually_line = ((line >> 8)*bitmapper_index_params.acctuall_bwt_gap) >> 5;



	ans += (bwt[actually_line] >> ((bitmapper_index_params.single_occ_bwt - 1 - delta) << 4))&bitmapper_index_params.mode_16;





	/**
	if (line%bitmapper_index_params.compress_occ == 0)
		return ans;
	**/

	if ((line & 255) == 0)
		return ans;


	/**
	unsigned int need_line = line % bitmapper_index_params.compress_occ;

	actually_line = actually_line + 4 / bitmapper_index_params.single_occ_bwt;
	**/

	unsigned int need_line = (line & 255);

	actually_line = actually_line + 1;



	bwt_string_type tmp_bwt;

	bwt_string_type mode_low_2 = (bwt_string_type)3;


	bwt_string_type P_A, P_B;


	i = 0;

	while (i + bitmapper_index_params.bwt_warp_number <= need_line)
	{
		P_A = bwt[actually_line++];
		P_B = P_A^bitmapper_index_params.pop_count_mode[delta];
		P_A = P_B >> 1;
		P_A = P_A & P_B;
		P_A = P_A & bitmapper_index_params.mode_low;

		ans = ans + __builtin_popcountll(P_A);

		i = i + bitmapper_index_params.bwt_warp_number;

	}


	need_line = need_line - i;


	if (need_line != 0)
	{



		P_A = bwt[actually_line++];

		P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line) * 2);

		P_A = P_A & P_B;


		if (delta == 0)
		{
			P_B = ~P_B;

			P_A = P_A | P_B;
		}


		P_B = P_A^bitmapper_index_params.pop_count_mode[delta];
		P_A = P_B >> 1;
		P_A = P_A & P_B;
		P_A = P_A & bitmapper_index_params.mode_low;


		ans = ans + __builtin_popcountll(P_A);


	}

	/**
	if ((delta == 1) && (bitmapper_index_params.shapline >= (line / bitmapper_index_params.compress_occ)
		*bitmapper_index_params.compress_occ) && (bitmapper_index_params.shapline<line))
		ans--;
	**/

	if ((delta == 1) && (bitmapper_index_params.shapline >= ((line >> 8) << 8)) && (bitmapper_index_params.shapline<line))
		ans--;

	return ans;
}



inline unsigned int count(char* pattern, unsigned int length,
	unsigned int* sp, unsigned int* ep, unsigned int* sp1, unsigned int* ep1)
{

	long long j = length - 1;

	bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];


	if (bitmapper_index_params.delta > 3)
	{
		return 0;
	}

	unsigned int top = bitmapper_index_params.nacgt[bitmapper_index_params.delta];
	unsigned int bot = bitmapper_index_params.nacgt[bitmapper_index_params.delta + 1];


	for (j = j - 1; j >= 1; j--)
	{

		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];

		if (bitmapper_index_params.delta > 3)
		{
			return 0;
		}


		top = find_occ_fm_index(top, bitmapper_index_params.delta,
			bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table);
		bot = find_occ_fm_index(bot, bitmapper_index_params.delta,
			bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table);

		if (bot <= top)
		{
			break;
		}

	}

	(*sp1) = top;
	(*ep1) = bot;


	for (; j >= 0; j--)
	{

		if (bot <= top)
		{
			break;
		}


		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];


		if (bitmapper_index_params.delta > 3)
		{
			return 0;
		}



		top = find_occ_fm_index(top, bitmapper_index_params.delta, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table);
		bot = find_occ_fm_index(bot, bitmapper_index_params.delta, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table);


	}










	(*sp) = top;
	(*ep) = bot;


	if (bot <= top)
	{
		return 0;
	}
	else
	{


		return bot - top;
	}



}


inline unsigned int count_less_than_4(char* pattern, unsigned int length,
	unsigned int* sp, unsigned int* ep)
{


	long long j = length - 1;

	bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];

	if (bitmapper_index_params.delta > 3)
	{
		return 0;
	}

	unsigned int top = bitmapper_index_params.nacgt[bitmapper_index_params.delta];
	unsigned int bot = bitmapper_index_params.nacgt[bitmapper_index_params.delta + 1];



	j--;



	for (; j >= 0; j--)
	{
		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[j]];


		if (bitmapper_index_params.delta > 3)
		{
			return 0;
		}



		top = find_occ_fm_index(top, bitmapper_index_params.delta, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table);
		bot = find_occ_fm_index(bot, bitmapper_index_params.delta, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table);

		if (bot <= top)
		{
			return 0;
		}

	}




	(*sp) = top;
	(*ep) = bot;


	return bot - top;

}





void debug_information();


#endif

