#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <assert.h>
#include "kthread.h" // multi-threading models: pipeline and multi-threaded for loop
#include "yak-priv.h"
#include "ketopt.h"
#include "bseq.h"
#include <map>
#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "kseq.h" // FASTA/Q parser

using namespace std;

KSEQ_INIT(gzFile, gzread)
#define CHUNK_SIZE 200000000



typedef struct {
	int n, m;
	uint64_t n_ins;
	uint64_t *a;
	uint32_t *r;
	uint32_t *pos;
	bool *f;
} ch_buf_t;

static inline void ch_insert_buf(ch_buf_t *buf, int p, uint64_t y, uint32_t seq_id, bool f, uint32_t pos) // insert a k-mer $y to a linear buffer
{
	int pre = y & ((1<<p) - 1);
	ch_buf_t *b = &buf[pre];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->a, b->m);
		REALLOC(b->r, b->m);
		REALLOC(b->pos, b->m);
		REALLOC(b->f, b->m);
	}
	b->a[b->n] = y;
	b->f[b->n] = f;
	b->pos[b->n] = pos;
	b->r[b->n++] = seq_id;
}

static void count_seq_buf(ch_buf_t *buf, int k, int p, int len, const char *seq, uint32_t seq_id) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	bool f;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				f = x[0] < x[1];
				uint64_t y = f ? x[0] : x[1];
				ch_insert_buf(buf, p, yak_hash64(y, mask), seq_id, f, i);
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

static void count_seq_buf_long(ch_buf_t *buf, int k, int p, int len, const char *seq, uint32_t seq_id) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	bool f;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 1 | (c&1))  & mask;
			x[1] = (x[1] << 1 | (c>>1)) & mask;
			x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
			x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
			if (++l >= k)
				ch_insert_buf(buf, p, yak_hash_long(x), seq_id, f, i);
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
}

typedef struct { // global data structure for kt_pipeline()
	const yak_copt_t *opt;
	int create_new;
	kseq_t *ks;
	yak_ch_t *h;
	yak_ch_pos_t *h_pos;
	uint64_t global_counter;
	int seg_n;
	std::vector<char*>* names;
	std::vector<uint64_t>* lengths;
} pldat_t;

typedef struct { // data structure for each step in kt_pipeline()
	pldat_t *p;
	int n, m, sum_len, nk;
	uint64_t global_bias;
	int *len;
	char **seq;
	ch_buf_t *buf;
} stepdat_t;

static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
	stepdat_t *s = (stepdat_t*)data;
	yak_ch_t *h = s->p->h;
	yak_ch_pos_t *h_pos = s->p->h_pos;
	ch_buf_t *b = &s->buf[i];
	b->n_ins += yak_ch_insert_list_kmer_pos_haplo(h, h_pos, s->p->create_new, b->n, b->a, b->f, b->r, b->pos, i);
}

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
	pldat_t *p = (pldat_t*)data;
	if (step == 0) { // step 1: read a block of sequences
		int ret;
		stepdat_t *s;
		CALLOC(s, 1);
		s->p = p;
		s->global_bias = p->global_counter;
		while ((ret = kseq_read(p->ks)) >= 0) {
			p->global_counter++;
			int l = p->ks->seq.l;
			if (l < p->opt->k) continue;
			if (s->n == s->m) {
				s->m = s->m < 16? 16 : s->m + (s->n>>1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
			}
			MALLOC(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);

			p->names->push_back(strdup(p->ks->name.s));
			p->lengths->push_back((uint64_t)l);
			
			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += l - p->opt->k + 1;
			if (s->sum_len >= p->opt->chunk_size)
				break;
		}
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->opt->pre, m;
		CALLOC(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i) {
			s->buf[i].m = m;
			MALLOC(s->buf[i].a, m);
			MALLOC(s->buf[i].r, m);
			MALLOC(s->buf[i].pos, m);
			MALLOC(s->buf[i].f, m);
		}
		for (i = 0; i < s->n; ++i) {
			if (p->opt->k < 32)
				count_seq_buf(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i],s->global_bias+i);
			else
				count_seq_buf_long(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i],s->global_bias+i);
			free(s->seq[i]);
		}
		free(s->seq); free(s->len);
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->opt->pre;
		uint64_t n_ins = 0;

		kt_for(p->opt->n_thread, worker_for, s, n);
		for (i = 0; i < n; i++) {
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
			free(s->buf[i].r);
			free(s->buf[i].pos);
			free(s->buf[i].f);
		}
		p->h->tot += n_ins;
		free(s->buf);
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences; %ld distinct k-mers in the hash table\n", __func__,
				yak_realtime(), yak_cputime() / yak_realtime(), s->n, (long)p->h->tot);
		fprintf(stderr, "Total Processed %ld sequences\n", p->global_counter);
		free(s);
	}
	return 0;
}

pldat_t *extract_kmers_from_haplotypes(const char *fn, const yak_copt_t *opt, yak_ch_t *h0)
{
	pldat_t* pl;
	CALLOC(pl,1);
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	pl->ks = kseq_init(fp);
	pl->opt = opt;
	pl->create_new = 1;
	pl->h = yak_ch_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
	pl->h_pos = yak_ch_pos_init(opt->k, YAK_COUNTER_BITS_LONG, opt->bf_n_hash, opt->bf_shift);
	pl->names = new std::vector<char*>();
	pl->lengths = new std::vector<uint64_t>();
	kt_pipeline(3, worker_pipeline, pl, 3);
	kseq_destroy(pl->ks);
	gzclose(fp);
	return pl;
}




// Map unitigs to hic data.


typedef struct {
	int c[16];
	int sc[2];
	int nk;
} tb_cnt_t;

typedef struct {
	int max;
	uint32_t *s;
} tb_buf_t;

typedef struct {
	char* name;
	uint32_t unit_id;
	uint64_t pos;
	// bool forward;
	bool haplo;
} mapping_res_haplo_t;

typedef struct {
	int k, n_threads, print_diff;
	double ratio_thres;
	bseq_file_t *fp;
	yak_ch_t *ch;
	yak_ch_pos_t *ch_pos;
	tb_buf_t *buf;
	std::vector<mapping_res_haplo_t>* mappings;
	int record_num;
} tb_shared_t;

typedef struct {
	int n_seq;
	tb_shared_t *aux;
	bseq1_t *seq;
	uint32_t *mappings;
	uint32_t *mappings_pos;
	// bool *mappings_forward;
	bool *mappings_haplo;
} tb_step_t;

static void tb_worker(void *_data, long k, int tid)
{
	tb_step_t *t = (tb_step_t*)_data;
	tb_shared_t *aux = t->aux;
	bseq1_t *s = &t->seq[k];
	tb_buf_t *b = &aux->buf[tid];
	uint64_t x[4], mask;
	int i, l, shift;
	bool is_forward;
	uint32_t s_pos;
	std::map<uint32_t,int> counts;
	// std::map<uint32_t,bool> forward;
	std::map<uint32_t,int> pos;
	std::map<uint32_t,bool> haplo;

	if (aux->ch->k < 32) {
		mask = (1ULL<<2*aux->ch->k) - 1;
		shift = 2 * (aux->ch->k - 1);
	} else {
		mask = (1ULL<<aux->ch->k) - 1;
		shift = aux->ch->k - 1;
	}
	if (s->l_seq > b->max) {
		b->max = s->l_seq;
		kroundup32(b->max);
		b->s = (uint32_t*)realloc(b->s, b->max * sizeof(uint32_t));
	}
	memset(b->s, 0, s->l_seq * sizeof(uint32_t));
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->l_seq; ++i) {
		int c = seq_nt4_table[(uint8_t)s->seq[i]];
		uint32_t res;
		if (c < 4) {
			if (aux->ch->k < 32) {
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
			} else {
				x[0] = (x[0] << 1 | (c&1))  & mask;
				x[1] = (x[1] << 1 | (c>>1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
			}
			if (++l >= aux->k) {
				// int type = 0, c1, c2;
				uint64_t y;
				if (aux->ch->k < 32){
					is_forward = x[0] < x[1];
					y = yak_hash64(is_forward ? x[0] : x[1], mask);
				}else{
					y = yak_hash_long(x);
				}
				s_pos = 0;
				res = yak_ch_get_pos_haplo(aux->ch, aux->ch_pos, y, &s_pos);
				if(res != -1){
					counts[res&YAK_KEY_MASK_MID]++;
					// forward[res&YAK_KEY_MASK_MID] = !(((res&YAK_FORWARD_MASK_MID)!=0) ^ is_forward);
					pos[res&YAK_KEY_MASK_MID] = s_pos;
					haplo[res&YAK_KEY_MASK_MID] = (res&YAK_HAPLO_MASK_MID)!=0;
				}
			}
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0;
	}
	uint32_t max_idx = -1;
	int max_val = 0;
	for(auto idx: counts) {
		if(max_val < idx.second){
			max_val = idx.second;
			max_idx = idx.first;
		}
	}
	t->mappings[k] = max_idx;
	t->mappings_pos[k] = pos[max_idx];
	t->mappings_haplo[k] = haplo[max_idx];
	// t->mappings_forward[k] = forward[max_idx];
}

static void *tb_pipeline(void *shared, int step, void *_data)
{
	tb_shared_t *aux = (tb_shared_t*)shared;
	if (step == 0) {
		tb_step_t *s;
		s = (tb_step_t*)calloc(1, sizeof(tb_step_t));
		s->seq = bseq_read(aux->fp, CHUNK_SIZE, 0, &s->n_seq);
		s->aux = aux;
		if (s->n_seq) {
			s->mappings = (uint32_t*)calloc(s->n_seq, sizeof(uint32_t));
			s->mappings_pos = (uint32_t*)calloc(s->n_seq, sizeof(uint32_t));
			// s->mappings_forward = (bool*)calloc(s->n_seq, sizeof(bool));
			s->mappings_haplo = (bool*)calloc(s->n_seq, sizeof(bool));
			// fprintf(stderr, "[M::%s] read %d sequences\n", __func__, s->n_seq);
			return s;
		} else free(s);
	} else if (step == 1) {
		int i;
		tb_step_t *s = (tb_step_t*)_data;
		kt_for(aux->n_threads, tb_worker, s, s->n_seq);
			mapping_res_haplo_t res;
		for (i = 0; i < s->n_seq; ++i) {
			res.name = strdup(s->seq[i].name);
			res.unit_id = (s->mappings)[i];
			// res.forward = s->mappings_forward[i];
			res.pos = (s->mappings_pos)[i];

			res.haplo = (s->mappings_haplo)[i];
			aux->mappings->push_back(res);
			free(s->seq[i].name); free(s->seq[i].seq); free(s->seq[i].qual); free(s->seq[i].comment);
		}
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences;\n", __func__,
			yak_realtime(), yak_cputime() / yak_realtime(), s->n_seq);
		free(s->seq); free(s->mappings); free(s->mappings_pos); free(s->mappings_haplo); free(s);
	}
	return 0;
}

void map_hic_to_haplo_kmers(pldat_t* pl, bseq_file_t* hic_fn1, bseq_file_t* hic_fn2, char* out_fn)
{

	std::vector<char*> names = *(pl->names);
	std::vector<uint64_t> lengths = *(pl->lengths);

	if(hic_fn1 == 0 || hic_fn2 == 0){
		fprintf(stderr, "ERROR: Please give two hic files\n");
		exit(1);
	}

	// ketopt_t o = KETOPT_INIT;
	// int i, c, min_cnt = 2, mid_cnt = 5;
	tb_shared_t aux;

	memset(&aux, 0, sizeof(tb_shared_t));
	aux.n_threads = pl->opt->n_thread, aux.print_diff = 0;
	aux.ratio_thres = 0.33;
	aux.k = pl->h->k;
	aux.ch = pl->h;
	aux.ch_pos = pl->h_pos;
	aux.record_num = pl->global_counter;

	aux.fp = hic_fn1;
	std::vector<mapping_res_haplo_t>* map1 = new std::vector<mapping_res_haplo_t>();
	aux.mappings = map1;
	aux.buf = (tb_buf_t*)calloc(aux.n_threads, sizeof(tb_buf_t));
	kt_pipeline(2, tb_pipeline, &aux, 2);
	bseq_close(aux.fp);
	for (int i = 0; i < aux.n_threads; ++i){
		free(aux.buf[i].s);
	}
	free(aux.buf);



	aux.fp = hic_fn2;
	std::vector<mapping_res_haplo_t>* map2 = new std::vector<mapping_res_haplo_t>();
	aux.mappings = map2;
	aux.buf = (tb_buf_t*)calloc(aux.n_threads, sizeof(tb_buf_t));
	kt_pipeline(2, tb_pipeline, &aux, 2);
	bseq_close(aux.fp);
	for (int i = 0; i < aux.n_threads; ++i){
		free(aux.buf[i].s);
	}
	free(aux.buf);

	printf("Map1 size: %ld;\n", map1->size());
	printf("Map2 size: %ld;\n", map2->size());

	uint64_t success_counter_result = 0;
	uint32_t** connections;

	CALLOC(connections,pl->global_counter*2);
	for(uint i = 0; i< pl->global_counter*2; i++){
		CALLOC(connections[i], pl->global_counter*2);
		memset(connections[i], 0, sizeof(connections[i]));
	}
	std::vector<uint64_t> failed_matches;
	// const uint32_t max_count = -1;
	
	uint32_t counter_1 = 0;
	uint32_t counter_2 = 0;
	uint32_t counter_3 = 0;
	uint32_t counter_4 = 0;
	std::set<uint32_t> to_connect_from_1;
	std::set<uint32_t> to_connect_to_1;
	std::set<uint32_t> to_connect_from_2;
	std::set<uint32_t> to_connect_to_2;
	for(uint64_t j = 0; j < std::min(map1->size(), map2->size()); j++){

		if((*map1)[j].unit_id != -1 && (*map2)[j].unit_id != -1 ){
			success_counter_result++;
			mapping_res_haplo_t map_1 = (*map1)[j];
			mapping_res_haplo_t map_2 = (*map2)[j];
			uint32_t map_1_id = map_1.unit_id << 1;
			uint32_t map_2_id = map_2.unit_id << 1;
			bool map_1_front = (*map1)[j].pos < (lengths[(*map1)[j].unit_id]/2);
			bool map_2_front = (*map2)[j].pos < (lengths[(*map2)[j].unit_id]/2);
			// printf("%u, %u\n",(*map2)[j].pos, lengths[(*map2)[j].unit_id]);
			if( (map_1.unit_id>>1) != (map_2.unit_id>>1) ){
				if(map_1_front && map_2_front){
					counter_1++;
					to_connect_from_1.insert((map_1_id^1));
					to_connect_to_1.insert(map_2_id);

					to_connect_from_2.insert((map_2_id^1));
					to_connect_to_2.insert(map_1_id);

					if(map_1.haplo){
						to_connect_from_1.insert((map_1_id^1)^2);
						to_connect_to_2.insert(map_1_id^2);
					}
					if(map_2.haplo){
						to_connect_to_1.insert(map_2_id^2);
						to_connect_from_2.insert((map_2_id^1)^2);
					}
				}else if(map_1_front && !map_2_front){
					counter_2++;
					to_connect_from_1.insert((map_1_id^1));
					to_connect_to_1.insert((map_2_id^1));

					to_connect_from_2.insert(map_2_id);
					to_connect_to_2.insert(map_1_id);

					if(map_1.haplo){
						to_connect_from_1.insert((map_1_id^1)^2);
						to_connect_to_2.insert(map_1_id^2);
					}
					if(map_2.haplo){
						to_connect_to_1.insert((map_2_id^1)^2);
						to_connect_from_2.insert(map_2_id^2);
					}
				}else if(!map_1_front && map_2_front){
					counter_3++;
					to_connect_from_1.insert(map_1_id);
					to_connect_to_1.insert(map_2_id);

					to_connect_from_2.insert((map_2_id^1));
					to_connect_to_2.insert((map_1_id^1));

					if(map_1.haplo){
						to_connect_from_1.insert(map_1_id^2);
						to_connect_to_2.insert((map_1_id^1)^2);
					}
					if(map_2.haplo){
						to_connect_to_1.insert(map_2_id^2);
						to_connect_from_2.insert((map_2_id^1)^2);
					}
				}else if(!map_1_front && !map_2_front){
					counter_4++;
					to_connect_from_1.insert(map_1_id);
					to_connect_to_1.insert((map_2_id^1));

					to_connect_from_2.insert(map_2_id);
					to_connect_to_2.insert((map_1_id^1));

					if(map_1.haplo){
						to_connect_from_1.insert(map_1_id^2);
						to_connect_to_2.insert((map_1_id^1)^2);
					}
					if(map_2.haplo){
						to_connect_to_1.insert((map_2_id^1)^2);
						to_connect_from_2.insert(map_2_id^2);
					}
				}
				for(auto i: to_connect_from_1){
					for(auto j: to_connect_to_1){
						connections[i][j]++;
					}
				}
				to_connect_from_1.clear();
				to_connect_to_1.clear();
				for(auto i: to_connect_from_2){
					for(auto j: to_connect_to_2){
						connections[i][j]++;
					}
				}
				to_connect_from_2.clear();
				to_connect_to_2.clear();
			}
		}
	}
	// printf("%d,%d,%d,%d\n",counter_1,counter_2,counter_3,counter_4);
	for(uint i = 0; i < pl->global_counter*2; i++){
		double ratio_i = 1;
		if(lengths[i>>1]<20000000){
			ratio_i = 20000000.0/((double)lengths[i>>1]);
		}
		for(uint j = 0; j < pl->global_counter*2; j++){
			double ratio_j = 1;
			if(lengths[j>>1]<20000000){
				ratio_j = 20000000.0/((double)lengths[j>>1]);
			}
			connections[i][j] = (uint32_t)(connections[i][j] * ratio_i * ratio_j);
		}
	}
	printf("Useful Records: %ld\n" , success_counter_result);
	printf("Not matched due to Wrong names: %ld\n" , failed_matches.size());


	// {
	// 	FILE* fp = fopen ("coverage_report.output","w");
	// 	uint32_t total_no_coverage = 0;
	// 	uint32_t total_have_coverage = 0;
	// 	for(int i = 0; i < pl->global_counter; i++){
	// 		if(coverage[i]==0){
	// 			total_no_coverage++;
	// 		}else{
	// 			total_have_coverage++;
	// 		}
	// 	}
	// 	fprintf(fp, "Total No Coverage: %d\n", total_no_coverage);
	// 	for(int i = 0; i < pl->global_counter; i++){
	// 		if(coverage[i]==0){
	// 			fprintf(fp, "%s\n",(*pl->names)[i]);
	// 		}
	// 	}
	// 	fprintf(fp, "Total Have Coverage: %d\n",total_have_coverage);
	// 	for(int i = 0; i < pl->global_counter; i++){
	// 		if(coverage[i]!=0){
	// 			fprintf(fp, "%s, %d\n",(*pl->names)[i],coverage[i]);
	// 		}
	// 	}
	// 	fclose(fp);
	// }

	if(out_fn){
		ofstream outFileScaffold;
		outFileScaffold.open(string(out_fn), ofstream::out | ofstream::trunc);
		
		// for(int i = 0; i < pl->global_counter; i++){
		// 	uint32_t is_forward = forward[i] > backward[i] ? 1 : 0;
		// 	fprintf(fp, "%d\t%d\t%d\n", is_forward, backward[i], forward[i]);
		// }
		for(uint i = 0; i < pl->global_counter*2; i++){
			for(uint j = 0; j < pl->global_counter*2; j++){
				if((i&2)==0){
					connections[i][j] += connections[i^2][j^2];
				}else{
					connections[i][j] = 0;
				}
				if(connections[i][j]>0){
					outFileScaffold << string(names[i>>1]) << (i%2==0?"+":"-") << "\t" << string(names[j>>1]) << (j%2==0?"+":"-") << "\t" << connections[i][j] << endl;
				}
			}
		}
   	 	outFileScaffold.close();
	}

	// return connections;
}

int main_hic_map_haplo(int argc, char *argv[])
{
	pldat_t *h;
	int c;
	char *fn_out = 0;
	yak_copt_t opt;
	ketopt_t o = KETOPT_INIT;
	yak_copt_init(&opt);
	opt.pre = YAK_COUNTER_BITS_MID;
	opt.n_thread = 32;
	while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:b:H:o:", 0)) >= 0) {
		if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 'K') opt.chunk_size = atoi(o.arg);
		else if (c == 't') opt.n_thread = atoi(o.arg);
		else if (c == 'b') opt.bf_shift = atoi(o.arg);
		else if (c == 'o') fn_out = o.arg;
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: pstools hic_mapping_haplo [options] <pred_haplotypes.fa> <hic.R1.fastq.gz> <hic.R2.fastq.gz>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -b INT     set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", opt.bf_shift);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
		fprintf(stderr, "  -K INT     chunk size [100m]\n");
		fprintf(stderr, "  -o FILE    save mapping relationship to FILE []\n");
		fprintf(stderr, "Note: -b37 is recommended for human reads\n");
		return 1;
	}
	if (opt.k >= 64) {
		fprintf(stderr, "ERROR: -k must be smaller than 64\n");
		return 1;
	} else if (opt.k >= 32) {
		fprintf(stderr, "WARNING: counts are inexact if -k is greater than 31\n");
	}

	if(argv[o.ind+1] == 0 || argv[o.ind+2] == 0){
		fprintf(stderr, "ERROR: Please give two hic files\n");
		exit(1);
	}

	bseq_file_t* hic_fn1 = bseq_open(argv[o.ind+1]);
	bseq_file_t* hic_fn2 = bseq_open(argv[o.ind+2]);
	if (hic_fn1 == 0 || hic_fn2 == 0) {
		fprintf(stderr, "ERROR: fail to open hic files\n");
		exit(1);
	}

	h = extract_kmers_from_haplotypes(argv[o.ind], &opt, 0);

	
	map_hic_to_haplo_kmers(h, hic_fn1, hic_fn2, fn_out);


	free(h->h);
	free(h);
	return 0;
}
