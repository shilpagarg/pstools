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
#include <algorithm> 
#include <string>
#include <tuple>
#include <iostream>
#include "kseq.h" // FASTA/Q parser

KSEQ_INIT(gzFile, gzread)
#define CHUNK_SIZE 200000000

typedef struct {
	int n, m;
	uint64_t n_ins;
	uint64_t *a;
	uint16_t *r;
	uint32_t *pos;
} ch_buf_t;

static inline void ch_insert_buf(ch_buf_t *buf, int p, uint64_t y, uint16_t seq_id, uint32_t pos) // insert a k-mer $y to a linear buffer
{
	int pre = y & ((1<<p) - 1);
	ch_buf_t *b = &buf[pre];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->a, b->m);
		REALLOC(b->r, b->m);
		REALLOC(b->pos, b->m);
	}
	b->a[b->n] = y;
	b->pos[b->n] = pos;
	b->r[b->n++] = seq_id;
}

static void count_seq_buf(ch_buf_t *buf, int k, int p, int len, const char *seq, uint16_t seq_id) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				ch_insert_buf(buf, p, yak_hash64(y, mask), seq_id, i);
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

static void count_seq_buf_long(ch_buf_t *buf, int k, int p, int len, const char *seq, uint16_t seq_id) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 1 | (c&1))  & mask;
			x[1] = (x[1] << 1 | (c>>1)) & mask;
			x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
			x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
			if (++l >= k)
				ch_insert_buf(buf, p, yak_hash_long(x), seq_id, i);
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
}

typedef struct { // global data structure for kt_pipeline()
	const yak_copt_t *opt;
	int create_new;
	kseq_t *ks;
	yak_ch_t *h;
	uint64_t global_counter;
	std::vector<std::string> name;
	int seg_n;
	uint64_t tot_len;
	recordset_ps_t *record_set;
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
	ch_buf_t *b = &s->buf[i];
	recordset_ps_t *r = &s->p->record_set[i];
	b->n_ins += yak_ch_insert_list_kmer_pos(h, s->p->create_new, b->n, b->a, r, b->r, b->pos, i, NULL);
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
			p->name.push_back(std::string(strdup(p->ks->name.s)));
			if (l < p->opt->k) continue;
			p->tot_len+=l;
			if (s->n == s->m) {
				s->m = s->m < 16? 16 : s->m + (s->n>>1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
			}
			MALLOC(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);
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
			// printf("%d\n",s->p->record_set[i].current_length);
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
			free(s->buf[i].r);
			free(s->buf[i].pos);
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

pldat_t *get_haplo_kmers(const char *fn, const yak_copt_t *opt)
{
	pldat_t* pl;
	CALLOC(pl,1);
	gzFile fp;
	pl->opt = opt;
	pl->create_new = 1;
	pl->h = yak_ch_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
	CALLOC(pl->record_set, 1<<opt->pre);
	for(int i=0; i<(1<<opt->pre);i++){
		pl->record_set[i].current_length = 16;
		CALLOC(pl->record_set[i].records, 16);
	}
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	pl->ks = kseq_init(fp);
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
	uint16_t unit_id;
	uint64_t pos;
} mapping_res_t;

typedef struct {
	int k, n_threads, print_diff;
	double ratio_thres;
	bseq_file_t *fp;
	yak_ch_t *ch;
	tb_buf_t *buf;
	std::vector<mapping_res_t>* mappings;
	int record_num;
	recordset_ps_t *record_set;
} tb_shared_t;

typedef struct {
	int n_seq;
	tb_shared_t *aux;
	bseq1_t *seq;
	uint16_t *mappings;
	uint64_t *map_pos;
} tb_step_t;

static void tb_worker(void *_data, long k, int tid)
{
	uint64_t max_64 = -1;
	tb_step_t *t = (tb_step_t*)_data;
	tb_shared_t *aux = t->aux;
	bseq1_t *s = &t->seq[k];
	tb_buf_t *b = &aux->buf[tid];
	uint64_t x[4], mask;
	int i, l, shift;
	std::map<uint16_t,int> counts;
	std::map<uint16_t,uint64_t> pos;
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
		record_ps_t* res;
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
				int type = 0, c1, c2;
				uint64_t y;
				if (aux->ch->k < 32)
					y = yak_hash64(x[0] < x[1]? x[0] : x[1], mask);
				else
					y = yak_hash_long(x);
				res = yak_ch_get_pos(aux->ch, y, aux->record_set);
				if(res != NULL){
					counts[res->uni_id]++;
					pos[res->uni_id] = res->pos;
				}
			}
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0;
	}
	uint16_t max_idx = -1;
	int max_val = 0;
	for(auto idx: counts) {
		if(max_val < idx.second){
			max_val = idx.second;
			max_idx = idx.first;
		}
	}
	t->mappings[k] = max_idx;
	t->map_pos[k] = pos[max_idx];
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
			s->mappings = (uint16_t*)calloc(s->n_seq, sizeof(uint16_t));
			s->map_pos = (uint64_t*)calloc(s->n_seq, sizeof(uint64_t));
			// fprintf(stderr, "[M::%s] read %d sequences\n", __func__, s->n_seq);
			return s;
		} else free(s);
	} else if (step == 1) {
		int i;
		tb_step_t *s = (tb_step_t*)_data;
		kt_for(aux->n_threads, tb_worker, s, s->n_seq);
		for (i = 0; i < s->n_seq; ++i) {
			mapping_res_t res;
			res.unit_id = s->mappings[i];
			res.pos = s->map_pos[i];
			aux->mappings->push_back(res);
			free(s->seq[i].name); free(s->seq[i].seq); free(s->seq[i].qual); free(s->seq[i].comment);
		}
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences;\n", __func__,
			yak_realtime(), yak_cputime() / yak_realtime(), s->n_seq);
		free(s->seq); free(s->mappings); free(s->map_pos); free(s);
	}
	return 0;
}

void update_best_buddy(bool** seen, float** best_buddy, uint32_t** connections,uint32_t len){
	for(int i = 0; i < len; i++){
		for(int j = 0; j < len; j++){
			uint32_t num = connections[i*2][j*2] + connections[i*2+1][j*2] + connections[i*2][j*2+1] + connections[i*2+1][j*2+1];
			if(!seen[i][j] && num > 0){
				uint32_t max_other = 0;
				for(int p = 0; p < len; p++){
					if(p!=j){
						uint32_t cur_num = connections[i*2][p*2] + connections[i*2+1][p*2] + connections[i*2][p*2+1] + connections[i*2+1][p*2+1];
						if(!seen[i][p] && cur_num > max_other){
							max_other = cur_num;
						}
					}
				}

				for(int p = 0; p < len; p++){
					if(p!=i){
						uint32_t cur_num = connections[p*2][j*2] + connections[p*2+1][j*2] + connections[p*2][j*2+1] + connections[p*2+1][j*2+1];
						if(!seen[p][j] && cur_num > max_other){
							max_other = cur_num;
						}
					}
				}
				if(max_other>0){
					best_buddy[i][j] = ((float)num)/max_other;
				}else{
					best_buddy[i][j] = 2;
				}
			}else{
				best_buddy[i][j] = 0;
			}
		}
	}
}

void map_hic_pairs(pldat_t* pl, bseq_file_t* hic_fn1, bseq_file_t* hic_fn2)
{
	if(hic_fn1 == 0 || hic_fn2 == 0){
		fprintf(stderr, "ERROR: Please give two hic files\n");
		exit(1);
	}

	ketopt_t o = KETOPT_INIT;
	int i, c, min_cnt = 2, mid_cnt = 5;
	tb_shared_t aux;

	memset(&aux, 0, sizeof(tb_shared_t));
	aux.n_threads = pl->opt->n_thread, aux.print_diff = 0;
	aux.ratio_thres = 0.33;
	aux.k = pl->h->k;
	aux.ch = pl->h;
	aux.record_set = pl->record_set;
	aux.record_num = pl->global_counter;
	
	aux.fp = hic_fn1;
	std::vector<mapping_res_t>* map1 = new std::vector<mapping_res_t>();
	aux.mappings = map1;
	aux.buf = (tb_buf_t*)calloc(aux.n_threads, sizeof(tb_buf_t));
	kt_pipeline(2, tb_pipeline, &aux, 2);
	bseq_close(aux.fp);
	for (i = 0; i < aux.n_threads; ++i){
		free(aux.buf[i].s);
	}
	free(aux.buf);



	aux.fp = hic_fn2;
	std::vector<mapping_res_t>* map2 = new std::vector<mapping_res_t>();
	aux.mappings = map2;
	aux.buf = (tb_buf_t*)calloc(aux.n_threads, sizeof(tb_buf_t));
	kt_pipeline(2, tb_pipeline, &aux, 2);
	bseq_close(aux.fp);
	for (i = 0; i < aux.n_threads; ++i){
		free(aux.buf[i].s);
	}
	free(aux.buf);

	// std::vector<uint64_t> connections[pl->global_counter][pl->global_counter];
	uint32_t** connect_num;
	float **best_buddy;
	bool **seen;
	CALLOC(connect_num, pl->global_counter);
	CALLOC(best_buddy, pl->global_counter/2);
	CALLOC(seen, pl->global_counter/2);
	for(int i = 0; i < pl->global_counter; ++i){
		CALLOC(connect_num[i], pl->global_counter);
		memset(connect_num[i], 0, pl->global_counter*sizeof(uint32_t));
	}
	for(int i = 0; i < pl->global_counter/2; ++i){
		CALLOC(seen[i], pl->global_counter/2);
		memset(seen[i], 0, pl->global_counter/2*sizeof(bool));
		CALLOC(best_buddy[i], pl->global_counter/2);
		memset(best_buddy[i], 0, pl->global_counter/2*sizeof(float));
	}

	for(uint64_t j = 0; j < std::min(map1->size(), map2->size()); j++){
		if((*map1)[j].unit_id != 65535 && (*map2)[j].unit_id != 65535){
			uint16_t id_1 = (*map1)[j].unit_id;
			uint16_t id_2 = (*map2)[j].unit_id;
			uint64_t pos_1 = (*map2)[j].pos;
			uint64_t pos_2 = (*map2)[j].pos;
			
			if(id_1>>1 != id_2){
				// connections[id_1][id_2].push_back(pos_1);
				// connections[id_2][id_1].push_back(pos_2);
				connect_num[id_1][id_2]++;
				connect_num[id_2][id_1]++;
			}
		}
	}
	// for(int i = 0; i < pl->global_counter/2; i++){
	// 	for(int j = 0; j < pl->global_counter/2; j++){
	// 		printf("%d\t%d\t%d\n",i,j,connect_num[i<<1][j<<1]+connect_num[(i<<1)+1][j<<1]+connect_num[(i<<1)][(j<<1)+1]+connect_num[(i<<1)+1][(j<<1)+1]);
	// 	}
	// }
	// std::vector<std::vector<uint32_t>> seen_sets(23);

	uint32_t len = pl->global_counter/2;
	for(int i = 0; i < len; i++){
		for(int j = 0; j < len; j++){
			if(connect_num[i*2][j*2] + connect_num[i*2+1][j*2] + connect_num[i*2][j*2+1] + connect_num[i*2+1][j*2+1] ==0 || i == j){
				seen[i][j] = true;
			}
		}
	}
	// for(int s = 0; s < 23; s++){
	// 	bool next = true;
	// 	while(next){
	// 		std::vector<std::pair<float, std::pair<uint32_t,uint32_t>>> result;
	// 		next = false;
	// 		update_best_buddy(seen, best_buddy,connect_num,len);
	// 		for(int i = 0; i < len; i++){
	// 			for(int j = 0; j < len; j++){
	// 				if(i<j && !seen[i][j] && best_buddy[i][j] > 1){
	// 					result.push_back(std::make_pair(best_buddy[i][j], std::make_pair(i,j)));
	// 				}
	// 			}
	// 		}
	// 		std::sort(result.begin(), result.end(),
	// 		[](const std::pair<float, std::pair<uint32_t,uint32_t>>& r1, const std::pair<float, std::pair<uint32_t,uint32_t>>& r2){
	// 			return r1.first > r2.first;
	// 		});
			
	// 		for(auto res: result){
	// 			bool found = false;
	// 			for(int q = 0; q < s; q++){
	// 				if(( std::find(seen_sets[q].begin(), seen_sets[q].end(), res.second.first) != seen_sets[q].end() ) || ( std::find(seen_sets[q].begin(), seen_sets[q].end(), res.second.second) != seen_sets[q].end() )){
	// 					found = true;
	// 					next = true;
	// 					seen[res.second.first][res.second.second] = true;
	// 					seen[res.second.second][res.second.first] = true;
	// 					break;
	// 				}
	// 			}
	// 			if(!found){
	// 				if(seen_sets[s].size()==0){
	// 					seen_sets[s].push_back(res.second.first);
	// 					seen_sets[s].push_back(res.second.second);
	// 					next = true;
	// 					seen[res.second.first][res.second.second] = true;
	// 					seen[res.second.second][res.second.first] = true;
	// 				}else{
	// 					if(( std::find(seen_sets[s].begin(), seen_sets[s].end(), res.second.first) != seen_sets[s].end() ) && ( std::find(seen_sets[s].begin(), seen_sets[s].end(), res.second.second) != seen_sets[s].end() )){
	// 						next = true;
	// 						seen[res.second.first][res.second.second] = true;
	// 						seen[res.second.second][res.second.first] = true;
	// 					}else if(( std::find(seen_sets[s].begin(), seen_sets[s].end(), res.second.first) != seen_sets[s].end() ) ){
	// 						seen_sets[s].insert(seen_sets[s].begin(), res.second.second);
	// 						next = true;
	// 						seen[res.second.first][res.second.second] = true;
	// 						seen[res.second.second][res.second.first] = true;
	// 					}else if(( std::find(seen_sets[s].begin(), seen_sets[s].end(), res.second.second) != seen_sets[s].end() ) ){
	// 						seen_sets[s].insert(seen_sets[s].begin(), res.second.first);
	// 						next = true;
	// 						seen[res.second.first][res.second.second] = true;
	// 						seen[res.second.second][res.second.first] = true;
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	// for(int i = 0; i < len; i++){
	// 	for(int j = 0; j < len; j++){
	// 		if(connect_num[i][j]==0 || i == j){
	// 			seen[i][j] = true;
	// 		}else{
	// 			seen[i][j] = false;
	// 		}
	// 	}
	// }
	// bool next = true;
	// while(next){
	// 	std::vector<std::pair<float, std::pair<uint32_t,uint32_t>>> result;
	// 	next = false;
	// 	update_best_buddy(seen, best_buddy,connect_num,len);
	// 	for(int i = 0; i < len; i++){
	// 		for(int j = 0; j < len; j++){
	// 			if(i<j && !seen[i][j] && best_buddy[i][j] > 1){
	// 				result.push_back(std::make_pair(best_buddy[i][j], std::make_pair(i,j)));
	// 			}
	// 		}
	// 	}
	// 	std::sort(result.begin(), result.end(),
	// 	[](const std::pair<float, std::pair<uint32_t,uint32_t>>& r1, const std::pair<float, std::pair<uint32_t,uint32_t>>& r2){
	// 		return r1.first > r2.first;
	// 	});
		
	// 	for(auto res: result){
	// 		bool found = false;
	// 		bool once = true;
	// 		uint32_t id = 0;
	// 		for(int q = 0; q < 23; q++){
	// 			if(( std::find(seen_sets[q].begin(), seen_sets[q].end(), res.second.first) != seen_sets[q].end() ) || ( std::find(seen_sets[q].begin(), seen_sets[q].end(), res.second.second) != seen_sets[q].end() )){
	// 				if(!found){
	// 					found = true;
	// 				}else{
	// 					once = false;
	// 				}
	// 				// found = true;
	// 				next = true;
	// 				seen[res.second.first][res.second.second] = true;
	// 				seen[res.second.second][res.second.first] = true;
	// 				id = q;
	// 				// break;
	// 			}
	// 		}
	// 		if(found&&once){
	// 			if(( std::find(seen_sets[id].begin(), seen_sets[id].end(), res.second.first) != seen_sets[id].end() ) && ( std::find(seen_sets[id].begin(), seen_sets[id].end(), res.second.second) != seen_sets[id].end() )){
	// 				next = true;
	// 				seen[res.second.first][res.second.second] = true;
	// 				seen[res.second.second][res.second.first] = true;
	// 			}else if(( std::find(seen_sets[id].begin(), seen_sets[id].end(), res.second.first) != seen_sets[id].end() ) ){
	// 				seen_sets[id].insert(seen_sets[id].begin(), res.second.second);
	// 				next = true;
	// 				seen[res.second.first][res.second.second] = true;
	// 				seen[res.second.second][res.second.first] = true;
	// 			}else if(( std::find(seen_sets[id].begin(), seen_sets[id].end(), res.second.second) != seen_sets[id].end() ) ){
	// 				seen_sets[id].insert(seen_sets[id].begin(), res.second.first);
	// 				next = true;
	// 				seen[res.second.first][res.second.second] = true;
	// 				seen[res.second.second][res.second.first] = true;
	// 			}
	// 		}
	// 	}
	// }
	
	std::vector<std::vector<uint32_t>> graph_connection;
	// bool next = true;
	std::set<uint32_t> seen_node;
	// while(next){
	while(seen_node.size()!=pl->global_counter/2){
		// for(auto n: seen_node){
		// 	printf("%d, ",n);
		// }
		// printf("\n");
		std::vector<std::pair<float, std::pair<uint32_t,uint32_t>>> result;
		// next = false;
		update_best_buddy(seen, best_buddy,connect_num,len);
		for(int i = 0; i < len; i++){
			for(int j = 0; j < len; j++){
				if(i<j && !seen[i][j] && best_buddy[i][j] > 1){
					result.push_back(std::make_pair(best_buddy[i][j], std::make_pair(i,j)));
				}
			}
		}
		
		std::sort(result.begin(), result.end(),
		[](const std::pair<float, std::pair<uint32_t,uint32_t>>& r1, const std::pair<float, std::pair<uint32_t,uint32_t>>& r2){
			return r1.first > r2.first;
		});
	
		for(auto res: result){
			// printf("%d, %d\n",res.second.first,res.second.second);
			// next = true;
			bool found = false;
			uint32_t idx_i = 65535;
			uint32_t idx_j = 65535;
			uint32_t idx_i_sub = 65535;
			uint32_t idx_j_sub = 65535;
			for(int i = 0; i < graph_connection.size(); i++){
				if(find(graph_connection[i].begin(), graph_connection[i].end(), res.second.first)!=graph_connection[i].end()){
					for(int idx=0; idx<graph_connection[i].size(); idx++){
						if(graph_connection[i][idx]==res.second.first){
							idx_i_sub = idx;
						}
					}
					idx_i = i;
					found = true;
				}
				if(find(graph_connection[i].begin(), graph_connection[i].end(), res.second.second)!=graph_connection[i].end()){
					for(int idx=0; idx<graph_connection[i].size(); idx++){
						if(graph_connection[i][idx]==res.second.second){
							idx_j_sub = idx;
						}
					}
					idx_j = i;
					found = true;
				}
			}
			if(!found){
				std::vector<uint32_t> buf;
				buf.push_back(std::min(res.second.first, res.second.second));
				buf.push_back(std::max(res.second.first, res.second.second));
				graph_connection.push_back(buf);
				seen[res.second.first][res.second.second] = true;
				seen[res.second.second][res.second.first] = true;
				seen_node.insert(res.second.first);
				seen_node.insert(res.second.second);
			}else{
				std::vector<uint32_t> connect_i;
				std::vector<uint32_t> connect_j;
				if(idx_i!=65535){
					std::vector<uint32_t> connect_i = graph_connection[idx_i];
				}
				if(idx_j!=65535){
					std::vector<uint32_t> connect_j = graph_connection[idx_j];
				}
				if(idx_i == idx_j){
					seen[res.second.first][res.second.second] = true;
					seen[res.second.second][res.second.first] = true;
				}else{
					if(idx_i == 65535){
						if(idx_j_sub >= graph_connection[idx_j].size()/2){
							graph_connection[idx_j].push_back(res.second.first);
						}else{
							graph_connection[idx_j].insert(graph_connection[idx_j].begin(), res.second.first);
						}
						seen[res.second.first][res.second.second] = true;
						seen[res.second.second][res.second.first] = true;
						seen_node.insert(res.second.first);
						seen_node.insert(res.second.second);
					}else if(idx_j == 65535){
						if(idx_i_sub >= graph_connection[idx_i].size()/2){
							graph_connection[idx_i].push_back(res.second.second);
						}else{
							graph_connection[idx_i].insert(graph_connection[idx_i].begin(), res.second.second);
						}
						seen[res.second.first][res.second.second] = true;
						seen[res.second.second][res.second.first] = true;
						seen_node.insert(res.second.first);
						seen_node.insert(res.second.second);
					}else{
						// graph_connection.erase(graph_connection.begin()+idx_i);
						// graph_connection.erase(graph_connection.begin()+idx_j);
						if(idx_i_sub >= graph_connection[idx_i].size()/2 && idx_j_sub < graph_connection[idx_j].size()/2){
							for(auto n: graph_connection[idx_j]){
								graph_connection[idx_i].push_back(n);
							}
						}else if(idx_i_sub >= graph_connection[idx_i].size()/2 && idx_j_sub >= graph_connection[idx_j].size()/2){
							std::reverse(graph_connection[idx_j].begin(), graph_connection[idx_j].end());
							for(auto n: graph_connection[idx_j]){
								graph_connection[idx_i].push_back(n);
							}
						}else if(idx_i_sub < graph_connection[idx_i].size()/2 && idx_j_sub < graph_connection[idx_j].size()/2){
							std::reverse(graph_connection[idx_i].begin(), graph_connection[idx_i].end());
							for(auto n: graph_connection[idx_j]){
								graph_connection[idx_i].push_back(n);
							}
						}else if(idx_i_sub < graph_connection[idx_i].size()/2 && idx_j_sub >= graph_connection[idx_j].size()/2){
							std::reverse(graph_connection[idx_i].begin(), graph_connection[idx_i].end());
							std::reverse(graph_connection[idx_j].begin(), graph_connection[idx_j].end());
							for(auto n: graph_connection[idx_j]){
								graph_connection[idx_i].push_back(n);
							}
						}
						graph_connection.erase(graph_connection.begin()+idx_j);
						seen[res.second.first][res.second.second] = true;
						seen[res.second.second][res.second.first] = true;
						// connect_i.insert( connect_i.end(), connect_j.begin(), connect_j.end() );
						// graph_connection.push_back(connect_i);
						seen_node.insert(res.second.first);
						seen_node.insert(res.second.second);
					}
				}

			}
		}
		printf("%d\n",seen_node.size());
	}
	
	printf("Result:\n");
	for(int i = 0; i < graph_connection.size(); i++){
		std::vector<uint32_t> cur_set = graph_connection[i];
		for(auto p:cur_set){
			printf("%d, ", p);
			// printf("%s, ", pl->name[p<<1].c_str());
		}
		printf("\n");
	}
}

int main_haplo_map(yak_copt_t opt, char* hic_file1,char* hic_file2, char* pred_file)
{
	pldat_t *h;
	bseq_file_t* hic_fn1 = bseq_open(hic_file1);
	bseq_file_t* hic_fn2 = bseq_open(hic_file2);
	if (hic_fn1 == 0 || hic_fn2 == 0) {
		fprintf(stderr, "ERROR: fail to open hic files\n");
		exit(1);
	}


	h = get_haplo_kmers(pred_file, &opt);
	map_hic_pairs(h, hic_fn1, hic_fn2);
	
	// bseq_close(hic_fn1);
	// bseq_close(hic_fn2);
	yak_ch_destroy(h->h);
	free(h);
	return 0;
}
