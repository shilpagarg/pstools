#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <assert.h>
#include "kthread.h" // multi-threading models: pipeline and multi-threaded for loop
#include "yak-priv.h"
// #include "ketopt.h"
#include "bseq.h"
#include <map>
#include <vector>
#include <set>
#include <algorithm> 
#include <string>
#include <iostream>
#include "kseq.h" // FASTA/Q parser
#include <tuple>
#include <inttypes.h>
#include <fstream>

using namespace std;
KSEQ_INIT(gzFile, gzread)
#define CHUNK_SIZE 200000000

typedef struct {
	int n, m;
	uint64_t n_ins;
	uint64_t *a;
	uint32_t *r;
	uint32_t *pos;
} ch_buf_t;

static inline void ch_insert_buf(ch_buf_t *buf, int p, uint64_t y, uint32_t seq_id, uint32_t pos) // insert a k-mer $y to a linear buffer
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

static void count_seq_buf(ch_buf_t *buf, int k, int p, int len, const char *seq, uint32_t seq_id) // insert k-mers in $seq to linear buffer $buf
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

static void count_seq_buf_long(ch_buf_t *buf, int k, int p, int len, const char *seq, uint32_t seq_id) // insert k-mers in $seq to linear buffer $buf
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
	yak_ch_pos_t *h_pos;
	uint64_t global_counter;
	uint64_t hap1_final_id;
	int seg_n;
	uint64_t tot_len;
	std::vector<std::string>* names;
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
	b->n_ins += yak_ch_insert_list_kmer_pos(h, h_pos, s->p->create_new, b->n, b->a, b->r, b->pos, i);
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
			p->names->push_back(std::string(p->ks->name.s));
			p->lengths->push_back((uint64_t)l);
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

pldat_t *yak_counting_switch_error(string fn1, string fn2, const yak_copt_t *opt)
{
	pldat_t* pl;
	CALLOC(pl,1);
	gzFile fp;
	pl->names = new std::vector<std::string>();
	pl->lengths = new std::vector<uint64_t>();
	pl->opt = opt;
	pl->create_new = 1;
	pl->h = yak_ch_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
	pl->h_pos = yak_ch_pos_init(opt->k, YAK_COUNTER_BITS_LONG, opt->bf_n_hash, opt->bf_shift);
	if ((fp = gzopen(fn1.c_str(), "r")) == 0) return 0;
	pl->ks = kseq_init(fp);
	kt_pipeline(3, worker_pipeline, pl, 3);
	kseq_destroy(pl->ks);
	gzclose(fp);
	pl->hap1_final_id = pl->global_counter;
	if ((fp = gzopen(fn2.c_str(), "r")) == 0) return 0;
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
	uint32_t unit_id;
	uint64_t pos;
} mapping_res_t;

typedef struct {
	int k, n_threads, print_diff;
	double ratio_thres;
	bseq_file_t *fp;
	yak_ch_t *ch;
	yak_ch_pos_t *ch_pos;
	tb_buf_t *buf;
	std::vector<mapping_res_t>* mappings;
	int record_num;
} tb_shared_t;

typedef struct {
	int n_seq;
	tb_shared_t *aux;
	bseq1_t *seq;
	uint32_t *mappings;
	uint32_t *map_pos;
} tb_step_t;

static void tb_worker(void *_data, long k, int tid)
{
	tb_step_t *t = (tb_step_t*)_data;
	tb_shared_t *aux = t->aux;
	bseq1_t *s = &t->seq[k];
	tb_buf_t *b = &aux->buf[tid];
	uint64_t x[4], mask;
	int i, l, shift;
	uint32_t s_pos;
	std::map<uint32_t,int> counts;
	std::map<uint32_t,int> pos;
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
				if (aux->ch->k < 32)
					y = yak_hash64(x[0] < x[1]? x[0] : x[1], mask);
				else
					y = yak_hash_long(x);
				uint32_t res = yak_ch_get_pos(aux->ch, aux->ch_pos, y, &s_pos);
				if(res != -1){
					counts[res]++;
					pos[res] = s_pos;
					// printf("%d\t%d\t%ld\n",res,s_pos,k);
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
	t->map_pos[k] = pos[max_idx];
	// printf("%d\t%d\t%ld\n",max_idx,pos[max_idx],k);
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
			s->map_pos = (uint32_t*)calloc(s->n_seq, sizeof(uint32_t));
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

void get_switch_error(pldat_t* pl, bseq_file_t* hic_fn1, bseq_file_t* hic_fn2, std::map<std::string,std::string> contig2_1_map)
{

	std::vector<std::string> names = *pl->names;
	std::vector<uint64_t> lengths = *pl->lengths;

	if(hic_fn1 == 0 || hic_fn2 == 0){
		fprintf(stderr, "ERROR: Please give two hic files\n");
		exit(1);
	}

	// ketopt_t o = KETOPT_INIT;
	int i;
	// int c;
	tb_shared_t aux;

	memset(&aux, 0, sizeof(tb_shared_t));
	aux.n_threads = pl->opt->n_thread, aux.print_diff = 0;
	aux.ratio_thres = 0.33;
	aux.k = pl->h->k;
	aux.ch = pl->h;
	aux.ch_pos = pl->h_pos;
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

	map<uint32_t, int> unid_id;
	map<int, string> id_names;
	map<int, uint64_t> id_lengths;
	int ids = 0;
	for(auto i: contig2_1_map){
		for(uint idx = 0; idx < names.size()/2; idx++){
			if(names[idx]==i.first){
				unid_id[idx] = ids;
				id_names[ids] = names[idx];
				id_lengths[ids] = lengths[idx];
				break;
			}
		}
		for(uint idx = names.size()/2; idx < names.size(); idx++){
			if(names[idx]==i.second){
				unid_id[idx] = ids+1;
				id_names[ids+1] = names[idx];
				id_lengths[ids+1] = lengths[idx];
				break;
			}
		}
		ids+=2;
	}
	// for(auto i : unid_id){
	// 	printf("%d\t%d\n",i.first,i.second);
	// 	printf("%d\t%s\n",i.first,id_names[i.second].c_str());
	// 	printf("%d\t%d\n",i.first,id_lengths[i.second]);

	// }
	// uint64_t overall_total_length = pl->tot_len;
	// uint64_t support_count[unid_id.size()] = {0};
	// memset(support_count,0,sizeof(support_count));
	// uint64_t unsupport_count[unid_id.size()] = {0};
	// memset(unsupport_count,0,sizeof(unsupport_count));
	std::map<uint32_t, std::vector<uint32_t>> unsupported_position;
	std::map<uint32_t, std::vector<uint32_t>> supported_position;


    ofstream outFileConnections;
	outFileConnections.open(string("./") + string("hic_connection_in_haps.txt"), ofstream::out | ofstream::trunc);
	for(uint64_t j = 0; j < std::min(map1->size(), map2->size()); j++){
		

		// if((*map1)[j].unit_id == 65535 && (*map2)[j].unit_id != 65535 && unid_id.find((*map2)[j].unit_id)!=unid_id.end()) {
		// // 		supported_hic++;
		// 		supported_position[unid_id[(*map2)[j].unit_id]].push_back((*map2)[j].pos);
		// 		support_count[unid_id[(*map2)[j].unit_id]]++;
		// }else if((*map2)[j].unit_id == 65535 && (*map1)[j].unit_id != 65535 && unid_id.find((*map1)[j].unit_id)!=unid_id.end()) {
		// // 		supported_hic++;
		// 		supported_position[unid_id[(*map1)[j].unit_id]].push_back((*map1)[j].pos);
		// 		support_count[unid_id[(*map1)[j].unit_id]]++;
		// }
		
		
		if((*map1)[j].unit_id != 65535 && (*map2)[j].unit_id != 65535 && unid_id.find((*map1)[j].unit_id)!=unid_id.end() && unid_id.find((*map2)[j].unit_id)!=unid_id.end()
		&& (unid_id[(*map1)[j].unit_id]>>1) == (unid_id[(*map2)[j].unit_id]>>1) ){
		// if((*map1)[j].unit_id != 65535 && (*map2)[j].unit_id != 65535){
		// printf("%d, %d\n", (*map1)[j].unit_id, (*map2)[j].unit_id);
		// printf("%d, %d\n", (*map1)[j].pos, (*map2)[j].pos);
			outFileConnections << id_names[unid_id[(*map1)[j].unit_id]] << "\t" << (*map1)[j].pos << "\t";
			outFileConnections << id_names[unid_id[(*map2)[j].unit_id]] << "\t" << (*map2)[j].pos << "\n";
			if(unid_id[(*map1)[j].unit_id] == unid_id[(*map2)[j].unit_id]){
				// if((*map1)[j].unit_id >= aux.record_num/2){
					// support_count[(*map1)[j].unit_id - aux.record_num/2]++;
				// }else{
				// support_count[unid_id[(*map1)[j].unit_id]]++;
				// support_count[unid_id[(*map2)[j].unit_id]]++;
				// }
				supported_position[unid_id[(*map1)[j].unit_id]].push_back((*map1)[j].pos);
				supported_position[unid_id[(*map2)[j].unit_id]].push_back((*map2)[j].pos);
			}else{
				// unsupport_count[unid_id[(*map1)[j].unit_id]]++;
				// unsupport_count[unid_id[(*map2)[j].unit_id]]++;
				unsupported_position[unid_id[(*map1)[j].unit_id]].push_back((*map1)[j].pos);
				unsupported_position[unid_id[(*map2)[j].unit_id]].push_back((*map2)[j].pos);
			}
		}
	}
	outFileConnections.close();
	
	// uint32_t switch_error_count = 0;
	// uint64_t switch_error_length = 0;
	// uint32_t hamming_distance = 0;

	uint32_t overall_total_vars = 0;
	uint32_t overall_hamming_distance = 0;
	uint32_t overall_switch_error_count = 0;
	// uint64_t overall_switch_error_length = 0;
	uint32_t overall_potiential_switch_count = 0;

	for(auto i: id_names){
		string name = i.second;
		// uint64_t total_length = id_lengths[i.first];
		uint32_t total_vars = 0;
		uint32_t switch_error_count = 0;
		// uint64_t switch_error_length = 0;
		uint32_t hamming_distance = 0;
		// uint32_t potiential_switch_count = 1;
		vector<uint32_t> sup_pos;
		if(supported_position.find(i.first)!=supported_position.end()){
			sup_pos = supported_position[i.first];
		}
		vector<uint32_t> unsup_pos;
		if(unsupported_position.find(i.first)!=unsupported_position.end()){
			unsup_pos = unsupported_position[i.first];
		}

		vector<pair<uint32_t,bool>> all_counting;
		for(auto p: unsup_pos){
			all_counting.push_back(make_pair(p, true));
		}
		for(auto p: sup_pos){
			all_counting.push_back(make_pair(p, false));
		}
		sort(all_counting.begin(), all_counting.end(), [ ]( const auto& lhs, const auto& rhs )
        {
        return lhs.first < rhs.first;
        });
		vector<pair<uint32_t, int>> unsupport_counting;
		uint32_t cur_pos = 0;
		int cur_count = 0;
		for(auto p: all_counting){
			if(cur_pos != p.first){
				unsupport_counting.push_back(make_pair(cur_pos,cur_count));
				cur_pos = p.first;
				cur_count = 0;
			}
			if(p.second){
				cur_count++;
			}else{
				cur_count--;
			}
		}
		unsupport_counting.push_back(make_pair(cur_pos,cur_count));
		uint32_t unsup_start = 0;
		int kmer_count = 0;
		bool init = false;
		vector<pair<uint32_t,int>> condensed_counting;
		uint idx = 0;
		while(idx < unsupport_counting.size()){
			if(!init){
					init = true;
					unsup_start = unsupport_counting[idx].first;
					kmer_count = unsupport_counting[idx].second;
					idx++;
			}else{
				while(idx < unsupport_counting.size() && unsup_start >= unsupport_counting[idx].first - 150){
					kmer_count += unsupport_counting[idx].second;
					idx++;
				}
				condensed_counting.push_back(make_pair(unsup_start, kmer_count));
				init = false;
			}
		}
		if(condensed_counting.size() > 0){
			total_vars = condensed_counting.size() - 1;
		}
		// vector<uint32_t> switch_pos_vec;
		idx = 0;
		while(idx < condensed_counting.size()){
			if(condensed_counting[idx].second>0){
				uint32_t first_start = idx;
				uint32_t cur_start = idx;
				while(idx < condensed_counting.size() && cur_start >= idx-5){
					if(condensed_counting[idx].second>=0){
						cur_start = idx;
					}
					idx++;
				}
				if(first_start < cur_start-1){
					for(uint i = first_start; i <= cur_start; i++){
						condensed_counting[i].second = 1;
					}
				}else{
					for(uint i = first_start; i <= cur_start; i++){
						condensed_counting[i].second = -1;
					}
				}
			}else{
				idx++;
			}
		}

		if(condensed_counting.size()>0){
			for(uint p = 0; p < condensed_counting.size()-1; p++){
				if(condensed_counting[p].second > 0){
					// printf("Switch at %"PRIu64", %d\n", condensed_counting[p].first, condensed_counting[p].second);
					hamming_distance += 1;
				}else{
					// printf("Unswitch at %"PRIu64", %d\n", condensed_counting[p].first, condensed_counting[p].second);
				}
				if( (condensed_counting[p].second>0) ^ (condensed_counting[p+1].second>0)){
					if( condensed_counting[p].first){
						switch_error_count += 1;
					}
					// switch_pos_vec.push_back(condensed_counting[p+1].first);
				}
			// printf("%d\n", condensed_counting[p].first);
			}
				if(condensed_counting[condensed_counting.size()-1].second > 0){
					// printf("Switch at %"PRIu64", %d\n", condensed_counting[condensed_counting.size()-1].first, condensed_counting[condensed_counting.size()-1].second);
					hamming_distance += 1;
				}else{
					// printf("Unswitch at %"PRIu64", %d\n", condensed_counting[condensed_counting.size()-1].first, condensed_counting[condensed_counting.size()-1].second);
				}
		}
		// printf("%d\n", condensed_counting[condensed_counting.size()-1].first);
		// for(auto i:condensed_counting){
		// 	if(i.second>0){
		// 		printf("%"PRIu64"\n", i.first);
		// 	}
		// }
		// if(switch_pos_vec.size()>=2){
		// 	for(int p = 0; p < switch_pos_vec.size()-1; p+=2){
		// 		switch_error_length+=switch_pos_vec[p+1] - switch_pos_vec[p];	
		// 	}
		// }
		// printf("\n%s\t%d\n",name.c_str(), id_lengths[i.first]);
		// printf("Unswitched: %"PRIu64"\n", total_vars - hamming_distance);
		// printf("Switched: %"PRIu64"\n", hamming_distance);
		// printf("Hamming distance rate: %.4f%%\n", ((double) hamming_distance)/ (total_vars) * 100);
		// printf("Switch error number: %"PRIu64", rate: %.4f%%\n", switch_error_count, ((double) switch_error_count)/ (condensed_counting.size()-1) * 100);
		// printf("Switch error length: %"PRIu64", rate: %.4f%%\n", switch_error_length, ((double) switch_error_length)/ total_length * 100);
		overall_total_vars += total_vars;
		overall_hamming_distance += hamming_distance;
		overall_switch_error_count += switch_error_count;
		// overall_switch_error_length += switch_error_length;
		overall_potiential_switch_count += (condensed_counting.size()-1);
		// if(unsup_pos_len.size()>0){
		// 	vector<pair<uint32_t,uint32_t>> merged;
		// 	uint32_t cur_starting = unsup_pos_len[0].first;
		// 	uint32_t cur_ending = unsup_pos_len[0].first + unsup_pos_len[0].second;
		// 	for(int idx = 1; idx < unsup_pos_len.size(); idx++){
		// 		if(cur_ending+150 > unsup_pos_len[idx].first){
		// 			cur_ending = unsup_pos_len[idx].first+unsup_pos_len[idx].second;
		// 		}else{
		// 			merged.push_back(make_pair(cur_starting, cur_ending-cur_starting));
		// 			cur_starting = unsup_pos_len[idx].first;
		// 			cur_ending = unsup_pos_len[idx].first + unsup_pos_len[idx].second;
		// 		}
		// 	}
		// 	merged.push_back(make_pair(cur_starting, cur_ending-cur_starting));
		// 	for(auto switches: merged){
		// 		switch_error_count++;
		// 		switch_error_length+=switches.second;
		// 		printf("%s : switch at: %"PRIu32", length: %"PRIu32" \n", id_names[i.first].c_str(), switches.first, switches.second);
		// 	}
		// }
	}
	// printf("Overall total length: %"PRIu64"\n",overall_total_length);
	fprintf(stderr, "[M::Result]\n");
	printf("%s\t%s\n","Hamming Distance", "Switch Error Rate");
	printf("%.4f%%\t%.4f%%\n", ((double) overall_hamming_distance)/ (overall_total_vars) * 100, ((double) overall_switch_error_count)/ overall_potiential_switch_count * 100);
	// printf("\nOverall Unswitched: %"PRIu64"\n", overall_total_vars - overall_hamming_distance);
	// printf("Overall Switched: %"PRIu64"\n", overall_hamming_distance);
	// printf("Overall Hamming distance rate: %.4f%%\n", ((double) overall_hamming_distance)/ (overall_total_vars) * 100);
	// printf("Overall Switch error number: %"PRIu64", rate: %.4f%%\n", overall_switch_error_count, ((double) overall_switch_error_count)/ overall_potiential_switch_count * 100);
	// printf("Overall Switch error length: %"PRIu64", rate: %.4f%%\n", overall_switch_error_length, ((double) overall_switch_error_length)/ overall_total_length * 100);
	// for(int i = 0; i < ids; i++){
	// 	printf("%s : sup: %"PRIu64", unsup: %"PRIu64", rate: %.4f%%\n", id_names[i].c_str(), support_count[i], unsupport_count[i],
	// 	                                      ((double) support_count[i])/ (support_count[i] + unsupport_count[i]) * 100 );
	// }
	// printf("Total support pairs: %"PRIu64"\n", supported_hic);
	// printf("Total unsupport pairs: %"PRIu64"\n", unsupported_hic);
	// printf("Total hamming distance rate: %.4f%%\n", ((double) unsupported_hic)/ (supported_hic + unsupported_hic) * 100);
	// printf("Switch error number: %"PRIu64", length: %"PRIu64", length rate: %.4f%%\n", switch_error_count, switch_error_length, (((double) switch_error_length)/ ((double) total_length)) * 100);
}

int main_switch_error(yak_copt_t opt, string hic_file1,string hic_file2, string hap_file1, string hap_file2, std::map<std::string,std::string> contig2_1_map)
{
	pldat_t *h;
	bseq_file_t* hic_fn1 = bseq_open(hic_file1.c_str());
	bseq_file_t* hic_fn2 = bseq_open(hic_file2.c_str());
	if (hic_fn1 == 0 || hic_fn2 == 0) {
		fprintf(stderr, "ERROR: fail to open hic files\n");
		exit(1);
	}
	h = yak_counting_switch_error(hap_file1,hap_file2, &opt);

	get_switch_error(h, hic_fn1, hic_fn2, contig2_1_map);
	
	yak_ch_destroy(h->h);
	free(h);
	return 0;
}
