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
#include <iostream>
#include "kseq.h" // FASTA/Q parser

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define CHUNK_SIZE 200000000
using namespace std;

typedef struct {
	int appear_counting;
	int total_counting;
} tb_cnt_t;

typedef struct {
	int max;
	uint32_t *s;
} tb_buf_t;

typedef struct {
	int k, n_threads, print_diff;
	double ratio_thres;
	bseq_file_t *fp;
	const yak_ch_t *ch;
	tb_buf_t *buf;
	uint64_t appear_counting;
	uint64_t total_counting;
} tb_shared_t;

typedef struct {
	int n_seq;
	tb_shared_t *aux;
	bseq1_t *seq;
	tb_cnt_t *cnt;
} tb_step_t;

static void tb_worker(void *_data, long k, int tid)
{
	tb_step_t *t = (tb_step_t*)_data;
	tb_shared_t *aux = t->aux;
	bseq1_t *s = &t->seq[k];
	tb_buf_t *b = &aux->buf[tid];
	uint64_t x[4], mask;
	int i, l, shift;
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
		int flag, c = seq_nt4_table[(uint8_t)s->seq[i]];
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
				if(yak_ch_get(aux->ch, y)!=-1){
					t->cnt[k].appear_counting++;
				}
				t->cnt[k].total_counting++;
			}
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0;
	}
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
			s->cnt = (tb_cnt_t*)calloc(s->n_seq, sizeof(tb_cnt_t));
			fprintf(stderr, "[M::%s] read %d sequences\n", __func__, s->n_seq);
			return s;
		} else free(s);
	} else if (step == 1) {
		int i;
		tb_step_t *s = (tb_step_t*)_data;
		kt_for(aux->n_threads, tb_worker, s, s->n_seq);
		for (i = 0; i < s->n_seq; ++i) {
			aux->appear_counting += s->cnt[i].appear_counting;
			aux->total_counting += s->cnt[i].total_counting;
			free(s->seq[i].name); free(s->seq[i].seq); free(s->seq[i].qual); free(s->seq[i].comment);
		}
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences, current completeness: %.4f%%;\n", __func__,
			yak_realtime(), yak_cputime() / yak_realtime(), s->n_seq, ((double)aux->appear_counting) / aux->total_counting * 100);
		free(s->seq); free(s->cnt); free(s);
	}
	return 0;
}

double main_completeness(yak_ch_t *ch, int n_threads, string file1, string file2)
{
	ketopt_t o = KETOPT_INIT;
	int i, c, min_cnt = 2, mid_cnt = 5;
	tb_shared_t aux;

	memset(&aux, 0, sizeof(tb_shared_t));
	aux.n_threads = n_threads, aux.print_diff = 0;
	aux.ratio_thres = 0.33;



	aux.k = ch->k;
	aux.ch = ch;

	aux.fp = bseq_open(file1.c_str());
	if (aux.fp == 0) {
		fprintf(stderr, "ERROR: fail to open file '%s'\n", file1.c_str());
		exit(1);
	}
	aux.buf = (tb_buf_t*)calloc(aux.n_threads, sizeof(tb_buf_t));
	kt_pipeline(2, tb_pipeline, &aux, 2);
	bseq_close(aux.fp);
	for (i = 0; i < aux.n_threads; ++i) free(aux.buf[i].s);
	free(aux.buf);

	aux.fp = bseq_open(file2.c_str());
	if (aux.fp == 0) {
		fprintf(stderr, "ERROR: fail to open file '%s'\n", file2.c_str());
		exit(1);
	}
	aux.buf = (tb_buf_t*)calloc(aux.n_threads, sizeof(tb_buf_t));
	kt_pipeline(2, tb_pipeline, &aux, 2);
	bseq_close(aux.fp);
	for (i = 0; i < aux.n_threads; ++i) free(aux.buf[i].s);
	free(aux.buf);

	return ((double)aux.appear_counting) / aux.total_counting * 100;
}