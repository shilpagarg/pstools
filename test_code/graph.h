#ifndef __GRAPH__
#define __GRAPH__

#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "kalloc.h"
#include "ketopt.h"
#include "khash.h"
#include "kseq.h"
#include "ksort.h"
#include "kstring.h"
#include "kvec.h"
#define GFA_MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define Get_NAME_LENGTH(R_INF, ID) ((R_INF).name_index[(ID)+1] - (R_INF).name_index[(ID)])
typedef khint32_t khint_t;
KHASH_MAP_INIT_STR(seg, uint32_t)
typedef khash_t(seg) seghash_t;
KSTREAM_INIT(gzFile, gzread, 65536)
#define VERBOSE 0
// #define gfa_arc_key(a) ((a).ul)
// KRADIX_SORT_INIT(arc, asg_arc_t, gfa_arc_key, 8)
// #define generic_key(x) (x)
// KRADIX_SORT_INIT(gfa64, uint64_t, generic_key, 8)

#define asg_arc_len(arc) ((uint32_t)(arc).ul)
#define asg_arc_n(g, v) ((uint32_t)(g)->idx[(v)])
#define asg_arc_a(g, v) (&(g)->arc[(g)->idx[(v)]>>32])
#define asg_n_vtx(g) ((g)->n_seq << 1)
#define asg_arc_key(a) ((a).ul)
#define gfa_arc_head(a) ((uint32_t)((a).ul>>32))
#define gfa_arc_tail(a) ((a).v)

int gfa_verbose;

typedef struct {
	uint64_t ul;
	uint32_t v;
	uint32_t ol:31, del:1;
	uint8_t strong;
	uint8_t el;
	uint8_t no_l_indel;

	// uint64_t v_lv; // higher 32 bits: vertex_id; lower 32 bits: lv; packed together for sorting
	// uint32_t w;
	// int32_t rank;
	// int32_t ov, ow;
	// uint64_t link_id:61, strong:1, del:1, comp:1; // link_id: a pair of dual arcs are supposed to have the same link_id
} asg_arc_t;
KRADIX_SORT_INIT(asg, asg_arc_t, asg_arc_key, 8)

typedef struct {
	uint32_t len:31, circ:1; // len: length of the unitig; circ: circular if non-zero
	uint32_t start, end; // start: starting vertex in the string graph; end: ending vertex
	uint32_t m, n; // number of reads
	uint64_t *a; // list of reads
	char *s; // unitig sequence is not null
} ma_utg_t;

typedef struct {
	uint32_t len:31, del:1;
	char *name, *seq;
	uint8_t c;
} asg_seq_t;
typedef struct {
	uint32_t m_arc, n_arc:31, is_srt:1;  // m_arc is the size of edges; n_arc is the size of edges; 
	asg_arc_t *arc;  // gives us arcs coming out of the node
	uint32_t m_seq, n_seq:31, is_symm:1;  // m_seq is the size of nodes; n_seq is the number of nodes; 
	uint32_t r_seq;
	void *h_names;
	void *h_snames;

	asg_seq_t *seq;
	uint64_t *idx;  // "utg028945l	-" is the 12th line. 2*(12-1) + 1(sign == -) -> 23
					// idx(23) gives us the index into arc for "utg028945l	-"

// 					p ((asg_arc_t *)&g->arc[g->idx[2]>>32])[1]
// $77 = {ul = 12884901888, v = 21, ol = 14899, del = 0, strong = 0 '\000', 
//   el = 0 '\000', no_l_indel = 0 '\000'}

	uint8_t* seq_vis;

	uint32_t n_F_seq;
	ma_utg_t* F_seq;
} asg_t;

// void gfa_destory(asg_t* g){
// 	// for(int i = 0; i < g->n_arc; i++){
// 	// 	free(&g->arc[i]);
// 	// }
// 	free(g->arc);
// 	free(g->idx);
// 	for(int i = 0; i < g->n_seq; i++){
// 		free(g->seq[i].name);
// 		free(g->seq[i].seq);
// 		// free(&g->seq[i]);
// 	}
// 	free(g->seq);
// 	// kh_destory(g->h_names);
// 	// kh_destory(g->h_snames);
// 	free(g->h_names);
// 	free(g->h_snames);
// 	free(g);
// }

asg_t *gfa_init(void)
{
	asg_t *g;
	g = (asg_t*)calloc(1, sizeof(asg_t));
	g->h_names = kh_init(seg);
	g->h_snames = kh_init(seg);
	return g;
}

char *gfa_strdup(const char *src)
{
	int32_t len;
	char *dst;
	len = strlen(src);
	GFA_MALLOC(dst, len + 1);
	memcpy(dst, src, len + 1);
	return dst;//e.g. "utg000001l"
}

int32_t gfa_add_seg(asg_t *g, const char *name)
{
	khint_t k;
	int absent;
	seghash_t *h = (seghash_t*)g->h_names;//p ((seghash_t*)g->h_names)->keys[k]
											//$18 = (kh_cstr_t) 0x55555558b3c0 "utg000001l"

	k = kh_put(seg, h, name, &absent);
	if (absent) {
		asg_seq_t *s;
		if (g->n_seq == g->m_seq) {
			uint32_t old_m = g->m_seq;
			g->m_seq = g->m_seq? g->m_seq<<1 : 16;
			g->seq = (asg_seq_t*)realloc(g->seq, g->m_seq * sizeof(asg_seq_t));
			memset(&g->seq[old_m], 0, (g->m_seq - old_m) * sizeof(asg_seq_t));
		}
		s = &g->seq[g->n_seq++];
		kh_key(h, k) = s->name = gfa_strdup(name);//kh_key(h, k) == h->keys[k]
		s->del = s->len = 0;
		// s->snid = s->soff = s->rank = -1;
		kh_val(h, k) = g->n_seq - 1;
	}
	return kh_val(h, k);
}

int gfa_parse_S(asg_t *g, char *s)
{
	int i, is_ok = 0;
	char *p, *q, *seg = 0, *seq = 0, *rest = 0;
	uint32_t sid, len = 0;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) seg = q;
			else if (i == 1) {
				seq = q[0] == '*'? 0 : gfa_strdup(q);
				is_ok = 1, rest = c? p + 1 : 0;
				break;
			}
			++i, q = p + 1;
			if (c == 0) break;
		}
	}
    if (is_ok) { // all mandatory fields read
		int l_aux, m_aux = 0, LN = -1;
		uint8_t *aux = 0, *s_LN = 0;
		asg_seq_t *s;
		if (s_LN && s_LN[0] == 'i') {
			LN = *(int32_t*)(s_LN + 1);
		}
		if (seq == 0) {
			if (LN >= 0) len = LN;
		} else len = strlen(seq);
		if (LN >= 0 && len != LN && gfa_verbose >= 2)
			fprintf(stderr, "[W] for segment '%s', LN:i:%d tag is different from sequence length %d\n", seg, LN, len);
		sid = gfa_add_seg(g, seg);
		s = &g->seq[sid];
		s->len = len, s->seq = seq;
	} else return -1;
	return 0;
}

asg_arc_t *gfa_add_arc1(asg_t *g, uint32_t v, uint32_t w, int32_t ov, int32_t ow, int64_t link_id, int comp)
{
	asg_arc_t *a;
	if (g->m_arc == g->n_arc) {
		uint64_t old_m = g->m_arc;
		g->m_arc = g->m_arc? g->m_arc<<1 : 16;
		g->arc = (asg_arc_t*)realloc(g->arc, g->m_arc * sizeof(asg_arc_t));
		memset(&g->arc[old_m], 0, (g->m_arc - old_m) * sizeof(asg_arc_t));
	}
	a = &g->arc[g->n_arc++];
	a->ul = (uint64_t)v << 32;
	a->v = w, a->ol = ov, a->ol = ow;
	a->del = a->strong = 0;
	return a;
}

int gfa_parse_L(asg_t *g, char *s)
{
	int i, oriv = -1, oriw = -1, is_ok = 0;
	char *p, *q, *segv = 0, *segw = 0, *rest = 0;
	int32_t ov = INT32_MAX, ow = INT32_MAX;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) {
				segv = q;
			} else if (i == 1) {
				if (*q != '+' && *q != '-') return -2;
				oriv = (*q != '+');
			} else if (i == 2) {
				segw = q;
			} else if (i == 3) {
				if (*q != '+' && *q != '-') return -2;
				oriw = (*q != '+');
			} else if (i == 4) {
				if (*q == ':') {
					ov = INT32_MAX;
					ow = isdigit(*(q+1))? strtol(q+1, &q, 10) : INT32_MAX;
				} else if (isdigit(*q)) {
					char *r;
					ov = strtol(q, &r, 10);
					if (isupper(*r)) { // CIGAR
						ov = ow = 0;
						do {
							long l;
							l = strtol(q, &q, 10);
							if (*q == 'M' || *q == 'D' || *q == 'N') ov += l;
							if (*q == 'M' || *q == 'I' || *q == 'S') ow += l;
							++q;
						} while (isdigit(*q));
					} else if (*r == ':') { // overlap lengths
						ow = isdigit(*(r+1))? strtol(r+1, &r, 10) : INT32_MAX;
					} else break;
				} else break;
				is_ok = 1, rest = c? p + 1 : 0;
				break;
			}
			++i, q = p + 1;
			if (c == 0) break;
		}
	}
	if (is_ok) {
		uint32_t v, w;
		int l_aux, m_aux = 0;
		uint8_t *aux = 0;
		asg_arc_t *arc;
		v = gfa_add_seg(g, segv) << 1 | oriv;
		w = gfa_add_seg(g, segw) << 1 | oriw;
		arc = gfa_add_arc1(g, v, w, ov, ow, -1, 0);
	} else return -1;
	return 0;
}

uint32_t gfa_fix_no_seg(asg_t *g)
{
	uint32_t i, n_err = 0;
	for (i = 0; i < g->n_seq; ++i) {
		asg_seq_t *s = &g->seq[i];
		if (s->len == 0) {
			++n_err, s->del = 1;
			if (gfa_verbose >= 2)
				fprintf(stderr, "[W] segment '%s' is used on an L-line but not defined on an S-line\n", s->name);
		}
	}
	return n_err;
}

void asg_arc_sort(asg_t *g)
{
	radix_sort_asg(g->arc, g->arc + g->n_arc);
}

uint64_t *gfa_arc_index_core(size_t max_seq, size_t n, const asg_arc_t *a)
{
	size_t i, last;
	uint64_t *idx;
	idx = (uint64_t*)calloc(max_seq * 2, 8);
	/**
	 * ul: |____________31__________|__________1___________|______________32_____________|
						 qns            direction of overlap       length of this node (not overlap length)
	**/
    ///so if we use high 32-bit, we store the index of each qn with two direction
	for (i = 1, last = 0; i <= n; ++i)
		if (i == n || gfa_arc_head(a[i-1]) != gfa_arc_head(a[i]))
			idx[gfa_arc_head(a[i-1])] = (uint64_t)last<<32 | (i - last), last = i;
	return idx;
}

void gfa_arc_index(asg_t *g)
{
	if (g->idx) free(g->idx);
	g->idx = gfa_arc_index_core(g->n_seq, g->n_arc, g->arc);
}

uint32_t gfa_fix_semi_arc(asg_t *g)
{
	uint32_t n_err = 0, v, n_vtx = asg_n_vtx(g);
	int i, j;
	for (v = 0; v < n_vtx; ++v) {
		int nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			if (!av[i].del && (av[i].ol == INT32_MAX || av[i].ol == INT32_MAX)) { // overlap length is missing
				uint32_t w = av[i].v^1;
				int is_multi = 0, c, jv = -1, nw = asg_arc_n(g, w);
				asg_arc_t *aw = asg_arc_a(g, w);
				for (j = 0, c = 0; j < nw; ++j)
					if (!aw[j].del && aw[j].v == (v^1)) ++c, jv = j;
				if (c == 1) {
					if (av[i].ol != INT32_MAX && aw[jv].ol != INT32_MAX && av[i].ol != aw[jv].ol) is_multi = 1;
					if (av[i].ol != INT32_MAX && aw[jv].ol != INT32_MAX && av[i].ol != aw[jv].ol) is_multi = 1;
				}
				if (c == 1 && !is_multi) {
					if (aw[jv].ol != INT32_MAX) av[i].ol = aw[jv].ol;
					if (aw[jv].ol != INT32_MAX) av[i].ol = aw[jv].ol;
				} else {
					if (gfa_verbose >= 2)
						fprintf(stderr, "[W] can't infer overlap length for %s%c -> %s%c\n",
								g->seq[v>>1].name, "+-"[v&1], g->seq[w>>1].name, "+-"[(w^1)&1]);
					++n_err;
					av[i].del = 1;
				}
			}
		}
	}
	return n_err;
}

void gfa_arc_rm(asg_t *g)
{
	uint32_t e, n;
	for (e = n = 0; e < g->n_arc; ++e) {
		uint32_t u = g->arc[e].ul>>32, v = g->arc[e].v;
		if (!g->arc[e].del && !g->seq[u>>1].del && !g->seq[v>>1].del)
			g->arc[n++] = g->arc[e];
		else {
		}
	}
	if (n < g->n_arc) { // arc index is out of sync
		if (g->idx) free(g->idx);
		g->idx = 0;
	}
	g->n_arc = n;
}

int gfa_arc_is_sorted(const asg_t *g)
{
	uint64_t e;
	for (e = 1; e < g->n_arc; ++e)
		if (g->arc[e-1].ul > g->arc[e].ul)
			break;
	return (e == g->n_arc);
}

void gfa_cleanup(asg_t *g)
{
	gfa_arc_rm(g);
	if (!gfa_arc_is_sorted(g)) {
		asg_arc_sort(g);
		if (g->idx) free(g->idx);
		g->idx = 0;
	}
	if (g->idx == 0) gfa_arc_index(g);
}

uint32_t gfa_fix_symm_add(asg_t *g)
{
	uint32_t n_err = 0, v, n_vtx = asg_n_vtx(g);
	int i;
	for (v = 0; v < n_vtx; ++v) {
		int nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			int j, nw;
			asg_arc_t *aw, *avi = &av[i];
			// if (avi->del || avi->comp) continue;
			if (avi->del) continue;
			nw = asg_arc_n(g, avi->v^1);
			aw = asg_arc_a(g, avi->v^1);
			for (j = 0; j < nw; ++j) {
				asg_arc_t *awj = &aw[j];
				// if (awj->del || awj->comp) continue;
				if (awj->del) continue;
				if (awj->v == (v^1) && awj->ol == avi->ol && awj->ol == avi->ol) { // complement found
					// awj->comp = 1;
					// awj->link_id = avi->link_id;
					break;
				}
			}
			if (j == nw) {
				asg_arc_t *arc_old = g->arc, *arc_new;
				// arc_new = gfa_add_arc1(g, avi->v^1, v^1, avi->ol, avi->ol, avi->link_id, 1);
				arc_new = gfa_add_arc1(g, avi->v^1, v^1, avi->ol, avi->ol, '*', 1);
				if (arc_old != g->arc) av = asg_arc_a(g, v); // g->arc may be reallocated
				// arc_new->rank = av[i].rank;
			}
		}
	}
	if (n_vtx < asg_n_vtx(g)) {
		asg_arc_sort(g);
		gfa_arc_index(g);
	}
	return n_err;
}

void gfa_fix_arc_len(asg_t *g)
{
	uint64_t k;
	for (k = 0; k < g->n_arc; ++k) {
		asg_arc_t *a = &g->arc[k];
		uint32_t v = gfa_arc_head(*a), w = gfa_arc_tail(*a);
		const asg_seq_t *sv = &g->seq[v>>1];
		if (!sv->del && sv->len < a->ol) {
			if (gfa_verbose >= 2)
				fprintf(stderr, "[W] overlap length longer than segment length for '%s': %d > %d\n", sv->name, a->ol, sv->len);
			a->ol = sv->len;
		}
		if (sv->del || g->seq[w>>1].del) {
			a->del = 1;
		} else {
			a->ul |= sv->len - a->ol;
		}
	}
}

void gfa_finalize(asg_t *g)  // TODO
{
	// gfa_fix_no_seg(g);
	asg_arc_sort(g);
	gfa_arc_index(g);
	// gfa_fix_semi_arc(g);
	gfa_fix_symm_add(g);  // doubles the arcs
	gfa_fix_arc_len(g);
	// gfa_fix_utg(g);
	gfa_cleanup(g);  // uses the doubled arcs to create idx
}

// delete multi-arcs
/**
 * remove edges like:   v has two out-edges to w
**/
int asg_arc_del_multi(asg_t *g)
{
	//the number of nodes are number of read times 2
	uint32_t *cnt, n_vtx = g->n_seq * 2, n_multi = 0, v;
	cnt = (uint32_t*)calloc(n_vtx, 4);
	for (v = 0; v < n_vtx; ++v) {
		///out-nodes of v
		asg_arc_t *av = asg_arc_a(g, v);
		int32_t i, nv = asg_arc_n(g, v);
		///if v just have one out-node, there is no muti-edge
		if (nv < 2) continue;
		for (i = nv - 1; i >= 0; --i) ++cnt[av[i].v];
		for (i = nv - 1; i >= 0; --i)
			if (--cnt[av[i].v] != 0)
				av[i].del = 1, ++n_multi;
	}
	free(cnt);
	if (n_multi) gfa_cleanup(g);

    if(VERBOSE >= 1)
    {
        fprintf(stderr, "[M::%s] removed %d multi-arcs\n", __func__, n_multi);
    }
	
	return n_multi;
}

// remove asymmetric arcs: u->v is present, but v'->u' not
int asg_arc_del_asymm(asg_t *g)
{
	uint32_t e, n_asymm = 0;
	///g->n_arc is the number of overlaps 
	for (e = 0; e < g->n_arc; ++e) {
		uint32_t v = g->arc[e].v^1, u = g->arc[e].ul>>32^1;
		uint32_t i, nv = asg_arc_n(g, v);
		asg_arc_t *av = asg_arc_a(g, v);
		for (i = 0; i < nv; ++i)
			if (av[i].v == u) break;
		if (i == nv) g->arc[e].del = 1, ++n_asymm;
	}
	if (n_asymm) gfa_cleanup(g);
    if(VERBOSE >= 1)
    {
	    fprintf(stderr, "[M::%s] removed %d asymmetric arcs\n", __func__, n_asymm);
    }
	return n_asymm;
}

void asg_symm(asg_t *g)
{
	asg_arc_del_multi(g);
	asg_arc_del_asymm(g);
	g->is_symm = 1;
}

asg_t *gfa_read(const char *fn)
{
	gzFile fp;
	asg_t *g;
	kstring_t s = {0,0,0}, fa_seq = {0,0,0};
	kstream_t *ks;
	int dret, is_fa = 0;
	asg_seq_t *fa_seg = 0;
	uint64_t lineno = 0;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	g = gfa_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
		int ret = 0;
		++lineno;
		// if (s.l > 0 && s.s[0] == '>') { // FASTA header
		// 	is_fa = 1;
		// 	if (fa_seg) gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s);
		// 	fa_seg = gfa_parse_fa_hdr(g, s.s);
		// 	fa_seq.l = 0;
		// } else if (is_fa) { // FASTA mode
		// 	if (s.l >= 3 && s.s[1] == '\t') { // likely a GFA line
		// 		gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s); // finalize fa_seg
		// 		fa_seg = 0;
		// 		is_fa = 0;
		// 	} else kputsn(s.s, s.l, &fa_seq); // likely a FASTA sequence line
		// }
		if (is_fa) continue;
		if (s.l < 3 || s.s[1] != '\t') continue; // empty line
		if (s.s[0] == 'S') ret = gfa_parse_S(g, s.s);
		else if (s.s[0] == 'L') ret = gfa_parse_L(g, s.s);
		// else if (s.s[0] == 'A') ret = gfa_parse_A(g, s.s);
		if (ret < 0 && gfa_verbose >= 1)
			fprintf(stderr, "[E] invalid %c-line at line %ld (error code %d)\n", s.s[0], (long)lineno, ret);
	}
	// if (is_fa && fa_seg) gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s);
	free(fa_seq.s);
	free(s.s);
	gfa_finalize(g);
	ks_destroy(ks);
	gzclose(fp);
	return g;
}

#endif