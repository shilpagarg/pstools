int main_switch_error(int argc, char *argv[])
{
	pldat_t *h;
	int c;
	yak_copt_t opt;
	ketopt_t o = KETOPT_INIT;
	yak_copt_init(&opt);
	opt.pre = YAK_POS_ID_BITS;
	opt.n_thread = 32;
	while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:b:H:o:", 0)) >= 0) {
		if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 'p') opt.pre = atoi(o.arg);
		else if (c == 'K') opt.chunk_size = atoi(o.arg);
		else if (c == 't') opt.n_thread = atoi(o.arg);
		else if (c == 'b') opt.bf_shift = atoi(o.arg);
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: pstools eval switch_error [options] <hap1.fa> <hap2.fa> <hic.R1.fastq.gz> <hic.R2.fastq.gz> \n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", opt.pre);
		fprintf(stderr, "  -b INT     set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", opt.bf_shift);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
		fprintf(stderr, "  -K INT     chunk size [100m]\n");
		fprintf(stderr, "Note: -b37 is recommended for human reads\n");
		return 1;
	}
	if (opt.pre < YAK_POS_ID_BITS) {
		fprintf(stderr, "ERROR: -p should be at least %d\n", YAK_POS_ID_BITS);
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

	bseq_file_t* hic_fn1 = bseq_open(argv[o.ind+2]);
	bseq_file_t* hic_fn2 = bseq_open(argv[o.ind+3]);
	if (hic_fn1 == 0 || hic_fn2 == 0) {
		fprintf(stderr, "ERROR: fail to open hic files\n");
		exit(1);
	}


	h = yak_counting_switch_error(argv[o.ind],argv[o.ind+1], &opt, 0);

	get_switch_error(h, hic_fn1, hic_fn2);


	free(h->h);
	free(h);
	return 0;
}



int main_completeness(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int i, c, min_cnt = 2, mid_cnt = 5;
	yak_ch_t *ch;
	tb_shared_t aux;

	memset(&aux, 0, sizeof(tb_shared_t));
	aux.n_threads = 8, aux.print_diff = 0;
	aux.ratio_thres = 0.33;
	while ((c = ketopt(&o, argc, argv, 1, "t:", 0)) >= 0) {
		if (c == 't') aux.n_threads = atoi(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: pstools eval completeness [options] <kmer.yak> <prediction.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t INT     number of threads [%d]\n", aux.n_threads);
		return 1;
	}

	ch = yak_ch_restore(argv[o.ind]);

	aux.k = ch->k;
	aux.fp = bseq_open(argv[o.ind+1]);
	if (aux.fp == 0) {
		fprintf(stderr, "ERROR: fail to open file '%s'\n", argv[o.ind+1]);
		exit(1);
	}
	aux.ch = ch;
	aux.buf = (tb_buf_t*)calloc(aux.n_threads, sizeof(tb_buf_t));
	kt_pipeline(2, tb_pipeline, &aux, 2);
	bseq_close(aux.fp);
	yak_ch_destroy(ch);
	for (i = 0; i < aux.n_threads; ++i) free(aux.buf[i].s);
	free(aux.buf);
	return 0;
}


int main_qv_check(int argc, char *argv[])
{
	yak_qopt_t opt;
	yak_ch_t *ch = 0;
	ketopt_t o = KETOPT_INIT;
	int64_t cnt[YAK_N_COUNTS], hist[YAK_N_COUNTS];
	int c, i, kmer;
	yak_qstat_t qs;

	yak_qopt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "K:t:l:f:pe:", 0)) >= 0) {
		if (c == 'K') opt.chunk_size = mm_parse_num(o.arg);
		else if (c == 'l') opt.min_len = mm_parse_num(o.arg);
		else if (c == 'f') opt.min_frac = atof(o.arg);
		else if (c == 't') opt.n_threads = atoi(o.arg);
		else if (c == 'p') opt.print_each = 1;
		else if (c == 'e') opt.fpr = atof(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: pstools eval qv [options] <kmers.yak> <prediction.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l NUM      min sequence length [%d]\n", opt.min_len);
		fprintf(stderr, "  -f FLOAT    min k-mer fraction [%g]\n", opt.min_frac);
		fprintf(stderr, "  -e FLOAT    false positive rate [%g]\n", opt.fpr);
		fprintf(stderr, "  -p          print QV for each sequence\n");
		fprintf(stderr, "  -t INT      number of threads [%d]\n", opt.n_threads);
		fprintf(stderr, "  -K NUM      batch size [1g]\n");
		return 1;
	}

	ch = yak_ch_restore(argv[o.ind]);
	assert(ch);
	kmer = ch->k;
	yak_ch_hist(ch, hist, opt.n_threads);
	printf("CC\tCT  kmer_occurrence    short_read_kmer_count  raw_input_kmer_count  adjusted_input_kmer_count\n");
	printf("CC\tFR  fpr_lower_bound    fpr_upper_bound\n");
	printf("CC\tER  total_input_kmers  adjusted_error_kmers\n");
	printf("CC\tCV  coverage\n");
	printf("CC\tQV  raw_quality_value  adjusted_quality_value\n");
	printf("CC\n");
	yak_qv(&opt, argv[o.ind+1], ch, cnt);
	yak_qv_solve(hist, cnt, kmer, opt.fpr, &qs);
	for (i = (1<<YAK_COUNTER_BITS) - 1; i >= 0; --i)
		printf("CT\t%d\t%ld\t%ld\t%.3f\n", i, (long)hist[i], (long)cnt[i], qs.adj_cnt[i]);
	printf("FR\t%.3g\t%.3g\n", qs.fpr_lower, qs.fpr_upper);
	printf("ER\t%ld\t%.3f\n", (long)qs.tot, qs.err);
	printf("CV\t%.3f\n", qs.cov);
	printf("QV\t%.3f\t%.3f\n", qs.qv_raw, qs.qv);
	yak_ch_destroy(ch);
	return 0;
}







int main_eval(int argc, char *argv[]){
	int ret = 0;
	yak_copt_t opt;
	yak_qopt_t qopt;
	yak_ch_t *haplo_h;
	yak_ch_t *hic_h;
	ketopt_t o = KETOPT_INIT;
	int64_t cnt[YAK_N_COUNTS], hist[YAK_N_COUNTS];
	int c, i, kmer;
	yak_qstat_t qs;
	double completeness;

	opt.pre = YAK_COUNTER_BITS;
	yak_copt_init(&opt);
	yak_qopt_init(&qopt);
	while ((c = ketopt(&o, argc, argv, 1, "k:K:t:", 0)) >= 0) {
		if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 'K') opt.chunk_size = atoi(o.arg);
		else if (c == 't') opt.n_thread = atoi(o.arg);
	}
	qopt.n_threads = opt.n_threads;
	if (argc - o.ind < 4) {
		fprintf(stderr, "Usage: pstools eval [options] <hap1.fa> <hap2.fa> <hic.R1.fastq.gz> <hic.R2.fastq.gz>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
		fprintf(stderr, "  -K INT     chunk size [100m]\n");
		return 1;
	}
	if (opt.k >= 64) {
		fprintf(stderr, "ERROR: -k must be smaller than 64\n");
		return 1;
	} else if (opt.k >= 32) {
		fprintf(stderr, "WARNING: counts are inexact if -k is greater than 31\n");
	}
	
	char* hap1_file = argv[o.ind];
	char* hap2_file = argv[o.ind+1];
	char* hic1_file = argv[o.ind+2];
	char* hic2_file = argv[o.ind+3];
	// bseq_file_t* hic_fn1 = bseq_open(hic1_file);
	// bseq_file_t* hic_fn2 = bseq_open(hic2_file);

	// if (hic_fn1 == 0 || hic_fn2 == 0) {
	// 	fprintf(stderr, "ERROR: fail to open hic files\n");
	// 	exit(1);
	// }else{
	// 	bseq_close(hic_fn1);
	// 	bseq_close(hic_fn2);
	// }

	hic_h = yak_count(hic1_file, &opt, 0);
	hic_h = yak_count(hic2_file, &opt, hic_h);
	kmer = hic_h->k;
	yak_ch_hist(hic_h, hist, opt.n_threads);

	yak_qv(&qopt, argv[o.ind+1], hic_h, cnt);
	yak_qv_solve(hist, cnt, kmer, qopt.fpr, &qs);
	yak_ch_destroy(hic_h);

	haplo_h = yak_count(hap1_file, &opt, 0);
	haplo_h = yak_count(hap2_file, &opt, haplo_h);
	completeness = main_completeness(haplo_h, opt.n_thread, hic1_file, hic2_file);
	yak_ch_destroy(haplo_h);

	main_switch_error(opt, hic1_file, hic2_file, hap1_file, hap2_file);




	printf("CC\tCT  kmer_occurrence    short_read_kmer_count  raw_input_kmer_count  adjusted_input_kmer_count\n");
	printf("CC\tFR  fpr_lower_bound    fpr_upper_bound\n");
	printf("CC\tER  total_input_kmers  adjusted_error_kmers\n");
	printf("CC\tCV  coverage\n");
	printf("CC\tQV  raw_quality_value  adjusted_quality_value\n");
	printf("CC\n");
	for (i = (1<<YAK_COUNTER_BITS) - 1; i >= 0; --i)
		printf("CT\t%d\t%ld\t%ld\t%.3f\n", i, (long)hist[i], (long)cnt[i], qs.adj_cnt[i]);
	printf("FR\t%.3g\t%.3g\n", qs.fpr_lower, qs.fpr_upper);
	printf("ER\t%ld\t%.3f\n", (long)qs.tot, qs.err);
	printf("CV\t%.3f\n", qs.cov);
	printf("QV\t%.3f\t%.3f\n", qs.qv_raw, qs.qv);
	printf("Completeness\t%.4f\n",completeness);


	
    return ret;
}

int main_count(int argc, char *argv[])
{
	yak_ch_t *h;
	int c;
	char *fn_out = 0;
	yak_copt_t opt;
	ketopt_t o = KETOPT_INIT;
	yak_copt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:b:H:o:", 0)) >= 0) {
		if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 'p') opt.pre = atoi(o.arg);
		else if (c == 'K') opt.chunk_size = atoi(o.arg);
		else if (c == 't') opt.n_thread = atoi(o.arg);
		else if (c == 'b') opt.bf_shift = atoi(o.arg);
		else if (c == 'H') opt.bf_n_hash = mm_parse_num(o.arg);
		else if (c == 'o') fn_out = o.arg;
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: pstools count [options] <in.fa> [in.fa]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", opt.pre);
		fprintf(stderr, "  -b INT     set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", opt.bf_shift);
		fprintf(stderr, "  -H INT     use INT hash functions for Bloom filter [%d]\n", opt.bf_n_hash);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
		fprintf(stderr, "  -o FILE    dump the count hash table to FILE []\n");
		fprintf(stderr, "  -K INT     chunk size [100m]\n");
		fprintf(stderr, "Note: -b37 is recommended for human reads\n");
		return 1;
	}
	if (opt.pre < YAK_COUNTER_BITS) {
		fprintf(stderr, "ERROR: -p should be at least %d\n", YAK_COUNTER_BITS);
		return 1;
	}
	if (opt.k >= 64) {
		fprintf(stderr, "ERROR: -k must be smaller than 64\n");
		return 1;
	} else if (opt.k >= 32) {
		fprintf(stderr, "WARNING: counts are inexact if -k is greater than 31\n");
	}
	h = yak_count(argv[o.ind], &opt, 0);
	if (opt.bf_shift > 0) {
		yak_ch_destroy_bf(h);
		yak_ch_clear(h, opt.n_thread);
		h = yak_count(argc - o.ind >= 2? argv[o.ind+1] : argv[o.ind], &opt, h);
		yak_ch_shrink(h, 2, YAK_MAX_COUNT, opt.n_thread);
		fprintf(stderr, "[M::%s] %ld distinct k-mers after shrinking\n", __func__, (long)h->tot);
	}

	if (fn_out) yak_ch_dump(h, fn_out);
	yak_ch_destroy(h);
	return 0;
}
