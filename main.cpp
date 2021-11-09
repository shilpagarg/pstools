extern "C"{
#include "paf.h"
}
#include <string.h>
#include <iostream>
#include "graph.h"
#include "sys.cpp"
#include <fstream>
#include "bubble_chain.h"
#include "paf_intersect.h"
#include "resolve_repeat_haplotype.h"
#include "hic_mapping.h"
#define PSTOOLS_VERSION "0.1"


typedef struct {
    int num_reads;
    char** read_file_names;
    char* output_file_name;
    int thread_num;
} ps_opt_t;

ps_opt_t asm_opt;

void init_opt(ps_opt_t* asm_opt)
{
	memset(asm_opt, 0, sizeof(ps_opt_t));
    asm_opt->num_reads = 0;
    asm_opt->read_file_names = NULL;
    asm_opt->thread_num = 1;
}

void destory_opt(ps_opt_t* asm_opt)
{
    if(asm_opt->read_file_names != NULL)
    {
        free(asm_opt->read_file_names);
    }
}

static ko_longopt_t long_options[] = {
	{ "version",       ko_no_argument, 300 },
	{ "write-paf",     ko_no_argument, 302 },
	{ 0, 0, 0 }
};

void Print_H(ps_opt_t* asm_opt)
{
    fprintf(stderr, "Usage: pstools [options] <input> <output> <...>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  Assembly:\n");
    fprintf(stderr, "    -o FILE       prefix of output files [%s]\n", asm_opt->output_file_name);
    fprintf(stderr, "    -t INT        number of threads [%d]\n", asm_opt->thread_num);
    fprintf(stderr, "    --version     show version number\n");
    fprintf(stderr, "    -h            show help information\n");


    fprintf(stderr, "Example: ./pstools bubble_chain myfile.gfa -o []\n");
    fprintf(stderr, "Example: ./pstools resolve_repeats myfile.gfa [].paf.gz -o []\n");
    fprintf(stderr, "Example: ./pstools asm_component \n");
    fprintf(stderr, "Example: ./pstools filter_paf -t[n] myfile.gfa [].paf.gz -o [].paf \n");
    fprintf(stderr, "See `man ./pstools.1' for detailed description of these command-line options.\n");
}

void get_queries(int argc, char *argv[], ketopt_t* opt, ps_opt_t* asm_opt)
{
    if(opt->ind == argc)
    {
        return;
    }

    asm_opt->num_reads = argc - opt->ind;
    asm_opt->read_file_names = (char**)malloc(sizeof(char*)*asm_opt->num_reads);

    long long i;
    gzFile dfp;
    for (i = 0; i < asm_opt->num_reads; i++)
    {
        asm_opt->read_file_names[i] = argv[i + opt->ind];
        dfp = gzopen(asm_opt->read_file_names[i], "r");
        if (dfp == 0)
        {
            fprintf(stderr, "[ERROR] Cannot find the input read file: %s\n",
                    asm_opt->read_file_names[i]);
		    exit(0);
        }
        gzclose(dfp);
    }
}

int check_option(ps_opt_t* asm_opt)
{
    if(asm_opt->read_file_names == NULL || asm_opt->num_reads == 0)
    {
        fprintf(stderr, "[ERROR] missing input: please specify a read file\n");
        return 0;
    }

    if(asm_opt->output_file_name == NULL)
    {
        fprintf(stderr, "[ERROR] missing output: please specify the output name (-o)\n");
        return 0;
    }

    if(asm_opt->thread_num < 1)
    {
        fprintf(stderr, "[ERROR] the number of threads must be > 0 (-t)\n");
        return 0;
    }

    return 1;
}

int CommandLine_process(int argc, char *argv[], ps_opt_t* asm_opt)
{
    ketopt_t opt = KETOPT_INIT;

    int c;

    while ((c = ketopt(&opt, argc, argv, 1, "hvt:o", long_options)) >= 0) {
        if (c == 'h')
        {
            Print_H(asm_opt);
            return 0;
        }
        else if (c == 't') asm_opt->thread_num = atoi(opt.arg);
        else if (c == 'o') asm_opt->output_file_name = opt.arg;
        else if (c == ':')
        {
			fprintf(stderr, "[ERROR] missing option argument in \"%s\"\n", argv[opt.i - 1]);
			return 1;
		}
        else if (c == '?')
        {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[opt.i - 1]);
			return 1;
		}
    }

    if (argc == opt.ind)
    {
        Print_H(asm_opt);
        return 0;
    }

    get_queries(argc, argv, &opt, asm_opt);

    return check_option(asm_opt);
}

static inline int64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm') x *= 1e6;
	else if (*p == 'K' || *p == 'k') x *= 1e3;
	return (int64_t)(x + .499);
}


asg_t* gfa_read_cerr(char* gfa_file) {
    asg_t* GFAobj = gfa_read(gfa_file);
    if (GFAobj == nullptr) {
        std::cerr << "ERROR: failed to read the graph." << std::endl;
  	} else {
        std::cout << "GFA file read." << std::endl;
    }
    return GFAobj;
}

// ./gfatools view -l utg028545l -r 50 ../../HG002-34X-r244.r_utg.x.gfa > ../../HG002-34X-r244.sub2.gfa
// ./gfatools view -l utg001423l -r 50 ../../HG002-34X-r244.r_utg.x.gfa > ../../HG002-34X-r244.sub4.gfa
// ./gfatools view -l utg036600l -r 100 ../../HG002-34X-r244.r_utg.x.gfa > ../../HG002-34X-r244.sub5.gfa
// ./bc intersect myfile0.gfa ../test1E4.paf ../test_sub.paf
// ./bc intersect ../HG002-34X-r244.sub3.gfa ../test1E4.paf ../test_sub2.paf
// utg000163l
// utg000609l
int main_intersect(int argc, char* argv[]) {
    paf_reader reader;
    printf("start main\n");
    cout << argc <<endl;
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " [gfa-file] [paf-file]" << std::endl;
        return 1;
    }
	// char *sub_gfa_name = argv[1];  // myfile0.gfa
	char *paf_name = argv[2];  // ../test1E4.paf  // I modified this so that some of the nodes would actually match
	char *paf_output_filename = argv[3];  // ../test_sub.paf

    asg_t *GFAobj = gfa_read_cerr(argv[1]);
	paf_file_t *paf_file = reader.paf_open(paf_name);
	output_paf_intersect(GFAobj, paf_file, paf_output_filename, reader);

    printf("finish main\n");
    return 0;

}

int main_clean_graph(int argc, char* argv[]) {
    printf("start main\n");
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input-gfa-file>  <output-gfa-file>" << std::endl;
        return 1;
    }

    asg_t* g;
	int ret = 1;
	int counter = 0;
	string infile = string(argv[1]);
	string outfile = string("pstools.clean_graph.temp.")+to_string(counter)+string(".out");
	string file_out = string(argv[2]);
	while(ret!=0){
		g = gfa_read_cerr(strdup(infile.c_str()));
    	ret = clean_graph(g,infile,outfile,NULL);
		counter++;
		infile = outfile;
		outfile = string("pstools.clean_graph.temp.")+to_string(counter)+string(".out");
		if(counter>=2){
			string rm_cmd = string("rm ") + string("pstools.clean_graph.temp.")+to_string(counter-2)+string(".out");
			system(rm_cmd.c_str());
		}
	}
	string line;
    ifstream finalOutput(infile);
    ofstream outputFile(file_out);
    while (getline(finalOutput, line))
    {
		outputFile << line << endl;
	}
	outputFile.close();
	finalOutput.close();
	string rm_cmd = string("rm ") + infile;
	system(rm_cmd.c_str());
	printf("finish main\n");
    return 0;
}
int main_bubble_chain(int argc, char* argv[]) {
    // printf("start main\n");
    // if (argc != 3) {
    //     std::cerr << "Usage: " << argv[0] << " <gfa-file> <output-directory>" << std::endl;
    //     return 1;
    // }

    // asg_t* g = gfa_read_cerr(argv[1]);
    // // print_gfa(GFAobj);
    // vector<uint32_t> sources = get_sources(g);
	// // get_bubble_chain_graph(g, sources);
    // get_bubbles(g,string(argv[2]));
    // printf("finish main\n");
    return 0;
}

int main_resolve_repeat(int argc, char* argv[]) {
    // printf("start main\n");
    // asg_t* graph;
    // // paf_file_t *paf_file;
    // // struct access *paf_index;
    // if(argc!=4){
    //     cout << "Usage: resolve_repeat <paf file> <gfa file> <output_directory>" << endl;
    //     return -1;
    // }
	// repeat_resolver resolver;
    // char* paf_filename = argv[1];
    // char* gfa_filename = argv[2];
    // graph = gfa_read(gfa_filename);

    // if (graph == nullptr) {
    //     cerr << "ERROR: failed to read the graph." << endl;
  	// } else {
    //     cout << "GFA file read." << endl;
    // }

    // resolver.name_index_mapping = name2idx(graph);
    // paf_reader reader;
    // vector<vector<paf_rec_str_t>>* ordered_records = resolver.get_records_from_paf_file(reader,paf_filename);

    // map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_graph = get_bubbles(graph,string(argv[3]));
    // // for(auto a : count_begin_end){
    // //     cout << graph->seq[a.first[0]>>1].name << " to " << graph->seq[a.first[1]>>1].name << ": " << a.second << endl;
    // // }
    // map<set<uint32_t>,vector<vector<paf_rec_str_t>>>* repeating_nodes = resolver.check_for_same_pos( ordered_records, graph, bubble_chain_graph);
    // resolver.checkForMatch(repeating_nodes, graph);
    // vector<vector<uint32_t>>* supported_pathes = resolver.count_support_reads_for_branches(ordered_records, bubble_chain_graph, graph);
    // cout << "finish getting supported pathes" << endl;
    // set<vector<uint32_t>>* after_covered = resolver.cover_gaps_in_long_path(supported_pathes);
    // cout << "finish trying to cover gaps in supported pathes" << endl;
    // free(supported_pathes);
    // get_seperate_haplotype(after_covered, ordered_records, bubble_chain_graph, graph);
    // free(ordered_records);
    // resolver.save_pathes_to_file(after_covered,bubble_chain_graph,graph,string(argv[3]));
    // cout << "finish saving pathes to files" << endl;
    // free(after_covered);
    // free(bubble_chain_graph);

    // printf("finish main\n");
    return 0;
}


int main_resolve_haplotypes(int argc, char* argv[]) {
	ketopt_t o = KETOPT_INIT;
    int c,n_threads = 8;
	while ((c = ketopt(&o, argc, argv, 1, "t:", 0)) >= 0) {
		if (c == 't') n_threads = atoi(o.arg);
	}
	if (argc - o.ind < 3) {
		fprintf(stderr, "Usage: pstools resolve_haplotypes [options] <hic_mapping file> <gfa file> <output_directory>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t INT     number of threads [%d]\n", n_threads);
		return 1;
	}



    char* connectionFile = argv[o.ind];
    char* gfa_filename = argv[o.ind+1];
    char* output_directory = argv[o.ind+2];
    printf("start main\n");
    asg_t* graph = gfa_read(gfa_filename);
	asg_t** graph_ptr = (&graph);
	// for(int i = 0; i<graph->n_seq; i++){
	    // cout << graph->seq[0].len << "\t" << string(graph->seq[0].seq).length() << endl; 
	    // cout << graph->seq[asg_arc_a(graph,0)[0].v>>1].len << "\t" << string(graph->seq[asg_arc_a(graph,0)[0].v>>1].seq).length() << endl; 
	    // cout << asg_arc_a(graph,0)[0].ol << endl; 
	// }
	// graph->seq[asg_arc_a(graph,0)[0].v>>1].len
    uint32_t **connections_foward;
	CALLOC(connections_foward,graph->n_seq);
	for(int i = 0; i< graph->n_seq; i++){
		CALLOC(connections_foward[i], graph->n_seq);
		memset(connections_foward[i], 0, sizeof(*connections_foward[i]));
	}
	uint32_t **connections_backward;
	CALLOC(connections_backward,graph->n_seq);
	for(int i = 0; i< graph->n_seq; i++){
		CALLOC(connections_backward[i], graph->n_seq);
		memset(connections_backward[i], 0, sizeof(*connections_backward[i]));
	}
    ifstream infile(connectionFile);
    uint32_t i,j,count_forward,count_backward;
	string tag, name, line;
	map<string, uint32_t> node_name_id;
	map<uint32_t, uint32_t> mapID_nodeID;
	for(int i = 0; i < graph->n_seq; i++){
		node_name_id[string(graph->seq[i].name)] = i;
	}
	while (getline(infile, line))
    {
        istringstream iss(line);
        iss >> tag;
        if(tag==string("M")){
			iss >> i >> name;
			mapID_nodeID[i] = node_name_id[name];
        }else if(tag==string("C")){
			iss >> i >> j >> count_backward >> count_forward;
			if(mapID_nodeID.find(i)!=mapID_nodeID.end() && mapID_nodeID.find(j)!=mapID_nodeID.end()){
				i = mapID_nodeID[i];
				j = mapID_nodeID[j];
				connections_backward[i][j] = count_backward;
				connections_backward[j][i] = count_backward;
				connections_foward[i][j] = count_forward;
				connections_foward[j][i] = count_forward;
			}
        }
        // process pair (a,b)
    }
	// gfa_destory(graph);
    map<string,string> excluded_nodes;
    map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_graph = get_bubbles(string(gfa_filename),graph_ptr,string(output_directory),connections_foward,connections_backward,(&excluded_nodes));
    get_haplotype_path(connections_foward, connections_backward, *graph_ptr, bubble_chain_graph, output_directory,n_threads,&excluded_nodes);

    return 0;
}

int main_haplotype_scaffold(int argc, char* argv[]) {
	ketopt_t o = KETOPT_INIT;
    int c,n_threads = 8;
	// string enzymes_unsplit;
	string identityFile;
	bool check_identity = false;
	while ((c = ketopt(&o, argc, argv, 1, "t:e:i:f:", 0)) >= 0) {
		if (c == 't') n_threads = atoi(o.arg);
		// else if (c == 'e') enzymes_unsplit = string(o.arg);
		else if (c == 'i') check_identity = strcmp(o.arg, "false");
		else if (c == 'f') identityFile = string(o.arg);
	}
	if (argc - o.ind < 3) {
		fprintf(stderr, "Usage: pstools haplotype_scaffold [options] <haplotype connection file> <pred_haplotypes.fa> <output_directory>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t INT     number of threads [%d]\n", n_threads);
		// fprintf(stderr, "  -e STR     enzymes separated by comma, optional [%s]\n", enzymes_unsplit.c_str());
		fprintf(stderr, "  -f STR     identity file path, identity check will be ran if it is enabled and no file is given [%s]\n", identityFile.c_str());
		fprintf(stderr, "  -i BOOL    enable identity check on contigs, exclude identical ones while scaffolding [%s]\n", (check_identity ? "true":"false"));
		return 1;
	}



    char* connectionFile = argv[o.ind];
    char* haploFile = argv[o.ind+1];
    char* output_directory = argv[o.ind+2];
    printf("start main\n");

	map<string, string> contig_hap1s;
    map<string, string> contig_hap2s;
    map<string, uint32_t> contig_lengths;
    map<string, uint32_t> contig_id;
    map<uint32_t, string> id_contig;

    ifstream inHaploFile(haploFile);
	string line;
	string name;
	string sequence;
	uint32_t id = 0;
	while (getline(inHaploFile, line)){
        istringstream iss1(line);
		iss1 >> name;
		name = name.substr(1,name.size()-6);
		getline(inHaploFile, line);
        istringstream iss2(line);
		iss2 >> sequence;
		if(id%2 == 0){
			contig_id[name] = (id>>1);
			id_contig[id>>1] = name;
			contig_hap1s[name] = sequence;
			contig_lengths[name] = sequence.size();
		}else{
			contig_hap2s[name] = sequence;
		}
		id++;
	}




    uint32_t **contig_connection;
	CALLOC(contig_connection,contig_lengths.size()*4);
	for(int i = 0; i< contig_lengths.size()*4; i++){
		CALLOC(contig_connection[i], contig_lengths.size()*4);
		memset(contig_connection[i], 0, contig_lengths.size()*4*sizeof(uint32_t));
	}

    ifstream infile(connectionFile);
    string outContig;
    string inContig;
	uint32_t count;
    while(infile >> outContig >> inContig >> count){
		bool inForward = inContig[inContig.size()-1]=='+';
		bool outForward = outContig[outContig.size()-1]=='+';
		bool inHap1 = inContig[inContig.size()-2]=='1';
		uint32_t inID = contig_id[inContig.substr(0,inContig.size()-6)]<<2;		
		uint32_t outID = contig_id[outContig.substr(0,outContig.size()-6)]<<2;
		if(!inForward){
			inID = inID^1;
		}
		if(!outForward){
			outID = outID^1;
		}
		if(!inHap1){
			inID = inID^2;
		}
		contig_connection[outID][inID] = count;
    }
	haplotypes_scaffold(contig_connection, contig_hap1s, contig_hap2s, contig_lengths, contig_id, id_contig, output_directory, identityFile, check_identity, n_threads);
    return 0;
}

int main_obtain_graph_sequence(int argc, char* argv[]){
    printf("start main\n");
    cout << argc <<endl;
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " [full-gfa-file] [gfa-file] [output-gfa-file]" << std::endl;
        return 1;
    }

    ifstream GFA_FULL(argv[1]);
    ofstream result(argv[3]);
    asg_t *GFA_PART = gfa_read_cerr(argv[2]);

    set<string> nodes;
    for(int i  = 0; i < GFA_PART->n_seq; i++){
        nodes.insert(string(GFA_PART->seq[i].name));
    }
    string line;
    while (getline(GFA_FULL, line))
    {
        istringstream iss(line);
        string a,b,c,d,e,f,g;
        iss >> a;
        if(a[0]=='S'){
            iss >> b;
            if(nodes.find(b)!=nodes.end()){
                result << line << endl;
            }
        }else if(a[0]=='L'){
            iss >> b;
            iss >> c;
            iss >> d;
            if(nodes.find(b)!=nodes.end() && nodes.find(d)!=nodes.end()){
                result << line << endl;
            }
        }
        // process pair (a,b)
    }
    GFA_FULL.close();
    result.close();
    printf("finish main\n");
    return 0;

}

int main_hic_mapping(int argc, char* argv[]){
    return main_hic_map(argc, argv);
}

int main_hic_mapping_haplo(int argc, char* argv[]){
    return main_hic_map_haplo(argc, argv);
}
// int main_completeness_check(int argc, char* argv[]){
//     return main_completeness(argc, argv);
// }

// int main_switch_error_check(int argc, char *argv[]){
//     return main_switch_error(argc, argv);
// }

// int main_qv_check(int argc, char *argv[])
// {
// 	yak_qopt_t opt;
// 	yak_ch_t *ch = 0;
// 	ketopt_t o = KETOPT_INIT;
// 	int64_t cnt[YAK_N_COUNTS], hist[YAK_N_COUNTS];
// 	int c, i, kmer;
// 	yak_qstat_t qs;

// 	yak_qopt_init(&opt);
// 	while ((c = ketopt(&o, argc, argv, 1, "K:t:l:f:pe:", 0)) >= 0) {
// 		if (c == 'K') opt.chunk_size = mm_parse_num(o.arg);
// 		else if (c == 'l') opt.min_len = mm_parse_num(o.arg);
// 		else if (c == 'f') opt.min_frac = atof(o.arg);
// 		else if (c == 't') opt.n_threads = atoi(o.arg);
// 		else if (c == 'p') opt.print_each = 1;
// 		else if (c == 'e') opt.fpr = atof(o.arg);
// 	}
// 	if (argc - o.ind < 2) {
// 		fprintf(stderr, "Usage: pstools eval qv [options] <kmers.yak> <prediction.fa>\n");
// 		fprintf(stderr, "Options:\n");
// 		fprintf(stderr, "  -l NUM      min sequence length [%d]\n", opt.min_len);
// 		fprintf(stderr, "  -f FLOAT    min k-mer fraction [%g]\n", opt.min_frac);
// 		fprintf(stderr, "  -e FLOAT    false positive rate [%g]\n", opt.fpr);
// 		fprintf(stderr, "  -p          print QV for each sequence\n");
// 		fprintf(stderr, "  -t INT      number of threads [%d]\n", opt.n_threads);
// 		fprintf(stderr, "  -K NUM      batch size [1g]\n");
// 		return 1;
// 	}

// 	ch = yak_ch_restore(argv[o.ind]);
// 	assert(ch);
// 	kmer = ch->k;
// 	yak_ch_hist(ch, hist, opt.n_threads);
// 	printf("CC\tCT  kmer_occurrence    short_read_kmer_count  raw_input_kmer_count  adjusted_input_kmer_count\n");
// 	printf("CC\tFR  fpr_lower_bound    fpr_upper_bound\n");
// 	printf("CC\tER  total_input_kmers  adjusted_error_kmers\n");
// 	printf("CC\tCV  coverage\n");
// 	printf("CC\tQV  raw_quality_value  adjusted_quality_value\n");
// 	printf("CC\n");
// 	yak_qv(&opt, argv[o.ind+1], ch, cnt);
// 	yak_qv_solve(hist, cnt, kmer, opt.fpr, &qs);
// 	for (i = (1<<YAK_COUNTER_BITS) - 1; i >= 0; --i)
// 		printf("CT\t%d\t%ld\t%ld\t%.3f\n", i, (long)hist[i], (long)cnt[i], qs.adj_cnt[i]);
// 	printf("FR\t%.3g\t%.3g\n", qs.fpr_lower, qs.fpr_upper);
// 	printf("ER\t%ld\t%.3f\n", (long)qs.tot, qs.err);
// 	printf("CV\t%.3f\n", qs.cov);
// 	printf("QV\t%.3f\t%.3f\n", qs.qv_raw, qs.qv);
// 	yak_ch_destroy(ch);
// 	return 0;
// }

int main_completeness_check(int argc, char *argv[]){
	int ret = 0;
	yak_copt_t opt;
	ketopt_t o = KETOPT_INIT;
	yak_ch_t *haplo_h;
	int c, i, kmer;
	double completeness;
	yak_copt_init(&opt);
	opt.pre = YAK_COUNTER_BITS;
	while ((c = ketopt(&o, argc, argv, 1, "k:K:t:", 0)) >= 0) {
		if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 'K') opt.chunk_size = atoi(o.arg);
		else if (c == 't') opt.n_thread = atoi(o.arg);
	}
	if (argc - o.ind < 3) {
		fprintf(stderr, "Usage: pstools completeness [options] <seq.fa> <hic.R1.fastq.gz> <hic.R2.fastq.gz>\n");
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
	char* hic1_file = argv[o.ind+1];
	char* hic2_file = argv[o.ind+2];
	haplo_h = yak_count_create_new(hap1_file, &opt, 0);
	completeness = main_completeness(haplo_h, opt.n_thread, string(hic1_file), string(hic2_file));
	yak_ch_destroy(haplo_h);
	fprintf(stderr, "\n[M::Result] ");
	printf("Completeness\n");
	printf("%.3f%%\n",completeness);
	fprintf(stderr, "\n");
}

int main_qv_check(int argc, char *argv[]){
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

	yak_copt_init(&opt);
	yak_qopt_init(&qopt);
	opt.pre = YAK_COUNTER_BITS;
	while ((c = ketopt(&o, argc, argv, 1, "k:K:t:", 0)) >= 0) {
		if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 'K') opt.chunk_size = atoi(o.arg);
		else if (c == 't') opt.n_thread = atoi(o.arg);
	}
	if (argc - o.ind < 3) {
		fprintf(stderr, "Usage: pstools qv [options] <seq.fa> <hic.R1.fastq.gz> <hic.R2.fastq.gz>\n");
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
	qopt.n_threads = opt.n_thread;

	char* hap1_file = argv[o.ind];
	char* hic1_file = argv[o.ind+1];
	char* hic2_file = argv[o.ind+2];
	hic_h = yak_count_create_new(hic1_file, &opt, 0);
	hic_h = yak_count_create_new(hic2_file, &opt, hic_h);
	kmer = hic_h->k;
	yak_ch_hist(hic_h, hist, opt.n_thread);

	yak_qv(&qopt, hap1_file, hic_h, cnt);
	yak_qv_solve(hist, cnt, kmer, qopt.fpr, &qs);
	yak_ch_destroy(hic_h);
	fprintf(stderr, "\n[M::Result] ");
	printf("QV_RAW\tQV\n");
	printf("%.3f\t%.3f\n", qs.qv_raw, qs.qv);
	fprintf(stderr, "\n");
}

int main_switch_error_check(int argc, char *argv[]){
	int ret = 0;
	yak_copt_t opt;
	ketopt_t o = KETOPT_INIT;
	int c, i, kmer;

	yak_copt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "k:K:t:", 0)) >= 0) {
		if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 'K') opt.chunk_size = atoi(o.arg);
		else if (c == 't') opt.n_thread = atoi(o.arg);
	}
	if (argc - o.ind < 4) {
		fprintf(stderr, "Usage: pstools phasing_error [options] <hap1.fa> <hap2.fa> <hic.R1.fastq.gz> <hic.R2.fastq.gz>\n");
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

	string minimap2_cmd = string("minimap2 -ax asm5 -t");
	minimap2_cmd += to_string(opt.n_thread) + string(" -o temp.sam ");
	minimap2_cmd += string(hap1_file) + " " + string(hap2_file);

	ret = system(minimap2_cmd.c_str());
	
	map<string, map<string, uint32_t>> contig_name_map;
	map<string, map<string, uint32_t>> contig_name_len;
	map<string, map<string, uint32_t>> contig_name_mq;
	map<string, uint32_t> contig_name_max_mq;
	map<string, uint32_t> contig_name_max_len;

    ifstream infile("temp.sam");
    string hap2_name, hap1_name, dump;
	uint32_t mapq, mlen;
	string line;
	while (std::getline(infile, line)) {
    	std::istringstream is(line);
		is >> hap1_name;
		if(hap1_name[0]=='@') continue;
		is >> dump >> hap2_name;
		is >> dump >> mapq;
		is >> dump >> dump >> dump >> dump >> dump;
		mlen = dump.size();
		if(mlen == 1 && dump[0]=='='){
			mlen = 0xFFFFFFFF;
		}
		// cout << hap1_name << " " << hap2_name << " " << mapq << " " << mlen << endl;
		if(strcmp(hap1_name.c_str() ,"*") !=0 && (strcmp(hap2_name.c_str() ,"*") !=0)){
			if(contig_name_map.find(hap2_name) == contig_name_map.end()){
				contig_name_map[hap2_name] = map<string, uint32_t>();
				contig_name_len[hap2_name] = map<string, uint32_t>();
				contig_name_mq[hap2_name] = map<string, uint32_t>();
				contig_name_max_mq[hap2_name] = 0;
				contig_name_max_len[hap2_name] = 0;
			}
			if(contig_name_len[hap2_name].find(hap1_name)==contig_name_len[hap2_name].end()){
				contig_name_len[hap2_name][hap1_name] = 0;
			}
			if(contig_name_mq[hap2_name].find(hap1_name)==contig_name_mq[hap2_name].end()){
				contig_name_mq[hap2_name][hap1_name] = 0;
			}
			contig_name_max_mq[hap2_name] = max(mapq, contig_name_max_mq[hap2_name]);
			contig_name_max_len[hap2_name] = max(mlen, contig_name_max_len[hap2_name]);
			contig_name_len[hap2_name][hap1_name] = max(mlen, contig_name_len[hap2_name][hap1_name]);
			contig_name_mq[hap2_name][hap1_name] = max(mapq, contig_name_mq[hap2_name][hap1_name]);
			contig_name_map[hap2_name][hap1_name]++;
		}
	}
	map<string,string> contig_map;
	set<string> matched_contigs;
	for(auto i : contig_name_map){
		uint32_t mq_thred = (contig_name_max_mq[i.first] == 60 ? 50 : 0);
		uint32_t max_count = 0;
		string max_matched;
		for(auto j : i.second){
			if(contig_name_mq[i.first][j.first] > mq_thred && contig_name_len[i.first][j.first] > contig_name_max_len[i.first]*0.8 && max_count < j.second && matched_contigs.find(j.first)==matched_contigs.end()){
				max_count = j.second;
				max_matched = j.first;
			}
		}
		if(max_count>0){
			contig_map[i.first] = max_matched;
			matched_contigs.insert(max_matched);
		}
	}
	for(auto i : contig_map){
		fprintf(stderr, "[Mapping relationship] %s and %s \n", i.first.c_str(), i.second.c_str());
		// cout << i.first << " and " << i.second << endl;
	}
	// cout << contig_map.size() << endl;
	system("rm temp.sam");
	opt.pre = YAK_COUNTER_BITS;
	main_switch_error(opt, string(hic1_file), string(hic2_file), string(hap1_file), string(hap2_file), contig_map);

}


// int main_identity_check(int argc, char *argv[]){
// 	stringstream minimap2_cmd;
// 	int c;
// 	yak_copt_t opt;
// 	ketopt_t o = KETOPT_INIT;
// 	yak_copt_init(&opt);
// 	while ((c = ketopt(&o, argc, argv, 1, "t:", 0)) >= 0) {
// 		if (c == 't') opt.n_thread = atoi(o.arg);
// 	}
// 	if (argc - o.ind < 1) {
// 		fprintf(stderr, "Usage: pstools identity_check [options] [in.fa]\n");
// 		fprintf(stderr, "Options:\n");
// 		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
// 		return 1;
// 	}
// 	minimap2_cmd << "minimap2 -I40G -x asm20 -Y -a --eqx -t" << opt.n_thread << " " << argv[o.ind] << " "<< argv[o.ind]; 
// 	minimap2_cmd << " | samtools view -F 4 -u - | samtools sort - > haplotype_identity.bam";

// 	int ret = system(minimap2_cmd.str().c_str());
// 	ret = system("python ./samIdentity.py --header haplotype_identity.bam | awk '$1 != $5' > output.tbl");
// 	return ret;
// }

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

int main(int argc, char *argv[])
{
	extern double yak_realtime(void);
	extern double yak_cputime(void);
	extern void yak_reset_realtime(void);
	double t_start;
	int ret = 0, i;

	if (argc == 1) {
		fprintf(stderr, "Usage: pstools <command> <arguments> <inputs>\n");
		fprintf(stderr, "Commands:\n");
		// fprintf(stderr, "  bubble_chain          returns bubble chains from hifi graph using ONT data\n");
		fprintf(stderr, "  clean_graph           cleaning the given graph\n");
		fprintf(stderr, "  intersect             extract corresponding paf records for sub graph\n");
		// fprintf(stderr, "  resolve_repeat        use ONT data to resolve repeats in hifi graph\n");
		fprintf(stderr, "  resolve_haplotypes    use hic data to resolve haplotypes\n");
		fprintf(stderr, "  haplotype_scaffold    scaffolding the predicted haplotypes\n");
		fprintf(stderr, "  hic_mapping_unitig    map hic data to sequences in the graph\n");
		fprintf(stderr, "  hic_mapping_haplo     map hic data to predicted haplotypes for scaffolding\n");
		// fprintf(stderr, "  count                 count k-mers\n");
		// fprintf(stderr, "  identity_check        generate identity table for contigs\n");
        fprintf(stderr, "  qv                    calculate qv score for prediction\n");
        fprintf(stderr, "  completeness          calculate completeness for prediction\n");
        fprintf(stderr, "  phasing_error         calculate switch error\n");
		fprintf(stderr, "  version               print version number\n");
		return 1;
	}
	yak_reset_realtime();
	t_start = yak_realtime();
	if (strcmp(argv[1], "intersect") == 0) ret = main_intersect(argc-1, argv+1);
	// else if (strcmp(argv[1], "bubble_chain") == 0) ret = main_bubble_chain(argc-1, argv+1);
	else if (strcmp(argv[1], "clean_graph") == 0) ret = main_clean_graph(argc-1, argv+1);
	// else if (strcmp(argv[1], "resolve_repeat") == 0) ret = main_resolve_repeat(argc-1, argv+1);
	else if (strcmp(argv[1], "obtain_graph_sequence") == 0) ret = main_obtain_graph_sequence(argc-1, argv+1);
	else if (strcmp(argv[1], "hic_mapping_unitig") == 0) ret = main_hic_mapping(argc-1, argv+1);
	else if (strcmp(argv[1], "hic_mapping_haplo") == 0) ret = main_hic_mapping_haplo(argc-1, argv+1);
	else if (strcmp(argv[1], "haplotype_scaffold") == 0) ret = main_haplotype_scaffold(argc-1, argv+1);
	else if (strcmp(argv[1], "resolve_haplotypes") == 0) ret = main_resolve_haplotypes(argc-1, argv+1);
	// else if (strcmp(argv[1], "count") == 0) ret = main_count(argc-1, argv+1);
	else if (strcmp(argv[1], "phasing_error") == 0) ret = main_switch_error_check(argc-1, argv+1);
	else if (strcmp(argv[1], "qv") == 0) ret = main_qv_check(argc-1, argv+1);
	else if (strcmp(argv[1], "completeness") == 0) ret = main_completeness_check(argc-1, argv+1);
	// else if (strcmp(argv[1], "identity_check") == 0) ret = main_identity_check(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		printf("gfa.h: %s\nps: %s\n", PSTOOLS_VERSION, PSTOOLS_VERSION);
		return 0;
	} else {
		fprintf(stderr, "[E::%s] unknown command\n", __func__);
		return 1;
	}
	if (ret == 0) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, PSTOOLS_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, yak_realtime() - t_start, yak_cputime());
	}
	return ret;
}
