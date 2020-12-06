int main_hic_map(int argc, char *argv[]);

double main_completeness(yak_ch_t *ch, int n_threads, char* file1, char* file2);

int main_switch_error(yak_copt_t opt, char* hic_file1,char* hic_file2, char* hap_file1, char* hap_file2, std::map<std::string, std::string> contig_map);

int main_haplo_map(yak_copt_t opt, char* hic_file1,char* hic_file2, char* pred_file);
