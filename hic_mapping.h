int main_hic_map_haplo(int argc, char *argv[]);
int main_hic_map(int argc, char *argv[]);
double main_completeness(yak_ch_t *ch, int n_threads, string file1, string file2);
int main_switch_error(yak_copt_t opt, string hic_r1,string hic_r2, string hap_file1, string hap_file2, std::map<std::string, std::string> contig_map);
int main_haplo_map(yak_copt_t opt, char* hic_file1,char* hic_file2, char* pred_file);