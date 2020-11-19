#include "graph.h"
#include <iostream>
#include <vector>
#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include "paf.h"
#include <map>
#include "bubble_chain.h"
#include <queue>
#include <assert.h>
#include <sstream>
#include "minimizer_sketch.h"
#include "minimap.h"
#include "modmap.c"
#include "kthread.h"

using namespace std;

string complement(string unitig){
    stringstream result;
    int c = 0;
    while(c < unitig.size()){
        char current = unitig[c];
        if(current == 'A'){
            result << 'T';
        }else if(current == 'T'){
            result << 'A';
        }else if(current == 'C'){
            result << 'G';
        }else if(current == 'G'){
            result << 'C';
        }else{
            result << 'N';
        }
        c++;
    }
    string result_str = result.str();
    reverse(result_str.begin(),result_str.end());
    return result_str;
}
string translate(string unitig, uint32_t start, string cs_tag){
    stringstream result;
    char cur_operator;
    char current_char;
    int c = 0;
    int u = start;

    while(c < cs_tag.size()){
        current_char = cs_tag[c];
        if(current_char == ':' || current_char == '*' || current_char == '=' || current_char == '+' || current_char == '-'){
            cur_operator = current_char;
            c++;
        }else{
            if(cur_operator == '='){
                while(cs_tag[c] == 'a' || cs_tag[c] == 'c' || cs_tag[c] == 'g' || cs_tag[c] == 't' || cs_tag[c] == 'n'){
                    result << (char)(cs_tag[c]-32);
                    c++;
                    u++;
                }
            }else if(cur_operator == ':'){
                stringstream current_num;
                while(isdigit(cs_tag[c])){
                    current_num << cs_tag[c];
                    c++;
                }
                int counter = stoi(current_num.str());
                for(int i = 0; i < counter; i++){
                    result << unitig[u];
                    u++;
                }
            }else if(cur_operator == '*'){
                c+=2;
                u++;
                result << (char)(cs_tag[c-1]-32);
            }else if(cur_operator == '+'){
                while(cs_tag[c] == 'a' || cs_tag[c] == 'c' || cs_tag[c] == 'g' || cs_tag[c] == 't' || cs_tag[c] == 'n'){
                    result << (char)(cs_tag[c]-32);
                    c++;
                }
            }else if(cur_operator == '-'){
                while(cs_tag[c] == 'a' || cs_tag[c] == 'c' || cs_tag[c] == 'g' || cs_tag[c] == 't' || cs_tag[c] == 'n'){
                    c++;
                    u++;
                }
            }
        }
    }
    return result.str();
}


set<string> find_variant(string unitig, uint32_t start, string cs_tag, int kmer_l){
    int pre = (kmer_l/2);
    int post = pre;
    if(kmer_l%2 == 1){
        post++;
    }
    char cur_operator;
    char current_char;
    int c = 0;
    int u = start;
    set<string> results = set<string>();
    stringstream buf;

    while(c < cs_tag.size()){
        current_char = cs_tag[c];
        if(current_char == ':' || current_char == '*' || current_char == '=' || current_char == '+' || current_char == '-'){
            cur_operator = current_char;
            c++;
        }else{
            if(cur_operator == '='){
                while(cs_tag[c] == 'a' || cs_tag[c] == 'c' || cs_tag[c] == 'g' || cs_tag[c] == 't' || cs_tag[c] == 'n'){
                    c++;
                    u++;
                }
            }else if(cur_operator == ':'){
                stringstream current_num;
                while(isdigit(cs_tag[c])){
                    current_num << cs_tag[c];
                    c++;
                }
                int counter = stoi(current_num.str());
                u+=counter;
            }else if(cur_operator == '*'){
                c+=2;
                if(u>=pre && u+post <= unitig.length()){
                    results.insert(unitig.substr(u-pre,kmer_l));
                }else if(u >= pre){
                    results.insert(unitig.substr(unitig.length() - kmer_l,kmer_l));
                }else{
                    results.insert(unitig.substr(0,kmer_l));
                }
                u++;
            }else if(cur_operator == '-'){
                if(u>=pre && u+post <= unitig.length()){
                    results.insert(unitig.substr(u-pre,kmer_l));
                }else if(u >= pre){
                    results.insert(unitig.substr(unitig.length() - kmer_l,kmer_l));
                }else{
                    results.insert(unitig.substr(0,kmer_l));
                }
                while(cs_tag[c] == 'a' || cs_tag[c] == 'c' || cs_tag[c] == 'g' || cs_tag[c] == 't' || cs_tag[c] == 'n'){
                    c++;
                }
            }else if(cur_operator == '+'){
                buf = stringstream();
                while(cs_tag[c] == 'a' || cs_tag[c] == 'c' || cs_tag[c] == 'g' || cs_tag[c] == 't' || cs_tag[c] == 'n'){
                    c++;
                    u++;
                    buf << cs_tag[c];
                }
                string buffer = buf.str();
                if(buffer.length()>kmer_l){
                    int buf_u = u - buffer.length();
                    while(buf_u+pre <= u && u+post < unitig.length()){
                        results.insert(unitig.substr(buf_u,kmer_l));
                        buf_u += kmer_l;
                    }

                    if(u>buffer.length()+post){
                        results.insert(unitig.substr(u-buffer.length()-post,kmer_l));
                    }else{
                        results.insert(unitig.substr(0,kmer_l));
                    }
                }else{
                    if(u>buffer.length()+1+(kmer_l-buffer.length())/2){
                        results.insert(unitig.substr(u-buffer.length()-1-(kmer_l-buffer.length())/2,kmer_l));
                    }else{
                        results.insert(unitig.substr(0,kmer_l));
                    }
                }

            }
        }
    }
    return results;
}








typedef struct {
    string qn;
	uint32_t tn;
	uint32_t qs;
    uint32_t ql, qe, tl, ts, te;
	uint32_t ml, bl;
    string seq;
} paf_rec_str_t;

paf_rec_str_t remove_redundent_info(paf_rec_t r, map<string, uint32_t>* name_index_mapping){
    paf_rec_str_t x;
    x.tn = (*name_index_mapping)[r.tn] << 1;
    x.tn = r.rev? x.tn^1 : x.tn;
    x.qs = r.qs;
    x.ml = r.ml;
    x.bl = r.bl;
    x.ql = r.ql;
    x.qe = r.qe;
    x.tl = r.tl;
    x.ts = r.ts;
    x.te = r.te;
    x.qn = string(r.qn);
    if(r.seq[0] == ':' || r.seq[0] == '+' || r.seq[0] == '-' || r.seq[0] == '*' || r.seq[0] == '='){
        x.seq = string(r.seq);
    }
    return x;
}

map<string, uint32_t>* name2idx(asg_t* graph){
    map<string,uint32_t>* result = new map<string,uint32_t>();
    for(uint32_t i = 0; i < graph->n_seq; i++){
        (*result)[string(graph->seq[i].name)] = i;
    }
    return result;
}


class repeat_resolver{
    public:
        // asg_t* graph;
        // FILE *paf_file;
        // struct access *paf_index;
        // mutex* paf_read_mutex = new mutex();
        map<string, uint32_t>* name_index_mapping;
        vector<vector<paf_rec_str_t>>* get_records_from_paf_file(paf_reader reader, char* filename);
        vector<vector<uint32_t>>* count_support_reads_for_branches(vector<vector<paf_rec_str_t>>* ordered_records  \
                                                                        ,map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_graph \
                                                                        ,asg_t* graph);
        map<set<uint32_t>,vector<vector<paf_rec_str_t>>>* check_for_same_pos(vector<vector<paf_rec_str_t>>* ordered_records, asg_t* graph, map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_graph);
        map<uint32_t, vector<uint32_t>> checkForMatch(map<set<uint32_t>,vector<vector<paf_rec_str_t>>>* repeating_nodes, asg_t* graph);
        map<uint32_t, vector<uint32_t>> aligning_for_variantion(map<set<uint32_t>,map<map<uint32_t,string>,string>>* repeating_nodes, asg_t* graph);
        set<vector<uint32_t>>* cover_gaps_in_long_path(vector<vector<uint32_t>>* supported_pathes);
        set<vector<uint32_t>>* merge_long_pathes(set<vector<uint32_t>>* path_set);
        set<vector<uint32_t>>* recover_unlinearized_branches(set<vector<uint32_t>>* path_set, map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_graph);
        map<string, vector<vector<paf_rec_str_t>*>*> get_queries_with_duplicated_matches(map<string, vector<paf_rec_str_t>*> unordered_records_by_query_name);
        vector<vector<uint32_t>> detect_cycles_in_graph(vector<uint32_t> duplicated_nodes, asg_t* graph);
        void save_pathes_to_file(set<vector<uint32_t>>* set_of_pathes, map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_graph, asg_t* graph, string output_directory);


};

vector<vector<paf_rec_str_t>>* repeat_resolver::get_records_from_paf_file(paf_reader reader, char* filename){
    paf_file_t *paf_file;
    paf_file = reader.paf_open(filename);

    if (paf_file == NULL) {
        cout << "Error: failed to read PAF file." << endl;
    }else{
        cout << "PAF file read." << endl;
    }
    paf_rec_t r;
    map<string, vector<paf_rec_str_t>*> unordered_records_by_query_name;
    vector<vector<paf_rec_str_t>>* result = new vector<vector<paf_rec_str_t>>();
    cout << "start read paf file." << endl;
    bool first_map = true;
    int total_count = 0;
    int long_count = 0;
	while (reader.paf_read(paf_file, &r) >= 0) {
        total_count++;
        if(r.ql >= 10000){
            long_count++;
            // cout << r.qn << "\t";
            // cout << r.ql << "\t";
            // cout << r.qs << "\t";
            // cout << r.qe << "\t";
            // cout << (r.rev ? "-" : "+") << "\t";  // TODO
            // cout << r.tn << "\t";
            // cout << r.tl << "\t";
            // cout << r.ts << "\t";
            // cout << r.te << "\t";
            // cout << r.ml << "\t";
            // cout << r.bl << "\t";
            // // Mapping Quality is not parsed
            // cout << endl;
            if(unordered_records_by_query_name[r.qn]==nullptr){
                unordered_records_by_query_name[r.qn] = new vector<paf_rec_str_t>;
            }
            unordered_records_by_query_name[r.qn]->push_back(remove_redundent_info(r, this->name_index_mapping));
        }
	}
    for(auto a:unordered_records_by_query_name){
        sort(a.second->begin(), a.second->end(), [ ]( const auto& lhs, const auto& rhs )
        {
        return lhs.qs < rhs.qs;
        });
        result->push_back(*a.second);
        free(a.second);
    }
    cout << "Total records: " << total_count <<endl;
    cout <<"Long records: " << long_count << endl;
    cout << "Percentage: " << ((double)long_count / (double)total_count) * 100 << "%" << endl;
	reader.paf_close(paf_file);
    return result;
}


map<set<uint32_t>,vector<vector<paf_rec_str_t>>>* repeat_resolver::check_for_same_pos(vector<vector<paf_rec_str_t>>* ordered_records, asg_t* graph, map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_graph){
    set<uint32_t>* nodes_to_begin_end = new set<uint32_t>[graph->n_seq*2];
    for(auto loop_begin: *bubble_chain_graph){
        for(auto loop_end: loop_begin.second){
            if(loop_end.second.size()>1){
                for(auto node: loop_end.second){
                    // cout << node;
                    nodes_to_begin_end[node].insert(loop_begin.first>>1);
                    nodes_to_begin_end[node].insert(loop_end.first>>1);
                    nodes_to_begin_end[node^1].insert(loop_begin.first>>1);
                    nodes_to_begin_end[node^1].insert(loop_end.first>>1);
                }
            }
        }
    }
    for(int i = 0; i<graph->n_seq*2; i++){
        if(nodes_to_begin_end[i].size() == 0){
            nodes_to_begin_end[i].insert(i);
        }
    }
    map<set<uint32_t>,vector<vector<paf_rec_str_t>>>* repeating_nodes = new map<set<uint32_t>, vector<vector<paf_rec_str_t>>>();
    map<set<uint32_t>,vector<vector<paf_rec_str_t>>> buf_repeating_nodes = map<set<uint32_t>, vector<vector<paf_rec_str_t>>>();
    for(auto record: *ordered_records){
        bool begin = true;
        int start = 0;
        set<uint32_t> current_result;
        set<set<uint32_t>> begin_end_set;
        for(auto query: record){
            if(begin){
                start = query.qs;
                current_result.insert(query.tn);
                begin_end_set.insert(nodes_to_begin_end[query.tn]);
                begin = false;
            }else{
                if(query.qs == start){
                    current_result.insert(query.tn);
                    begin_end_set.insert(nodes_to_begin_end[query.tn]);
                }else{
                    if(current_result.size() > 1 && begin_end_set.size() > 1){ //begin_end_set.size() > 1 ?

                        // bool flag = true;
                        // for(auto res: current_result){
                        //     if((strcmp(string(graph->seq[res>>1].name).c_str(),string("utg000163l").c_str())==0)||(strcmp(string(graph->seq[res>>1].name).c_str(),string("utg002848l").c_str())==0)){
                        //     }else{
                        //         flag = false;
                        //         break;
                        //     }
                        // }

                        // if(flag){
                        //     for(auto bufff: begin_end_set){
                        //         for(auto buffff: bufff){
                        //             cout << graph->seq[buffff>>1].name << "\t";
                        //         }
                        //         cout << endl;
                        //     }
                        //     cout << endl;
                        // }


                        if(buf_repeating_nodes.find(current_result)==buf_repeating_nodes.end()){
                            buf_repeating_nodes[current_result] = vector<vector<paf_rec_str_t>>();
                        }
                        vector<paf_rec_str_t> buf;
                        for(auto buffer: record){
                            if(buffer.qs == start){
                                buf.push_back(buffer);
                            }
                        }
                        buf_repeating_nodes[current_result].push_back(buf);
                    }
                    current_result = set<uint32_t>();
                    begin_end_set = set<set<uint32_t>>();
                    begin_end_set.insert(nodes_to_begin_end[query.tn]);
                    current_result.insert(query.tn);
                    start = query.qs;
                }
            }
        }
    }
    // for(auto x: *repeating_nodes){
    //     for(auto q: x.second[0]){
    //         if(x.first.find(q.tn)!=x.first.end()){
    //             cout << graph->seq[q.tn>>1].name << " ";
    //             // ha_sketch(graph->seq[q.tn>>1].seq, graph->seq[q.tn>>1].len, 30, 100, 0, 0, new ha_mz1_v(), nullptr);

    //             if(q.tn%2==1){
    //                 cout << "-" << endl;
    //                 int start = q.tl - q.te;
    //                 cout << complement(string(graph->seq[q.tn>>1].seq)).substr(start,120) << endl;
    //                 cout << (q.seq).substr(0,120) << endl;
    //                 cout << translate(complement(string(graph->seq[q.tn>>1].seq)), start , q.seq).substr(0,120) << endl;
    //             }else{
    //                 cout << "+" << endl;
    //                 cout << string(graph->seq[q.tn>>1].seq).substr(q.ts,120) << endl;
    //                 cout << (q.seq).substr(0,120) << endl;
    //                 cout << translate(graph->seq[q.tn>>1].seq, q.ts , q.seq).substr(0,120) << endl;
    //             }
    //             // break;
    //         }
    //         // break;
    //     }
    //     // for(auto y: x.first){
    //     //     cout << graph->seq[y>>1].name << "\t";
    //     // }
    //     // cout << ": " << x.second.size();
    //     // cout << endl;
    // }

    // for(auto x: *repeating_nodes){
    //     for(auto q: x.second[0]){
    //         if(x.first.find(q.tn)!=x.first.end()){
    //             cout << graph->seq[q.tn>>1].name << " ";
    //             // ha_sketch(graph->seq[q.tn>>1].seq, graph->seq[q.tn>>1].len, 30, 100, 0, 0, new ha_mz1_v(), nullptr);
    //             seq_mini_t* sequence_minimizer;
    //             if(q.tn%2==1){
    //                 cout << "-" << endl;
    //                 sequence_minimizer = get_minimizer(complement(string(graph->seq[q.tn>>1].seq)),100, 7);
    //             }else{
    //                 cout << "+" << endl;
    //                 sequence_minimizer = get_minimizer(string(graph->seq[q.tn>>1].seq),100, 7);
    //             }

    //             // for(int i = 0; i < current_pos->size(); i++){
    //             // }
    //             // break;
    //         }
    //         // break;
    //     }
    //     // for(auto y: x.first){
    //     //     cout << graph->seq[y>>1].name << "\t";
    //     // }
    //     // cout << ": " << x.second.size();
    //     // cout << endl;
    // }
    // free(nodes_to_begin_end);
    for(auto x: buf_repeating_nodes){
        if(x.second.size()>2){
            (*repeating_nodes)[x.first] = x.second;
        }
    }
    return repeating_nodes;
}


map<uint32_t, vector<uint32_t>> repeat_resolver::checkForMatch(map<set<uint32_t>,vector<vector<paf_rec_str_t>>>* repeating_nodes, asg_t* graph){
    map<uint32_t, map<string, double>> sequence_minimizers;
    map<set<uint32_t>, map<uint32_t, uint32_t>> result;
    map<uint32_t, vector<uint32_t>> substitution_relation;
    map<uint32_t, uint32_t> substitution_relation_reversed;
    map<set<uint32_t>,map<map<uint32_t,string>,string>>* to_align_new = new map<set<uint32_t>,map<map<uint32_t,string>,string>>();
    for(auto x: *repeating_nodes){
        (*to_align_new)[x.first] = map<map<uint32_t,string>,string>();
        result[x.first] = map<uint32_t, uint32_t>();
        for(auto y: x.first){
            sequence_minimizers[y] = map<string, double>();
            result[x.first][y] = 0;
        }
    }
    for(auto x: sequence_minimizers){
        string current_seq = string(graph->seq[x.first>>1].seq);
        if(x.first%2==1){
            current_seq = complement(current_seq);
        }
        seq_mini_t* buf_minis = get_minimizer(current_seq, 31 ,17);
        sequence_minimizers[x.first] = profile_minimizer(buf_minis);
        free(buf_minis);
        // vector<minimizer_profile_t> minimizer_profile = profile_minimizer(sequence_minimizers[x.first]);
        // string max_mini;
        // uint32_t max_count = 0;
        // uint64_t total = 0;
        // uint32_t  num_bin = 0;
        // for(auto b: minimizer_profile){
        //     total += b.count;
        //     num_bin += 1;
        //     if(max_count < b.count){
        //         max_mini = b.minimizer;
        //         max_count = b.count;
        //     }
        // }
        // cout << graph->seq[x.first>>1].name << " max mini: " << max_mini << "," <<" percentage: " << ((double)max_count / (double)total) * 100 << " against "<< ((double)1 / (double)num_bin) * 100 << endl;
    }
    for(auto x: *repeating_nodes){
        map<uint32_t, map<string, double>> full_sequence_minis;
        for(auto node: x.first){
            full_sequence_minis[node] = sequence_minimizers[node];
            // cout << graph->seq[node>>1].name << "\t";
        }
        // cout << endl;
        for(auto query: x.second){
            map<uint32_t,string> current_sub_map = map<uint32_t,string>();
            vector<string> read_sequences;
            string read_sequence_max;
            map<uint32_t, map<string, double>> sub_sequence_minis;
            uint32_t end = UINT32_MAX;
            for(auto record: query){
                if(record.qe < end){
                    end = record.qe;
                }
            }

            for(auto record: query){
                uint32_t buf_end = record.te;
                uint32_t buf_start = record.ts;
                string current_unitig_seq = graph->seq[record.tn>>1].seq;
                if(record.tn%2==1){
                    uint32_t buf_length = record.tl;
                    uint32_t buffer_end = buf_length - buf_start;
                    buf_start = buf_length - buf_end;
                    buf_end = buffer_end;
                    current_unitig_seq = complement(current_unitig_seq);
                }
                sub_sequence_minis[record.tn] = full_sequence_minis[record.tn];
                // sub_sequence_minis[record.tn] = sub_minimizer(full_sequence_minis[record.tn], buf_start, buf_end);
                // cout << "Total length " << record.tl << ": " << "start at " << record.ts << ", end at: " << record.te << endl;

                current_sub_map[record.tn] = string(current_unitig_seq,buf_start,buf_end);
                // current_sub_map[record.tn] = current_unitig_seq;
                if(record.tn % 2 == 0){
                    read_sequences.push_back(translate(graph->seq[record.tn>>1].seq, buf_start , record.seq));
                }else{
                    read_sequences.push_back(translate(complement(graph->seq[record.tn>>1].seq), buf_start , record.seq));
                }
                // cout << buf_start << " to " << buf_end << ": " << buf_end - buf_start << endl;
                // cout << sub_sequence_minis[record.tn]->pos.front()[0] << " to " << sub_sequence_minis[record.tn]->pos.back()[1] << ": " << sub_sequence_minis[record.tn]->pos.back()[1] - sub_sequence_minis[record.tn]->pos.front()[0] << endl;
            }
            for(auto seq: read_sequences){
                if(seq.size()>read_sequence_max.size()){
                    read_sequence_max = seq;
                }
            }
            (*to_align_new)[x.first][current_sub_map] = read_sequence_max;
            seq_mini_t* read_minimizer = get_minimizer(read_sequence_max, 31 ,17);
            map<string, double> read_minimizer_profile = profile_minimizer(read_minimizer);
            map<uint32_t, double> minimizer_match_result;
            for(auto xx: sub_sequence_minis){
                minimizer_match_result[xx.first] = minimizer_match(read_minimizer_profile,xx.second);
            }
            // map<uint32_t, uint64_t> match_result = match_minimizers(read_minimizer, sub_sequence_minis);
            double max_val = 0;
            for(auto res: minimizer_match_result){
                // cout << graph->seq[res.first>>1].name << ": " << res.second << endl;
                if(res.second > max_val){
                    max_val = res.second;
                }
            }
            for(auto res: minimizer_match_result){
                result[x.first][res.first] += 1;
                if(res.second <= max_val * 0.9){
                    result[x.first][res.first] -= 1;
                }
            }
            // for(auto res: minimizer_match_result){
            //     cout << graph->seq[res.first>>1].name << ": " << res.second << endl;
            // }

            // map<uint32_t, seq_mini_t*> repeat_variant = get_repeat_variant(sub_sequence_minis);
            // // for(auto variant: repeat_variant){
            // //     cout << graph->seq[variant.first>>1].name << ":" << endl;
            // //     for(int i = 0; i<variant.second->minimizers.size(); i++){
            // //         cout << variant.second->pos[i][0] << " to " << variant.second->pos[i][1] << ": " << variant.second->minimizers[i] << endl;
            // //     }
            // // }
            // map<uint32_t, vector<vector<uint32_t>>> variant_region;
            // for(auto x: repeat_variant){
            //     variant_region[x.first] = vector<vector<uint32_t>>();
            //     for(int i = 0; i < x.second->pos.size(); i++){
            //         if(i==0){
            //             variant_region[x.first].push_back(vector<uint32_t>());
            //             variant_region[x.first].back().push_back(x.second->pos.front().front());
            //         }else if(repeat_variant[x.first]->pos[i-1][1] != repeat_variant[x.first]->pos[i][0]){
            //             variant_region[x.first].back().push_back(x.second->pos[i].front());
            //             variant_region[x.first].push_back(vector<uint32_t>());
            //             variant_region[x.first].back().push_back(x.second->pos[i].front());
            //         }else if(i==x.second->pos.size()-1 && variant_region[x.first].back().size()==1){
            //             variant_region[x.first].back().push_back(x.second->pos[i].back());
            //         }
            //     }
            //     if(variant_region[x.first].back().size()==1){
            //         variant_region[x.first].back().push_back(x.second->pos.back().back());
            //     }
            // }
            // for(auto x: variant_region){
            //     cout << graph->seq[x.first>>1].name << ": " <<endl;
            //     for(auto y: x.second){
            //         cout << y[0] << " to " << y[1] << endl;
            //     }
            // }
        }
    }
    // map<uint32_t, vector<uint32_t>> buf_substitution;
    // for(auto res: result){
    //     uint32_t max_val = 0;
    //     uint32_t max_node = 0;
    //     for(auto x: res.second){
    //         if(max_val < x.second){
    //             max_val = x.second;
    //             max_node = x.first;
    //         }
    //     }
    //     for(auto x: res.second){
    //         if(max_val > 1 && x.first != max_node){
    //             if(substitution_relation.find(x.first)!=substitution_relation.end()){
    //                 if(substitution_relation[x.first].back()<max_val){
    //                     substitution_relation[x.first] = vector<uint32_t>();
    //                     substitution_relation[x.first].push_back(max_node);
    //                     substitution_relation[x.first].push_back(max_val);
    //                 }
    //             }else{
    //                 substitution_relation[x.first] = vector<uint32_t>();
    //                 substitution_relation[x.first].push_back(max_node);
    //                 substitution_relation[x.first].push_back(max_val);
    //             }
    //         }
    //         cout << graph->seq[x.first>>1].name << ": " << x.second <<"\t";
    //     }
    //     cout << endl;
    // }

    for(auto res: result){
        bool solvable = true;
        uint32_t max_val = 0;
        uint32_t max_result;
        for(auto buf: res.first){
            cout << graph->seq[buf>>1].name << "\t";
        }
        // cout << endl;
        for(auto x: res.second){
            if(x.second > max_val){
                max_val = x.second;
                max_result = x.first;
            }
            // cout << graph->seq[x.first>>1].name << ": " << x.second <<"\t" << endl;
        }
        for(auto x: res.second){
            if(x.first!=max_result && (double)x.second>=((double)max_val)*0.8){
                solvable = false;
            }
            // cout << graph->seq[x.first>>1].name << ": " << x.second <<"\t" << endl;
        }
        if(!solvable){
            cout << "Unsolvable" << endl;
        }else{
            cout << "Solveable: " << graph->seq[max_result>>1].name << endl;
            to_align_new->erase(res.first);
        }
        // for(auto x:res.second){
        //     cout << graph->seq[x.first>>1].name << ": " << x.second <<"\t" << endl;
        // }
    }
    // aligning_for_variantion(to_align_new, graph);
    return substitution_relation;
}

map<uint32_t, vector<uint32_t>> repeat_resolver::aligning_for_variantion(map<set<uint32_t>,map<map<uint32_t,string>,string>>* repeating_nodes, asg_t* graph){
    mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	mopt.min_cnt = 500;
	int n_threads = 3;
	mm_verbose = 2; // disable message output to stderr
	mm_set_opt(0, &iopt, &mopt);
	mopt.flag |= MM_F_CIGAR; // perform alignment
	mopt.flag |= MM_F_OUT_CS; // perform alignment
    params.k = 31 ;
    params.w = 31 ;
    params.s = 7 ;
    params.B = 28 ;
    for(auto repeats: *repeating_nodes){
        set<uint32_t> nodes = repeats.first;
        map<map<uint32_t,string>,string> queries = repeats.second;
        for(auto query: queries){
            // for(auto buf: queries){
            //     cout << buf.second.length() << "\t";
            // }
            // cout << endl;


            Seqhash *hasher = seqhashCreate (params.k, params.w, params.s) ;
            // cout << "factor1 " << hasher->factor1 << endl;
            // cout << "shift1 " << hasher->shift1 << endl;
            // cout << "factor2 " << hasher->factor2 << endl;
            // cout << "shift2 " << hasher->shift2 << endl;
            Modset *ms = modsetCreate (hasher, params.B, 0) ;
            Reference *ref = referenceCreate (ms, 1 << 26) ;
            string totLen = 0 ;
            int id ;
            dictAdd (ref->dict, strdup("query_sequence"), &id);
            array (ref->len, id, int) = query.second.length() ;
            totLen += query.second.length() ;
            char* seq_str = strdup(query.second.c_str());
            SeqhashRCiterator *mi = modRCiterator (ref->ms->hasher, seq_str, query.second.length()) ;
            U64 hash ; int pos ;
            while (modRCnext (mi, &hash, &pos, 0))
            { U32 index = modsetIndexFind (ref->ms, hash, TRUE) ;
            if (index)
                { if (ref->max+1 >= ref->size) die ("reference size overflow") ;
                ref->index[ref->max] = index ;
                ++ref->depth[index] ;
                ref->offset[ref->max] = pos ;
                ref->id[ref->max] = id ;
                ++ref->max ;
                }
            }
            seqhashRCiteratorDestroy (mi) ;

            // fprintf (stdout, "  %d hashes from %d reference sequences, total length %lld\n",
            //     ref->max, dictMax(ref->dict), totLen) ;
            int i ; U32 *d = &ref->depth[1] ; U32 n1 = 0, n2 = 0, nM = 0 ;
            for (i = 1 ; i <= ref->ms->max ; ++i, ++d)
                if (*d == 1) { msSetCopy1 (ref->ms, i) ; ++n1 ; }
                else if (*d == 2) { msSetCopy2 (ref->ms, i) ; ++n2 ; }
                else { msSetCopyM (ref->ms, i) ; ++nM ; }
            modsetPack (ref->ms) ;
            referencePack (ref) ;
            // cout << ref->dict->size << endl;



            map<uint32_t, set<string>> all_unique_kmers = map<uint32_t, set<string>>();


            for(auto node : query.first){
                const char* buf = node.second.c_str();
                const char** buf_seq = const_cast<const char**>(&buf);
                const char** buf_name = const_cast<const char**>(&graph->seq[node.first>>1].name);
                mm_idx_t *this_node = mm_idx_str(iopt.w, iopt.k, iopt.flag, iopt.bucket_bits, 1, buf_seq, buf_name);

                mm_tbuf_t *tbuf = mm_tbuf_init();
                mm_mapopt_update(&mopt, this_node);
                for(auto other_node: query.first){
                    if(other_node != node){
                        int j, i, n_reg;
                        // cout << graph->seq[other_node.first>>1].name << ((other_node.first%2==1) ? "-" : "+") << "\t";
                        // cout << graph->seq[node.first>>1].name << ((node.first%2==1) ? "-" : "+") << endl;
                        mm_reg1_t *reg = mm_map(this_node, other_node.second.length(), other_node.second.c_str(), &n_reg, tbuf, &mopt, graph->seq[other_node.first>>1].name);
                        for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
                            mm_reg1_t *r = &reg[j];
                            assert(r->p); // with MM_F_CIGAR, this should not be NULL
                            // printf("%s\t%d\t%d\t%d\t%c\t", graph->seq[other_node.first>>1].name, other_node.second.length(), r->qs, r->qe, "+-"[r->rev]);
                            // printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcs:Z:", this_node->seq[r->rid].name, this_node->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);

                            char *buf = NULL;
                            void *km = NULL;
                            int max_len = 0;
                            int result = mm_gen_cs(km,&buf,&max_len,const_cast<const mm_idx_t*>(this_node),const_cast<const mm_reg1_t*>(r),other_node.second.c_str(),1);
                            stringstream res = stringstream();
                            for(int i = 0 ; i < result; i++){
                                res << *(buf+i);
                            }
                            // cout << res.str() << endl;
                            set<string> unique_kmers = find_variant(other_node.second,  r->qs,  res.str(), 31);
                            for(auto k: unique_kmers){
                                // cout << k << endl;
                                if(all_unique_kmers.find(other_node.first)==all_unique_kmers.end()){
                                    all_unique_kmers[other_node.first] = set<string>();
                                }
                                all_unique_kmers[other_node.first].insert(k);
                            }
                            // for(auto k_debug: unique_kmers_debug){
                            //     cout << k_debug.first << " generate " << k_debug.second << endl;
                            // }

                            free(r->p);
                        }
                        // for(auto k: other_node_unique_kmers){
                        //     cout <<k << endl;
                        // }
                    }
                }
            }

            // Array seeds = 0 ;
            // map<uint32_t, float> shared_kmer_count = map<uint32_t,float>();
            // for(auto counting: all_unique_kmers){
            //     // cout << query.first.find(counting.first)->second << endl;
            //     SeqhashRCiterator *mi_buf = modRCiterator (ref->ms->hasher, strdup(query.first.find(counting.first)->second.c_str()), query.first.find(counting.first)->second.length());
            //     string hash ; int pos ;
            //     int buf_counter = 0;
            //     seeds = arrayReCreate (seeds, 1024, Seed) ;
            //     int missed = 0, copy[4] ; copy[1] = copy[2] = copy[3] = 0 ;
            //     while(modRCnext (mi_buf, &hash, &pos, 0)){
            //         U32 index = modsetIndexFind (ref->ms, hash, FALSE) ;
            //         Seed *s = arrayp(seeds,arrayMax(seeds),Seed) ;
            //         s->index = index ; s->pos = pos ;
            //         if (index) ++copy[msCopy(ref->ms,index)] ;
            //         else ++missed ;
            //     }
            //     seqhashRCiteratorDestroy (mi_buf) ;
            //     fprintf (stdout, "%d miss, %d copy1, %d copy2, %d multi, %.2f hit\n",
            //         missed, copy[1], copy[2], copy[3],
            //         (arrayMax(seeds)-missed)/(double)(arrayMax(seeds))) ;
            // }

            map<uint32_t, float> kmer_appearance_count = map<uint32_t,float>();
            for(auto counting: all_unique_kmers){
                uint32_t counter = 0;
                for(auto k_mer: counting.second){
                    SeqhashRCiterator *mi_buf = modRCiterator (ref->ms->hasher, strdup(k_mer.c_str()), 31);
                    U64 hash ; int pos ;
                    // int buf_counter = 0;
                    while(modRCnext (mi_buf, &hash, &pos, 0)){
                        // buf_counter++;
                        U32 index = modsetIndexFind (ref->ms, hash, FALSE) ;
                        // cout << index << endl;
                        // assert(index);
                        if(index){
                            counter++;
                        }
                    }
                    // cout << k_mer << endl;
                    // cout << k_mer.length()<< endl;
                    // cout << buf_counter << endl;
                    // assert(buf_counter<=1);

                }
                kmer_appearance_count[counting.first] = ((float)counter)/counting.second.size();
            }

            float buff_max_num = 0;
            for(auto buff: kmer_appearance_count){
                cout << graph->seq[buff.first>>1].name << "\t";
                if(buff.second>buff_max_num){
                    buff_max_num = buff.second;
                }
            }
            cout << endl;

            for(auto buff: kmer_appearance_count){
                cout << graph->seq[buff.first>>1].name << "\t" << (float)((int)(buff.second*10000))/100 << "% ";
                if(buff_max_num >= buff.second-0.000001 && buff_max_num <= buff.second+0.00001){
                    cout << "\tMax!!!";
                }
                cout << endl;
            }
        }
    }
    return map<uint32_t, vector<uint32_t>>();
}

vector<vector<uint32_t>>* repeat_resolver::count_support_reads_for_branches(vector<vector<paf_rec_str_t>>* ordered_records, map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_graph, asg_t* graph){
    vector<uint32_t>* nodes_to_begin_end = new vector<uint32_t>[graph->n_seq*2];
    map<vector<uint32_t>, int> edge_counting;
    vector<vector<uint32_t>>* edge_counting_filtered = new vector<vector<uint32_t>>();
    for(auto loop_begin: *bubble_chain_graph){
        for(auto loop_end: loop_begin.second){
            if(loop_end.second.size()>1){
                for(auto node: loop_end.second){
                    // cout << node;
                    nodes_to_begin_end[node].push_back(loop_begin.first);
                    nodes_to_begin_end[node].push_back(loop_end.first);
                }
            }
        }
    }
    for(auto i : *ordered_records){
        queue<uint32_t> supporting_branches;
        // set<uint32_t> observed_chains;
        uint32_t prev_begin = 0;
        uint32_t prev_end = 0;
        uint32_t prev_start = 0;
        uint32_t this_begin = 0;
        uint32_t this_end = 0;
        bool invalid = false;
        for(int a = 0; a < i.size(); a ++ ){
            if(prev_start == i[a].qs){
                continue;
            }
            prev_start = i[a].qs;
            if(invalid){
                break;
            }
            if(nodes_to_begin_end[i[a].tn].size() != 0){
                this_begin = nodes_to_begin_end[i[a].tn][0];
                this_end = nodes_to_begin_end[i[a].tn][1];

                            // cout << graph->seq[this_begin>>1].name << " " << graph->seq[this_end>>1].name << endl;
                                // if(this_begin!=this_end && !observed_chains.insert(this_begin>>1).second ){
                                // }else
                if (a==0||((this_begin) == (prev_end))){
                    // if(supporting_branches.size() == 0){
                        supporting_branches.push(this_begin);
                    // }
                    // if(this_end != this_begin){
                        supporting_branches.push(this_end);
                    // }
                }
                prev_begin = this_begin;
                prev_end = this_end;
            }
            // assert(found);
        }
        vector<uint32_t> current_result;
        while(!supporting_branches.empty()){
            current_result.push_back(supporting_branches.front());
            supporting_branches.pop();
        }
        if(current_result.size() > 0){
            vector<uint32_t> current_result_reversed;
            for(auto a: current_result){
                current_result_reversed.push_back(a^1);
            }
            reverse(current_result_reversed.begin(), current_result_reversed.end());
            if(edge_counting.find(current_result)==edge_counting.end()){
                edge_counting[current_result] = 0;
            }
            edge_counting[current_result]++;
            if(edge_counting.find(current_result_reversed)==edge_counting.end()){
                edge_counting[current_result_reversed] = 0;
            }
            edge_counting[current_result_reversed]++;
        }
    }
    for(auto a: edge_counting){
        // if(a.second>4 || (a.second>1 && a.first.size() > 2)){
        if(a.second>=1 || (a.second>=1 && a.first.size() >= 1)){
            edge_counting_filtered->push_back(a.first);
        }
    }
    sort(edge_counting_filtered->begin(), edge_counting_filtered->end(), [ ]( const auto& lhs, const auto& rhs )
    {
        return lhs.size() < rhs.size();
    });
    return edge_counting_filtered;
}

set<vector<uint32_t>>* repeat_resolver::cover_gaps_in_long_path(vector<vector<uint32_t>>* supported_pathes){
    set<vector<uint32_t>> short_edges;
    set<vector<uint32_t>> used_edges;
    set<vector<uint32_t>>* result = new set<vector<uint32_t>>();

    for(auto b : *supported_pathes){
        if(b.size()==2){
            short_edges.insert(b);
        }else{
            uint32_t prev_beg;
            uint32_t prev_end;
            uint32_t this_beg;
            uint32_t this_end;
            vector<uint32_t> current_result;
            bool valid = true;
            for(int a = 0; a < b.size()/2; a++){
                this_beg = b[a*2];
                this_end = b[a*2+1];
                if(a != 0 && this_beg != prev_end){
                    vector<uint32_t> wrapper;
                    wrapper.push_back(prev_end);
                    wrapper.push_back(this_beg);
                    if(short_edges.find(wrapper) != short_edges.end()){
                        current_result.push_back(prev_end);
                        current_result.push_back(this_beg);
                    }else{
                        valid = false;
                        break;
                    }
                }
                prev_beg = this_beg;
                prev_end = this_end;
                current_result.push_back(this_beg);
                current_result.push_back(this_end);
            }
            if(valid){
                for(int a = 0; a < current_result.size()/2; a++){
                    vector<uint32_t> wrapper;
                    wrapper.push_back(current_result[a*2]);
                    wrapper.push_back(current_result[a*2+1]);
                    used_edges.insert(wrapper);
                    wrapper[0] = current_result[a*2+1]^1;
                    wrapper[1] = current_result[a*2]^1;
                    used_edges.insert(wrapper);
                }
                vector<uint32_t> current_result_reversed = current_result;
                for(int a = 0; a < current_result_reversed.size(); a++){
                    current_result_reversed[a] = current_result_reversed[a] ^ 1;
                }
                reverse(current_result_reversed.begin(), current_result_reversed.end());
                if(result->find(current_result_reversed)==result->end() && result->find(current_result) == result->end()){
                    result->insert(current_result);
                    result->insert(current_result_reversed);
                }
            }
        }
    }
    for(auto a: short_edges){
        if(used_edges.find(a) == used_edges.end()){
            vector<uint32_t> current_result = a;
            vector<uint32_t> current_result_reversed = current_result;
            for(int b = 0; b < current_result_reversed.size(); b++){
                    current_result_reversed[b] = current_result_reversed[b] ^ 1;
                }
            reverse(current_result_reversed.begin(), current_result_reversed.end());
            if(result->find(current_result_reversed)==result->end() && result->find(current_result) == result->end()){
                result->insert(current_result);
            }
        }
    }
    return result;
}

set<vector<uint32_t>>* repeat_resolver::merge_long_pathes(set<vector<uint32_t>>* path_set){
    map<vector<uint32_t>,vector<vector<uint32_t>>> start_map;
    map<vector<uint32_t>,vector<vector<uint32_t>>> end_map;
    set<vector<uint32_t>> extended_path;
    set<set<uint32_t>> buffer_result;
    set<vector<uint32_t>>* result = new set<vector<uint32_t>>();
    set<vector<uint32_t>>* final_result = new set<vector<uint32_t>>();
    for(auto b : *path_set){
        if(b.size()>2){
            vector<uint32_t> start_wrapper;
            vector<uint32_t> end_wrapper;
            start_wrapper.push_back(b[0]);
            start_wrapper.push_back(b[1]);
            end_wrapper.push_back(b[b.size()-2]);
            end_wrapper.push_back(b[b.size()-1]);
            if(start_map.find(start_wrapper)==start_map.end()){
                start_map[start_wrapper] = vector<vector<uint32_t>>();
            }
            if(end_map.find(end_wrapper)==end_map.end()){
                end_map[end_wrapper] = vector<vector<uint32_t>>();
            }
            start_map[start_wrapper].push_back(b);
            end_map[end_wrapper].push_back(b);
        }
    }
    for(auto iterate_starts: start_map){
        if(end_map.find(iterate_starts.first)!=end_map.end()){
            vector<vector<uint32_t>> beginnings = iterate_starts.second;
            vector<vector<uint32_t>> endings = end_map[iterate_starts.first];
            for(auto begin: beginnings){
                extended_path.insert(begin);
                for(auto end: endings){
                    extended_path.insert(end);
                    vector<uint32_t> extended_result = end;
                    for(int i=2; i<begin.size(); i++){
                        extended_result.push_back(begin[i]);
                    }
                    vector<uint32_t> extended_result_reversed = extended_result;
                    for(int a = 0; a < extended_result_reversed.size(); a++){
                        extended_result_reversed[a] = extended_result_reversed[a] ^ 1;
                    }
                    reverse(extended_result_reversed.begin(), extended_result_reversed.end());
                    if(result->find(extended_result_reversed)==result->end() && result->find(extended_result) == result->end()){
                        result->insert(extended_result);
                    }
                }
            }
        }
    }
    for(auto b : *path_set){
        if(b.size()>2 && extended_path.find(b)==extended_path.end()){
            vector<uint32_t> extended_result_reversed = b;
            for(int a = 0; a < extended_result_reversed.size(); a++){
                extended_result_reversed[a] = extended_result_reversed[a] ^ 1;
            }
            reverse(extended_result_reversed.begin(), extended_result_reversed.end());
            if(result->find(extended_result_reversed)==result->end() && result->find(b) == result->end()){
                result->insert(b);
            }
        }
    }
    for(auto b: *result){
        set<uint32_t> current_set;
        for(auto a: b){
            current_set.insert(a>>1);
        }
        if(buffer_result.insert(current_set).second){
            final_result->insert(b);
        }
    }
    free(result);
    return final_result;
}

set<vector<uint32_t>>* repeat_resolver::recover_unlinearized_branches(set<vector<uint32_t>>* path_set, map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_graph){
    set<vector<uint32_t>>* result = new set<vector<uint32_t>>();
    map<uint32_t,map<uint32_t,set<uint32_t>>> begin_end_map = *bubble_chain_graph;
    map<uint32_t,map<uint32_t,set<uint32_t>>> end_begin_map;
    for(auto begin: begin_end_map){
        for(auto end: begin_end_map){
            if(end_begin_map.find(end.first)==end_begin_map.end()){
                end_begin_map[end.first] = map<uint32_t, set<uint32_t>>();
            }
            end_begin_map[end.first][begin.first] = begin_end_map[begin.first][end.first];
        }
    }
    for(auto p: *path_set){
        vector<uint32_t> current_result = p;
        // uint32_t counting = current_result.size();
        set<uint32_t> begin_expand_set;
        set<uint32_t> end_expand_set;
        set<uint32_t> expanded_begin_set;
        set<uint32_t> expanded_end_set;
        begin_expand_set.insert(p[0]);
        end_expand_set.insert(p[p.size()-1]);
        while(begin_expand_set.size()!=0){
            set<uint32_t> new_begin_expand_set;
            for(auto begin: begin_expand_set){
                if(end_begin_map[begin].size()!=0){
                    for(auto x: end_begin_map[begin]){
                        if(x.second.size() > 0 && expanded_begin_set.insert(x.first).second){
                            new_begin_expand_set.insert(x.first);
                            current_result.push_back(x.first);
                            current_result.push_back(begin);
                        }
                    }
                }
            }
            begin_expand_set = new_begin_expand_set;
        }
        while(end_expand_set.size()!=0){
            set<uint32_t> new_end_expand_set;
            for(auto end: end_expand_set){
                if(begin_end_map[end].size()!=0){
                    for(auto x: begin_end_map[end]){
                        if(x.second.size() > 0 && expanded_end_set.insert(x.first).second){
                            new_end_expand_set.insert(x.first);
                            current_result.push_back(end);
                            current_result.push_back(x.first);
                        }
                    }
                }
            }
            end_expand_set = new_end_expand_set;
        }
        // current_result.push_back(counting);
        result->insert(current_result);
    }
    return result;
}

void get_seperate_haplotype(set<vector<uint32_t>>* set_of_pathes, vector<vector<paf_rec_str_t>>* ordered_records, map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_graph, asg_t* graph){
    mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	mopt.min_cnt = 500;
	int n_threads = 3;
	mm_verbose = 2; // disable message output to stderr
	mm_set_opt(0, &iopt, &mopt);
	mopt.flag |= MM_F_CIGAR; // perform alignment
	mopt.flag |= MM_F_OUT_CS; // perform alignment
    params.k = 31 ;
    params.w = 31 ;
    params.s = 7 ;
    params.B = 28 ;
    Seqhash *hasher = seqhashCreate (params.k, params.w, params.s) ;

    // set<vector<uint32_t>> bubbles = set<vector<uint32_t>>();
    // set<set<uint32_t>> bubble_nodes = set<set<uint32_t>>();
    // map<uint32_t, set<string>> node_support_reads = map<uint32_t, set<string>>();
    map<uint32_t, map<U64,uint32_t>> node_discrimitive_kmers = map<uint32_t, map<U64,uint32_t>>();
    map<string, vector<paf_rec_str_t>> name_rec = map<string, vector<paf_rec_str_t>>();
    map<set<uint32_t>, set<string>> path_suppport_reads = map<set<uint32_t>, set<string>>();
    map<uint32_t, set<set<uint32_t>>> node_path_map = map<uint32_t, set<set<uint32_t>>>();

    for(auto rec: (*ordered_records)){
        name_rec[rec[0].qn] = rec;
    }

	uint32_t n_vtx = graph->n_seq * 2;
    int node_type[n_vtx];

    vector<bubble_t*> bubbles_bubbles;

    for(uint32_t i=0; i<n_vtx; i++){
        node_type[i] = 0;
        bubble_t* result = detect_bubble(graph, i);
        if(result!=nullptr){
            result->id = i;
            bubbles_bubbles.push_back(result);
        }
    }

    for (int b=0; b<bubbles_bubbles.size(); b++) {
        bubble_t* bubble = bubbles_bubbles[b];
        uint32_t bubble_beginning = bubble->begNode;
        uint32_t bubble_end = bubble->endNode;

        // cout << "start get bubble paths from " << g->seq[bubble_beginning /2].name<< " to " << g->seq[bubble_end /2].name << endl;

        vector<asg_arc_t*> arc_stack;
        vector<uint32_t> node_stack;  // DFS
        vector<uint32_t> node_vi_stack;
        node_stack.push_back(bubble_beginning);
        node_vi_stack.push_back(0);
        int stack_count = 0;
        uint32_t lastNode = bubble_beginning;
        while(!node_stack.empty()) {
            uint32_t u = node_stack.back();
            uint32_t vi = node_vi_stack.back();
            if (u==bubble_end) {
                for (int ui=0; ui<node_stack.size(); ui++) {
                    // cout << g->seq[node_stack[ui]/2].name << " ";
                    if(ui==0 || ui == node_stack.size()-1){
                        if(node_type[node_stack[ui]] != 2){
                            node_type[node_stack[ui]] = BUBBLE_END_BEGIN;
                        }
                    }else{
                        node_type[node_stack[ui]] = BUBBLE_INSIDE;
                    }
                }
                stack_count ++;
                // cout << endl;
                set<uint32_t> buf_path_node = set<uint32_t>();
                for(auto arc_buf: arc_stack){
                    // if(arc_buf->v!=bubble->endNode){
                        buf_path_node.insert(arc_buf->v);
                        node_discrimitive_kmers[arc_buf->v] = map<U64, uint32_t>();
                    // }
                }
                buf_path_node.insert(bubble->begNode);
                buf_path_node.insert(bubble->endNode);
                node_discrimitive_kmers[bubble->begNode] = map<U64, uint32_t>();
                node_discrimitive_kmers[bubble->endNode] = map<U64, uint32_t>();
                for(auto buf_node: buf_path_node){
                    if(node_path_map.find(buf_node)==node_path_map.end()){
                        node_path_map[buf_node] = set<set<uint32_t>>();
                    }
                    node_path_map[buf_node].insert(buf_path_node);
                }
                // bubble_nodes.insert(buf_path_node);
                path_suppport_reads[buf_path_node] = set<string>();
                bubble->paths_nodes.push_back(buf_path_node);
                bubble->paths.push_back(arc_stack);
                bubble->starting_arcs.insert(arc_stack.front()->v);
                bubble->ending_arcs.insert(lastNode);
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
                continue;
            }

            uint32_t num_outgoing_arcs = asg_arc_n(graph, u);
            if (vi < num_outgoing_arcs) {
                asg_arc_t *outgoing_arcs = asg_arc_a(graph, u);
                uint32_t v = outgoing_arcs[vi].v;
                node_vi_stack.back()++;
                arc_stack.push_back(outgoing_arcs + vi);
                node_stack.push_back(v);
                node_vi_stack.push_back(0);
            } else {
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
            }
            lastNode = u;
        }
        // cout << stack_count << endl;
        // cout << "finish get bubble paths from " << bubble_beginning << " to " << bubble_end << endl;
    }

    vector<bubble_t*> bubbles_pure;
    vector<bubble_t*> complex_bubbles;


    for(auto bubble: bubbles_bubbles){
        // if(node_type[bubble->begNode]!=2 && node_type[bubble->endNode]!=2 && bubble->paths_nodes.size()>20){
        if(node_type[bubble->begNode]!=2 && node_type[bubble->endNode]!=2){
            set<uint32_t> nodes;
            for(auto path:bubble->paths_nodes){
                for(auto node: path){
                    nodes.insert(node);
                }
            }
            if(nodes.size()>12){
                complex_bubbles.push_back(bubble);
            }else{
                bubbles_pure.push_back(bubble);
            }
            // cout << endl;
            // for(auto path: bubble->paths_nodes){
            //     for(auto node: path){
            //         cout << graph->seq[node>>1].name << "\t";
            //     }
            //     cout << endl;
            // }
        }else{
            delete bubble;
        }
    }

    map<uint32_t,uint32_t> nodes_count;
    for(auto bubble: complex_bubbles){
        set<uint32_t> nodes;
        for(auto path:bubble->paths_nodes){
            for(auto node: path){
                nodes.insert(node>>1);
                nodes_count[node>>1] = 0;
            }
        }
    }
    for(auto record: *ordered_records){
        for(auto single_record: record){
            if(nodes_count.find(single_record.tn>>1)!=nodes_count.end()){
                nodes_count[single_record.tn>>1]++;
            }
        }
    }
    for(auto bubble: complex_bubbles){
        cout << endl;
        set<uint32_t> nodes;
        set<uint32_t> double_nodes;
        map<uint32_t, set<string>> node_query_set;
        map<vector<uint32_t>, vector<uint32_t>> bubbles_counting;
        map<vector<uint32_t>, vector<set<U64>>> node_discrimitive_kmers = map<vector<uint32_t>, vector<set<U64>>>();
        map<string, vector<uint32_t>> reads_support_path = map<string, vector<uint32_t>>();
        map<uint32_t, set<string>> node_support_reads = map<uint32_t, set<string>>();
        for(auto path:bubble->paths_nodes){
            for(auto node: path){
                nodes.insert(node>>1);
                double_nodes.insert(node);
                double_nodes.insert(node^1);
                // nodes_count[node>>1] = 0;
            }
        }
        for(auto node: nodes){
            cout << graph->seq[node].name << "\t" << nodes_count[node] << endl;
        }
        cout << endl;

        for(auto record: *ordered_records){
            for(auto single_record: record){
                if(double_nodes.find(single_record.tn)!=double_nodes.end()){
                    if(node_query_set.find(single_record.tn)==node_query_set.end()){
                        node_query_set[single_record.tn] = set<string>();
                    }
                    node_query_set[single_record.tn].insert(single_record.qn);
                }
            }
        }

        for(auto x: double_nodes){
            // cout << graph->seq[x>>1].name << " " << (x%2==0 ? "+" : "-") << endl;
            for(auto y: double_nodes){
                if(x!=y && x!=(y^1) && x>y){
                    string x_seq = string(graph->seq[x>>1].seq);
                    uint32_t x_len = graph->seq[x>>1].len;
                    if(x%2==1){
                        x_seq = complement(string(x_seq));
                    }
                    string y_seq = string(graph->seq[y>>1].seq);
                    uint32_t y_len = graph->seq[y>>1].len;
                    if(y%2==1){
                        y_seq = complement(string(y_seq));
                    }

                    // x_seq = string("CCCCGCCCCGAAGCGCCCGCGTTGCGGTGGTTCATTTCTTCGCTTGCCCACTGGGCCTGGCAGCCTTCCGGCCCGTGGTCGTGCCTTGGCAGTCCCGCAC");
                    // y_seq = string("CCCCGCCCCGATGCCCCCGCGTTGCGGTGGTTCATTTCTTCGCTTGCCCACTGGGCCTGGCAGCCTTCCGGCCCGTGGTCGTGCCTTGGCAGTCCCGCAC");


                    // const char** buf_seq = const_cast<const char**>(&(x_seq.c_str()));
                    // const char** buf_name = const_cast<const char**>(&graph->seq[node.first>>1].name);
                    // mm_idx_t *this_node = mm_idx_str(iopt.w, iopt.k, iopt.flag, iopt.bucket_bits, 1, buf_seq, buf_name);
                    // mm_tbuf_t *tbuf = mm_tbuf_init();
                    // mm_mapopt_update(&mopt, this_node);


                    vector<uint32_t> bufff = vector<uint32_t>();
                    vector<uint32_t> bufff_count = vector<uint32_t>();
                    bufff.push_back(x);
                    bufff.push_back(y);
                    bufff_count.push_back(0);
                    bufff_count.push_back(0);
                    bubbles_counting[bufff] = bufff_count;


                    Modset *ms_x = modsetCreate (hasher, params.B, 0) ;
                    Reference *ref_x = referenceCreate (ms_x, 1 << 26) ;
                    int id_x ;
                    dictAdd (ref_x->dict, strdup("query_sequence"), &id_x);
                    array (ref_x->len, id_x, int) = x_len ;
                    SeqhashRCiterator *mi_x = modRCiterator (hasher, strdup(x_seq.c_str()), x_len) ;
                    U64 hash_x ; int pos_x ;
                    while (modRCnext (mi_x, &hash_x, &pos_x, 0))
                    { U32 index = modsetIndexFind (ref_x->ms, hash_x, TRUE) ;
                    if (index)
                        { if (ref_x->max+1 >= ref_x->size) die ("reference size overflow") ;
                        ref_x->index[ref_x->max] = index ;
                        ++ref_x->depth[index] ;
                        ref_x->offset[ref_x->max] = pos_x ;
                        ref_x->id[ref_x->max] = id_x ;
                        ++ref_x->max ;
                        }
                    }
                    seqhashRCiteratorDestroy (mi_x) ;
                    // int i_x ; U32 *d_x = &ref_x->depth[1] ;
                    // for (i_x = 1 ; i_x <= ref_x->ms->max ; ++i_x, ++d_x)
                    //     if (*d_x == 1) { msSetCopy1 (ref_x->ms, i_x) ; }
                    //     else if (*d_x == 2) { msSetCopy2 (ref_x->ms, i_x) ; }
                    //     else { msSetCopyM (ref_x->ms, i_x) ; }



                    Modset *ms_y = modsetCreate (hasher, params.B, 0) ;
                    Reference *ref_y = referenceCreate (ms_y, 1 << 26) ;
                    int id_y ;
                    dictAdd (ref_y->dict, strdup("query_sequence"), &id_y);
                    array (ref_y->len, id_y, int) = y_len ;
                    SeqhashRCiterator *mi_y = modRCiterator (hasher, strdup(y_seq.c_str()), y_len) ;
                    U64 hash_y ; int pos_y ;
                    while (modRCnext (mi_y, &hash_y, &pos_y, 0))
                    { U32 index = modsetIndexFind (ref_y->ms, hash_y, TRUE) ;
                    if (index)
                        { if (ref_y->max+1 >= ref_y->size) die ("reference size overflow") ;
                        ref_y->index[ref_y->max] = index ;
                        ++ref_y->depth[index] ;
                        ref_y->offset[ref_y->max] = pos_y ;
                        ref_y->id[ref_y->max] = id_y ;
                        ++ref_y->max ;
                        }
                    }
                    seqhashRCiteratorDestroy (mi_y) ;
                    // int i_y ; U32 *d_y = &ref_y->depth[1] ;
                    // for (i_y = 1 ; i_y <= ref_y->ms->max ; ++i_y, ++d_y)
                    //     if (*d_y == 1) { msSetCopy1 (ref_y->ms, i_y) ; }
                    //     else if (*d_y == 2) { msSetCopy2 (ref_y->ms, i_y) ; }
                    //     else { msSetCopyM (ref_y->ms, i_y) ; }


                    set<U64> y_discrimitive_kmers = set<U64>();
                    set<U64> x_discrimitive_kmers = set<U64>();

                    SeqhashRCiterator *mi_query_y = modRCiterator (hasher, strdup(y_seq.c_str()), y_len) ;
                    SeqhashRCiterator *mi_query_x = modRCiterator (hasher, strdup(x_seq.c_str()), x_len) ;

                    U64 hash_query_x ; int pos_query_x ;
                    while (modRCnext (mi_query_x, &hash_query_x, &pos_query_x, 0))
                    { U32 index = modsetIndexFind (ref_y->ms, hash_query_x, FALSE) ;
                    if (!index)
                        {
                            x_discrimitive_kmers.insert(hash_query_x);
                        }
                    }
                    seqhashRCiteratorDestroy (mi_query_x) ;

                    U64 hash_query_y ; int pos_query_y ;
                    while (modRCnext (mi_query_y, &hash_query_y, &pos_query_y, 0))
                    { U32 index = modsetIndexFind (ref_x->ms, hash_query_y, FALSE) ;
                    if (!index)
                        {
                            y_discrimitive_kmers.insert(hash_query_y);
                        }
                    }
                    seqhashRCiteratorDestroy (mi_query_y) ;
                    cout << graph->seq[x>>1].name << " " << x_discrimitive_kmers.size() << " against "  << graph->seq[y>>1].name << " " << y_discrimitive_kmers.size() << endl;


                    modsetDestroy(ms_x);
                    referenceDestroy(ref_x);
                    modsetDestroy(ms_y);
                    referenceDestroy(ref_y);
                    node_discrimitive_kmers[bufff] = vector<set<U64>>();
                    node_discrimitive_kmers[bufff].push_back(x_discrimitive_kmers);
                    node_discrimitive_kmers[bufff].push_back(y_discrimitive_kmers);
                }
            }
        }

        vector<vector<paf_rec_str_t>> filtered_ordered_records;
        for(auto query_path: *ordered_records){
            vector<paf_rec_str_t> fake_query_path = query_path;
            paf_rec_str_t fake_final;
            fake_final.qs = query_path[-1].qs+100000;
            fake_query_path.push_back(fake_final);
            uint32_t current_pos = query_path[0].qs;
            vector<paf_rec_str_t> current_result;
            vector<paf_rec_str_t> whole_result;
            for(auto single_query: query_path){
                if(double_nodes.find(single_query.tn)!=double_nodes.end()){
                    if(single_query.qs>current_pos+5000){
                        current_pos = single_query.qs;
                        if(current_result.size()>1){
                            //TODO:: Implement here!
                            map<uint32_t, set<U64>> current_node_kmers;
                            map<uint32_t, set<U64>> current_full_node_kmers;
                            paf_rec_str_t max_query = current_result[0];
                            for(auto buf_query: current_result){
                                current_full_node_kmers[buf_query.tn] = set<U64>();
                                if(buf_query.qe - buf_query.qs > max_query.qe - buf_query.qs){
                                    max_query = buf_query;
                                }
                            }

                            for(auto node_1: current_full_node_kmers){
                                for(auto node_2: current_full_node_kmers){
                                    vector<uint32_t> buf_pair;
                                    buf_pair.push_back(max(node_1.first,node_2.first));
                                    buf_pair.push_back(min(node_1.first,node_2.first));
                                    if(node_discrimitive_kmers.find(buf_pair)!=node_discrimitive_kmers.end()){
                                        for(auto kmer: node_discrimitive_kmers[buf_pair][0]){
                                            current_full_node_kmers[max(node_1.first,node_2.first)].insert(kmer);
                                        }
                                        for(auto kmer: node_discrimitive_kmers[buf_pair][1]){
                                            current_full_node_kmers[min(node_1.first,node_2.first)].insert(kmer);
                                        }
                                    }
                                }
                            }
                            for(auto node_1: current_full_node_kmers){
                                set<U64> filtered_kmers;
                                for(auto kmer: node_1.second){
                                    bool appear = false;
                                    for(auto node_2: current_full_node_kmers){
                                        if(node_1.first!=node_2.first){
                                            appear = appear || node_2.second.find(kmer) != node_2.second.end();
                                        }
                                        if(appear){
                                            break;
                                        }
                                    }
                                    if(!appear){
                                        filtered_kmers.insert(kmer);
                                    }
                                }
                                current_node_kmers[node_1.first] = filtered_kmers;
                            }

                            string node_seq_str = string(graph->seq[max_query.tn>>1].seq);
                            uint32_t buf_start = max_query.ts;
                            uint32_t buf_end = max_query.te;
                            if(max_query.tn%2==1){
                                node_seq_str = complement(node_seq_str);
                                uint32_t buf_length = max_query.tl;
                                uint32_t buffer_end = buf_length - buf_start;
                                buf_start = buf_length - buf_end;
                                buf_end = buffer_end;
                            }

                            string query_sequence = translate(node_seq_str,buf_start,max_query.seq);
                            Modset *ms_query = modsetCreate (hasher, params.B, 0) ;
                            Reference *ref_query = referenceCreate (ms_query, 1 << 26) ;
                            int id_query ;
                            dictAdd (ref_query->dict, strdup("query_sequence"), &id_query);
                            array (ref_query->len, id_query, int) = query_sequence.length() ;
                            SeqhashRCiterator *mi_query = modRCiterator (hasher, strdup(query_sequence.c_str()), query_sequence.length()) ;
                            U64 hash_query ; int pos_query ;
                            while (modRCnext (mi_query, &hash_query, &pos_query, 0))
                            { U32 index = modsetIndexFind (ref_query->ms, hash_query, TRUE) ;
                            if (index)
                                { if (ref_query->max+1 >= ref_query->size) die ("reference size overflow") ;
                                ref_query->index[ref_query->max] = index ;
                                ++ref_query->depth[index] ;
                                ref_query->offset[ref_query->max] = pos_query ;
                                ref_query->id[ref_query->max] = id_query ;
                                ++ref_query->max ;
                                }
                            }
                            seqhashRCiteratorDestroy (mi_query) ;
                            // int i_query ; U32 *d_query = &ref_query->depth[1] ;
                            // for (i_query = 1 ; i_query <= ref_query->ms->max ; ++i_query, ++d_query)
                            //     if (*d_query == 1) { msSetCopy1 (ref_query->ms, i_query) ; }
                            //     else if (*d_query == 2) { msSetCopy2 (ref_query->ms, i_query) ; }
                            //     else { msSetCopyM (ref_query->ms, i_query) ; }
                            map<uint32_t, float> node_support_ratio;

                            for(auto current_node: current_node_kmers){
                                uint32_t current_kmer_appear_count = 0;
                                for(auto kmer: current_node.second){
                                    if(modsetIndexFind (ref_query->ms, kmer, FALSE)!=0){
                                        current_kmer_appear_count++;
                                    }
                                }
                                if(current_kmer_appear_count==0){
                                    node_support_ratio[current_node.first] = 0;
                                }else{
                                    node_support_ratio[current_node.first] = ((float)current_kmer_appear_count)/((float)current_node.second.size());
                                }
                                cout << graph->seq[current_node.first>>1].name << "\t" <<  current_kmer_appear_count << "\t" <<  current_node.second.size() << ";\t";
                            }
                            cout << endl;

                            float max_ratio = 0;
                            for(auto buf_node: node_support_ratio){
                                if(buf_node.second>max_ratio){
                                    max_ratio = buf_node.second;
                                }
                            }
                            set<uint32_t> valid_nodes;
                            for(auto buf_node: node_support_ratio){
                                if(buf_node.second>=0.99*max_ratio && max_ratio!=0){
                                    valid_nodes.insert(buf_node.first);
                                }
                            }
                            for(auto buf_result: current_result){
                                if(valid_nodes.find(buf_result.tn)!=valid_nodes.end()){
                                    whole_result.push_back(buf_result);
                                }
                            }
                            modsetDestroy(ms_query);
                            referenceDestroy(ref_query);



                        }else if(current_result.size()==1){
                            whole_result.push_back(current_result[0]);
                        }
                        current_result = vector<paf_rec_str_t>();
                        current_result.push_back(single_query);
                    }else{
                        current_result.push_back(single_query);
                    }
                }
            }

            if(whole_result.size()>1){
                for(auto query: query_path){
                    cout << graph->seq[query.tn>>1].name << ",\t";
                }
                cout << endl;
                for(auto query: whole_result){
                    cout << graph->seq[query.tn>>1].name << ",\t";
                }
                cout << endl;

                filtered_ordered_records.push_back(whole_result);
            }
            for(auto buf_single_record: whole_result){
                if(double_nodes.find(buf_single_record.tn)!=double_nodes.end()){
                    if(node_support_reads.find(buf_single_record.tn)==node_support_reads.end()){
                        node_support_reads[buf_single_record.tn] = set<string>();
                    }
                    if(reads_support_path.find(buf_single_record.qn)==reads_support_path.end()){
                        reads_support_path[buf_single_record.qn] = vector<uint32_t>();
                    }
                    node_support_reads[buf_single_record.tn].insert(buf_single_record.qn);
                    reads_support_path[buf_single_record.qn].push_back(buf_single_record.tn);
                }
            }
        }
        // cout << endl;
        // cout << endl;
        // cout << endl;
        // for(auto read: reads_support_path){
        //     if(read.second.size()>1){
        //         for(auto node: read.second){
        //             cout << graph->seq[node>>1].name<< ",\t";
        //         }
        //     }
        //     cout << endl;
        // }

        cout << endl;
        cout << endl;
        for(auto node: node_support_reads){
            if(node.second.size()>1){
                cout << graph->seq[node.first>>1].name << "\t" << node.second.size() << endl;
            }
        }

        for(auto read: reads_support_path){
            if(read.second.size()>1){
                for(auto node: read.second){
                    cout << graph->seq[node>>1].name << ", ";
                }
                cout << endl;
            }
        }

        // set<set<uint32_t>> pathes;
        // for(auto path: bubble->paths_nodes){
        //     bool appears = true;
        //     for(auto node: path){
        //         if(node_support_reads.find(node)==node_support_reads.end() && node_support_reads.find(node^1)==node_support_reads.end()){
        //             appears = false;
        //             break;
        //         }
        //     }
        //     if(appears){
        //         pathes.insert(path);
        //     }
        // }
        // for(auto path: pathes){
        //     for(auto node: path){
        //         cout << graph->seq[node>>1].name << ", ";
        //     }
        //     cout << endl;
        // }

        map<uint32_t, map<uint32_t, uint32_t>> node_connections;
        map<vector<uint32_t>, uint32_t> connection_count;
        for(auto node: node_support_reads){
            node_connections[node.first] = map<uint32_t,uint32_t>();
            node_connections[node.first^1] = map<uint32_t,uint32_t>();
            for(auto node_2: node_support_reads){
                if(node.first!=node_2.first && node.first!=(node_2.first^1)){
                    node_connections[node.first][node_2.first] = 0;
                    vector<uint32_t> buf_to_insert;
                    buf_to_insert.push_back(max(node.first,node_2.first));
                    buf_to_insert.push_back(min(node.first,node_2.first));
                    connection_count[buf_to_insert]=0;
                    node_connections[node.first^1][node_2.first^1] = 0;
                    vector<uint32_t> buf_to_insert_2;
                    buf_to_insert_2.push_back(max(node.first,node_2.first)^1);
                    buf_to_insert_2.push_back(min(node.first,node_2.first)^1);
                    connection_count[buf_to_insert_2]=0;
                }
            }
        }
        for(auto read: reads_support_path){
            for(auto node_1: read.second){
                for(auto node_2: read.second){
                    if(node_1!=node_2 && node_1!=(node_2^1)){
                        bool linear = false;
                        for(auto path: bubble->paths_nodes){
                            bool found1 = false;
                            bool found1_c = false;
                            bool found2 = false;
                            bool found2_c = false;
                            for(auto node: path){
                                if(node_1==node){
                                    found1 = true;
                                }else if(node_2==node){
                                    found2 = true;
                                }else if((node_1^1)==node){
                                    found1_c = true;
                                }else if((node_2^1)==node){
                                    found2_c = true;
                                }
                                if((found1 && found2)||(found1_c && found2_c)){
                                    linear = true;
                                    break;
                                }
                            }
                            if(linear){
                                break;
                            }
                        }
                        if(linear){
                            node_connections[node_1][node_2]++;
                            vector<uint32_t> buf_to_insert;
                            buf_to_insert.push_back(max(node_1,node_2));
                            buf_to_insert.push_back(min(node_1,node_2));
                            connection_count[buf_to_insert]++;
                            node_connections[node_1^1][node_2^1]++;
                            vector<uint32_t> buf_to_insert_2;
                            buf_to_insert_2.push_back(max(node_1,node_2)^1);
                            buf_to_insert_2.push_back(min(node_1,node_2)^1);
                            connection_count[buf_to_insert_2]++;
                        }
                    }
                }
            }
        }
        // map<uint32_t, set<vector<uint32_t>>> count_connection;
        // uint32_t max_count = 0;
        // for(auto connection: connection_count){
        //     if(connection.second>0){
        //         if(count_connection.find(connection.second)==count_connection.end()){
        //             count_connection[connection.second] = set<vector<uint32_t>>();
        //         }
        //         count_connection[connection.second].insert(connection.first);
        //         max_count = max(max_count, connection.second);
        //         cout << graph->seq[connection.first[0]>>1].name << ", " << graph->seq[connection.first[1]>>1].name << "\t" << connection.second << endl;
        //     }
        // }

        set<set<uint32_t>> set_of_path;
        uint32_t round = 10;
        while(round--){
            map<uint32_t, set<vector<uint32_t>>> count_connection;
            uint32_t max_count = 0;
            for(auto connection: connection_count){
                if(connection.second>0){
                    if(count_connection.find(connection.second)==count_connection.end()){
                        count_connection[connection.second] = set<vector<uint32_t>>();
                    }
                    count_connection[connection.second].insert(connection.first);
                    max_count = max(max_count, connection.second);
                    cout << graph->seq[connection.first[0]>>1].name << ", " << graph->seq[connection.first[1]>>1].name << "\t" << connection.second << endl;
                }
            }
            uint32_t min_flow = 10000;
            set<uint32_t> current_path;
            for(uint32_t i = max_count; i > 2; i--){
                if(count_connection.find(i)!=count_connection.end()){
                    for(auto connection: count_connection[i]){
                        if(current_path.size()==0){
                            current_path.insert(connection[0]>>1);
                            current_path.insert(connection[1]>>1);
                            min_flow = min(min_flow,i);
                        }else{
                            set<uint32_t> buf_new_path = current_path;
                            // if(buf_new_path.find(connection[0]>>1)!=buf_new_path.end() || buf_new_path.find(connection[1]>>1)!=buf_new_path.end()){
                                buf_new_path.insert(connection[0]>>1);
                                buf_new_path.insert(connection[1]>>1);
                            // }
                            bool valid = false;
                            for(auto path: bubble->paths_nodes){
                                bool linear = true;
                                for(auto node: buf_new_path){
                                    if(path.find(node<<1)==path.end() && path.find((node<<1)^1)==path.end()){
                                        linear = false;
                                    }
                                }
                                if(linear){
                                    valid = true;
                                    break;
                                }
                            }

                            if(valid){
                                min_flow = min(min_flow,i);
                                current_path = buf_new_path;
                            }
                        }
                    }
                }
            }
            set_of_path.insert(current_path);
            for(auto node_1: current_path){
                for(auto node_2: current_path){
                    vector<uint32_t> buf_to_insert_1;
                    buf_to_insert_1.push_back(max(node_1<<1,node_2<<1));
                    buf_to_insert_1.push_back(min(node_1<<1,node_2<<1));
                    if(connection_count.find(buf_to_insert_1)!=connection_count.end()){
                        if(connection_count[buf_to_insert_1] > min_flow){
                            connection_count[buf_to_insert_1]-=min_flow;
                        }else{
                            connection_count[buf_to_insert_1] = 0;
                        }
                    }
                    vector<uint32_t> buf_to_insert_2;
                    buf_to_insert_2.push_back(max((node_1<<1)^1,node_2<<1));
                    buf_to_insert_2.push_back(min((node_1<<1)^1,node_2<<1));
                    if(connection_count.find(buf_to_insert_2)!=connection_count.end()){
                        if(connection_count[buf_to_insert_2] > min_flow){
                            connection_count[buf_to_insert_2]-=min_flow;
                        }else{
                            connection_count[buf_to_insert_2] = 0;
                        }
                    }
                    vector<uint32_t> buf_to_insert_3;
                    buf_to_insert_3.push_back(max(node_1<<1,(node_2<<1)^1));
                    buf_to_insert_3.push_back(min(node_1<<1,(node_2<<1)^1));
                    if(connection_count.find(buf_to_insert_3)!=connection_count.end()){
                        if(connection_count[buf_to_insert_3] > min_flow){
                            connection_count[buf_to_insert_3]-=min_flow;
                        }else{
                            connection_count[buf_to_insert_3] = 0;
                        }
                    }
                    vector<uint32_t> buf_to_insert_4;
                    buf_to_insert_4.push_back(max((node_1<<1)^1,(node_2<<1)^1));
                    buf_to_insert_4.push_back(min((node_1<<1)^1,(node_2<<1)^1));
                    if(connection_count.find(buf_to_insert_4)!=connection_count.end()){
                        if(connection_count[buf_to_insert_4] > min_flow){
                            connection_count[buf_to_insert_4]-=min_flow;
                        }else{
                            connection_count[buf_to_insert_4] = 0;
                        }
                    }
                }
            }
        }

        set<set<uint32_t>> clean_set_path;
        for(auto to_delete_path: set_of_path){
            bool needed = true;
            for(auto path: set_of_path){
                if(to_delete_path.size() < path.size()){
                    bool not_contained = false;
                    for(auto node: to_delete_path){
                        not_contained = (path.find(node)==path.end());
                        if(not_contained){
                            break;
                        }
                    }
                    needed = not_contained && needed;
                    if(!needed){
                        break;
                    }
                }
            }
            if(needed){
                clean_set_path.insert(to_delete_path);
            }
        }

        for(auto path: clean_set_path){
            for(auto node: path){
                cout << graph->seq[node].name << ", ";
            }
            cout << endl;
        }
    }
                    // uint32_t num_outgoing_arcs_x = asg_arc_n(graph, x);
                    // asg_arc_t *outgoing_arcs_x = asg_arc_a(graph, x);
                    // uint32_t num_outgoing_arcs_y = asg_arc_n(graph, y);
                    // asg_arc_t *outgoing_arcs_y = asg_arc_a(graph, y);
                    // uint32_t num_incoming_arcs_x = asg_arc_n(graph, x^1);
                    // asg_arc_t *incoming_arcs_x = asg_arc_a(graph, x^1);
                    // uint32_t num_incoming_arcs_y = asg_arc_n(graph, y^1);
                    // asg_arc_t *incoming_arcs_y = asg_arc_a(graph, y^1);
                    // bool flag = false;
                    // for(int i=0; i<num_outgoing_arcs_x; i++){
                    //     for(int j=0; j<num_outgoing_arcs_y; j++){
                    //         if(outgoing_arcs_x[i].v == outgoing_arcs_y[j].v){
                    //             flag = true;
                    //             break;
                    //         }
                    //     }
                    //     if(flag){
                    //         break;
                    //     }
                    // }
                    // bool flag2 = false;
                    // for(int i=0; i<num_incoming_arcs_x; i++){
                    //     for(int j=0; j<num_incoming_arcs_y; j++){
                    //         if(incoming_arcs_x[i].v == incoming_arcs_y[j].v){
                    //             flag2 = true;
                    //             break;
                    //         }
                    //     }
                    //     if(flag2){
                    //         break;
                    //     }
                    // }

        //             if(flag2||flag){
        //                 string x_seq = string(graph->seq[x>>1].seq);
        //                 uint32_t x_len = graph->seq[x>>1].len;
        //                 if(x%2==1){
        //                     x_seq = complement(string(x_seq));
        //                 }
        //                 string y_seq = string(graph->seq[y>>1].seq);
        //                 uint32_t y_len = graph->seq[y>>1].len;
        //                 if(y%2==1){
        //                     y_seq = complement(string(y_seq));
        //                 }

        //                 // x_seq = string("CCCCGCCCCGAAGCGCCCGCGTTGCGGTGGTTCATTTCTTCGCTTGCCCACTGGGCCTGGCAGCCTTCCGGCCCGTGGTCGTGCCTTGGCAGTCCCGCAC");
        //                 // y_seq = string("CCCCGCCCCGATGCCCCCGCGTTGCGGTGGTTCATTTCTTCGCTTGCCCACTGGGCCTGGCAGCCTTCCGGCCCGTGGTCGTGCCTTGGCAGTCCCGCAC");


        //                 // const char** buf_seq = const_cast<const char**>(&(x_seq.c_str()));
        //                 // const char** buf_name = const_cast<const char**>(&graph->seq[node.first>>1].name);
        //                 // mm_idx_t *this_node = mm_idx_str(iopt.w, iopt.k, iopt.flag, iopt.bucket_bits, 1, buf_seq, buf_name);
        //                 // mm_tbuf_t *tbuf = mm_tbuf_init();
        //                 // mm_mapopt_update(&mopt, this_node);


        //                 vector<uint32_t> bufff = vector<uint32_t>();
        //                 vector<uint32_t> bufff_count = vector<uint32_t>();
        //                 bufff.push_back(max(x,y));
        //                 bufff.push_back(min(x,y));
        //                 bufff_count.push_back(0);
        //                 bufff_count.push_back(0);
        //                 bubbles_counting[bufff] = bufff_count;


        //                 Modset *ms_x = modsetCreate (hasher, params.B, 0) ;
        //                 Reference *ref_x = referenceCreate (ms_x, 1 << 26) ;
        //                 int id_x ;
        //                 dictAdd (ref_x->dict, strdup("query_sequence"), &id_x);
        //                 array (ref_x->len, id_x, int) = x_len ;
        //                 SeqhashRCiterator *mi_x = modRCiterator (hasher, strdup(x_seq.c_str()), x_len) ;
        //                 U64 hash_x ; int pos_x ;
        //                 while (modRCnext (mi_x, &hash_x, &pos_x, 0))
        //                 { U32 index = modsetIndexFind (ref_x->ms, hash_x, TRUE) ;
        //                 if (index)
        //                     { if (ref_x->max+1 >= ref_x->size) die ("reference size overflow") ;
        //                     ref_x->index[ref_x->max] = index ;
        //                     ++ref_x->depth[index] ;
        //                     ref_x->offset[ref_x->max] = pos_x ;
        //                     ref_x->id[ref_x->max] = id_x ;
        //                     ++ref_x->max ;
        //                     }
        //                 }
        //                 seqhashRCiteratorDestroy (mi_x) ;
        //                 // int i_x ; U32 *d_x = &ref_x->depth[1] ;
        //                 // for (i_x = 1 ; i_x <= ref_x->ms->max ; ++i_x, ++d_x)
        //                 //     if (*d_x == 1) { msSetCopy1 (ref_x->ms, i_x) ; }
        //                 //     else if (*d_x == 2) { msSetCopy2 (ref_x->ms, i_x) ; }
        //                 //     else { msSetCopyM (ref_x->ms, i_x) ; }



        //                 Modset *ms_y = modsetCreate (hasher, params.B, 0) ;
        //                 Reference *ref_y = referenceCreate (ms_y, 1 << 26) ;
        //                 int id_y ;
        //                 dictAdd (ref_y->dict, strdup("query_sequence"), &id_y);
        //                 array (ref_y->len, id_y, int) = y_len ;
        //                 SeqhashRCiterator *mi_y = modRCiterator (hasher, strdup(y_seq.c_str()), y_len) ;
        //                 U64 hash_y ; int pos_y ;
        //                 while (modRCnext (mi_y, &hash_y, &pos_y, 0))
        //                 { U32 index = modsetIndexFind (ref_y->ms, hash_y, TRUE) ;
        //                 if (index)
        //                     { if (ref_y->max+1 >= ref_y->size) die ("reference size overflow") ;
        //                     ref_y->index[ref_y->max] = index ;
        //                     ++ref_y->depth[index] ;
        //                     ref_y->offset[ref_y->max] = pos_y ;
        //                     ref_y->id[ref_y->max] = id_y ;
        //                     ++ref_y->max ;
        //                     }
        //                 }
        //                 seqhashRCiteratorDestroy (mi_y) ;
        //                 // int i_y ; U32 *d_y = &ref_y->depth[1] ;
        //                 // for (i_y = 1 ; i_y <= ref_y->ms->max ; ++i_y, ++d_y)
        //                 //     if (*d_y == 1) { msSetCopy1 (ref_y->ms, i_y) ; }
        //                 //     else if (*d_y == 2) { msSetCopy2 (ref_y->ms, i_y) ; }
        //                 //     else { msSetCopyM (ref_y->ms, i_y) ; }


        //                 map<U64,uint32_t> y_discrimitive_kmers = map<U64,uint32_t>();
        //                 map<U64,uint32_t> x_discrimitive_kmers = map<U64,uint32_t>();

        //                 SeqhashRCiterator *mi_query_y = modRCiterator (hasher, strdup(y_seq.c_str()), y_len) ;
        //                 SeqhashRCiterator *mi_query_x = modRCiterator (hasher, strdup(x_seq.c_str()), x_len) ;

        //                 U64 hash_query_x ; int pos_query_x ;
        //                 while (modRCnext (mi_query_x, &hash_query_x, &pos_query_x, 0))
        //                 { U32 index = modsetIndexFind (ref_y->ms, hash_query_x, FALSE) ;
        //                 if (!index)
        //                     {
        //                         x_discrimitive_kmers[hash_query_x] = pos_query_x;
        //                     }
        //                 }
        //                 seqhashRCiteratorDestroy (mi_query_x) ;

        //                 U64 hash_query_y ; int pos_query_y ;
        //                 while (modRCnext (mi_query_y, &hash_query_y, &pos_query_y, 0))
        //                 { U32 index = modsetIndexFind (ref_x->ms, hash_query_y, FALSE) ;
        //                 if (!index)
        //                     {
        //                         y_discrimitive_kmers[hash_query_y] = pos_query_y;
        //                     }
        //                 }
        //                 seqhashRCiteratorDestroy (mi_query_y) ;


        //                 modsetDestroy(ms_x);
        //                 referenceDestroy(ref_x);
        //                 modsetDestroy(ms_y);
        //                 referenceDestroy(ref_y);
        //                 // set<string> y_discrimitive_kmers = set<string>();
        //                 // set<string> x_discrimitive_kmers = set<string>();

        //                 // set<string> x_kmers = set<string>();
        //                 // set<string> y_kmers = set<string>();

        //                 // for(int i = 0; i < x_seq.length()-31; i++){
        //                 //     x_kmers.insert(x_seq.substr(i,31));
        //                 // }
        //                 // for(int q = 0; q < y_seq.length()-31; q++){
        //                 //     y_kmers.insert(y_seq.substr(q,31));
        //                 // }
        //                 // for(auto k_mer: x_kmers){
        //                 //     if(y_kmers.find(k_mer)==y_kmers.end()){
        //                 //         x_discrimitive_kmers.insert(k_mer);

        //                 //     }
        //                 // }
        //                 // for(auto k_mer: y_kmers){
        //                 //     if(x_kmers.find(k_mer)==x_kmers.end()){
        //                 //         y_discrimitive_kmers.insert(k_mer);
        //                 //     }
        //                 // }

        //                 // cout << x_seq << endl;
        //                 // cout << y_seq << endl;
        //                 cout << graph->seq[x>>1].name << " kmers: " << x_discrimitive_kmers.size() << endl;
        //                 // for(auto k_mer : x_discrimitive_kmers){
        //                 //     cout << k_mer << endl;
        //                 // }
        //                 // cout << counter_valid << " against " << valid_kmers.size() << endl;
        //                 cout << graph->seq[y>>1].name << " kmers: " << y_discrimitive_kmers.size() << endl;
        //                 // for(auto k_mer : y_discrimitive_kmers){
        //                 //     cout << k_mer << endl;
        //                 // }
        //                 node_discrimitive_kmers[x] = x_discrimitive_kmers;
        //                 node_discrimitive_kmers[y] = y_discrimitive_kmers;

        //                 set<string> joint_set;
        //                 for(auto n: node_query_set[x]){
        //                     joint_set.insert(n);
        //                 }
        //                 for(auto n: node_query_set[y]){
        //                     joint_set.insert(n);
        //                 }
        //                 uint32_t node_x = x;
        //                 uint32_t node_y = y;
        //                 for(auto query: joint_set){
        //                     vector<paf_rec_str_t> record = name_rec[query];
        //                     float node_x_support_record_ratio = 0;
        //                     float node_y_support_record_ratio = 0;
        //                     bool have_x = false;
        //                     bool have_y = false;
        //                     for(auto buf_record: record){
        //                         string query_sequence = string();
        //                         if(buf_record.tn==node_x && buf_record.tl >= graph->seq[node_x].len*0.2){
        //                             have_x = true;
        //                             string x_seq_str = string(graph->seq[node_x>>1].seq);
        //                             if(node_x%2==1){
        //                                 x_seq_str = complement(x_seq_str);
        //                             }
        //                             query_sequence = translate(x_seq_str,buf_record.ts,buf_record.seq);
        //                         }else if(buf_record.tn == node_y && buf_record.tl >= graph->seq[node_y].len*0.2){
        //                             have_y = true;
        //                             string y_seq_str = string(graph->seq[node_y>>1].seq);
        //                             if(node_y%2==1){
        //                                 y_seq_str = complement(y_seq_str);
        //                             }
        //                             query_sequence = translate(y_seq_str,buf_record.ts,buf_record.seq);
        //                         }
        //                         if(query_sequence.length()>0){

        //                             // set<string> query_kmers = set<string>();

        //                             // for(int i = 0; i < query_sequence.length()-31; i++){
        //                             //     query_kmers.insert(query_sequence.substr(i,31));
        //                             // }

        //                             // // int i ; U32 *d = &ref_query->depth[1] ; U32 n1 = 0, n2 = 0, nM = 0 ;
        //                             // // for (i = 1 ; i <= ref_query->ms->max ; ++i, ++d)
        //                             // //     if (*d == 1) { msSetCopy1 (ref_query->ms, i) ; ++n1 ; }
        //                             // //     else if (*d == 2) { msSetCopy2 (ref_query->ms, i) ; ++n2 ; }
        //                             // //     else { msSetCopyM (ref_query->ms, i) ; ++nM ; }
        //                             // // modsetPack (ref_query->ms) ;
        //                             // // referencePack (ref_query) ;


        //                             // for(auto kmer_x: node_discrimitive_kmers[node_x]){
        //                             //     if(query_kmers.find(kmer_x)!=query_kmers.end()){
        //                             //         node_x_support_record+=1;
        //                             //     }
        //                             // }
        //                             // for(auto kmer_y: node_discrimitive_kmers[node_y]){
        //                             //     if(query_kmers.find(kmer_y)!=query_kmers.end()){
        //                             //         node_y_support_record+=1;
        //                             //     }
        //                             // }

        //                             Modset *ms_query = modsetCreate (hasher, params.B, 0) ;
        //                             Reference *ref_query = referenceCreate (ms_query, 1 << 26) ;
        //                             int id_query ;
        //                             dictAdd (ref_query->dict, strdup("query_sequence"), &id_query);
        //                             array (ref_query->len, id_query, int) = query_sequence.length() ;
        //                             SeqhashRCiterator *mi_query = modRCiterator (hasher, strdup(query_sequence.c_str()), query_sequence.length()) ;
        //                             U64 hash_query ; int pos_query ;
        //                             while (modRCnext (mi_query, &hash_query, &pos_query, 0))
        //                             { U32 index = modsetIndexFind (ref_query->ms, hash_query, TRUE) ;
        //                             if (index)
        //                                 { if (ref_query->max+1 >= ref_query->size) die ("reference size overflow") ;
        //                                 ref_query->index[ref_query->max] = index ;
        //                                 ++ref_query->depth[index] ;
        //                                 ref_query->offset[ref_query->max] = pos_query ;
        //                                 ref_query->id[ref_query->max] = id_query ;
        //                                 ++ref_query->max ;
        //                                 }
        //                             }
        //                             seqhashRCiteratorDestroy (mi_query) ;
        //                             // int i_query ; U32 *d_query = &ref_query->depth[1] ;
        //                             // for (i_query = 1 ; i_query <= ref_query->ms->max ; ++i_query, ++d_query)
        //                             //     if (*d_query == 1) { msSetCopy1 (ref_query->ms, i_query) ; }
        //                             //     else if (*d_query == 2) { msSetCopy2 (ref_query->ms, i_query) ; }
        //                             //     else { msSetCopyM (ref_query->ms, i_query) ; }
        //                             uint32_t node_x_support_record = 0;
        //                             uint32_t node_y_support_record = 0;
        //                             for(auto kmer_x: node_discrimitive_kmers[node_x]){
        //                                 if(modsetIndexFind (ref_query->ms, kmer_x.first, FALSE)!=0){
        //                                     node_x_support_record+=1;
        //                                 }
        //                             }
        //                             for(auto kmer_y: node_discrimitive_kmers[node_y]){
        //                                 if(modsetIndexFind (ref_query->ms, kmer_y.first, FALSE)!=0){
        //                                     node_y_support_record+=1;
        //                                 }
        //                             }

        //                             modsetDestroy(ms_query);
        //                             referenceDestroy(ref_query);


        //                             if(buf_record.tn==node_x && node_discrimitive_kmers[node_x].size()!=0){
        //                                 node_x_support_record_ratio += (float)node_x_support_record/node_discrimitive_kmers[node_x].size();
        //                             }
        //                             if(buf_record.tn==node_y &&node_discrimitive_kmers[node_y].size()!=0){
        //                                 node_y_support_record_ratio += (float)node_y_support_record/node_discrimitive_kmers[node_y].size();
        //                             }
        //                             // cout << node_x_support_record_ratio << " against " << node_y_support_record_ratio << endl;
        //                             // cout << graph->seq[node_x].name << " against " << graph->seq[node_y].name << endl;

        //                         }
        //                     }

        //                     if(node_x_support_record_ratio > node_y_support_record_ratio*1.5){
        //                         if(node_support_reads.find(node_x)==node_support_reads.end()){
        //                             node_support_reads[node_x] = set<string>();
        //                         }
        //                         if(have_x && have_y){
        //                             bubbles_counting[bufff][0] += 1;
        //                         }
        //                         node_support_reads[node_x].insert(string(record[0].qn));
        //                     }else if(node_y_support_record_ratio > node_x_support_record_ratio*1.5){
        //                         if(node_support_reads.find(node_y)==node_support_reads.end()){
        //                             node_support_reads[node_y] = set<string>();
        //                         }
        //                         if(have_x && have_y){
        //                             bubbles_counting[bufff][0] += 1;
        //                         }
        //                         bubbles_counting[bufff][1] += 1;
        //                         node_support_reads[node_y].insert(string(record[0].qn));
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }
        // for(auto pair: bubbles_counting){
        //     cout << graph->seq[pair.first[0]>>1].name << " " << pair.second[0] << " against " << graph->seq[pair.first[1]>>1].name << " " << pair.second[1] << endl;
        // }
        // map<string, set<uint32_t>> read_node_map;
        // for(auto node: node_support_reads){
        //     for(auto query_name: node.second){
        //         if(read_node_map.find(query_name)==read_node_map.end()){
        //             read_node_map[query_name] = set<uint32_t>();
        //         }
        //         read_node_map[query_name].insert(node.first);
        //     }
        // }
        // for(auto read: read_node_map){
        //     for(auto node: read.second){
        //         cout << graph->seq[node>>1].name << ",\t";
        //     }
        //     cout << endl;
    //     }
    // }



    set<set<uint32_t>> seen;
    for(auto bubble: bubbles_pure){
        for(int i = 0; i < bubble->paths_nodes.size(); i++){
            for(int j = 0; j < bubble->paths_nodes.size(); j++){
                if(i!=j){
                    set<uint32_t> path_nodes_1 = bubble->paths_nodes[i];
                    set<uint32_t> path_nodes_2 = bubble->paths_nodes[j];
                    cout << i << "\t" << j << endl;
                    for(auto x: path_nodes_1){
                        for(auto y: path_nodes_2){
                            set<uint32_t> buff_seen;
                            buff_seen.insert(x);
                            buff_seen.insert(y);

                            bool comparable = (x>>1!= y>>1) && seen.insert(buff_seen).second;
                            bool flag = false;
                            bool flag2 = false;
                            if(comparable){
                                uint32_t num_outgoing_arcs_x = asg_arc_n(graph, x);
                                asg_arc_t *outgoing_arcs_x = asg_arc_a(graph, x);
                                uint32_t num_outgoing_arcs_y = asg_arc_n(graph, y);
                                asg_arc_t *outgoing_arcs_y = asg_arc_a(graph, y);
                                uint32_t num_incoming_arcs_x = asg_arc_n(graph, x^1);
                                asg_arc_t *incoming_arcs_x = asg_arc_a(graph, x^1);
                                uint32_t num_incoming_arcs_y = asg_arc_n(graph, y^1);
                                asg_arc_t *incoming_arcs_y = asg_arc_a(graph, y^1);
                                for(int i=0; i<num_outgoing_arcs_x; i++){
                                    for(int j=0; j<num_outgoing_arcs_y; j++){
                                        if(outgoing_arcs_x[i].v == outgoing_arcs_y[j].v){
                                            flag = true;
                                            break;
                                        }
                                    }
                                    if(flag){
                                        break;
                                    }
                                }
                                for(int i=0; i<num_incoming_arcs_x; i++){
                                    for(int j=0; j<num_incoming_arcs_y; j++){
                                        if(incoming_arcs_x[i].v == incoming_arcs_y[j].v){
                                            flag2 = true;
                                            break;
                                        }
                                    }
                                    if(flag2){
                                        break;
                                    }
                                }
                            }
                            comparable = comparable && (flag||flag2);
                            if(comparable){
                                string x_seq = string(graph->seq[x>>1].seq);
                                uint32_t x_len = graph->seq[x>>1].len;
                                if(x%2==1){
                                    x_seq = complement(string(x_seq));
                                }
                                string y_seq = string(graph->seq[y>>1].seq);
                                uint32_t y_len = graph->seq[y>>1].len;
                                if(y%2==1){
                                    y_seq = complement(string(y_seq));
                                }


                                Modset *ms_x = modsetCreate (hasher, params.B, 0) ;
                                Reference *ref_x = referenceCreate (ms_x, 1 << 26) ;
                                int id_x ;
                                dictAdd (ref_x->dict, strdup("query_sequence"), &id_x);
                                array (ref_x->len, id_x, int) = x_len ;
                                SeqhashRCiterator *mi_x = modRCiterator (hasher, strdup(x_seq.c_str()), x_len) ;
                                U64 hash_x ; int pos_x ;
                                while (modRCnext (mi_x, &hash_x, &pos_x, 0))
                                { U32 index = modsetIndexFind (ref_x->ms, hash_x, TRUE) ;
                                if (index)
                                    { if (ref_x->max+1 >= ref_x->size) die ("reference size overflow") ;
                                    ref_x->index[ref_x->max] = index ;
                                    ++ref_x->depth[index] ;
                                    ref_x->offset[ref_x->max] = pos_x ;
                                    ref_x->id[ref_x->max] = id_x ;
                                    ++ref_x->max ;
                                    }
                                }
                                seqhashRCiteratorDestroy (mi_x) ;
                                // int i_x ; U32 *d_x = &ref_x->depth[1] ;
                                // for (i_x = 1 ; i_x <= ref_x->ms->max ; ++i_x, ++d_x)
                                //     if (*d_x == 1) { msSetCopy1 (ref_x->ms, i_x) ; }
                                //     else if (*d_x == 2) { msSetCopy2 (ref_x->ms, i_x) ; }
                                //     else { msSetCopyM (ref_x->ms, i_x) ; }



                                Modset *ms_y = modsetCreate (hasher, params.B, 0) ;
                                Reference *ref_y = referenceCreate (ms_y, 1 << 26) ;
                                int id_y ;
                                dictAdd (ref_y->dict, strdup("query_sequence"), &id_y);
                                array (ref_y->len, id_y, int) = y_len ;
                                SeqhashRCiterator *mi_y = modRCiterator (hasher, strdup(y_seq.c_str()), y_len) ;
                                U64 hash_y ; int pos_y ;
                                while (modRCnext (mi_y, &hash_y, &pos_y, 0))
                                { U32 index = modsetIndexFind (ref_y->ms, hash_y, TRUE) ;
                                if (index)
                                    { if (ref_y->max+1 >= ref_y->size) die ("reference size overflow") ;
                                    ref_y->index[ref_y->max] = index ;
                                    ++ref_y->depth[index] ;
                                    ref_y->offset[ref_y->max] = pos_y ;
                                    ref_y->id[ref_y->max] = id_y ;
                                    ++ref_y->max ;
                                    }
                                }
                                seqhashRCiteratorDestroy (mi_y) ;
                                // int i_y ; U32 *d_y = &ref_y->depth[1] ;
                                // for (i_y = 1 ; i_y <= ref_y->ms->max ; ++i_y, ++d_y)
                                //     if (*d_y == 1) { msSetCopy1 (ref_y->ms, i_y) ; }
                                //     else if (*d_y == 2) { msSetCopy2 (ref_y->ms, i_y) ; }
                                //     else { msSetCopyM (ref_y->ms, i_y) ; }


                                map<U64,uint32_t> y_discrimitive_kmers = map<U64,uint32_t>();
                                map<U64,uint32_t> x_discrimitive_kmers = map<U64,uint32_t>();

                                SeqhashRCiterator *mi_query_y = modRCiterator (hasher, strdup(y_seq.c_str()), y_len) ;
                                SeqhashRCiterator *mi_query_x = modRCiterator (hasher, strdup(x_seq.c_str()), x_len) ;

                                U64 hash_query_x ; int pos_query_x ;
                                while (modRCnext (mi_query_x, &hash_query_x, &pos_query_x, 0))
                                { U32 index = modsetIndexFind (ref_y->ms, hash_query_x, FALSE) ;
                                if (!index)
                                    {
                                        x_discrimitive_kmers[hash_query_x] = pos_query_x;
                                    }
                                }
                                seqhashRCiteratorDestroy (mi_query_x) ;

                                U64 hash_query_y ; int pos_query_y ;
                                while (modRCnext (mi_query_y, &hash_query_y, &pos_query_y, 0))
                                { U32 index = modsetIndexFind (ref_x->ms, hash_query_y, FALSE) ;
                                if (!index)
                                    {
                                        y_discrimitive_kmers[hash_query_y] = pos_query_y;
                                    }
                                }
                                seqhashRCiteratorDestroy (mi_query_y) ;


                                modsetDestroy(ms_x);
                                referenceDestroy(ref_x);
                                modsetDestroy(ms_y);
                                referenceDestroy(ref_y);


                                // cout << graph->seq[x>>1].name << " kmers: " << x_discrimitive_kmers.size() << endl;
                                for(auto k_mer : x_discrimitive_kmers){
                                    node_discrimitive_kmers[x].insert(k_mer);
                                }
                                // cout << counter_valid << " against " << valid_kmers.size() << endl;
                                // cout << graph->seq[y>>1].name << " kmers: " << y_discrimitive_kmers.size() << endl;
                                for(auto k_mer : y_discrimitive_kmers){
                                    node_discrimitive_kmers[y].insert(k_mer);
                                }
                            }
                        }
                    }
                }
            }
        }
    }



    // set<uint32_t> nodes = set<uint32_t>();
    // for(int b = 0; b < graph->n_seq * 2; b++){
    //     nodes.insert(b);
    // }
    // cout << nodes.size() << endl;
    // for(auto x: nodes){
    //     // cout << graph->seq[x>>1].name << " " << (x%2==0 ? "+" : "-") << endl;
    //     for(auto y: nodes){
    //         if(x!=y && x!=(y^1)){
    //             uint32_t num_outgoing_arcs_x = asg_arc_n(graph, x);
    //             asg_arc_t *outgoing_arcs_x = asg_arc_a(graph, x);
    //             uint32_t num_outgoing_arcs_y = asg_arc_n(graph, y);
    //             asg_arc_t *outgoing_arcs_y = asg_arc_a(graph, y);
    //             uint32_t num_incoming_arcs_x = asg_arc_n(graph, x^1);
    //             asg_arc_t *incoming_arcs_x = asg_arc_a(graph, x^1);
    //             uint32_t num_incoming_arcs_y = asg_arc_n(graph, y^1);
    //             asg_arc_t *incoming_arcs_y = asg_arc_a(graph, y^1);
    //             bool flag = false;
    //             for(int i=0; i<num_outgoing_arcs_x; i++){
    //                 for(int j=0; j<num_outgoing_arcs_y; j++){
    //                     if(outgoing_arcs_x[i].v == outgoing_arcs_y[j].v){
    //                         flag = true;
    //                         break;
    //                     }
    //                 }
    //                 if(flag){
    //                     break;
    //                 }
    //             }
    //             bool flag2 = false;
    //             for(int i=0; i<num_incoming_arcs_x; i++){
    //                 for(int j=0; j<num_incoming_arcs_y; j++){
    //                     if(incoming_arcs_x[i].v == incoming_arcs_y[j].v){
    //                         flag2 = true;
    //                         break;
    //                     }
    //                 }
    //                 if(flag2){
    //                     break;
    //                 }
    //             }

    //             if(flag2||flag){
    //                 string x_seq = string(graph->seq[x>>1].seq);
    //                 uint32_t x_len = graph->seq[x>>1].len;
    //                 if(x%2==1){
    //                     x_seq = complement(string(x_seq));
    //                 }
    //                 string y_seq = string(graph->seq[y>>1].seq);
    //                 uint32_t y_len = graph->seq[y>>1].len;
    //                 if(y%2==1){
    //                     y_seq = complement(string(y_seq));
    //                 }

    //                 // x_seq = string("CCCCGCCCCGAAGCGCCCGCGTTGCGGTGGTTCATTTCTTCGCTTGCCCACTGGGCCTGGCAGCCTTCCGGCCCGTGGTCGTGCCTTGGCAGTCCCGCAC");
    //                 // y_seq = string("CCCCGCCCCGATGCCCCCGCGTTGCGGTGGTTCATTTCTTCGCTTGCCCACTGGGCCTGGCAGCCTTCCGGCCCGTGGTCGTGCCTTGGCAGTCCCGCAC");


    //                 // const char** buf_seq = const_cast<const char**>(&(x_seq.c_str()));
    //                 // const char** buf_name = const_cast<const char**>(&graph->seq[node.first>>1].name);
    //                 // mm_idx_t *this_node = mm_idx_str(iopt.w, iopt.k, iopt.flag, iopt.bucket_bits, 1, buf_seq, buf_name);
    //                 // mm_tbuf_t *tbuf = mm_tbuf_init();
    //                 // mm_mapopt_update(&mopt, this_node);


    //                 vector<uint32_t> bufff = vector<uint32_t>();
    //                 bufff.push_back(max(x,y));
    //                 bufff.push_back(min(x,y));
    //                 bubbles.insert(bufff);
    //                 bubble_nodes.insert(x);
    //                 bubble_nodes.insert(y);


    //                 Modset *ms_x = modsetCreate (hasher, params.B, 0) ;
    //                 Reference *ref_x = referenceCreate (ms_x, 1 << 26) ;
    //                 int id_x ;
    //                 dictAdd (ref_x->dict, strdup("query_sequence"), &id_x);
    //                 array (ref_x->len, id_x, int) = x_len ;
    //                 SeqhashRCiterator *mi_x = modRCiterator (hasher, strdup(x_seq.c_str()), x_len) ;
    //                 U64 hash_x ; int pos_x ;
    //                 while (modRCnext (mi_x, &hash_x, &pos_x, 0))
    //                 { U32 index = modsetIndexFind (ref_x->ms, hash_x, TRUE) ;
    //                 if (index)
    //                     { if (ref_x->max+1 >= ref_x->size) die ("reference size overflow") ;
    //                     ref_x->index[ref_x->max] = index ;
    //                     ++ref_x->depth[index] ;
    //                     ref_x->offset[ref_x->max] = pos_x ;
    //                     ref_x->id[ref_x->max] = id_x ;
    //                     ++ref_x->max ;
    //                     }
    //                 }
    //                 seqhashRCiteratorDestroy (mi_x) ;
    //                 // int i_x ; U32 *d_x = &ref_x->depth[1] ;
    //                 // for (i_x = 1 ; i_x <= ref_x->ms->max ; ++i_x, ++d_x)
    //                 //     if (*d_x == 1) { msSetCopy1 (ref_x->ms, i_x) ; }
    //                 //     else if (*d_x == 2) { msSetCopy2 (ref_x->ms, i_x) ; }
    //                 //     else { msSetCopyM (ref_x->ms, i_x) ; }



    //                 Modset *ms_y = modsetCreate (hasher, params.B, 0) ;
    //                 Reference *ref_y = referenceCreate (ms_y, 1 << 26) ;
    //                 int id_y ;
    //                 dictAdd (ref_y->dict, strdup("query_sequence"), &id_y);
    //                 array (ref_y->len, id_y, int) = y_len ;
    //                 SeqhashRCiterator *mi_y = modRCiterator (hasher, strdup(y_seq.c_str()), y_len) ;
    //                 U64 hash_y ; int pos_y ;
    //                 while (modRCnext (mi_y, &hash_y, &pos_y, 0))
    //                 { U32 index = modsetIndexFind (ref_y->ms, hash_y, TRUE) ;
    //                 if (index)
    //                     { if (ref_y->max+1 >= ref_y->size) die ("reference size overflow") ;
    //                     ref_y->index[ref_y->max] = index ;
    //                     ++ref_y->depth[index] ;
    //                     ref_y->offset[ref_y->max] = pos_y ;
    //                     ref_y->id[ref_y->max] = id_y ;
    //                     ++ref_y->max ;
    //                     }
    //                 }
    //                 seqhashRCiteratorDestroy (mi_y) ;
    //                 // int i_y ; U32 *d_y = &ref_y->depth[1] ;
    //                 // for (i_y = 1 ; i_y <= ref_y->ms->max ; ++i_y, ++d_y)
    //                 //     if (*d_y == 1) { msSetCopy1 (ref_y->ms, i_y) ; }
    //                 //     else if (*d_y == 2) { msSetCopy2 (ref_y->ms, i_y) ; }
    //                 //     else { msSetCopyM (ref_y->ms, i_y) ; }


    //                 map<U64,uint32_t> y_discrimitive_kmers = map<U64,uint32_t>();
    //                 map<U64,uint32_t> x_discrimitive_kmers = map<U64,uint32_t>();

    //                 SeqhashRCiterator *mi_query_y = modRCiterator (hasher, strdup(y_seq.c_str()), y_len) ;
    //                 SeqhashRCiterator *mi_query_x = modRCiterator (hasher, strdup(x_seq.c_str()), x_len) ;

    //                 U64 hash_query_x ; int pos_query_x ;
    //                 while (modRCnext (mi_query_x, &hash_query_x, &pos_query_x, 0))
    //                 { U32 index = modsetIndexFind (ref_y->ms, hash_query_x, FALSE) ;
    //                 if (!index)
    //                     {
    //                         x_discrimitive_kmers[hash_query_x] = pos_query_x;
    //                     }
    //                 }
    //                 seqhashRCiteratorDestroy (mi_query_x) ;

    //                 U64 hash_query_y ; int pos_query_y ;
    //                 while (modRCnext (mi_query_y, &hash_query_y, &pos_query_y, 0))
    //                 { U32 index = modsetIndexFind (ref_x->ms, hash_query_y, FALSE) ;
    //                 if (!index)
    //                     {
    //                         y_discrimitive_kmers[hash_query_y] = pos_query_y;
    //                     }
    //                 }
    //                 seqhashRCiteratorDestroy (mi_query_y) ;


    //                 modsetDestroy(ms_x);
    //                 referenceDestroy(ref_x);
    //                 modsetDestroy(ms_y);
    //                 referenceDestroy(ref_y);
    //                 // set<string> y_discrimitive_kmers = set<string>();
    //                 // set<string> x_discrimitive_kmers = set<string>();

    //                 // set<string> x_kmers = set<string>();
    //                 // set<string> y_kmers = set<string>();

    //                 // for(int i = 0; i < x_seq.length()-31; i++){
    //                 //     x_kmers.insert(x_seq.substr(i,31));
    //                 // }
    //                 // for(int q = 0; q < y_seq.length()-31; q++){
    //                 //     y_kmers.insert(y_seq.substr(q,31));
    //                 // }
    //                 // for(auto k_mer: x_kmers){
    //                 //     if(y_kmers.find(k_mer)==y_kmers.end()){
    //                 //         x_discrimitive_kmers.insert(k_mer);

    //                 //     }
    //                 // }
    //                 // for(auto k_mer: y_kmers){
    //                 //     if(x_kmers.find(k_mer)==x_kmers.end()){
    //                 //         y_discrimitive_kmers.insert(k_mer);
    //                 //     }
    //                 // }

    //                 // cout << x_seq << endl;
    //                 // cout << y_seq << endl;
    //                 cout << graph->seq[x>>1].name << " kmers: " << x_discrimitive_kmers.size() << endl;
    //                 // for(auto k_mer : x_discrimitive_kmers){
    //                 //     cout << k_mer << endl;
    //                 // }
    //                 // cout << counter_valid << " against " << valid_kmers.size() << endl;
    //                 cout << graph->seq[y>>1].name << " kmers: " << y_discrimitive_kmers.size() << endl;
    //                 // for(auto k_mer : y_discrimitive_kmers){
    //                 //     cout << k_mer << endl;
    //                 // }
    //                 node_discrimitive_kmers[x] = x_discrimitive_kmers;
    //                 node_discrimitive_kmers[y] = y_discrimitive_kmers;
    //             }
    //         }
    //     }
    // }



    map<set<uint32_t>, set<string>> path_query_map = map<set<uint32_t>, set<string>>();
    map<set<uint32_t>, set<string>> path_query_support_map = map<set<uint32_t>, set<string>>();


    for(auto record: *ordered_records){
        set<set<uint32_t>> inserted_path;
        for(auto single_record: record){
            if(node_path_map.find(single_record.tn)!=node_path_map.end() && record.size()>1){
                for(auto buf_path:node_path_map[single_record.tn]){
                    if(inserted_path.insert(buf_path).second){
                        if(path_query_map.find(buf_path)==path_query_map.end()){
                            path_query_map[buf_path] = set<string>();
                        }
                        path_query_map[buf_path].insert(single_record.qn);
                    }
                }
            }
        }
    }

    for(auto bubble: bubbles_pure){
        vector<vector<paf_rec_str_t>> joint_set = vector<vector<paf_rec_str_t>>();
        set<string> inserted_set = set<string>();
        for(auto path: bubble->paths_nodes){
            for(auto buf_query: path_query_map[path]){
                if(inserted_set.insert(buf_query).second){
                    joint_set.push_back(name_rec[buf_query]);
                }
            }
        }


        for(auto buf_query: joint_set){
            vector<uint32_t> path_support_count(bubble->paths_nodes.size(),0);
            vector<float> path_support_ratio(bubble->paths_nodes.size(),0);
            for(auto buf_record: buf_query){
                bool found = false;
                for(auto path: bubble->paths_nodes){
                    found = path.find(buf_record.tn)!=path.end();
                    if(found){
                        break;
                    }
                }
                if(found){
                    float current_ratio = 0;
                    uint32_t current_support_count = 0;
                    uint32_t node_x = buf_record.tn;
                    string x_seq_str = string(graph->seq[node_x>>1].seq);
                    if(node_x%2==1){
                        x_seq_str = complement(x_seq_str);
                    }
                    string query_sequence = translate(x_seq_str,buf_record.ts,buf_record.seq);
                    Modset *ms_query = modsetCreate (hasher, params.B, 0) ;
                    Reference *ref_query = referenceCreate (ms_query, 1 << 26) ;
                    int id_query ;
                    dictAdd (ref_query->dict, strdup("query_sequence"), &id_query);
                    array (ref_query->len, id_query, int) = query_sequence.length() ;
                    SeqhashRCiterator *mi_query = modRCiterator (hasher, strdup(query_sequence.c_str()), query_sequence.length()) ;
                    U64 hash_query ; int pos_query ;
                    while (modRCnext (mi_query, &hash_query, &pos_query, 0))
                    { U32 index = modsetIndexFind (ref_query->ms, hash_query, TRUE) ;
                    if (index)
                        { if (ref_query->max+1 >= ref_query->size) die ("reference size overflow") ;
                        ref_query->index[ref_query->max] = index ;
                        ++ref_query->depth[index] ;
                        ref_query->offset[ref_query->max] = pos_query ;
                        ref_query->id[ref_query->max] = id_query ;
                        ++ref_query->max ;
                        }
                    }
                    seqhashRCiteratorDestroy (mi_query) ;

                    for(auto kmer_x: node_discrimitive_kmers[node_x]){
                        if(modsetIndexFind (ref_query->ms, kmer_x.first, FALSE)!=0){
                            current_support_count+=1;
                        }
                    }
                    if(node_discrimitive_kmers[node_x].size()!=0){
                        current_ratio = (float)current_support_count/node_discrimitive_kmers[node_x].size();
                    }
                    modsetDestroy(ms_query);
                    referenceDestroy(ref_query);

                    for(int idx = 0; idx < bubble->paths_nodes.size(); idx++){
                        if(bubble->paths_nodes[idx].find(node_x)!=bubble->paths_nodes[idx].end()){
                            path_support_count[idx]++;
                            path_support_ratio[idx]+=current_ratio;
                        }
                    }
                }
            }

            uint32_t max_idx = 0;
            float max_ratio = 0;
            for(int idx = 0; idx < bubble->paths_nodes.size(); idx++){
                path_support_ratio[idx]/=path_support_count[idx];
                if(max_ratio < path_support_ratio[idx]){
                    max_idx = idx;
                    max_ratio = path_support_ratio[idx];
                }
            }
            path_query_support_map[bubble->paths_nodes[max_idx]].insert(string(buf_query[0].qn));
        }
    }


    for(auto path: path_query_support_map){
        cout << endl;
        for(auto node: path.first){
            cout << graph->seq[node>>1].name << ", ";
        }
        cout << endl;
        cout << path.second.size() << endl;
        cout << endl;
    }
    // for(auto bubble: bubbles){
    //     uint32_t node_x = bubble[0];
    //     uint32_t node_y = bubble[1];
    //     cout << graph->seq[node_x>>1].name << " against " << graph->seq[node_y>>1].name << endl;
    //     if(node_query_set.find(node_x)!=node_query_set.end() && node_query_set.find(node_y)!=node_query_set.end()){
    //         vector<vector<paf_rec_str_t>> joint_set = node_query_set[node_x];
    //         for(auto buf_query: node_query_set[node_y]){
    //             joint_set.push_back(buf_query);
    //         }
    //         for(auto buf_query: joint_set){
    //             uint32_t node_x_support_record = 0;
    //             uint32_t node_y_support_record = 0;
    //             for(auto buf_record: buf_query){
    //                 string query_sequence = string();
    //                 if(buf_record.tn==node_x){
    //                     string x_seq_str = string(graph->seq[node_x>>1].seq);
    //                     if(node_x%2==1){
    //                         x_seq_str = complement(x_seq_str);
    //                     }
    //                     query_sequence = translate(x_seq_str,buf_record.ts,buf_record.seq);
    //                 }else if(buf_record.tn == node_y){
    //                     string y_seq_str = string(graph->seq[node_y>>1].seq);
    //                     if(node_y%2==1){
    //                         y_seq_str = complement(y_seq_str);
    //                     }
    //                     query_sequence = translate(y_seq_str,buf_record.ts,buf_record.seq);
    //                 }
    //                 if(query_sequence.length()>0){

    //                     // set<string> query_kmers = set<string>();

    //                     // for(int i = 0; i < query_sequence.length()-31; i++){
    //                     //     query_kmers.insert(query_sequence.substr(i,31));
    //                     // }

    //                     // // int i ; U32 *d = &ref_query->depth[1] ; U32 n1 = 0, n2 = 0, nM = 0 ;
    //                     // // for (i = 1 ; i <= ref_query->ms->max ; ++i, ++d)
    //                     // //     if (*d == 1) { msSetCopy1 (ref_query->ms, i) ; ++n1 ; }
    //                     // //     else if (*d == 2) { msSetCopy2 (ref_query->ms, i) ; ++n2 ; }
    //                     // //     else { msSetCopyM (ref_query->ms, i) ; ++nM ; }
    //                     // // modsetPack (ref_query->ms) ;
    //                     // // referencePack (ref_query) ;


    //                     // for(auto kmer_x: node_discrimitive_kmers[node_x]){
    //                     //     if(query_kmers.find(kmer_x)!=query_kmers.end()){
    //                     //         node_x_support_record+=1;
    //                     //     }
    //                     // }
    //                     // for(auto kmer_y: node_discrimitive_kmers[node_y]){
    //                     //     if(query_kmers.find(kmer_y)!=query_kmers.end()){
    //                     //         node_y_support_record+=1;
    //                     //     }
    //                     // }

    //                     Modset *ms_query = modsetCreate (hasher, params.B, 0) ;
    //                     Reference *ref_query = referenceCreate (ms_query, 1 << 26) ;
    //                     int id_query ;
    //                     dictAdd (ref_query->dict, strdup("query_sequence"), &id_query);
    //                     array (ref_query->len, id_query, int) = query_sequence.length() ;
    //                     SeqhashRCiterator *mi_query = modRCiterator (hasher, strdup(query_sequence.c_str()), query_sequence.length()) ;
    //                     U64 hash_query ; int pos_query ;
    //                     while (modRCnext (mi_query, &hash_query, &pos_query, 0))
    //                     { U32 index = modsetIndexFind (ref_query->ms, hash_query, TRUE) ;
    //                     if (index)
    //                         { if (ref_query->max+1 >= ref_query->size) die ("reference size overflow") ;
    //                         ref_query->index[ref_query->max] = index ;
    //                         ++ref_query->depth[index] ;
    //                         ref_query->offset[ref_query->max] = pos_query ;
    //                         ref_query->id[ref_query->max] = id_query ;
    //                         ++ref_query->max ;
    //                         }
    //                     }
    //                     seqhashRCiteratorDestroy (mi_query) ;
    //                     // int i_query ; U32 *d_query = &ref_query->depth[1] ;
    //                     // for (i_query = 1 ; i_query <= ref_query->ms->max ; ++i_query, ++d_query)
    //                     //     if (*d_query == 1) { msSetCopy1 (ref_query->ms, i_query) ; }
    //                     //     else if (*d_query == 2) { msSetCopy2 (ref_query->ms, i_query) ; }
    //                     //     else { msSetCopyM (ref_query->ms, i_query) ; }

    //                     for(auto kmer_x: node_discrimitive_kmers[node_x]){
    //                         if(modsetIndexFind (ref_query->ms, kmer_x.first, FALSE)!=0){
    //                             node_x_support_record+=1;
    //                         }
    //                     }
    //                     for(auto kmer_y: node_discrimitive_kmers[node_y]){
    //                         if(modsetIndexFind (ref_query->ms, kmer_y.first, FALSE)!=0){
    //                             node_y_support_record+=1;
    //                         }
    //                     }

    //                     modsetDestroy(ms_query);
    //                     referenceDestroy(ref_query);

    //                     float node_x_support_record_ratio = 0;
    //                     float node_y_support_record_ratio = 0;
    //                     if(node_discrimitive_kmers[node_x].size()!=0){
    //                         node_x_support_record_ratio = (float)node_x_support_record/node_discrimitive_kmers[node_x].size();
    //                     }
    //                     if(node_discrimitive_kmers[node_y].size()!=0){
    //                         node_y_support_record_ratio = (float)node_y_support_record/node_discrimitive_kmers[node_y].size();
    //                     }
    //                     // cout << node_x_support_record_ratio << " against " << node_y_support_record_ratio << endl;
    //                     // cout << graph->seq[node_x].name << " against " << graph->seq[node_y].name << endl;
    //                     if(node_x_support_record_ratio > node_y_support_record_ratio*1.5){
    //                         if(node_support_reads.find(node_x)==node_support_reads.end()){
    //                             node_support_reads[node_x] = set<string>();
    //                         }
    //                         node_support_reads[node_x].insert(string(buf_record.qn));
    //                     }else if(node_y_support_record_ratio > node_x_support_record_ratio*1.5){
    //                         if(node_support_reads.find(node_y)==node_support_reads.end()){
    //                             node_support_reads[node_y] = set<string>();
    //                         }
    //                         node_support_reads[node_y].insert(string(buf_record.qn));
    //                     }
    //                     break;
    //                 }
    //             }
    //         }
    //     }
    // }

    vector<bubble_t*> pure_bubbles = vector<bubble_t*>();
    set<vector<uint32_t>> added_begin_end = set<vector<uint32_t>>();
    map<set<uint32_t>, set<string>> pure_path_query_support_map = map<set<uint32_t>, set<string>>();
    map<vector<uint32_t>,vector<set<uint32_t>>> begEnd_pathes = map<vector<uint32_t>,vector<set<uint32_t>>>();
    map<set<uint32_t>,vector<uint32_t>> path_begEnd = map<set<uint32_t>,vector<uint32_t>>();
    map<set<uint32_t>, map<set<uint32_t>, set<string>>> path_connection = map<set<uint32_t>, map<set<uint32_t>, set<string>>>();
    map<set<uint32_t>, map<set<uint32_t>, uint32_t>> path_connection_count = map<set<uint32_t>, map<set<uint32_t>, uint32_t>>();
    map<uint32_t, set<vector<set<uint32_t>>>> ordered_path_connection_count = map<uint32_t, set<vector<set<uint32_t>>>>();
    uint32_t max_connection_support = 0;
    map<vector<uint32_t>, map<uint32_t,set<uint32_t>>> begEnd_haplo_path = map<vector<uint32_t>, map<uint32_t,set<uint32_t>>>();

    for(auto bubble: bubbles_pure){
        vector<uint32_t> buf = vector<uint32_t>();
        buf.push_back(bubble->begNode>>1);
        buf.push_back(bubble->endNode>>1);
        if(added_begin_end.insert(buf).second){
            begEnd_haplo_path[buf] = map<uint32_t,set<uint32_t>>();
            begEnd_haplo_path[buf][1] = set<uint32_t>();
            begEnd_haplo_path[buf][2] = set<uint32_t>();
            vector<set<uint32_t>> all_pure_pathes = vector<set<uint32_t>>();
            for(auto path: bubble->paths_nodes){
                set<uint32_t> pure_path = set<uint32_t>();
                for(auto buf_node: path){
                    pure_path.insert(buf_node>>1);
                }
                all_pure_pathes.push_back(pure_path);
                path_begEnd[pure_path] = buf;
            }
            begEnd_pathes[buf] = all_pure_pathes;
            pure_bubbles.push_back(bubble);
        }else{
            delete bubble;
        }
    }

    for(auto path: path_query_support_map){
        set<uint32_t> pure_path = set<uint32_t>();
        for(auto buf_node: path.first){
            pure_path.insert(buf_node>>1);
        }
        if(pure_path_query_support_map.find(pure_path)==pure_path_query_support_map.end()){
            pure_path_query_support_map[pure_path] = set<string>();
        }
        for(auto read_name: path.second){
            pure_path_query_support_map[pure_path].insert(read_name);
        }
    }

    for(auto path_a: pure_path_query_support_map){
        path_connection[path_a.first] = map<set<uint32_t>, set<string>>();
        for(auto path_b: pure_path_query_support_map){
            path_connection[path_a.first][path_b.first] = set<string>();
            for(auto read_b: path_b.second){
                if(path_a.second.find(read_b)!=path_a.second.end()){
                    path_connection[path_a.first][path_b.first].insert(read_b);
                }
            }
        }
    }

    for(auto path_a: pure_path_query_support_map){
        cout << endl;
        for(auto node: path_a.first){
            cout << graph->seq[node].name << ", ";
        }
        cout << endl;
        for(auto name: path_a.second){
            vector<paf_rec_str_t> records = name_rec[name];
            for(auto record: records){
                if(path_a.first.find(record.tn>>1)!=path_a.first.end()){
                    cout << graph->seq[record.tn >>1].name << ", ";
                }
            }
            cout << endl;
        }
        cout << endl;
    }

    for(auto path_a:path_connection){
        path_connection_count[path_a.first] = map<set<uint32_t>, uint32_t>();
        for(auto path_b:path_a.second){
            if(path_b.second.size()>max_connection_support){
                max_connection_support = path_b.second.size();
            }
            if(ordered_path_connection_count.find(path_b.second.size())==ordered_path_connection_count.end()){
                ordered_path_connection_count[path_b.second.size()] = set<vector<set<uint32_t>>>();
            }
            vector<set<uint32_t>> buf_begEnd = vector<set<uint32_t>>();
            buf_begEnd.push_back(path_a.first);
            buf_begEnd.push_back(path_b.first);
            ordered_path_connection_count[path_b.second.size()].insert(buf_begEnd);
            path_connection_count[path_a.first][path_b.first] = path_b.second.size();
        }
    }

    set<set<uint32_t>> haplo_2 = set<set<uint32_t>>();
    set<set<uint32_t>> haplo_1 = set<set<uint32_t>>();
    bool first = true;
    bool finished = false;
    bool changed = true;
    uint32_t round = 1000;
    while(!finished && round-- && changed){
        changed = false;
        for(int i = max_connection_support; i>0; i--){
            if(ordered_path_connection_count.find(i)!=ordered_path_connection_count.end()){
                for(auto connection: ordered_path_connection_count[i]){
                    set<uint32_t> path_x = connection[0];
                    set<uint32_t> path_y = connection[1];
                    if(first){
                        haplo_1.insert(path_x);
                        haplo_1.insert(path_y);
                        begEnd_haplo_path[path_begEnd[path_x]][1] = path_x;
                        begEnd_haplo_path[path_begEnd[path_y]][1] = path_y;
                        first = false;
                        changed = true;
                        continue;
                    }

                    uint32_t seen_x = 0;
                    uint32_t seen_y = 0;

                    if(haplo_2.find(path_x)!=haplo_2.end()){
                        seen_x = 2;
                    }else if(haplo_1.find(path_x)!=haplo_1.end()){
                        seen_x = 1;
                    }

                    if(haplo_2.find(path_y)!=haplo_2.end()){
                        seen_y = 2;
                    }else if(haplo_1.find(path_y)!=haplo_1.end()){
                        seen_y = 1;
                    }

                    if(seen_x!=0 && seen_y==0){
                        if(begEnd_haplo_path[path_begEnd[path_y]][seen_x].size()==0){
                            changed = true;
                            begEnd_haplo_path[path_begEnd[path_y]][seen_x] = path_y;
                            if(seen_x == 1){
                                haplo_1.insert(path_y);
                            }else{
                                haplo_2.insert(path_y);
                            }
                        }
                    }else if(seen_y!=0 && seen_x==0){
                        if(begEnd_haplo_path[path_begEnd[path_x]][seen_y].size()==0){
                            changed = true;
                            begEnd_haplo_path[path_begEnd[path_x]][seen_y] = path_x;
                            if(seen_y == 1){
                                haplo_1.insert(path_x);
                            }else{
                                haplo_2.insert(path_x);
                            }
                        }
                    }else if(seen_x==0 && seen_y==0){
                        if(begEnd_pathes[path_begEnd[path_x]].size()==2){
                            if(begEnd_haplo_path[path_begEnd[path_x]][1].size()!=0 && begEnd_haplo_path[path_begEnd[path_x]][2].size()!=0){
                                seen_x = -1;
                            }else if(begEnd_haplo_path[path_begEnd[path_x]][1].size()!=0){
                                seen_x = 2;
                            }else if(begEnd_haplo_path[path_begEnd[path_x]][2].size()!=0){
                                seen_x = 1;
                            }
                        }
                        if(begEnd_pathes[path_begEnd[path_y]].size()==2){
                            if(begEnd_haplo_path[path_begEnd[path_y]][1].size()!=0 && begEnd_haplo_path[path_begEnd[path_y]][2].size()!=0){
                                seen_x = -1;
                            }else if(begEnd_haplo_path[path_begEnd[path_y]][1].size()!=0){
                                seen_y = 2;
                            }else if(begEnd_haplo_path[path_begEnd[path_y]][2].size()!=0){
                                seen_y = 1;
                            }
                        }
                        if(seen_x!=-1 && seen_y!=-1){
                            if((seen_y == seen_x || seen_y==0) && seen_x != 0){
                                changed = true;
                                begEnd_haplo_path[path_begEnd[path_y]][seen_x] = path_y;
                                begEnd_haplo_path[path_begEnd[path_x]][seen_x] = path_x;
                                if(seen_x == 1){
                                    haplo_1.insert(path_x);
                                    haplo_1.insert(path_y);
                                }else{
                                    haplo_2.insert(path_x);
                                    haplo_2.insert(path_y);
                                }
                            }else if(seen_x == 0 && seen_y != 0){
                                changed = true;
                                begEnd_haplo_path[path_begEnd[path_y]][seen_y] = path_y;
                                begEnd_haplo_path[path_begEnd[path_x]][seen_y] = path_x;
                                if(seen_y == 1){
                                    haplo_1.insert(path_x);
                                    haplo_1.insert(path_y);
                                }else{
                                    haplo_2.insert(path_x);
                                    haplo_2.insert(path_y);
                                }
                            }
                        }
                    }
                }
            }
        }
        finished = true;
        for(auto begEnd:begEnd_haplo_path){
            if(begEnd.second[1].size()==0 || begEnd.second[2].size()==0){
                finished = false;
                break;
            }
        }
    }

    cout << endl;
    for(auto path: haplo_1){
        for(auto node: path){
            cout << graph->seq[node].name << ", ";
        }
        cout << "\t";
    }
    cout << endl;

    for(auto path: haplo_2){
        for(auto node: path){
            cout << graph->seq[node].name << ", ";
        }
        cout << "\t";
    }
    cout << endl;
    // map<uint32_t,map<uint32_t,uint32_t>> node_connection = map<uint32_t,map<uint32_t,uint32_t>>();
    // for(auto bubble_a: pure_bubbles){
    //     if(pure_node_support_reads[bubble_a[0]>>1].size()==0 && pure_node_support_reads[bubble_a[1]>>1].size()==0){
    //         continue;
    //     }
    //     bubble_connection[bubble_a] = map<vector<uint32_t>, vector<set<string>>>();
    //     for(auto bubble_b: pure_bubbles){
    //         if(bubble_a[0]!=bubble_b[0] || bubble_a[1]!=bubble_b[1]){
    //             if(pure_node_support_reads[bubble_b[0]>>1].size()==0 && pure_node_support_reads[bubble_b[1]>>1].size()==0){
    //                 continue;
    //             }
    //             bubble_connection[bubble_a][bubble_b] = vector<set<string>>();
    //             bubble_connection[bubble_a][bubble_b].push_back(set<string>());
    //             bubble_connection[bubble_a][bubble_b].push_back(set<string>());
    //             bubble_connection[bubble_a][bubble_b].push_back(set<string>());
    //             bubble_connection[bubble_a][bubble_b].push_back(set<string>());
    //             set<string> set_record_node_b_0 = pure_node_support_reads[bubble_b[0]>>1];
    //             set<string> set_record_node_b_1 = pure_node_support_reads[bubble_b[1]>>1];
    //             uint32_t counter_0_0 = 0;
    //             uint32_t counter_0_1 = 0;
    //             uint32_t counter_1_0 = 0;
    //             uint32_t counter_1_1 = 0;
    //             for(auto record_node_a_0: pure_node_support_reads[bubble_a[0]>>1]){
    //                 if(set_record_node_b_0.find(record_node_a_0)!=set_record_node_b_0.end()){
    //                     bubble_connection[bubble_a][bubble_b][0].insert(record_node_a_0);
    //                     counter_0_0++;
    //                 }
    //                 if(set_record_node_b_1.find(record_node_a_0)!=set_record_node_b_1.end()){
    //                     bubble_connection[bubble_a][bubble_b][1].insert(record_node_a_0);
    //                     counter_0_1++;
    //                 }
    //             }
    //             for(auto record_node_a_1: pure_node_support_reads[bubble_a[1]>>1]){
    //                 if(set_record_node_b_0.find(record_node_a_1)!=set_record_node_b_0.end()){
    //                     bubble_connection[bubble_a][bubble_b][2].insert(record_node_a_1);
    //                     counter_1_0++;
    //                 }
    //                 if(set_record_node_b_1.find(record_node_a_1)!=set_record_node_b_1.end()){
    //                     bubble_connection[bubble_a][bubble_b][3].insert(record_node_a_1);
    //                     counter_1_1++;
    //                 }
    //             }
    //             if(counter_0_0 + counter_1_1 > (counter_1_0 + counter_0_1)*1.5){
    //                 if(node_connection.find(bubble_a[0]>>1)==node_connection.end()){
    //                     node_connection[bubble_a[0]>>1] = map<uint32_t,uint32_t>();
    //                 }
    //                 node_connection[bubble_a[0]>>1].insert({bubble_b[0]>>1,counter_0_0 + counter_1_1});
    //                 if(node_connection.find(bubble_a[1]>>1)==node_connection.end()){
    //                     node_connection[bubble_a[1]>>1] = map<uint32_t,uint32_t>();
    //                 }
    //                 node_connection[bubble_a[1]>>1].insert({bubble_b[1]>>1,counter_0_0 + counter_1_1});
    //             }else if(counter_1_0 + counter_0_1 > (counter_0_0 + counter_1_1)*1.5){
    //                 if(node_connection.find(bubble_a[1]>>1)==node_connection.end()){
    //                     node_connection[bubble_a[1]>>1] = map<uint32_t,uint32_t>();
    //                 }
    //                 node_connection[bubble_a[1]>>1].insert({bubble_b[0]>>1,counter_1_0 + counter_0_1});
    //                 if(node_connection.find(bubble_a[0]>>1)==node_connection.end()){
    //                     node_connection[bubble_a[0]>>1] = map<uint32_t,uint32_t>();
    //                 }
    //                 node_connection[bubble_a[0]>>1].insert({bubble_b[1]>>1,counter_1_0 + counter_0_1});
    //             }else{
    //                 if(counter_0_0 > 2){
    //                     if(node_connection.find(bubble_a[0]>>1)==node_connection.end()){
    //                         node_connection[bubble_a[0]>>1] = map<uint32_t,uint32_t>();
    //                     }
    //                     node_connection[bubble_a[0]>>1].insert({bubble_b[0]>>1,counter_0_0});
    //                 }
    //                 if(counter_0_1 > 2){
    //                     if(node_connection.find(bubble_a[0]>>1)==node_connection.end()){
    //                         node_connection[bubble_a[0]>>1] = map<uint32_t,uint32_t>();
    //                     }
    //                     node_connection[bubble_a[0]>>1].insert({bubble_b[1]>>1,counter_0_1});
    //                 }
    //                 if(counter_1_0 > 2){
    //                     if(node_connection.find(bubble_a[1]>>1)==node_connection.end()){
    //                         node_connection[bubble_a[1]>>1] = map<uint32_t,uint32_t>();
    //                     }
    //                     node_connection[bubble_a[1]>>1].insert({bubble_b[0]>>1,counter_1_0});
    //                 }
    //                 if(counter_1_1 > 2){
    //                     if(node_connection.find(bubble_a[1]>>1)==node_connection.end()){
    //                         node_connection[bubble_a[1]>>1] = map<uint32_t,uint32_t>();
    //                     }
    //                     node_connection[bubble_a[1]>>1].insert({bubble_b[1]>>1,counter_1_1});
    //                 }
    //             }
    //         }
    //     }
    // }


    // uint32_t max_support = 0;
    // map<uint32_t,vector<vector<uint32_t>>> support_count_connections = map<uint32_t,vector<vector<uint32_t>>>();
    // uint32_t path_num = 0;
    // map<uint32_t, map<uint32_t,uint32_t>> node_path_map = map<uint32_t, map<uint32_t,uint32_t>>();
    // vector<set<uint32_t>> path_set = vector<set<uint32_t>>();
    // path_set.push_back(set<uint32_t>());
    // path_set.push_back(set<uint32_t>());

    // for(auto node_1: node_connection){
    //     for(auto node_2: node_1.second){
    //         if(node_1.first<node_2.first){
    //             if(node_2.second>max_support){
    //                 max_support = node_2.second;
    //             }
    //             if(support_count_connections.find(node_2.second)==support_count_connections.end()){
    //                 support_count_connections[node_2.second] = vector<vector<uint32_t>>();
    //             }
    //             vector<uint32_t> buff;
    //             buff.push_back(node_1.first);
    //             buff.push_back(node_2.first);
    //             support_count_connections[node_2.second].push_back(buff);
    //         }
    //     }
    // }
    // // for(int i = max_support; i>0; i--){
    // //     if(support_count_connections.find(i)!=support_count_connections.end()){
    // //         for(auto connection: support_count_connections[i]){
    // //             uint32_t node_1 = connection[0];
    // //             uint32_t node_2 = connection[1];
    // //             bool node_1_seen = node_path_map.find(node_1)!=node_path_map.end();
    // //             bool node_2_seen = node_path_map.find(node_2)!=node_path_map.end();
    // //             if(!node_1_seen && !node_2_seen){
    // //                 node_path_map[node_1] = set<uint32_t>();
    // //                 node_path_map[node_2] = set<uint32_t>();
    // //                 node_path_map[node_1].insert(path_set.size());
    // //                 node_path_map[node_2].insert(path_set.size());
    // //                 path_set.push_back(set<uint32_t>());
    // //                 path_set[path_set.size()-1].insert(node_1);
    // //                 path_set[path_set.size()-1].insert(node_2);
    // //             }else if(node_1_seen && !node_2_seen){
    // //                 node_path_map[node_2] = set<uint32_t>();
    // //                 for(auto x: node_path_map[node_1]){
    // //                     node_path_map[node_2].insert(x);
    // //                     path_set[x].insert(node_2);
    // //                 }
    // //             }else if(node_2_seen && !node_1_seen){
    // //                 node_path_map[node_1] = set<uint32_t>();
    // //                 for(auto x: node_path_map[node_2]){
    // //                     node_path_map[node_1].insert(x);
    // //                     path_set[x].insert(node_1);
    // //                 }
    // //             }else{
    // //                 for(auto x: node_path_map[node_1]){
    // //                     node_path_map[node_2].insert(x);
    // //                     path_set[x].insert(node_2);
    // //                 }
    // //                 for(auto x: node_path_map[node_2]){
    // //                     node_path_map[node_1].insert(x);
    // //                     path_set[x].insert(node_1);
    // //                 }
    // //             }
    // //         }
    // //     }
    // // }
    // set<vector<uint32_t>> counted_connections = set<vector<uint32_t>>();
    // for(auto connection: support_count_connections[max_support]){
    //     counted_connections.insert(connection);
    //     vector<uint32_t> other_connection = vector<uint32_t>();
    //     other_connection.push_back(connection[1]);
    //     other_connection.push_back(connection[0]);
    //     counted_connections.insert(other_connection);
    //     uint32_t node_1 = connection[0];
    //     uint32_t node_2 = connection[1];
    //     bool node_1_seen = node_path_map.find(node_1)!=node_path_map.end();
    //     bool node_2_seen = node_path_map.find(node_2)!=node_path_map.end();
    //     if(!node_1_seen&&!node_2_seen){
    //         node_path_map[node_1] = map<uint32_t,uint32_t>();
    //         node_path_map[node_1][0] = 0;
    //         node_path_map[node_1][1] = 0;
    //         node_path_map[node_1][0] += max_support;
    //         node_path_map[node_2] = map<uint32_t,uint32_t>();
    //         node_path_map[node_2][0] = 0;
    //         node_path_map[node_2][1] = 0;
    //         node_path_map[node_2][0] += max_support;
    //     }
    //     set<uint32_t> grow_set = set<uint32_t>();
    //     for(auto next_node: node_connection[node_1]){
    //         if(node_path_map.find(next_node.first)==node_path_map.end()){
    //             node_path_map[next_node.first] = map<uint32_t,uint32_t>();
    //             node_path_map[next_node.first][0] = 0;
    //             node_path_map[next_node.first][1] = 0;
    //             grow_set.insert(next_node.first);
    //             vector<uint32_t> connection = vector<uint32_t>();
    //             connection.push_back(node_1);
    //             connection.push_back(next_node.first);
    //             vector<uint32_t> other_connection = vector<uint32_t>();
    //             other_connection.push_back(next_node.first);
    //             other_connection.push_back(node_1);
    //             counted_connections.insert(connection);
    //             counted_connections.insert(other_connection);
    //         }
    //         node_path_map[next_node.first][0] += next_node.second;
    //     }
    //     for(auto next_node: node_connection[node_2]){
    //         if(node_path_map.find(next_node.first)==node_path_map.end()){
    //             node_path_map[next_node.first] = map<uint32_t,uint32_t>();
    //             node_path_map[next_node.first][0] = 0;
    //             node_path_map[next_node.first][1] = 0;
    //             grow_set.insert(next_node.first);
    //             vector<uint32_t> connection = vector<uint32_t>();
    //             connection.push_back(node_2);
    //             connection.push_back(next_node.first);
    //             vector<uint32_t> other_connection = vector<uint32_t>();
    //             other_connection.push_back(next_node.first);
    //             other_connection.push_back(node_2);
    //             counted_connections.insert(connection);
    //             counted_connections.insert(other_connection);
    //         }
    //         node_path_map[next_node.first][0] += next_node.second;
    //     }

    //     while(grow_set.size()>0){
    //         set<uint32_t> new_grow_set = set<uint32_t>();
    //         for(auto node: grow_set){
    //             for(auto next_node: node_connection[node]){
    //                 vector<uint32_t> connection = vector<uint32_t>();
    //                 connection.push_back(node);
    //                 connection.push_back(next_node.first);
    //                 vector<uint32_t> other_connection = vector<uint32_t>();
    //                 other_connection.push_back(next_node.first);
    //                 other_connection.push_back(node);
    //                 bool seen = (!counted_connections.insert(connection).second) || (!counted_connections.insert(other_connection).second);
    //                 if(!seen){
    //                     if(node_path_map.find(next_node.first)==node_path_map.end()){
    //                         node_path_map[next_node.first] = map<uint32_t,uint32_t>();
    //                         node_path_map[next_node.first][0] = 0;
    //                         node_path_map[next_node.first][1] = 0;
    //                     }
    //                     if(node_path_map[node][0]>node_path_map[node][1]){
    //                         node_path_map[next_node.first][0] += next_node.second;
    //                     }else if(node_path_map[node][0]<node_path_map[node][1]){
    //                         node_path_map[next_node.first][1] += next_node.second;
    //                     }
    //                     new_grow_set.insert(next_node.first);
    //                 }
    //             }
    //         }
    //         grow_set = new_grow_set;
    //     }
    //     // break;
    // }

    // // set<uint32_t> path_0;
    // // set<uint32_t> path_1;
    // set<uint32_t> path_0;
    // set<uint32_t> path_1;


    // for(auto bubble:bubbles){
    //     uint32_t node_x = bubble[0]>>1;
    //     uint32_t node_y = bubble[1]>>1;
    //     uint32_t path_0_node_x_support = 0;
    //     uint32_t path_0_node_y_support = 0;
    //     if(node_path_map.find(node_x)!=node_path_map.end()){
    //         path_0_node_x_support = node_path_map[node_x][0];
    //     }
    //     if(node_path_map.find(node_y)!=node_path_map.end()){
    //         path_0_node_y_support = node_path_map[node_y][0];
    //     }
    //     if(path_0_node_x_support>path_0_node_y_support){
    //         path_0.insert(node_x);
    //         path_1.insert(node_y);
    //     }else if(path_0_node_y_support!=0){
    //         path_0.insert(node_y);
    //         path_1.insert(node_x);
    //     }

    //     cout << endl;
    //     cout << graph->seq[node_x].name << "\t" << graph->seq[node_y].name << endl;;
    //     cout << path_0_node_x_support << "\t" << path_0_node_y_support << endl;

    // }

    // for(auto node : path_0){
    //     cout << graph->seq[node].name << ",\t";
    // }
    // cout << endl;

    // for(auto node : path_1){
    //     cout << graph->seq[node].name << ",\t";
    // }
    // cout << endl;

    for(auto bubble: pure_bubbles){
        delete bubble;
    }
    // for(auto node_path :node_path_map){
    //     cout << graph->seq[node_path.first].name <<": ";
    //     for(auto path: node_path.second){
    //         if(path==1){
    //             path_1.insert(node_path.first);

    //         }else if(path==0||path==6||path==7){
    //             path_0.insert(node_path.first);
    //         }
    //         cout << path << "\t";
    //     }
    //     cout << endl;
    // }

    // for(auto node : path_0){
    //     cout << graph->seq[node].name << ";\t";
    // }
    // cout << endl;

    // for(auto node : path_1){
    //     cout << graph->seq[node].name << ";\t";
    // }
    // cout << endl;
}

void repeat_resolver::save_pathes_to_file(set<vector<uint32_t>>* set_of_pathes, map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_graph, asg_t* graph, string output_directory){
    stringstream ss;
    ss << "repeat_resolve_graph.gfa";
    ofstream current_file(output_directory+string("/")+ss.str());
    int counting = 0;
    // set<uint32_t> common_nodes;
    for(auto a : *set_of_pathes){
        // uint32_t num_of_nodes = a[a.size()-1];
        for(auto b: a){
            cout << graph->seq[b>>1].name << " to ";
        }
        set<uint32_t> nodes;
        for(int b = 0; b < a.size()/2; b++){
            if(nodes.insert(a[b*2]>>1).second){
                current_file << "S\t" << string(graph->seq[a[b*2]>>1].name) + string("_") + to_string(counting) << "\t" << "*" << endl;
            }
            if(nodes.insert(a[b*2+1]>>1).second){
                current_file << "S\t" << string(graph->seq[a[b*2+1]>>1].name + string("_")) + to_string(counting) << "\t" << "*" << endl;
            }
            for(auto q : (*bubble_chain_graph)[a[b*2]][a[b*2+1]]){
                if(nodes.insert(q>>1).second){
                    current_file << "S\t" << string(graph->seq[q>>1].name) + string("_") + to_string(counting) << "\t" << "*" << endl;
                }
            }
        }
        for(auto q : nodes){
            uint32_t num_outgoing_arcs = asg_arc_n(graph, (q<<1));
            asg_arc_t *outgoing_arcs = asg_arc_a(graph, (q<<1));
            uint32_t num_outgoing_reversed_arcs = asg_arc_n(graph, (q<<1)^1);
            asg_arc_t *outgoing_reversed_arcs = asg_arc_a(graph, (q<<1)^1);
            for(int i = 0; i < num_outgoing_arcs; i++){
                uint32_t other_node = outgoing_arcs[i].v>>1;
                if(nodes.find(other_node)!=nodes.end()){
                    current_file << "L\t" << string(graph->seq[q].name) + string("_") + to_string(counting) << "\t"<< "+" << "\t" << string(graph->seq[outgoing_arcs[i].v>>1].name) + string("_") + to_string(counting) << "\t" << (outgoing_arcs[i].v%2==0 ? "+" : "-") << "\t0M\t" <<endl;
                }
            }
            for(int i = 0; i < num_outgoing_reversed_arcs; i++){
                uint32_t other_node = outgoing_reversed_arcs[i].v>>1;
                if(nodes.find(other_node)!=nodes.end()){
                    current_file << "L\t" << string(graph->seq[q].name) + string("_") + to_string(counting) << "\t"<< "-" << "\t" << string(graph->seq[outgoing_reversed_arcs[i].v>>1].name) + string("_") + to_string(counting) << "\t" << (outgoing_reversed_arcs[i].v%2==0 ? "+" : "-") << "\t0M\t" <<endl;
                }
            }
        }
        cout << endl;
        counting += 1;
    }
    current_file.close();
}
// map<string, vector<vector<paf_rec_str_t>*>*> repeat_resolver::get_queries_with_duplicated_matches(map<string, vector<paf_rec_str_t>*> unordered_records_by_query_name){
//     map<string, vector<vector<paf_rec_str_t>*>*> queires_with_duplicated_node_matches_by_node_name;
// 	for(auto i : unordered_records_by_query_name){
//         vector<paf_rec_str_t>* ordered_record = i.second;
//         set<string> matched_nodes;
//         set<string> duplicated_nodes;
//         bool flag = false;
//         for(auto q: *ordered_record){
//                 flag = true;
//                 if(!matched_nodes.insert(q.tn).second && duplicated_nodes.insert(q.tn).second){
//                     if(queires_with_duplicated_node_matches_by_node_name.find(q.tn) == queires_with_duplicated_node_matches_by_node_name.end()){
//                         queires_with_duplicated_node_matches_by_node_name[q.tn] = new vector<vector<paf_rec_str_t>*>;
//                     }
//                     queires_with_duplicated_node_matches_by_node_name[q.tn]->push_back(ordered_record);

//                 }
//         }
//         sort(ordered_record->begin(), ordered_record->end(), [ ]( const auto& lhs, const auto& rhs )
//         {
//         return lhs.qs < rhs.qs;
//         });
//     }
//     return queires_with_duplicated_node_matches_by_node_name;
// }

// vector<vector<uint32_t>> repeat_resolver::detect_cycles_in_graph(vector<uint32_t> duplicated_nodes, asg_t* graph){
//     vector<vector<uint32_t>> result;
//     for(auto i: duplicated_nodes){
//         vector<uint32_t> node_vector;  // DFS
//         vector<uint32_t> node_vi_vector;
//         node_vector.push_back(i);
//         node_vi_vector.push_back(0);
//         while(!node_vector.empty()) {
//             uint32_t u = node_vector.back();
//             uint32_t vi = node_vi_vector.back();
//             if ( node_vector.size() > 1 && (u==i || u == i^1)) {
//                 result.push_back(node_vector);
//                 cout << "Cycle detected: ";
//                 for (int ui=0; ui<node_vector.size(); ui++) {
//                     cout << graph->seq[node_vector[ui]/2].name << " ";
//                 }
//                 cout << endl;
//                 break;
//             }
//             uint32_t num_outgoing_arcs = asg_arc_n(graph, u);
//             if (vi < num_outgoing_arcs) {
//                 asg_arc_t *outgoing_arcs = asg_arc_a(graph, u);
//                 uint32_t v = outgoing_arcs[vi].v;
//                 node_vi_vector.back()++;
//                 node_vector.push_back(v);
//                 node_vi_vector.push_back(0);
//             } else {
//                 node_vector.pop_back();
//                 node_vi_vector.pop_back();
//             }
//         }
//     }
//     return result;
// }


// void repeat_resolver::solve_repeat_hic(char* hic_filename_1, char* hic_filename_2, asg_t* graph)
// {

//     mm_idxopt_t iopt;
//     mm_mapopt_t mopt;
//     mopt.min_cnt = 100;
//     mm_verbose = 2; // disable message output to stderr
//     mm_set_opt(0, &iopt, &mopt);
//     mopt.flag |= MM_F_CIGAR; // perform alignment
//     mopt.flag |= MM_F_OUT_CS; // perform alignment

//     map<uint32_t, mm_idx_t*> kmers;

//     for(uint32_t j = 0; j < graph->n_seq*2; j++){
//         string query = string(graph->seq[j/2].seq);
//         if(j%2==1){
//             query = complement(query);
//         }
//         const char* buf = query.c_str();
//         const char** buf_seq = const_cast<const char**>(&buf);
//         const char** buf_name = const_cast<const char**>(&graph->seq[j>>1].name);
//         mm_idx_t *current_kmers = mm_idx_str(iopt.w, iopt.k, iopt.flag, iopt.bucket_bits, 1, buf_seq, buf_name);

//         kmers[j] = current_kmers;
//     }



//     hic_reader reader;
//     hic_file_t* file1 = reader.hic_open(hic_filename_1);
//     hic_file_t* file2 = reader.hic_open(hic_filename_2);
//     hic_rec_t r;

//     map<vector<uint32_t>, uint32_t> connection_count;
//     while (reader.hic_read(file1, file2, &r) >= 0) {
//         Reference* ref;
//         map<uint32_t, float> count_1;
//         map<uint32_t, float> count_2;
//         for(auto unitig: kmers){


//             mm_tbuf_t *tbuf = mm_tbuf_init();
//             mm_mapopt_update(&mopt, unitig.second);
//             int j, i, n_reg;
//             mm_reg1_t *reg = mm_map(unitig.second, r.len_1, r.seq_1, &n_reg, tbuf, &mopt, "query_name_1");
//             if(n_reg>0){
//                 assert(reg->p);
//                 count_1[unitig.first] = reg->mapq;
//                 // cout << reg->mapq << endl;
//                     printf("%s\t%d\t%d\t%d\t%c\t", "query_name_1", r.len_1, reg->qs, reg->qe, "+-"[reg->rev]);
//                     printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcs:Z:", graph->seq[unitig.first/2].name, graph->seq[unitig.first/2].len, reg->rs, reg->re, reg->mlen, reg->blen, reg->mapq);
//                     char *buf = NULL;
//                     void *km = NULL;
//                     int max_len = 0;
//                     int result = mm_gen_cs(km,&buf,&max_len,const_cast<const mm_idx_t*>(unitig.second),const_cast<const mm_reg1_t*>(reg),r.seq_1,1);
//                     stringstream res = stringstream();
//                     for(int i = 0 ; i < result; i++){
//                         res << *(buf+i);
//                     }
//                     cout << res.str() << endl;
//             }else{
//                 count_1[unitig.first] = 0;
//             }

//             tbuf = mm_tbuf_init();
//             mm_mapopt_update(&mopt, unitig.second);
//             j=0; i=0; n_reg=0;
//             reg = mm_map(unitig.second, r.len_2, r.seq_2, &n_reg, tbuf, &mopt, "query_name_2");
//             if(n_reg>0){
//                 assert(reg->p);
//                 count_2[unitig.first] = reg->mapq;
//                 // cout << reg->mapq << endl;
//                     // printf("%s\t%d\t%d\t%d\t%c\t", "query_name_2", r.len_2, reg->qs, reg->qe, "+-"[reg->rev]);
//                     // printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcs:Z:", graph->seq[unitig.first/2].name, graph->seq[unitig.first/2].len, reg->rs, reg->re, reg->mlen, reg->blen, reg->mapq);
//             }else{
//                 count_2[unitig.first] = 0;
//             }

//         }
//         uint32_t matched_1 = 0;
//         uint32_t matched_2 = 0;
//         float max_count_1 = 0;
//         float max_count_2 = 0;
//         for(auto u: count_1){
//             if(max_count_1 < u.second){
//                 max_count_1 = u.second;
//                 matched_1 = u.first;
//             }
//         }
//         for(auto u: count_2){
//             if(max_count_2 < u.second){
//                 max_count_2 = u.second;
//                 matched_2 = u.first;
//             }
//         }
//         uint32_t max_unitig_num = 0;
//         for(auto u: count_1){
//             if(max_count_1 == u.second){
//                 max_unitig_num++;
//             }
//         }
//         for(auto u: count_2){
//             if(max_count_2 == u.second){
//                 max_unitig_num++;
//             }
//         }
//         if(matched_1!=matched_2 && matched_2!=(matched_1^1)){
//             vector<uint32_t> buf_to_insert;
//             buf_to_insert.push_back(max(matched_1, matched_2));
//             buf_to_insert.push_back(min(matched_1, matched_2));
//             if(connection_count.find(buf_to_insert)==connection_count.end()){
//                 connection_count[buf_to_insert] = 0;
//             }
//             connection_count[buf_to_insert]++;
//         }
//     }

// }


typedef struct { // global data structure for kt_pipeline()
	asg_t* graph;
    map<uint32_t,bubble_t*>* node_bubble_map;
    map<uint16_t,map<uint16_t,uint16_t>>* connections_count;
    map<uint32_t,set<uint32_t>>* node_path_id_map;
} shared_data;

typedef struct { // data structure for each step in kt_pipeline()
	shared_data *p;
	vector<uint32_t>* beg_node;
    vector<uint32_t>* end_node;
    vector<set<uint32_t>>* current_nodes;
    string* results;
} step_data;


string get_haplotype_sequence(asg_t* graph, set<uint32_t> set_of_nodes, uint32_t begin, uint32_t end, int id){
    stringstream result;
    uint32_t node = begin;
    uint32_t num_outgoing_arcs = asg_arc_n(graph, node);
    asg_arc_t *outgoing_arcs = asg_arc_a(graph, node);
    set<uint32_t> visited_nodes;
    vector<uint32_t> visited_arcs;
    // visited_arcs.push_back(node);
    result << string(">")+string(graph->seq[begin>>1].name) + string("_") + string(string(graph->seq[end>>1].name))+string("_hap")+to_string(id);
    result << endl;
    string pre_seq = string(graph->seq[node>>1].seq);
    if(node&1){
        pre_seq = complement(pre_seq);
    }
    while(node!=-1){
        visited_nodes.insert(node);
        // bool found = false;
        if(node == end){
            result << pre_seq;
            break;
        }
        map<uint32_t, asg_arc_t> next_node_arc;
        for(int i = 0; i < num_outgoing_arcs; i++){
            next_node_arc[outgoing_arcs[i].v] = outgoing_arcs[i];
        }

        uint32_t next_node = -1;
        uint32_t max_counter = 0;
        for(auto node: next_node_arc){
            if(set_of_nodes.find(node.first>>1)!=set_of_nodes.end() && visited_nodes.find(node.first)==visited_nodes.end() ){
                uint32_t num_outgoing_arcs = asg_arc_n(graph, node.first);
                asg_arc_t *outgoing_arcs = asg_arc_a(graph, node.first);
                set<uint32_t> current_next_next;
                for(int i = 0; i < num_outgoing_arcs; i++){
                    if(set_of_nodes.find(outgoing_arcs[i].v>>1)!=set_of_nodes.end()){
                        current_next_next.insert(outgoing_arcs[i].v);
                    }
                }
                uint32_t score = current_next_next.size();
                if(visited_nodes.find(node.first^1)==visited_nodes.end()){
                    score+=100;
                }
                if(score>=max_counter){
                    next_node = node.first;
                    max_counter = score;
                }
            }
        }
        visited_arcs.push_back(node);
        if(next_node == -1){
            bool found = false;
            uint32_t num_outgoing_arcs_temp = asg_arc_n(graph, node);
            asg_arc_t *outgoing_arcs_temp = asg_arc_a(graph, node);
            for(int i = 0; i < num_outgoing_arcs_temp; i++){
                if(found){
                    break;
                }
                uint32_t num_outgoing_arcs_temp_sec = asg_arc_n(graph, outgoing_arcs_temp[i].v);
                asg_arc_t *outgoing_arcs_temp_sec = asg_arc_a(graph, outgoing_arcs_temp[i].v);
                for(int j = 0; j < num_outgoing_arcs_temp_sec; j++){
                    if(set_of_nodes.find(outgoing_arcs_temp_sec[j].v>>1)!=set_of_nodes.end() && visited_nodes.find(outgoing_arcs_temp_sec[j].v)==visited_nodes.end()){
                        found = true;
                        next_node = outgoing_arcs_temp[i].v;
                        break;
                    }
                }
            }
            if(next_node == -1){
                stringstream res;
                res << endl;
                res << "Start-End: " << string(graph->seq[begin>>1].name) << ", " << string(graph->seq[end>>1].name)+string(" hap")+to_string(id) << endl;
                res << "All nodes: ";
                for(auto node_buf: set_of_nodes){
                    res << graph->seq[node_buf].name <<", ";
                }   
                res << endl;
                res << "Visited: ";
                for(auto node_buf: visited_arcs){
                    res << graph->seq[node_buf>>1].name << (node_buf%2==1?'-':'+')<< ", ";
                }   
                res << endl;
                cout << res.str();
                cout << endl;
                break;
            }
        }
        assert(next_node!=-1);

        node = next_node;
        result << pre_seq.substr(0,pre_seq.length() - next_node_arc[node].ol);
        pre_seq = string(graph->seq[node>>1].seq);
        if(node&1){
            pre_seq = complement(pre_seq);
        }

        // if(!found) {
        //     break;
        // }
        num_outgoing_arcs = asg_arc_n(graph, node);
        outgoing_arcs = asg_arc_a(graph, node);
    }
    result << endl;
    return result.str();
}

static void worker_for_single_step(void *data, long i, int tid) // callback for kt_for()
{
	step_data *s = (step_data*)data;
    set<uint32_t> current_nodes = (*s->current_nodes)[i];
    asg_t* graph = s->p->graph;
    map<uint32_t,set<uint32_t>>* node_path_id_map = s->p->node_path_id_map;
    map<uint32_t,bubble_t*>* node_bubble_map = s->p->node_bubble_map;
    map<uint16_t,map<uint16_t,uint16_t>>* connections_count = s->p->connections_count;

    // if(current_nodes.size()>4){
        set<bubble_t*> cur_bubbles;
        for(auto node:current_nodes){
            if(node_bubble_map->find(node)!=node_bubble_map->end()){
                cur_bubbles.insert((*node_bubble_map)[node]);
            }
        }
        // cout << current_nodes.size() << endl;
        map<bubble_t*, uint32_t> bubble_path_id_1;
        map<bubble_t*, uint32_t> bubble_path_id_2;
        set<uint32_t> path_1_not_included;
        set<uint32_t> path_2_not_included;
        map<bubble_t*, map<bubble_t*, map<uint32_t, map<uint32_t,uint32_t>>>> bubble_connection;
        for(auto node_1: current_nodes){
            for(auto node_2: current_nodes){

                if(node_1 > node_2
                && node_bubble_map->find(node_1) != node_bubble_map->end()
                && node_bubble_map->find(node_2) != node_bubble_map->end()
                && connections_count->find(node_1) != connections_count->end()
                && (*connections_count)[node_1].find(node_2) != (*connections_count)[node_1].end()){
                    bubble_t* bub_1 = (*node_bubble_map)[node_1];
                    bubble_t* bub_2 = (*node_bubble_map)[node_2];
                    set<uint32_t> pids_1 = (*node_path_id_map)[node_1];
                    set<uint32_t> pids_2 = (*node_path_id_map)[node_2];
                    // vector<bubble_t*> to_insert_bub;
                    // to_insert_bub.push_back(bub_1);
                    // to_insert_bub.push_back(bub_2);
                    if(bubble_connection.find(bub_1)==bubble_connection.end()){
                        bubble_connection[bub_1] = map<bubble_t*, map<uint32_t, map<uint32_t,uint32_t>>>();
                    }
                    if(bubble_connection[bub_1].find(bub_2)==bubble_connection[bub_1].end()){
                        bubble_connection[bub_1][bub_2] = map<uint32_t, map<uint32_t,uint32_t>>();
                    }
                    for(auto pid_1: pids_1){
                        for(auto pid_2:pids_2){
                            if(bubble_connection[bub_1][bub_2].find(pid_1)==bubble_connection[bub_1][bub_2].end()){
                                bubble_connection[bub_1][bub_2][pid_1] = map<uint32_t,uint32_t>();
                            }
                            if(bubble_connection[bub_1][bub_2][pid_1].find(pid_2)==bubble_connection[bub_1][bub_2][pid_1].end()){
                                bubble_connection[bub_1][bub_2][pid_1][pid_2] = 0;
                            }
                            bubble_connection[bub_1][bub_2][pid_1][pid_2] += (*connections_count)[node_1][node_2];
                        }
                    }

                    if(bubble_connection.find(bub_2)==bubble_connection.end()){
                        bubble_connection[bub_2] = map<bubble_t*, map<uint32_t, map<uint32_t,uint32_t>>>();
                    }
                    if(bubble_connection[bub_2].find(bub_1)==bubble_connection[bub_2].end()){
                        bubble_connection[bub_2][bub_1] = map<uint32_t, map<uint32_t,uint32_t>>();
                    }
                    for(auto pid_1: pids_1){
                        for(auto pid_2:pids_2){
                            if(bubble_connection[bub_2][bub_1].find(pid_2)==bubble_connection[bub_2][bub_1].end()){
                                bubble_connection[bub_2][bub_1][pid_2] = map<uint32_t,uint32_t>();
                            }
                            if(bubble_connection[bub_2][bub_1][pid_2].find(pid_1)==bubble_connection[bub_2][bub_1][pid_2].end()){
                                bubble_connection[bub_2][bub_1][pid_2][pid_1] = 0;
                            }
                            bubble_connection[bub_2][bub_1][pid_2][pid_1] += (*connections_count)[node_1][node_2];
                        }
                    }
                }
            }
        }
        cout << "Get connections between bubble" << endl;
        bubble_t* max_bubble_1;
        bubble_t* max_bubble_2;
        uint32_t max_path_id_1;
        uint32_t max_path_id_2;
        uint32_t max_count = 0;
        for(auto bub_1: bubble_connection){
            for(auto bub_2: bub_1.second){
                for(auto path_1: bub_2.second){
                    for(auto path_2: path_1.second){
                        if(max_count < path_2.second){
                            max_count = path_2.second;
                            max_bubble_1 = bub_1.first;
                            max_bubble_2 = bub_2.first;
                            max_path_id_1 = path_1.first;
                            max_path_id_2 = path_2.first;
                        }
                    }
                }
            }
        }
        if(max_count > 0){
            bubble_path_id_1[max_bubble_1] = max_path_id_1;
            bubble_path_id_1[max_bubble_2] = max_path_id_2;
            if(max_bubble_1->paths_nodes.size()==2){
                bubble_path_id_2[max_bubble_1] = 1;
                if(max_path_id_1==1){
                    bubble_path_id_2[max_bubble_1] = 0;
                }
            }
            if(max_bubble_2->paths_nodes.size()==2){
                bubble_path_id_2[max_bubble_2] = 1;
                if(max_path_id_2==1){
                    bubble_path_id_2[max_bubble_2] = 0;
                }
            }
            while(bubble_connection.size()!= bubble_path_id_1.size()){
                bubble_t* max_bubble;
                uint32_t max_path_id;
                uint32_t max_count = 0;

                map<bubble_t*, map<uint32_t,uint32_t>> extension_count;
                for(auto bub_1: bubble_path_id_1){
                    for(auto new_bub: bubble_connection[bub_1.first]){
                        if(bubble_path_id_1.find(new_bub.first)==bubble_path_id_1.end()){
                            for(auto new_path_id: new_bub.second[bub_1.second]){
                                if(extension_count.find(new_bub.first)==extension_count.end()){
                                    extension_count[new_bub.first] = map<uint32_t,uint32_t>();
                                }
                                if(extension_count[new_bub.first].find(new_path_id.first)==extension_count[new_bub.first].end()){
                                    extension_count[new_bub.first][new_path_id.first] = 0;
                                }
                                extension_count[new_bub.first][new_path_id.first] += new_path_id.second;
                            }
                        }
                    }
                }

                for(auto new_bub: extension_count){
                    for(auto new_path_id: new_bub.second){
                        if(new_path_id.second>max_count){
                            max_bubble = new_bub.first;
                            max_path_id = new_path_id.first;
                            max_count = new_path_id.second;
                        }
                    }
                }
                if(max_count!= 0){
                    bubble_path_id_1[max_bubble] = max_path_id;
                    if(max_bubble->paths_nodes.size()==2){
                        bubble_path_id_2[max_bubble] = 1;
                        if(max_path_id==1){
                            bubble_path_id_2[max_bubble] = 0;
                        }
                    }
                }else{
                    break;
                }
            }
            cout << "Get connected path in hap1" << endl;
            max_count = 0;
            for(auto bub_1: bubble_connection){
                for(auto bub_2: bub_1.second){
                    for(auto path_1: bub_2.second){
                        for(auto path_2: path_1.second){
                            if(max_count < path_2.second){
                                if( (bubble_path_id_1.find(bub_1.first)==bubble_path_id_1.end() || bubble_path_id_1[bub_1.first] != path_1.first)
                                && (bubble_path_id_1.find(bub_2.first)==bubble_path_id_1.end() || bubble_path_id_1[bub_2.first] != path_2.first)
                                ){
                                    max_count = path_2.second;
                                    max_bubble_1 = bub_1.first;
                                    max_bubble_2 = bub_2.first;
                                    max_path_id_1 = path_1.first;
                                    max_path_id_2 = path_2.first;
                                }
                            }
                        }
                    }
                }
            }
            if(max_count > 0){
                bubble_path_id_2[max_bubble_1] = max_path_id_1;
                bubble_path_id_2[max_bubble_2] = max_path_id_2;
                if(bubble_path_id_1.find(max_bubble_1) == bubble_path_id_1.end()){
                    if(max_path_id_1>0){
                        bubble_path_id_1[max_bubble_1] = 0;
                    }else{
                        bubble_path_id_1[max_bubble_1] = 1;
                    }
                }
                if(bubble_path_id_1.find(max_bubble_2) == bubble_path_id_1.end()){
                    if(max_path_id_2>0){
                        bubble_path_id_1[max_bubble_2] = 0;
                    }else{
                        bubble_path_id_1[max_bubble_2] = 1;
                    }
                }


            }
            while(bubble_connection.size()!= bubble_path_id_2.size()){
                bubble_t* max_bubble;
                uint32_t max_path_id;
                uint32_t max_count = 0;

                map<bubble_t*, map<uint32_t,uint32_t>> extension_count;
                for(auto bub_1: bubble_path_id_2){
                    for(auto new_bub: bubble_connection[bub_1.first]){
                        if(bubble_path_id_2.find(new_bub.first)==bubble_path_id_2.end()){
                            for(auto new_path_id: new_bub.second[bub_1.second]){
                                if(extension_count.find(new_bub.first)==extension_count.end()){
                                    extension_count[new_bub.first] = map<uint32_t,uint32_t>();
                                }
                                if(extension_count[new_bub.first].find(new_path_id.first)==extension_count[new_bub.first].end()){
                                    extension_count[new_bub.first][new_path_id.first] = 0;
                                }
                                extension_count[new_bub.first][new_path_id.first] += new_path_id.second;
                            }
                        }
                    }
                }

                for(auto new_bub: extension_count){
                    for(auto new_path_id: new_bub.second){
                        if(new_path_id.second>max_count){
                            max_bubble = new_bub.first;
                            max_path_id = new_path_id.first;
                            max_count = new_path_id.second;
                        }
                    }
                }
                if(max_count!= 0){
                    bubble_path_id_2[max_bubble] = max_path_id;
                    if(bubble_path_id_1.find(max_bubble) == bubble_path_id_1.end()){
                        if(max_path_id>0){
                            bubble_path_id_1[max_bubble] = 0;
                        }else{
                            bubble_path_id_1[max_bubble] = 1;
                        }
                    }
                }else{
                    break;
                }
            }
            cout << "Get connected path in hap2" << endl;
        }

        for(auto bubble: cur_bubbles){
            if(bubble_path_id_1.find(bubble)== bubble_path_id_1.end()){
                bubble_path_id_1[bubble] = 0;
            }
        }
        for(auto bubble: cur_bubbles){
            if(bubble_path_id_2.find(bubble)== bubble_path_id_2.end()){
                bubble_path_id_2[bubble] = 0;
            }
        }

        set<uint32_t> path_1_included;
        set<uint32_t> path_2_included;

        for(auto bubble: bubble_path_id_1){
            set<uint32_t> cur_path = bubble.first->paths_nodes[bubble.second];
            set<uint32_t> current_path;
            for(auto node: cur_path){
                current_path.insert(node>>1);
                path_1_included.insert(node>>1);
            }
            for(auto path: bubble.first->paths_nodes){
                for(auto node: path){
                    if(current_path.find(node>>1)== current_path.end()){
                        path_1_not_included.insert(node>>1);
                    }
                }
            }
        }

        for(auto bubble: bubble_path_id_2){
            set<uint32_t> cur_path = bubble.first->paths_nodes[bubble.second];
            set<uint32_t> current_path;
            for(auto node: cur_path){
                current_path.insert(node>>1);
                path_2_included.insert(node>>1);
            }
            for(auto path: bubble.first->paths_nodes){
                for(auto node: path){
                    if(current_path.find(node>>1)== current_path.end()){
                        path_2_not_included.insert(node>>1);
                    }
                }
            }
        }

        cout << "Fill gaps in haps" << endl;
        cout << bubble_connection.size() << "\t";
        cout << bubble_path_id_1.size() << "\t";
        cout << bubble_path_id_2.size() << endl;
        // cout << "haplotype 1:" << endl;
        for(auto node: current_nodes){
            if(path_1_not_included.find(node)== path_1_not_included.end()){
                path_1_included.insert(node);
                // cout << graph->seq[node].name << ", ";
            }
        }
        cout << endl;
        // cout << "haplotype 2:" << endl;
        for(auto node: current_nodes){
            if(path_2_not_included.find(node)== path_2_not_included.end()){
                path_2_included.insert(node);
                // cout << graph->seq[node].name << ", ";
            }
        }
        cout << graph->seq[(*s->beg_node)[i]>>1].name << " to " << graph->seq[(*s->end_node)[i]>>1].name << endl;


        stringstream res;
        res << get_haplotype_sequence(graph, path_1_included, (*s->beg_node)[i], (*s->end_node)[i], 1);
        res << get_haplotype_sequence(graph, path_2_included, (*s->beg_node)[i], (*s->end_node)[i], 2);
        (s->results)[i] = res.str();
    // }
}

void get_haplotype_path(uint16_t** connection_count,asg_t *graph, map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_graph, char* output_directory, int n_threads){
    cout << "Start get haplotypes" << endl;
    map<uint16_t,map<uint16_t,uint16_t>> connections_count;
    for(int i = 0; i < graph->n_seq; i++){
        for(int j = 0; j < graph->n_seq; j++){
            if(i>j && connection_count[i][j]>0){
                if(connections_count.find(i)==connections_count.end()){
                    connections_count[i] = map<uint16_t,uint16_t>();
                }
                connections_count[i][j] = connection_count[i][j];
            }
        }
        free(connection_count[i]);
    }
    free(connection_count);

    uint32_t n_vtx = graph->n_seq * 2;
    int node_type[n_vtx];

    vector<bubble_t*> bubbles;
    vector<bubble_t*> pure_bubbles;
    vector<bubble_t*> complex_bubbles;
    map<uint32_t,bubble_t*> node_bubble_map;
    map<uint32_t,set<uint32_t>> node_path_id_map;

    for(uint32_t i=0; i<n_vtx; i++){
        node_type[i] = 0;
        bubble_t* result = detect_bubble(graph, i);
        if(result!=nullptr){
            result->id = i;
            bubbles.push_back(result);
        }
    }


    for (int b=0; b<bubbles.size(); b++) {
        bubble_t* bubble = bubbles[b];
        uint32_t bubble_beginning = bubble->begNode;
        uint32_t bubble_end = bubble->endNode;

        // cout << "start get bubble paths from " << g->seq[bubble_beginning /2].name<< " to " << g->seq[bubble_end /2].name << endl;

        vector<asg_arc_t*> arc_stack;
        vector<uint32_t> node_stack;  // DFS
        vector<uint32_t> node_vi_stack;
        node_stack.push_back(bubble_beginning);
        node_vi_stack.push_back(0);
        int stack_count = 0;
        uint32_t lastNode = bubble_beginning;
        while(!node_stack.empty()) {
            uint32_t u = node_stack.back();
            uint32_t vi = node_vi_stack.back();
            if (u==bubble_end) {
                for (int ui=0; ui<node_stack.size(); ui++) {
                    // cout << g->seq[node_stack[ui]/2].name << " ";
                    if(ui==0 || ui == node_stack.size()-1){
                        if(node_type[node_stack[ui]] != 2){
                            node_type[node_stack[ui]] = BUBBLE_END_BEGIN;
                        }
                    }else{
                        node_type[node_stack[ui]] = BUBBLE_INSIDE;
                    }
                }
                stack_count ++;
                // cout << endl;
                set<uint32_t> buf_path_node = set<uint32_t>();
                for(auto arc_buf: arc_stack){
                    // if(arc_buf->v!=bubble->endNode){
                        buf_path_node.insert(arc_buf->v);
                    // }
                }
                buf_path_node.insert(bubble->begNode);
                buf_path_node.insert(bubble->endNode);
                // bubble_nodes.insert(buf_path_node);
                bubble->paths_nodes.push_back(buf_path_node);
                bubble->paths.push_back(arc_stack);
                bubble->starting_arcs.insert(arc_stack.front()->v);
                bubble->ending_arcs.insert(lastNode);
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
                continue;
            }

            uint32_t num_outgoing_arcs = asg_arc_n(graph, u);
            if (vi < num_outgoing_arcs) {
                asg_arc_t *outgoing_arcs = asg_arc_a(graph, u);
                uint32_t v = outgoing_arcs[vi].v;
                node_vi_stack.back()++;
                arc_stack.push_back(outgoing_arcs + vi);
                node_stack.push_back(v);
                node_vi_stack.push_back(0);
            } else {
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
            }
            lastNode = u;
        }
        // cout << stack_count << endl;
        // cout << "finish get bubble paths from " << bubble_beginning << " to " << bubble_end << endl;
    }



    for(auto bubble: bubbles){
        // if(node_type[bubble->begNode]!=2 && node_type[bubble->endNode]!=2 && bubble->paths_nodes.size()>20){
        if(node_type[bubble->begNode]!=2 && node_type[bubble->endNode]!=2){
            set<uint32_t> nodes;
            for(auto path:bubble->paths_nodes){
                for(auto node: path){
                    nodes.insert(node);
                }
            }
            if(nodes.size()>12){
                complex_bubbles.push_back(bubble);
            }else{
                pure_bubbles.push_back(bubble);
            }
            // cout << endl;
            // for(auto path: bubble->paths_nodes){
            //     for(auto node: path){
            //         cout << graph->seq[node>>1].name << "\t";
            //     }
            //     cout << endl;
            // }
        }else{
            delete bubble;
        }
    }

    // set<uint32_t> nodes;
    // set<uint32_t> double_nodes;
    for(auto bubble: pure_bubbles){
        for(int i = 0; i < bubble->paths_nodes.size(); i++){
            for(auto node: bubble->paths_nodes[i]){
                if(node!=bubble->endNode && node!=bubble->begNode){
                    node_bubble_map[node>>1] = bubble;
                    if(node_path_id_map.find(node>>1)==node_path_id_map.end()){
                        node_path_id_map[node>>1] = set<uint32_t>();
                    }
                    node_path_id_map[node>>1].insert(i);
                }
            }
        }
    }

    cout << "Bubbles Retrieved" << endl;

    // map<uint32_t,map<uint32_t,uint32_t>> cnt_count_branches;
    map<uint32_t,map<uint32_t,set<uint32_t>>> bubble_chain_graph_full_buf;
    map<uint32_t,map<uint32_t,set<uint32_t>>> bubble_chain_graph_full;
    map<uint32_t,vector<vector<uint32_t>>> node_beg_end;
    for(auto begin: *bubble_chain_graph){
        for(auto end: begin.second){
            set<uint32_t> to_insert;
            for(auto node: end.second){
                to_insert.insert(node>>1);
            }


            if(bubble_chain_graph_full_buf.find(begin.first) != bubble_chain_graph_full_buf.end()){
                bubble_chain_graph_full_buf[begin.first] = map<uint32_t,set<uint32_t>>();
            }
            bubble_chain_graph_full_buf[begin.first][end.first] = to_insert;

            if(bubble_chain_graph_full_buf.find(end.first^1) != bubble_chain_graph_full_buf.end()){
                bubble_chain_graph_full_buf[end.first^1] = map<uint32_t,set<uint32_t>>();
            }
            bubble_chain_graph_full_buf[end.first^1][begin.first^1] = to_insert;

        }
    }
    delete bubble_chain_graph;
    
    
    for(auto begin: bubble_chain_graph_full_buf){
        uint64_t max_len = 0;
        for(auto buf_end: begin.second){
            uint64_t cur_len = 0;
            for(auto cur_node: buf_end.second){
                cur_len+=graph->seq[cur_node].len;
            }
            max_len = max(cur_len, max_len);
        }
        for(auto buf_end: (bubble_chain_graph_full_buf[begin.first^1])){
            uint64_t cur_len = 0;
            for(auto cur_node: buf_end.second){
                cur_len+=graph->seq[cur_node].len;
            }
            max_len = max(cur_len, max_len);
        }
        for(auto end: begin.second){
            uint64_t cur_len = 0;
            for(auto cur_node: end.second){
                cur_len+=graph->seq[cur_node].len;
            }
            if(cur_len>max_len*0.3){
                set<uint32_t> to_insert;
                vector<uint32_t> vec1;
                vec1.push_back(begin.first);
                vec1.push_back(end.first);
                vector<uint32_t> vec2;
                vec1.push_back(end.first^1);
                vec1.push_back(begin.first^1);
                for(auto node: end.second){
                    to_insert.insert(node);
                    node_beg_end[node] = vector<vector<uint32_t>>();
                    node_beg_end[node].push_back(vec1);
                    node_beg_end[node].push_back(vec2);
                }


                if(bubble_chain_graph_full.find(begin.first) != bubble_chain_graph_full.end()){
                    bubble_chain_graph_full[begin.first] = map<uint32_t,set<uint32_t>>();
                }
                bubble_chain_graph_full[begin.first][end.first] = to_insert;

                if(bubble_chain_graph_full.find(end.first^1) != bubble_chain_graph_full.end()){
                    bubble_chain_graph_full[end.first^1] = map<uint32_t,set<uint32_t>>();
                }
                bubble_chain_graph_full[end.first^1][begin.first^1] = to_insert;

            }
        }
    }

    map<vector<uint32_t>, map<vector<uint32_t>,uint32_t>> branch_connections;

    for(auto x: connections_count){
        if(node_beg_end.find(x.first)!=node_beg_end.end()){
            vector<vector<uint32_t>> vec1s = node_beg_end[x.first];
            for(auto y: x.second){
                if(node_beg_end.find(y.first)!=node_beg_end.end() && y.second>0){
                    vector<vector<uint32_t>> vec2s = node_beg_end[y.first];
                    for(auto vec1: vec1s){
                        for(auto vec2: vec2s){
                            bool valid = false;
                            for(auto node_1: vec1){
                                for(auto node_2: vec2){
                                    if(node_1 == node_2){
                                        valid = (!valid);
                                    }
                                }
                            }
                            if(valid){
                                if(branch_connections.find(vec1)==branch_connections.end()){
                                    branch_connections[vec1] = map<vector<uint32_t>,uint32_t>();
                                }
                                if(branch_connections[vec1].find(vec2)==branch_connections[vec1].end()){
                                    branch_connections[vec1][vec2] = 0;
                                }
                                branch_connections[vec1][vec2] += y.second;

                                if(branch_connections.find(vec2)==branch_connections.end()){
                                    branch_connections[vec2] = map<vector<uint32_t>,uint32_t>();
                                }
                                if(branch_connections[vec2].find(vec1)==branch_connections[vec2].end()){
                                    branch_connections[vec2][vec1] = 0;
                                }
                                branch_connections[vec2][vec1] += y.second;
                            }
                        }
                    }
                }
            }
        }
    }

    cout << "Get connections between branches" << endl;
    set<vector<vector<uint32_t>>> to_expand;
    for(auto begin: bubble_chain_graph_full){
        for(auto ending: begin.second){
            uint32_t beg = begin.first;
            uint32_t end = ending.first;
            vector<uint32_t> to_insert_single;
            to_insert_single.push_back(beg);
            to_insert_single.push_back(end);
            vector<vector<uint32_t>> to_insert;
            to_insert.push_back(to_insert_single);
            to_expand.insert(to_insert);
        }
    }
    set<vector<vector<uint32_t>>> all_pathes = to_expand;
    while(to_expand.size()>0){
        set<vector<vector<uint32_t>>> new_to_expand;
        set<vector<vector<uint32_t>>> not_touched;
        for(auto existing_path : to_expand){
            uint32_t path_beg = existing_path[0][0];
            uint32_t path_end = existing_path[existing_path.size()-1][1];
            set<vector<uint32_t>> beginning_expansion;
            set<vector<uint32_t>> ending_expansion;
            if(bubble_chain_graph_full.find(path_beg^1)!=bubble_chain_graph_full.end()){
                if(bubble_chain_graph_full[path_beg^1].size()==1){
                    uint32_t cur_end = path_beg;
                    for(auto ending: bubble_chain_graph_full[path_beg^1]){
                        uint32_t cur_beg = ending.first^1;
                        vector<uint32_t> to_insert_single;
                        to_insert_single.push_back(cur_beg);
                        to_insert_single.push_back(cur_end);
                        bool seen = false;
                        for(auto p: existing_path){
                            if(p==to_insert_single){
                                seen = true;
                                break;
                            }
                        }
                        if(!seen){
                            beginning_expansion.insert(to_insert_single);
                        }
                    }
                }else if(bubble_chain_graph_full[path_beg^1].size()>0){
                    uint32_t cur_end = path_beg;
                    for(auto ending: bubble_chain_graph_full[path_beg^1]){
                        uint32_t cur_beg = ending.first^1;
                        vector<uint32_t> to_insert_single;
                        to_insert_single.push_back(cur_beg);
                        to_insert_single.push_back(cur_end);
                        if( branch_connections.find(to_insert_single)!=branch_connections.end()
                            && branch_connections[to_insert_single].find(existing_path[0]) != branch_connections[to_insert_single].end()
                            && branch_connections[to_insert_single][existing_path[0]] > 5){
                            
                            bool seen = false;
                            for(auto p: existing_path){
                                if(p==to_insert_single){
                                    seen = true;
                                    break;
                                }
                            }
                            if(!seen){
                                beginning_expansion.insert(to_insert_single);
                            }
                        }
                    }
                }
            }

            if(bubble_chain_graph_full.find(path_end)!=bubble_chain_graph_full.end()){
                if(bubble_chain_graph_full[path_end].size()==1){
                    uint32_t cur_beg = path_end;
                    for(auto ending: bubble_chain_graph_full[path_end]){
                        uint32_t cur_end = ending.first;
                        vector<uint32_t> to_insert_single;
                        to_insert_single.push_back(cur_beg);
                        to_insert_single.push_back(cur_end);
                        bool seen = false;
                        for(auto p: existing_path){
                            if(p==to_insert_single){
                                seen = true;
                                break;
                            }
                        }
                        if(!seen){
                            ending_expansion.insert(to_insert_single);
                        }
                    }
                }else if(bubble_chain_graph_full[path_end].size()>0){
                    uint32_t cur_beg = path_beg;
                    for(auto ending: bubble_chain_graph_full[path_end]){
                        uint32_t cur_end = ending.first;
                        vector<uint32_t> to_insert_single;
                        to_insert_single.push_back(cur_beg);
                        to_insert_single.push_back(cur_end);
                        if( branch_connections.find(to_insert_single)!=branch_connections.end()
                            && branch_connections[to_insert_single].find(existing_path[existing_path.size()-1]) != branch_connections[to_insert_single].end()
                            && branch_connections[to_insert_single][existing_path[existing_path.size()-1]] > 5){
                            
                            bool seen = false;
                            for(auto p: existing_path){
                                if(p==to_insert_single){
                                    seen = true;
                                    break;
                                }
                            }
                            if(!seen){
                                ending_expansion.insert(to_insert_single);
                            }
                        }
                    }
                }
            }

            if(ending_expansion.size() > 0 || beginning_expansion.size() > 0){
                for(auto to_end: ending_expansion){
                    vector<vector<uint32_t>> to_insert_path = existing_path;
                    to_insert_path.push_back(to_end);
                    new_to_expand.insert(to_insert_path);
                }
                for(auto to_begin: beginning_expansion){
                    vector<vector<uint32_t>> to_insert_path = existing_path;
                    to_insert_path.insert(to_insert_path.begin(),to_begin);
                    new_to_expand.insert(to_insert_path);
                }
                for(auto to_end: ending_expansion){
                    for(auto to_begin: beginning_expansion){
                        vector<vector<uint32_t>> to_insert_path = existing_path;
                        to_insert_path.push_back(to_end);
                        to_insert_path.insert(to_insert_path.begin(),to_begin);
                        new_to_expand.insert(to_insert_path);
                    }
                }
            }else{
                not_touched.insert(existing_path);
            }
        }
        all_pathes.insert(not_touched.begin(), not_touched.end());
        to_expand = new_to_expand;
    }
    cout << "Get all possible pathes, total num: ";
    cout << all_pathes.size() << endl;

    for(auto x: all_pathes){
        uint32_t beg = -1;
        
        for(auto y: x){
            assert(beg==-1 || beg==y[0]);
            beg = y[1];
        }
    }

    set<vector<vector<uint32_t>>> expanded_pathes;
    for(auto path: all_pathes){
        if(expanded_pathes.find(path) == expanded_pathes.end()){
            set<uint32_t> path_set;
            for(auto sing_path: path){
                path_set.insert(bubble_chain_graph_full[sing_path[0]][sing_path[1]].begin(),bubble_chain_graph_full[sing_path[0]][sing_path[1]].end());
            }
            bool valid = true;
            for(auto to_check: all_pathes){
                set<uint32_t> check_path;
                for(auto sing_path: to_check){
                    check_path.insert(bubble_chain_graph_full[sing_path[0]][sing_path[1]].begin(),bubble_chain_graph_full[sing_path[0]][sing_path[1]].end());
                }
                bool all_found = true;
                uint32_t not_found_count = 0;
                for(auto node_buf: path_set){
                    if(check_path.find(node_buf)==check_path.end()){
                        all_found = false;
                        not_found_count++;
                        // break;
                    }
                }
                if( (all_found || (not_found_count <= ((double)(path_set.size()) * 0.1) && path.size() < to_check.size()) ) && path != to_check){
                    valid = false;
                    break;
                }
            }
            if(valid){
                expanded_pathes.insert(path);
            }
        }
    }

    cout << "Filtering short path, total num: ";
    cout << expanded_pathes.size() << endl;

    map<vector<uint32_t>,set<uint32_t>> expanded_pathes_nodes;
    map<uint32_t,set<uint32_t>> path_reach_out_set_beg;
    map<uint32_t,set<uint32_t>> path_reach_out_set_end;
    for(auto path: expanded_pathes){
        uint32_t beg = path[0][0];
        uint32_t end = path[path.size()-1][1];  
        if(path_reach_out_set_beg.find(beg)==path_reach_out_set_beg.end()){
            path_reach_out_set_beg[beg] = set<uint32_t>();
            uint32_t num_incoming_arcs_beg = asg_arc_n(graph, beg^1);
            asg_arc_t *incoming_arcs_beg = asg_arc_a(graph, beg^1);
            for(int i = 0; i < num_incoming_arcs_beg; i++){
                path_reach_out_set_beg[beg].insert(incoming_arcs_beg[i].v^1);
                uint32_t num_incoming_arcs_beg_sec = asg_arc_n(graph, incoming_arcs_beg[i].v^1);
                asg_arc_t *incoming_arcs_beg_sec = asg_arc_a(graph, incoming_arcs_beg[i].v^1);
                for(int j = 0; j < num_incoming_arcs_beg_sec; j++){
                    path_reach_out_set_beg[beg].insert(incoming_arcs_beg_sec[j].v^1);
                }
            }
        }
        if(path_reach_out_set_end.find(end)==path_reach_out_set_end.end()){
            path_reach_out_set_end[end] = set<uint32_t>();
            uint32_t num_outgoing_arcs_end = asg_arc_n(graph, end);
            asg_arc_t *outgoing_arcs_end = asg_arc_a(graph, end);
            for(int i = 0; i < num_outgoing_arcs_end; i++){
                path_reach_out_set_end[end].insert(outgoing_arcs_end[i].v);
                uint32_t num_outgoing_arcs_end_sec = asg_arc_n(graph, outgoing_arcs_end[i].v);
                asg_arc_t *outgoing_arcs_end_sec = asg_arc_a(graph, outgoing_arcs_end[i].v);
                for(int j = 0; j < num_outgoing_arcs_end_sec; j++){
                    path_reach_out_set_end[end].insert(outgoing_arcs_end_sec[j].v);
                }
            }
        }
    }

    for(auto path: expanded_pathes){
        set<uint32_t> nodes;
        for(auto single_path: path){
            nodes.insert(bubble_chain_graph_full[single_path[0]][single_path[1]].begin(),
            bubble_chain_graph_full[single_path[0]][single_path[1]].end());
            nodes.insert(single_path[0]>>1);
            nodes.insert(single_path[1]>>1);
        }
        uint32_t beg = path[0][0];
        uint32_t end = path[path.size()-1][1];
        bool inserted = false;
        for(auto to_reach_out: expanded_pathes){
            uint32_t r_beg = to_reach_out[0][0];
            uint32_t r_end = to_reach_out[to_reach_out.size()-1][1]; 
            set<uint32_t> cur_nodes;
            for(auto single_path: to_reach_out){
                cur_nodes.insert(bubble_chain_graph_full[single_path[0]][single_path[1]].begin(),
                bubble_chain_graph_full[single_path[0]][single_path[1]].end());
                cur_nodes.insert(single_path[0]>>1);
                cur_nodes.insert(single_path[1]>>1);
            } 
            if(path_reach_out_set_beg[r_beg].find(end)!=path_reach_out_set_beg[r_beg].end() && 
            ( (cur_nodes.size()<10 || nodes.size()<10) || (branch_connections.find(path[path.size()-1])!=branch_connections.end()
                && branch_connections[path[path.size()-1]].find(to_reach_out[0])!=branch_connections[path[path.size()-1]].end()
                && branch_connections[path[path.size()-1]][to_reach_out[0]]>5))){
                vector<uint32_t> vec;
                vec.push_back(beg);
                vec.push_back(to_reach_out[to_reach_out.size()-1][1]);
                inserted = true;
                cur_nodes.insert(nodes.begin(),nodes.end());
                bool found = false;
                uint32_t num_outgoing_arcs_end = asg_arc_n(graph, end);
                asg_arc_t *outgoing_arcs_end = asg_arc_a(graph, end);
                for(int i = 0; i < num_outgoing_arcs_end; i++){
                    if(found){
                        break;
                    }
                    if(outgoing_arcs_end[i].v==to_reach_out[0][0]){
                        found = true;
                        break;
                    }
                    uint32_t num_outgoing_arcs_end_sec = asg_arc_n(graph, outgoing_arcs_end[i].v);
                    asg_arc_t *outgoing_arcs_end_sec = asg_arc_a(graph, outgoing_arcs_end[i].v);
                    for(int j = 0; j < num_outgoing_arcs_end_sec; j++){
                        if(outgoing_arcs_end_sec[j].v==to_reach_out[0][0]){
                            found = true;
                            cur_nodes.insert(outgoing_arcs_end[i].v);
                            break;
                        }
                    }
                }
                cout << "Connected " << graph->seq[beg>>1].name << " , " << graph->seq[end>>1].name << " , ";
                cout << graph->seq[to_reach_out[0][0]>>1].name << " , " << graph->seq[to_reach_out[to_reach_out.size()-1][1]>>1].name << " , " << "Total Num: " << cur_nodes.size() << endl;
                expanded_pathes_nodes[vec] = cur_nodes;
            }else if(path_reach_out_set_end[r_end].find(beg)!=path_reach_out_set_end[r_end].end() && 
            ( (cur_nodes.size()<10 || nodes.size()<10) || (branch_connections.find(to_reach_out[to_reach_out.size()-1])!=branch_connections.end()
                && branch_connections[to_reach_out[to_reach_out.size()-1]].find(path[0])!=branch_connections[to_reach_out[to_reach_out.size()-1]].end()
                && branch_connections[to_reach_out[to_reach_out.size()-1]][path[0]]>5))){
                vector<uint32_t> vec;
                vec.push_back(to_reach_out[0][0]);
                vec.push_back(end);
                inserted = true;
                cur_nodes.insert(nodes.begin(),nodes.end());
                bool found = false;
                uint32_t num_incoming_arcs_beg = asg_arc_n(graph, beg^1);
                asg_arc_t *incoming_arcs_beg = asg_arc_a(graph, beg^1);
                for(int i = 0; i < num_incoming_arcs_beg; i++){
                    if(found){
                        break;
                    }
                    if(incoming_arcs_beg[i].v^1==to_reach_out[0][0]){
                        found = true;
                        break;
                    }
                    uint32_t num_incoming_arcs_beg_sec = asg_arc_n(graph, incoming_arcs_beg[i].v^1);
                    asg_arc_t *incoming_arcs_beg_sec = asg_arc_a(graph, incoming_arcs_beg[i].v^1);
                    for(int j = 0; j < num_incoming_arcs_beg_sec; j++){
                        if(incoming_arcs_beg_sec[j].v^1==to_reach_out[0][0]){
                            found = true;
                            cur_nodes.insert(incoming_arcs_beg[i].v^1);
                            break;
                        }
                    }
                }
                cout << "Connected " << graph->seq[to_reach_out[0][0]>>1].name << " , " << graph->seq[to_reach_out[to_reach_out.size()-1][1]>>1].name;
                cout  << " , "  << graph->seq[beg>>1].name << " , " << graph->seq[end>>1].name << " , " << "Total Num: " << cur_nodes.size() << endl;
                expanded_pathes_nodes[vec] = cur_nodes;
            }
        }
        if(!inserted){
            vector<uint32_t> vec;
            vec.push_back(beg);
            vec.push_back(end);
            expanded_pathes_nodes[vec] = nodes;
        }
    }

    uint32_t path_counting = 0;
    map<uint32_t,map<uint32_t,set<uint32_t>>> bubble_chain_graph_connected;
    for(auto path: expanded_pathes_nodes){
        uint32_t beg = path.first[0];
        uint32_t end = path.first[1];
        if((bubble_chain_graph_connected.find(end^1)!=bubble_chain_graph_connected.end()
        && bubble_chain_graph_connected[end^1].find(beg^1)!=bubble_chain_graph_connected[end^1].end()
        ) || (bubble_chain_graph_connected.find(beg)!=bubble_chain_graph_connected.end()
        && bubble_chain_graph_connected[beg].find(end)!=bubble_chain_graph_connected[beg].end()
        )){
            continue;
        }
        path_counting++;
        if(bubble_chain_graph_connected.find(beg)==bubble_chain_graph_connected.end()){
            bubble_chain_graph_connected[beg] = map<uint32_t,set<uint32_t>>();
        }
        bubble_chain_graph_connected[beg][end] = path.second;
    }

    cout << "Delete duplicated pathes, total unique pathes: ";
    cout << path_counting << endl;

    for(auto begin: bubble_chain_graph_connected){
        for(auto end: begin.second){
            cout << graph->seq[begin.first>>1].name <<" to " << graph->seq[end.first>>1].name << ": " << endl;
            for(auto node: end.second){
                cout << graph->seq[node].name << ", ";
            }
            cout << endl;
        }
    }

    shared_data shared;
    shared.graph = graph;
    shared.connections_count = &connections_count;
    shared.node_path_id_map = &node_path_id_map;
    shared.node_bubble_map = &node_bubble_map;
    step_data step;
    step.p = &shared;
    step.beg_node = new vector<uint32_t>();
    step.end_node = new vector<uint32_t>();
    step.current_nodes = new vector<set<uint32_t>>();
    // for(auto beg: bubble_chain_graph_full){
    for(auto beg: bubble_chain_graph_connected){
        for(auto end: beg.second){
            // if(end.second.size() < 100){
                step.beg_node->push_back(beg.first);
                step.end_node->push_back(end.first);
                step.current_nodes->push_back(end.second);
            // }
        }
    }
    step.results = (string*)calloc(step.beg_node->size(),sizeof(string));

    kt_for(n_threads, worker_for_single_step , &step, step.beg_node->size());

    ofstream outFile;
    outFile.open(string(output_directory)+string("/pred_haplotypes.fa"), ofstream::out | ofstream::trunc);
    for(int i = 0; i < step.beg_node->size(); i++){
        outFile << step.results[i];
    }
    outFile.close();
}
