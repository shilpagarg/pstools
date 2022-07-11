#ifndef PAF_INTERSECT_H
#define PAF_INTERSECT_H

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include "bubble_chain.h"
#include "paf.h"
using namespace std;

void output_paf_intersect(asg_t *g, paf_file_t *paf_file, char *paf_output_filename, paf_reader reader) {
    cout << "start get paf intersect" << endl;
    ofstream paf_output_file;
    paf_output_file.open(paf_output_filename);

    unordered_map<int, bool> graph_id_set;
    vector<uint32_t> topological_order = get_topological_order(g);
    for (uint u=0; u<topological_order.size(); u++) { // topological_order.size() is unsigned and u >= 0
        string name = g->seq[topological_order[u]/2].name;
        string subname = name.substr(3, 9);
        int subnameint = atoi(subname.c_str());
        graph_id_set[subnameint] = true;
    }

	paf_rec_t r;
	// int paf_line_count = 0;
	// while (paf_read(paf_file, &r) >= 0) {
	// 	cout << r.tn << "\t";
	// 	if (paf_line_count % 5 == 4) {
	// 		cout << endl;
	// 	}
	// 	paf_line_count++;
	// }
    cout << "start read paf file." << endl; 
    int count = 0;
	while (reader.paf_read(paf_file, &r) >= 0) {
        count++;
        string target_node_name = r.tn;
        string target_node_subname = target_node_name.substr(3, 9);
        int target_node_subnameint = atoi(target_node_subname.c_str());
        if (graph_id_set.find(target_node_subnameint) != graph_id_set.end()) {
            // see: https://github.com/lh3/miniasm/blob/master/PAF.md
            paf_output_file << r.qn << "\t";
            paf_output_file << r.ql << "\t";
            paf_output_file << r.qs << "\t";
            paf_output_file << r.qe << "\t";
            paf_output_file << (r.rev ? "-" : "+") << "\t";  // TODO
            paf_output_file << r.tn << "\t";
            paf_output_file << r.tl << "\t";
            paf_output_file << r.ts << "\t";
            paf_output_file << r.te << "\t";
            paf_output_file << r.ml << "\t";
            paf_output_file << r.bl << "\t";
            paf_output_file << "cs:Z:";
            paf_output_file << r.seq << "\t";
            // Mapping Quality is not parsed
            paf_output_file << endl;
            continue;  // TODO: continue the while not the for loop
        }
	}
	reader.paf_close(paf_file);
    paf_output_file.close();
    cout << count <<endl;
    cout << "finish get paf intersect" << endl;
    return;
}

#endif