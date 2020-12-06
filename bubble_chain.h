#ifndef _BUBBLE_CHAIN_H_
#define _BUBBLE_CHAIN_H_
#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <set>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "graph.h"
#include <map>
#include <sys/stat.h>
#include <string>
using namespace std;

typedef struct {
	uint32_t p; // the optimal parent vertex
	uint32_t d; // the shortest distance from the initial vertex
	uint32_t c; // max count of positive reads
	uint32_t m; // max count of negative reads
	uint32_t np; // max count of non-positive reads
	uint32_t nc; // max count of reads, no matter positive or negative
	uint32_t r:31, s:1; // r: the number of remaining incoming arc; s: state
	//s: state, s=0, this edge has not been visited, otherwise, s=1
} binfo_t;

typedef struct {
	///all information for each node
	binfo_t *a;
	kvec_t(uint32_t) S; // set of vertices without parents, nodes with all incoming edges visited
	kvec_t(uint32_t) T; // set of tips
	kvec_t(uint32_t) b; // visited vertices
	kvec_t(uint32_t) e; // visited edges/arcs
} buf_t;

typedef struct {
    uint32_t id;
    uint32_t begNode;
    uint32_t endNode;
    vector<vector<asg_arc_t*>> valid_paths;
    vector<vector<asg_arc_t*>> paths;
    vector<set<uint32_t>> paths_nodes;
    set<uint32_t> starting_arcs;
    set<uint32_t> ending_arcs;
} bubble_t;

typedef struct bubble_chain_t {
    uint32_t id;
    vector<bubble_t*> bubbles; //bubbles in order
    int *usl;  // unique sequence length of unique sequences between bubbles
};

typedef struct bubble_chain_graph_t {
    vector<uint32_t> normal_nodes;
    vector<bubble_chain_t*> condensed_bubble_chain;
};

#define NORMAL_NODE 0;
#define BUBBLE_END_BEGIN 1;
#define BUBBLE_INSIDE 2;
/*
 * Get the sources of a single component
 */
vector<uint32_t> get_sources(asg_t *g) {
    cout << "start source search" << endl;
    vector<uint32_t> sources;

	uint32_t n_vtx = g->n_seq * 2;
    bool frontier_added[n_vtx];
    for (int i=0; i<n_vtx; i++) {
        frontier_added[i] = false;
    }
    vector<uint32_t> frontier;
    uint32_t s0 = 0;  // arbitrary vertex
    frontier.push_back(s0);
    frontier_added[s0] = true;
    while(!frontier.empty()) {
        uint32_t u = frontier.back();
        frontier.pop_back();
        // cout << u << "\t" << g->seq[u/2].name << endl;

        uint32_t num_incoming_arcs = asg_arc_n(g, u^1);  // 1
        if (num_incoming_arcs == 0) {  // check if curr_v is a source by checking if no parents
            sources.push_back(u);
        }

        // add parents to frontier
        // uint32_t num_incoming_arcs = asg_arc_n(g, u^1);  // 1
        asg_arc_t *incoming_arcs_complement = asg_arc_a(g, u^1);  // incoming_arcs_complement[0].v^1 = 22;
        for (int vi=0; vi<num_incoming_arcs; vi++) {
            int p = incoming_arcs_complement[vi].v^1;
            if (!frontier_added[p]) {
                frontier.push_back(p);
                frontier_added[p] = true;
            }
        }

        // add children to frontier
        uint32_t num_outgoing_arcs = asg_arc_n(g, u);  // 2
        asg_arc_t *outgoing_arcs = asg_arc_a(g, u);  // p outgoing_arcs[0].v = 34; p outgoing_arcs[1].v = 37;
        for (int vi=0; vi<num_outgoing_arcs; vi++) {
            int c = outgoing_arcs[vi].v;
            if (!frontier_added[c]) {
                frontier.push_back(c);
                frontier_added[c] = true;
            }
        }
    }
    cout << "finish source search" << endl;
    return sources;
}

vector<uint32_t> get_topological_order(asg_t *g) {
    cout << "start get topological order" << endl;

    vector<uint32_t> reverse_topological_order;

	uint32_t n_vtx = g->n_seq * 2;
    bool frontier_added[n_vtx], child_frontier_added[n_vtx];
    for (int u=0; u<n_vtx; u++) {
        frontier_added[u] = false;
        child_frontier_added[u] = false;
    }
    vector<uint32_t> component, frontier, child_frontier, child_vi_frontier;
    uint32_t s0 = 0;  // arbitrary vertex
    frontier.push_back(s0);
    frontier_added[s0] = true;
    while(!frontier.empty()) {
        // for (int u=0; u<frontier.size(); u++) {
        //     cout << frontier[u] << ", ";
        // }
        // cout << endl;

        uint32_t u = frontier.back();
        frontier.pop_back();
        if (!child_frontier_added[u]) {
            child_frontier.push_back(u);
            child_vi_frontier.push_back(0);
            child_frontier_added[u] = true;
        }
        while(!child_frontier.empty()) {
            // for (int u=0; u<child_frontier.size(); u++) {
            //     cout << child_frontier[u] << ", ";
            // }
            // cout << endl;

            uint32_t u = child_frontier.back();
            uint32_t vi = child_vi_frontier.back();

            uint32_t num_outgoing_arcs = asg_arc_n(g, u);
            if (vi < num_outgoing_arcs) {
                asg_arc_t *outgoing_arcs = asg_arc_a(g, u);
                uint32_t v = outgoing_arcs[vi].v;
                child_vi_frontier.back()++;

                if (!child_frontier_added[v]) {
                    child_frontier.push_back(v);
                    child_vi_frontier.push_back(0);
                    frontier_added[v] = true;
                    child_frontier_added[v] = true;
                }
            } else {
                int num_incoming_arcs = asg_arc_n(g, u^1);
                asg_arc_t *incoming_arcs_complement = asg_arc_a(g, u^1);
                for (int uivj=0; uivj<num_incoming_arcs; uivj++) {
                    int vj = incoming_arcs_complement[uivj].v^1;
                    if (!frontier_added[vj]) {
                        frontier.push_back(vj);
                        frontier_added[vj] = true;
                    }
                }
                child_frontier.pop_back();
                child_vi_frontier.pop_back();
                reverse_topological_order.push_back(u);
            }
        }
    }

    reverse(reverse_topological_order.begin(), reverse_topological_order.end());
    vector<uint32_t> topological_order = reverse_topological_order;  // TODO
    // for (int u=0; u<topological_order.size(); u++) {
    //     cout << u << "\t" << g->seq[topological_order[u]/2].name << endl;
    // }

    cout << "finish get topological order" << endl;
    return topological_order;
}

vector<uint32_t> get_bubble_ends(asg_t *g, vector<uint32_t> sources, vector<uint32_t> *bubble_beginnings, vector<uint32_t> *bubble_ends) {
    cout << "start get bubble ends" << endl;

    uint32_t n_vtx = g->n_seq * 2;
    bool seen[n_vtx], visited[n_vtx];
    for (int u=0; u<n_vtx; u++) {
        seen[u] = false;
        visited[u] = false;
    }
    for(int ui = 0; ui<sources.size(); ui++){
        cout << "Source: " << sources[ui] << "\t" << g->seq[sources[ui]/2].name << endl;

    }
    vector<uint32_t> frontier = sources;
    int num_seen = frontier.size();  // seen = seen but not visited
    while(!frontier.empty()) {
        uint32_t u = frontier.back();
        frontier.pop_back();
        // cout << u << "\t" << g->seq[u/2].name << endl;  // topological order
            if(not visited[u]){
            visited[u] = true;
            seen[u] = false;
            num_seen--;

            uint32_t num_outgoing_arcs = asg_arc_n(g, u);
            uint32_t num_incoming_arcs = asg_arc_n(g, u^1);
            if (num_seen == 0) {
                if (num_outgoing_arcs > 1) {
                    bubble_beginnings->push_back(u);  // last one is meaningless/unpaired
                    cout << "Beginning:" << u << "\t" << g->seq[u/2].name<< endl;
                }
                bool sources_contains_u = find(sources.begin(), sources.end(), u) != sources.end();
                if (num_incoming_arcs > 1 ) {
                    bubble_ends->push_back(u);  // first one is meaningless/unpaired
                    cout << "Ending:" << u << "\t" << g->seq[u/2].name<< endl;
                }
            }

            // add children to frontier if they have no other unvisited parents
            // uint32_t num_outgoing_arcs = asg_arc_n(g, u);  // 2
            asg_arc_t *outgoing_arcs = asg_arc_a(g, u);  // p outgoing_arcs[0].v = 34; p outgoing_arcs[1].v = 37;
            for (int vi=0; vi<num_outgoing_arcs; vi++) {
                int c = outgoing_arcs[vi].v;

                uint32_t num_unvisited_wives = 0;
                // uint32_t num_incoming_arcs = asg_arc_n(g, c^1);
                num_incoming_arcs = asg_arc_n(g, c^1);
                asg_arc_t *incoming_arcs_complement = asg_arc_a(g, c^1);
                for (int viuj=0; viuj<num_incoming_arcs; viuj++) {
                    int w = incoming_arcs_complement[viuj].v^1;  // wives
                    if (!visited[w]) {
                        num_unvisited_wives++;
                    }
                }

                if (!seen[c]) {
                    num_seen++;
                }
                seen[c] = true;
                if (num_unvisited_wives == 0) {
                    frontier.push_back(c);
                }
            }
        }
    }
    cout << "finish get bubble ends: "<< num_seen << endl;
    for (int u=0; u<n_vtx; u++) {
        if(!visited[u]){
            cout <<"Unvisited: "<< u << endl;
        }
    }
    return *bubble_ends;
}

bubble_t* detect_bubble(asg_t *g, uint32_t source) {
    uint32_t num_source_outgoing_arcs = asg_arc_n(g, source);
	uint32_t n_vtx = g->n_seq * 2;
    uint32_t num_seen = 1;
    if(num_source_outgoing_arcs < 2){
        return nullptr;
    }

    bool seen[n_vtx], visited[n_vtx];
    for (int u=0; u<n_vtx; u++) {
        seen[u] = false;
        visited[u] = false;
    }
    vector<uint32_t> frontier;

    frontier.push_back(source);
    seen[source] = true;

    while(!frontier.empty()) {
        uint32_t v = frontier.back();
        frontier.pop_back();
        // cout << u << "\t" << g->seq[v/2].name << endl;  // topological order
        assert(seen[v]);
        assert(!visited[v]);
        visited[v] = true;
        seen[v] = false;
        num_seen--;
        uint32_t num_outgoing_arcs = asg_arc_n(g, v);

        if(num_outgoing_arcs == 0){
            return nullptr;
        }

        // add children to frontier if they have no other unvisited parents
        // uint32_t num_outgoing_arcs = asg_arc_n(g, v);  // 2
        asg_arc_t *outgoing_arcs = asg_arc_a(g, v);  // p outgoing_arcs[0].v = 34; p outgoing_arcs[1].v = 37;
        set<uint32_t> edge_set;
        for (int vi=0; vi<num_outgoing_arcs; vi++) {
            uint32_t u = outgoing_arcs[vi].v;
            if(edge_set.insert(u).second){
                if(u == v || u == (v^1) || visited[u^1] || u == source) {return nullptr;}

                assert(!visited[u]);
                if(!seen[u]){
                    num_seen++;
                }
                seen[u] = true;
                assert(asg_arc_n(g, u^1)>=1);
                bool has_unvisited_parents = false;
                uint32_t num_incoming_arcs = asg_arc_n(g, u^1);
                asg_arc_t *incoming_arcs_complement = asg_arc_a(g, u^1);
                for (int viuj=0; viuj<num_incoming_arcs; viuj++) {
                    uint32_t p = (incoming_arcs_complement[viuj].v)^1;  // parent
                    if (!visited[p]) has_unvisited_parents = true;
                }

                if (!has_unvisited_parents) {frontier.push_back(u);}
            }
        }

        if(frontier.size() == 1 && num_seen==1 && seen[frontier.front()]){
            uint32_t t = frontier.back();
            frontier.pop_back();
            uint32_t num_outgoing_arcs = asg_arc_n(g, t);
            if(num_outgoing_arcs > 0){
                asg_arc_t *outgoing_arcs = asg_arc_a(g, t);
                for(int vi=0; vi<num_outgoing_arcs; vi++){
                    if(outgoing_arcs[vi].v == source){
                        return nullptr;
                    }
                }
            }
            bubble_t *result = new bubble_t();
            result->begNode = source;
            result->endNode = t;
            return result;
        }
    }
    return nullptr;
}



map<uint32_t,map<uint32_t,set<uint32_t>>>* get_bubble_chain_graph(asg_t *g, set<uint32_t> bubble_chain_begin_end,string output_directory) {
    map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_begin_end_nodes = new map<uint32_t,map<uint32_t,set<uint32_t>>>();
    for(auto bubble_chain_begin: bubble_chain_begin_end){
        map<uint32_t,set<uint32_t>> current_map;
        vector<asg_arc_t*> arc_stack;
        vector<uint32_t> node_stack;  // DFS
        vector<uint32_t> node_vi_stack;
        node_stack.push_back(bubble_chain_begin);
        node_vi_stack.push_back(0);
        bool visited[g->n_seq * 2];
        for (int u=0; u<g->n_seq * 2; u++) {
            visited[u] = false;
        }
        while(!node_stack.empty()) {
            uint32_t u = node_stack.back();
            uint32_t vi = node_vi_stack.back();
            uint32_t num_outgoing_arcs = asg_arc_n(g, u);
            asg_arc_t *outgoing_arcs = asg_arc_a(g, u);
            visited[u] = true;
            int flag = 0;

            if (u != bubble_chain_begin && find(bubble_chain_begin_end.begin(),bubble_chain_begin_end.end(),u)!=bubble_chain_begin_end.end()) {
                flag = 1;
            }else{
                for(auto b:current_map){
                    if(flag){
                        break;
                    }
                    for(int i = 0; i < num_outgoing_arcs; i++){
                        if(b.second.find(outgoing_arcs[i].v)!=b.second.end()){
                            u = b.first;
                            flag = 2;
                            break;
                        }
                    }
                }
            }

            if(flag == 1){
                if(current_map.find(u)==current_map.end()){
                    set<uint32_t> new_set;
                    current_map.insert({u,new_set});
                }
                if(arc_stack.size() > 0 ){
                    for(int a = 0; a <arc_stack.size()-1; a++){
                        current_map[u].insert(arc_stack[a]->v);
                    }
                }
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
                continue;
            }else if(flag == 2){
                if(current_map.find(u)==current_map.end()){
                    set<uint32_t> new_set;
                    current_map.insert({u,new_set});
                }
                if(arc_stack.size() > 0 ){
                    for(int a = 0; a <arc_stack.size(); a++){
                        current_map[u].insert(arc_stack[a]->v);
                    }
                }
            }

            if (vi < num_outgoing_arcs) {
                while(vi < num_outgoing_arcs){
                    uint32_t v = outgoing_arcs[vi].v;
                    node_vi_stack.back()++;
                    if(!visited[v]){
                        arc_stack.push_back(outgoing_arcs + vi);
                        node_stack.push_back(v);
                        node_vi_stack.push_back(0);
                        break;
                    }else{
                        vi++;
                    }
                }
            }
            if (vi >= num_outgoing_arcs){
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
            }
        }
        bubble_chain_begin_end_nodes->insert({bubble_chain_begin,current_map});
    }

    for(auto a : *bubble_chain_begin_end_nodes){
        for(auto b : a.second){
            (*bubble_chain_begin_end_nodes)[a.first][b.first].insert(a.first);
            (*bubble_chain_begin_end_nodes)[a.first][b.first].insert(b.first);
        }
    }
    system((string("rm -r ")+output_directory).c_str());
    mkdir(output_directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    ofstream main_graph(output_directory+string("/")+string("main_graph.gfa"));
    for(auto a: *bubble_chain_begin_end_nodes){
            if(a.first%2==0){
                main_graph << "S\t" << g->seq[a.first>>1].name << "\t" << "*" <<endl;
            }
    }
    for(auto a: *bubble_chain_begin_end_nodes){
            for(auto b: a.second){
                main_graph << "L\t" << g->seq[a.first>>1].name << "\t"<< (a.first%2==0 ? "+" : "-") << "\t" << g->seq[b.first>>1].name << "\t" << (b.first%2==0 ? "+" : "-") << "\t0M\t" <<endl;
            }
    }
    main_graph.close();

    ofstream non_length_graph(output_directory+string("/")+string("no_longth_full_graph.gfa"));
    for(int i = 0; i<g->n_seq; i++){
        non_length_graph << "S\t" << g->seq[i].name << "\t" << "*" <<endl;
    }

    for(int i = 0; i<g->n_seq*2; i++){
        uint32_t num_outgoing_arcs = asg_arc_n(g, i);
        asg_arc_t *outgoing_arcs = asg_arc_a(g, i);
        for(int q = 0; q < num_outgoing_arcs; q++){
            non_length_graph << "L\t" << g->seq[i>>1].name << "\t"<< (i%2==0 ? "+" : "-") << "\t" << g->seq[outgoing_arcs[q].v>>1].name << "\t" << (outgoing_arcs[q].v%2==0 ? "+" : "-") << "\t0M\t" <<endl;
        }
    }
    non_length_graph.close();

    cout << "finish get bubble chain" << endl;

    return bubble_chain_begin_end_nodes;
}

map<uint32_t,map<uint32_t,set<uint32_t>>>* get_bubbles(asg_t *g, string output_directory, uint32_t** connection_forward, uint32_t** connection_backward, int detect_cycle_flag=false, int detect_branch_flag=false) {
    uint32_t **connections_count;
	CALLOC(connections_count,g->n_seq);
	for(int i = 0; i< g->n_seq; i++){
		CALLOC(connections_count[i], g->n_seq);
		memset(connections_count[i], 0, sizeof(*connections_count[i]));
	}
    for(int i = 0; i < g->n_seq; i++){
        for(int j = 0; j < g->n_seq; j++){
                connections_count[i][j] = connection_forward[i][j] + connection_backward[i][j];
                connections_count[i][j] += connection_forward[j][i] + connection_backward[j][i];
        }
    }
    cout << "start get bubbles" << endl;

	uint32_t n_vtx = g->n_seq * 2;
    int node_type[n_vtx];

    vector<bubble_t*> bubble_by_ending_begining[n_vtx];
    vector<bubble_t*> bubbles;
    for(uint32_t i=0; i<n_vtx; i++){
        node_type[i] = 0;
        bubble_t* result = detect_bubble(g, i);
        if(result!=nullptr){
            result->id = i;
            bubble_by_ending_begining[result->begNode].push_back(result);
            bubble_by_ending_begining[result->endNode].push_back(result);
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

                bubble->paths.push_back(arc_stack);
                bubble->starting_arcs.insert(arc_stack.front()->v);
                bubble->ending_arcs.insert(lastNode);
                arc_stack.pop_back();
                node_stack.pop_back();
                node_vi_stack.pop_back();
                continue;
            }

            uint32_t num_outgoing_arcs = asg_arc_n(g, u);
            if (vi < num_outgoing_arcs) {
                asg_arc_t *outgoing_arcs = asg_arc_a(g, u);
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

    set<uint32_t> bubble_chain_end_begin;
    map<uint32_t, uint32_t> pure_outgoing_num;
    map<uint32_t, uint32_t> pure_incoming_num;
    set<uint32_t> out_only;
    for(uint32_t i=0; i<n_vtx; i++){
        if(node_type[i] != 2){
            int num_outgoing_arcs = asg_arc_n(g, i);
            int num_incoming_arcs = asg_arc_n(g, i^1);
            asg_arc_t* all_outgoing_arcs = asg_arc_a(g, i);
            asg_arc_t* all_incoming_arcs = asg_arc_a(g, i^1);
            set<uint32_t> valid_outgoing_arcs, valid_incoming_arcs;
            for(int c = 0 ; c<num_outgoing_arcs; c++){
                valid_outgoing_arcs.insert(all_outgoing_arcs[c].v);
            }
            for(int c = 0 ; c<num_incoming_arcs; c++){
                valid_incoming_arcs.insert(all_incoming_arcs[c].v^1);
            }
            num_incoming_arcs = valid_incoming_arcs.size();
            num_outgoing_arcs = valid_outgoing_arcs.size();
            if(node_type[i] == 1){
                vector<bubble_t*> bubbles = bubble_by_ending_begining[i];
                set<uint32_t> outgoing_arcs_set,incoming_arcs_set;
                sort(bubbles.begin(), bubbles.end(), [ ]( const auto& lhs, const auto& rhs )
                {
                return lhs->starting_arcs.size() > rhs->starting_arcs.size();
                });
                for(int q = 0; q < bubbles.size(); q++){
                    bubble_t* bubble = bubbles[q];
                    if(bubble->begNode == i){
                        bool flag = false;
                        // cout << g->seq[bubble->begNode/2].name << " to " << g->seq[bubble->endNode/2].name << " starts: " <<bubble->starting_arcs.size() << " ends: " << bubble->ending_arcs.size() <<endl;
                        for(auto c: bubble->starting_arcs){
                            if(outgoing_arcs_set.insert(c).second){
                                num_outgoing_arcs--;
                                flag = true;
                            }
                        }
                        if(flag) num_outgoing_arcs ++;
                    }
                }
                sort(bubbles.begin(), bubbles.end(), [ ]( const auto& lhs, const auto& rhs )
                {
                return lhs->ending_arcs.size() > rhs->ending_arcs.size();
                });
                for(int q = 0; q < bubbles.size(); q++){
                    bubble_t* bubble = bubbles[q];
                    if(bubble->endNode == i){
                        bool flag = false;
                        for(auto c: bubble->ending_arcs){
                            if(incoming_arcs_set.insert(c).second){
                                num_incoming_arcs--;
                                flag = true;
                            }
                        }
                        if(flag) num_incoming_arcs ++;
                    }
                }
            }
            pure_outgoing_num[i] = num_outgoing_arcs;
            pure_incoming_num[i] = num_incoming_arcs;
            if(num_incoming_arcs==0 && num_outgoing_arcs>0){
                out_only.insert(i);
            }
            if(!(num_incoming_arcs==1 && num_outgoing_arcs==1)){
                // cout << g->seq[i/2].name << ": incoming " << num_incoming_arcs << "; outgoing " << num_outgoing_arcs <<"; bubbles:" << bubble_by_ending_begining[i].size() <<endl;
                bubble_chain_end_begin.insert(i);
            }
            //TODO:: implement

        }
    }

    
    // set<uint32_t> not_accessible;
    // for(auto beg: out_only){
    //     vector<asg_arc_t*> arc_stack;
    //     vector<uint32_t> node_stack;  // DFS
    //     vector<uint32_t> node_vi_stack;
    //     node_stack.push_back(beg);
    //     node_vi_stack.push_back(0);
    //     bool visited[g->n_seq * 2];
    //     for (int u=0; u<g->n_seq * 2; u++) {
    //         visited[u] = false;
    //     }
    //     uint32_t cur_node = beg;
    //     uint32_t linear_len = 1;
    //     while(pure_outgoing_num.find(cur_node) != pure_outgoing_num.end() && pure_outgoing_num[cur_node]!=1) {
    //         int num_outgoing_arcs = asg_arc_n(g, i);
    //         asg_arc_t* all_outgoing_arcs = asg_arc_a(g, i);

    //     }
    // }
    
    // set<string> names;
    // for(auto c : bubble_chain_end_begin){
    //     // cout << c << endl;
    //     names.insert(g->seq[c/2].name);
    // }

    // for(auto c : names){
    //     cout << c << endl;
    // }
    map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_begin_end_nodes_buf =  get_bubble_chain_graph(g,bubble_chain_end_begin,output_directory);
    for(auto a: *bubble_chain_begin_end_nodes_buf){
        for(auto b: a.second){
            for(auto bub: bubble_by_ending_begining[a.first]){
                if(bub->begNode == a.first){
                    for(auto node: bub->starting_arcs){
                        (*bubble_chain_begin_end_nodes_buf)[a.first][b.first].insert(node);
                    }
                }
            }
            for(auto bub: bubble_by_ending_begining[b.first]){
                if(bub->endNode == b.first){
                    for(auto node: bub->ending_arcs){
                        (*bubble_chain_begin_end_nodes_buf)[a.first][b.first].insert(node);
                    }
                }
            }
        }
    }
    map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_begin_end_nodes_buf2 = new map<uint32_t,map<uint32_t,set<uint32_t>>>();
    for(auto a: *bubble_chain_begin_end_nodes_buf){
        for(auto b: a.second){
            if((*bubble_chain_begin_end_nodes_buf2).find(a.first) == (*bubble_chain_begin_end_nodes_buf2).end()){
                (*bubble_chain_begin_end_nodes_buf2)[a.first] = map<uint32_t,set<uint32_t>>();
            }
            if((*bubble_chain_begin_end_nodes_buf2).find(b.first^1) == (*bubble_chain_begin_end_nodes_buf2).end()){
                (*bubble_chain_begin_end_nodes_buf2)[b.first^1] = map<uint32_t,set<uint32_t>>();
            }
            (*bubble_chain_begin_end_nodes_buf2)[a.first][b.first] = b.second;
            (*bubble_chain_begin_end_nodes_buf2)[b.first^1][a.first^1] = b.second;
        }
    }
    delete bubble_chain_begin_end_nodes_buf;
    set<uint32_t> unaccessable;
    set<uint32_t> not_filtered;
    for(auto a: *bubble_chain_begin_end_nodes_buf2){
        if(pure_outgoing_num[a.first]>1 && pure_incoming_num[a.first]>0){
            map<uint32_t,set<uint32_t>> non_outgoing;
            for(auto b: a.second){
                if(pure_outgoing_num[b.first]==0){
                    non_outgoing[b.first]=b.second;
                }
            }
            if(non_outgoing.size()==1){
                for(auto x: non_outgoing){
                    if(x.second.size() <= 8){
                        unaccessable.insert(x.first>>1);
                    }else{
                        not_filtered.insert(x.first>>1);
                    }
                }
            }else{
                uint32_t greater_count = 0;
                for(auto x: non_outgoing){
                    if(x.second.size() > 8){
                        greater_count++;
                    }
                }
                if(greater_count>0){
                    for(auto x: non_outgoing){
                        if(x.second.size() <= 8){
                            unaccessable.insert(x.first>>1);
                        }else{
                            not_filtered.insert(x.first>>1);
                        }
                    }   
                }else{
                    uint32_t max_idx = -1;
                    uint32_t max_count = 0;
                    for(auto x:non_outgoing){
                        uint32_t cur_len = 0;
                        for(auto i : x.second){
                            cur_len+=g->seq[i>>1].len;
                        }
                        if(cur_len >= max_count){
                            max_count = cur_len;
                            max_idx = x.first;
                        }
                    }
                    for(auto x:non_outgoing){
                        if(x.first!=max_idx){
                            unaccessable.insert(x.first>>1);
                        }else{
                            not_filtered.insert(x.first>>1);
                        }
                    }
                }
            }
        }
    }
    map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_begin_end_nodes_buf3 = new map<uint32_t,map<uint32_t,set<uint32_t>>>();
    map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_end_begin_nodes_buf3 = new map<uint32_t,map<uint32_t,set<uint32_t>>>();

    for(auto a: *bubble_chain_begin_end_nodes_buf2){
        if(unaccessable.find(a.first>>1)==unaccessable.end()){
            for(auto b: a.second){
                if(unaccessable.find(b.first>>1)==unaccessable.end() && b.second.size()>1){
                    if((*bubble_chain_begin_end_nodes_buf3).find(a.first) == (*bubble_chain_begin_end_nodes_buf3).end()){
                        (*bubble_chain_begin_end_nodes_buf3)[a.first] = map<uint32_t,set<uint32_t>>();
                    }
                    if((*bubble_chain_begin_end_nodes_buf3).find(b.first^1) == (*bubble_chain_begin_end_nodes_buf3).end()){
                        (*bubble_chain_begin_end_nodes_buf3)[b.first^1] = map<uint32_t,set<uint32_t>>();
                    }
                    (*bubble_chain_begin_end_nodes_buf3)[a.first][b.first] = b.second;
                    (*bubble_chain_begin_end_nodes_buf3)[b.first^1][a.first^1] = b.second;
                    if((*bubble_chain_end_begin_nodes_buf3).find(b.first) == (*bubble_chain_end_begin_nodes_buf3).end()){
                        (*bubble_chain_end_begin_nodes_buf3)[b.first] = map<uint32_t,set<uint32_t>>();
                    }
                    if((*bubble_chain_end_begin_nodes_buf3).find(a.first^1) == (*bubble_chain_end_begin_nodes_buf3).end()){
                        (*bubble_chain_end_begin_nodes_buf3)[a.first^1] = map<uint32_t,set<uint32_t>>();
                    }
                    (*bubble_chain_end_begin_nodes_buf3)[b.first][a.first] = b.second;
                    (*bubble_chain_end_begin_nodes_buf3)[a.first^1][b.first^1] = b.second;
                }
            }
        }
    }
    delete bubble_chain_begin_end_nodes_buf2;
    set<vector<uint32_t>> extensions;
    cout << "Start connect single branches" << endl;
    for(auto beg: (*bubble_chain_begin_end_nodes_buf3)){
        for(auto end: beg.second){
            vector<uint32_t> to_expand;
            set<uint32_t> visited;
            visited.insert(beg.first);
            visited.insert(end.first);
            to_expand.push_back(beg.first);
            to_expand.push_back(end.first);
            while((*bubble_chain_begin_end_nodes_buf3).find(to_expand[to_expand.size()-1])!=(*bubble_chain_begin_end_nodes_buf3).end() 
            && (*bubble_chain_begin_end_nodes_buf3)[to_expand[to_expand.size()-1]].size()==1
            && (*bubble_chain_end_begin_nodes_buf3).find(to_expand[to_expand.size()-1])!=(*bubble_chain_end_begin_nodes_buf3).end() 
            && (*bubble_chain_end_begin_nodes_buf3)[to_expand[to_expand.size()-1]].size()==1){
                bool seen = false;;
                for(auto i: (*bubble_chain_begin_end_nodes_buf3)[to_expand[to_expand.size()-1]]){
                    if(visited.find(i.first)!=visited.end()){
                        seen = true;
                        break;
                    }
                    // cout << i.first << endl;
                    to_expand.push_back(i.first);
                    visited.insert(i.first);
                }
                if(seen){
                    break;
                }
            }
            extensions.insert(to_expand);
        }
    }
    cout << "Start filter short connections:\t" << extensions.size() << endl;
    vector<vector<uint32_t>> pathes;
    for(auto path: extensions){
        set<uint32_t> path_nodes;
        for(auto n: path){
            path_nodes.insert(n);
        }
        bool valid = true;
        for(auto to_check: extensions){
            set<uint32_t> to_check_nodes;
            for(auto n: to_check){
                to_check_nodes.insert(n);
            }
            bool found = true;
            for(auto n: path_nodes){
                if(to_check_nodes.find(n)==to_check_nodes.end()){
                    found = false;
                    break;
                }
            }
            if(found && to_check_nodes.size()>path_nodes.size()){
                valid = false;
                break;
            }
        }
        if(valid){
            pathes.push_back(path);
        }
    }
    cout << "End filter short connections:\t" << pathes.size() << endl;
    map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_begin_end_nodes_buf4 = new map<uint32_t,map<uint32_t,set<uint32_t>>>();

    for(auto path: pathes){
        set<uint32_t> to_insert;
        for(int i = 0; i < path.size()-1; i++){
            to_insert.insert((*bubble_chain_begin_end_nodes_buf3)[path[i]][path[i+1]].begin(),(*bubble_chain_begin_end_nodes_buf3)[path[i]][path[i+1]].end());
        }
        if(to_insert.size()>=2){
            uint32_t beg = path[0];
            uint32_t end = path[path.size()-1];
            if((*bubble_chain_begin_end_nodes_buf4).find(beg) == (*bubble_chain_begin_end_nodes_buf4).end()){
                (*bubble_chain_begin_end_nodes_buf4)[beg] = map<uint32_t,set<uint32_t>>();
            }
            if((*bubble_chain_begin_end_nodes_buf4).find(end^1) == (*bubble_chain_begin_end_nodes_buf4).end()){
                (*bubble_chain_begin_end_nodes_buf4)[end^1] = map<uint32_t,set<uint32_t>>();
            }
            (*bubble_chain_begin_end_nodes_buf4)[beg][end] = to_insert;
            (*bubble_chain_begin_end_nodes_buf4)[end^1][beg^1] = to_insert;
        }

    }
    delete bubble_chain_begin_end_nodes_buf3;
    map<uint32_t,map<uint32_t,set<uint32_t>>>* bubble_chain_begin_end_nodes = new map<uint32_t,map<uint32_t,set<uint32_t>>>();
    
    for(auto beg: (*bubble_chain_begin_end_nodes_buf4)){
        (*bubble_chain_begin_end_nodes)[beg.first] = map<uint32_t,set<uint32_t>>();
        for(auto end: beg.second){
            set<uint32_t> to_insert;
            for(auto n: end.second){
                to_insert.insert(n>>1);
            }
            uint64_t all_len = 0;
            // if(to_insert.size()>8 || (pure_incoming_num[beg.first]==0 && pure_outgoing_num[end.first] == 0)){
                for(auto n: to_insert){
                    all_len += g->seq[n].len;
                }
            // }
            // if(to_insert.size()>8 || all_len > 10000000 || (pure_incoming_num[beg.first]==0 && pure_outgoing_num[end.first] == 0) ){
            // if(to_insert.size()>8 || all_len > 10000000 || (pure_incoming_num[beg.first]==0 && pure_outgoing_num[end.first] == 0) ){
            // if(all_len > 10000000 ){
            if(all_len > 1000000 ){
                (*bubble_chain_begin_end_nodes)[beg.first][end.first] = to_insert;
            }
        }
    }
    delete bubble_chain_begin_end_nodes_buf4;
    cout << "Filtered Nodes: " << endl;
    for(auto n: unaccessable){
        cout << g->seq[n].name << ", ";
    }
    cout << endl;


    cout << "Not Filtered Nodes: " << endl;
    for(auto n: not_filtered){
        cout << g->seq[n].name << ", ";
    }
    cout << endl;

    // for(auto a: *bubble_chain_begin_end_nodes){
    //     for(auto b: a.second){
    //         // all_nodes.insert(b.second.begin(),b.second.end());
    //         cout << g->seq[a.first>>1].name << ", " << g->seq[b.first>>1].name << endl;
    //     }
    // }

    uint32_t branch_counter = 0;
    set<uint32_t> all_nodes;
    for(auto a: *bubble_chain_begin_end_nodes){
        for(auto b: a.second){
            branch_counter++;
            all_nodes.insert(b.second.begin(),b.second.end());
        }
    }
    cout << "Nodes in original Graph: " << g->n_seq << endl;
    cout << "Total Nodes in Chain Graph: " << all_nodes.size() << endl;
    cout << "Branches: " << branch_counter << endl;
    

    for(auto bubble: bubbles){
        delete bubble;
    }

    for(int i = 0; i< g->n_seq; i++){
		free(connections_count[i]);
	}
    free(connections_count);
    return bubble_chain_begin_end_nodes;
    // cout << bubble_chain_end_begin.size() << endl;
    // cout << names.size() << endl;
    // cout << "finish get bubble chain" << endl;
}



#endif
