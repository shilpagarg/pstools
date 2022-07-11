// #include "zran.h"
// #include "graph.h"
#include <iostream>
#include <vector>
#include <zlib.h>
#include <stdio.h>
#include <string.h>
// #include "paf.h"
#include <map>
// #include "bubble_chain.h"
#include <queue>
#include <assert.h>
#include <sstream> 

using namespace std;


class sketch_minimizer{
    public:
        void get_minimizer();
};

class seq_mini_t{
    public:
        vector<vector<uint32_t>> pos;
        vector<string> minimizers;
};

seq_mini_t* get_minimizer(string sequence,int k, int m){
    seq_mini_t* result = new seq_mini_t();
    vector<uint32_t> buf;
    buf.push_back(0);
    result->pos.push_back(buf);
    string current_mini(m,'Z');
    string pre_mini;
    int endPos = sequence.size() - k - m;
    for(int i = 0; i < endPos; i++){
        if(i==0 || (sequence.substr(i-1,m).compare(pre_mini)==0 && sequence.substr(i+k-1,m).compare(pre_mini)>0)){
            current_mini = string(m, 'Z');
            for(int q = i; q < i+k; q++){
                if(sequence.substr(q,m).compare(current_mini) < 0){
                    current_mini = sequence.substr(q,m);
                }
            }
        }else if(sequence.substr(i+k-1,m).compare(pre_mini)<=0){
            current_mini = sequence.substr(i+k-1,m);
        }
        if(i!=0 && pre_mini.compare(current_mini)!=0){
            result->pos.back().push_back(i+k);
            vector<uint32_t> current_start;
            current_start.push_back(i+k);
            result->pos.push_back(current_start);
            result->minimizers.push_back(current_mini);
        }
        pre_mini = current_mini;
    }
    if(result->pos.size()!=result->minimizers.size()){
        result->pos.back().push_back(sequence.size());
        result->minimizers.push_back(current_mini);
    }
    return result;
}


seq_mini_t* sub_minimizer(seq_mini_t* seq_mini, uint32_t start, uint32_t end){
    seq_mini_t* result = new seq_mini_t();
    bool contain = false;
    for(uint i = 0; i<seq_mini->minimizers.size(); i++){
        if(!contain && start>=seq_mini->pos[i][0]){
            contain = true;
            result->pos.push_back(seq_mini->pos[i]);
            result->pos.back()[0] = start;
            result->minimizers.push_back(seq_mini->minimizers[i]);
        }else if(contain && seq_mini->pos[i][0]<=end && end<=seq_mini->pos[i][1] ){
            result->pos.push_back(seq_mini->pos[i]);
            result->pos.back()[1] = end;
            result->minimizers.push_back(seq_mini->minimizers[i]);
            break;
        }else if(contain){
            result->pos.push_back(seq_mini->pos[i]);
            result->minimizers.push_back(seq_mini->minimizers[i]);
        }
    }
    return result;
}



map<uint32_t, seq_mini_t*> get_repeat_variant(map<uint32_t, seq_mini_t*> seq_minis){
    map<uint32_t, seq_mini_t*> result;
    map<uint32_t, uint32_t> currentPos;
    map<uint32_t, string> currentMinimizer;
    map<uint32_t, string> prevMinimizer;
    map<uint32_t, uint32_t> start;
    uint32_t size = UINT32_MAX;
    bool first = true;
    // uint32_t prev_start = 0;
    for(auto x: seq_minis){
        if(x.second->pos.back().back() - x.second->pos.front().front() < size){
            size = x.second->pos.back().back() - x.second->pos.front().front();
        }
        start[x.first] = x.second->pos.front().front();
        result[x.first] = new seq_mini_t();
        currentPos[x.first] = 0;
        currentMinimizer[x.first] = string('Z',7);
        prevMinimizer[x.first] = string('Z',7);
    }
    string pre_minimizer = string('Z',7);
    for(uint32_t i = 0; i < size; i++){
        for(auto j: currentPos){
            if((seq_minis[j.first])->pos[j.second][1] <= start[j.first] + i){
                currentPos[j.first]++;
            }
        }
        set<string> minimizers;
        bool non_equal = true;
        for(auto x : currentPos){
            currentMinimizer[x.first] = seq_minis[x.first]->minimizers[x.second];
            if(!minimizers.insert(currentMinimizer[x.first]).second){
                non_equal = false;
            }
        }

        if(non_equal){
            if(first){
                for(auto x: result){
                    result[x.first]->pos.push_back(vector<uint32_t>());
                    result[x.first]->pos.back().push_back(i);
                    result[x.first]->minimizers.push_back(currentMinimizer[x.first]);
                }
                first = false;
            }else{
                bool equal_prev = true;
                for(auto buf: prevMinimizer){
                    equal_prev = (currentMinimizer[buf.first].compare(buf.second)==0);
                    if(!equal_prev){
                        break;
                    }
                }
                if(!equal_prev){
                    for(auto x: result){
                        result[x.first]->pos.back().push_back(i);
                        result[x.first]->pos.push_back(vector<uint32_t>());
                        result[x.first]->pos.back().push_back(i);
                        result[x.first]->minimizers.push_back(currentMinimizer[x.first]);
                    }
                }
            }
        }else{
            if(!first){
                for(auto x: result){
                    result[x.first]->pos.back().push_back(i);
                }
            first = true;
            }
        }
        prevMinimizer = currentMinimizer;
    }
    for(auto x: result){
        if(x.second->pos.back().size()==1){
            x.second->pos.back().push_back(size);
        }
    }
    return result;
}


map<uint32_t, uint64_t> match_minimizers(seq_mini_t* read_mini, map<uint32_t, seq_mini_t*> node_minis){
    map<uint32_t, uint64_t> result;
    map<string, uint64_t> read_mini_bin;
    for(uint i = 0; i<read_mini->pos.size(); i++){
        if(read_mini_bin.find(read_mini->minimizers[i])==read_mini_bin.end()){
            read_mini_bin[read_mini->minimizers[i]] = 0;
        }
        read_mini_bin[read_mini->minimizers[i]] += read_mini->pos[i][1] - read_mini->pos[i][0];
    }
    for(auto node_mini: node_minis){
        result[node_mini.first] = 0;
        map<string, uint64_t> buf_read_mini_bin = read_mini_bin;
        map<string, uint64_t> node_mini_bin;
        for(uint i = 0; i<node_mini.second->pos.size(); i++){
            if(node_mini_bin.find(node_mini.second->minimizers[i])==node_mini_bin.end()){
                node_mini_bin[node_mini.second->minimizers[i]] = 0;
            }
            node_mini_bin[node_mini.second->minimizers[i]] += node_mini.second->pos[i][1] - node_mini.second->pos[i][0];
        }
        for(auto pair: node_mini_bin){
            if(buf_read_mini_bin.find(pair.first)!=buf_read_mini_bin.end()){
                if(pair.second > buf_read_mini_bin[pair.first]){
                    result[node_mini.first] += buf_read_mini_bin[pair.first];
                    // buf_read_mini_bin[pair.first] = 0;
                }else{
                    result[node_mini.first] += pair.second;
                    // buf_read_mini_bin[pair.first] -= pair.second;
                }
            }
        }
    }
    return result;
}




map<string, double> profile_minimizer(seq_mini_t* to_profile){
    map<string, double> count;
    uint64_t total_num = 0;
    for(uint i = 0; i< to_profile->minimizers.size(); i++){
        if(count.find(to_profile->minimizers[i]) == count.end()){
            count[to_profile->minimizers[i]] = 0;
        }
        count[to_profile->minimizers[i]] += (to_profile->pos[i][1] - to_profile->pos[i][0]);
        total_num += (to_profile->pos[i][1] - to_profile->pos[i][0]);
    }
    for(auto i: count){
        count[i.first] = count[i.first] / (double)total_num * 100;
    }
    return count;
}

double minimizer_match(map<string, double> mini_seq_1, map<string, double> mini_seq_2){
    double result = 0;
    for(auto x: mini_seq_1){
        if(mini_seq_2.find(x.first)!=mini_seq_2.end()){
            result += mini_seq_2[x.first] * x.second;
        }
    }
    return result;
}





map<uint32_t, uint64_t> minimizers_match(seq_mini_t* read_mini, map<uint32_t, seq_mini_t*> node_minis){
    map<uint32_t, uint64_t> result;
    map<string, uint64_t> read_mini_bin;
    for(uint i = 0; i<read_mini->pos.size(); i++){
        if(read_mini_bin.find(read_mini->minimizers[i])==read_mini_bin.end()){
            read_mini_bin[read_mini->minimizers[i]] = 0;
        }
        read_mini_bin[read_mini->minimizers[i]] += read_mini->pos[i][1] - read_mini->pos[i][0];
    }
    for(auto node_mini: node_minis){
        result[node_mini.first] = 0;
        map<string, uint64_t> buf_read_mini_bin = read_mini_bin;
        map<string, uint64_t> node_mini_bin;
        for(uint i = 0; i<node_mini.second->pos.size(); i++){
            if(node_mini_bin.find(node_mini.second->minimizers[i])==node_mini_bin.end()){
                node_mini_bin[node_mini.second->minimizers[i]] = 0;
            }
            node_mini_bin[node_mini.second->minimizers[i]] += node_mini.second->pos[i][1] - node_mini.second->pos[i][0];
        }
        for(auto pair: node_mini_bin){
            if(buf_read_mini_bin.find(pair.first)!=buf_read_mini_bin.end()){
                if(pair.second > buf_read_mini_bin[pair.first]){
                    result[node_mini.first] += buf_read_mini_bin[pair.first];
                    // buf_read_mini_bin[pair.first] = 0;
                }else{
                    result[node_mini.first] += pair.second;
                    // buf_read_mini_bin[pair.first] -= pair.second;
                }
            }
        }
    }
    return result;
}

