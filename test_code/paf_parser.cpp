#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <future>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <experimental/filesystem>




using namespace std;

mutex mtx;

void readfile(string  filename, unordered_map<string, vector<pair<int,string> > >& target_map, int num_threads,unsigned long long int block_len, int index)
{
    // cout<< "thread" << index << endl;
    string target_seq_name, query_seq_name, temp,line;
    ifstream infile;
    infile.open(filename);

    int query_start;

    unsigned long long int read_len=0;

    unsigned long long int start=block_len*index;

    // cout <<"Start: " << start << endl;
    infile.seekg(start, ios::beg);


    char ch;

    if(start!=0)
    {
        infile.seekg(-2, infile.cur);
        infile >> noskipws >>ch;
        if(ch !='\n')
        {
          while(infile.peek() != '\n')
          {
            infile >> noskipws >>ch;
            // cout << ch;
            read_len++;
          }
          infile >> noskipws >>ch;
          read_len++; //for above \n
        }

    }


    while ((read_len<=block_len) && getline(infile, line))
    {
        read_len+=line.length()+1; //+1 for \n
        // cout << "Thread"<< index << "     read length: " << read_len << endl;

        istringstream record(line);
        if (!(record >> query_seq_name >> temp >> query_start >> temp >> temp >> target_seq_name))
        {
          cout<< "Indexing issue for chunks while reading" << endl;
          break;
        }

        mtx.lock();
        // cout << "thread" << index << " adds " << query_seq_name << " at " << query_start << " to " << target_seq_name << endl;
        // cout << "Read length: "<<read_len << endl;
        target_map[target_seq_name].push_back({query_start,query_seq_name});
        mtx.unlock();

    }


    infile.close();
}


void sort_query(unordered_map<string, vector<pair<int,string> > >& target_map, int index, unsigned int chunk_size, int num_threads)
{
    unordered_map<string, vector<pair<int,string> > >::iterator ptr = target_map.begin();

    advance(ptr,index*chunk_size);
    unsigned int i=0;

    while(i<chunk_size)
    {
        //cout << ptr->first << endl;
        // cout<< "Sorted " << ptr->first << endl;
        sort(ptr->second.begin(), ptr->second.end());
        i++;
        ptr++;
    }

    if(index==num_threads-1)
    {
        while(ptr!=target_map.end())
        {
            sort(ptr->second.begin(), ptr->second.end());
            ptr++;
        }
    }

}




int main(int argc, const char * argv[])
{
   string line;
   if(argc<4)
   {
    cout<< "Input & output filenames and number of threads needed" << endl;
    exit(1);
   }
   string filename=argv[1];
   int num_threads=stoi(argv[3]);

   size_t filesize = experimental::filesystem::file_size(filename);

   cout<< "Filesize: " << filesize << " bytes" << endl;

   int i;

   unordered_map<string, vector<pair<int,string> > > target_map;

    // future<unordered_map<string, vector<pair<int,string> > > > target_map_future = async(launch::async, readfile, filename);

   vector<thread> read_threads;

   unsigned long long int block_len=filesize/num_threads;

   cout<< "Block length: " << block_len << endl;

   auto start_read = chrono::high_resolution_clock::now();


   for(i=0;i<num_threads;i++)
   {

    read_threads.push_back(thread(readfile,filename,ref(target_map),num_threads,block_len,i));
   }

   for (thread &t : read_threads)
    {
        if (t.joinable())
        {
            t.join();
        }

    }

    auto end_read=chrono::high_resolution_clock::now();
    double time_taken_read=chrono::duration_cast<chrono::nanoseconds>(end_read - start_read).count();
    time_taken_read *= 1e-9;


    unsigned int target_map_size=target_map.size();

    unsigned int chunk_size=(target_map_size > num_threads) ? target_map_size/num_threads: 1;

    int num_threads_sort=(target_map_size > num_threads) ? num_threads : target_map_size;


    vector<thread> sort_threads;

    auto start_sort = chrono::high_resolution_clock::now();



    for(i=0;i<num_threads_sort;i++)
    {
        sort_threads.push_back(thread(sort_query,ref(target_map),i,chunk_size,num_threads_sort));
    }

   for (thread &t : sort_threads)
    {
        if (t.joinable())
        {
            t.join();
        }

    }


    auto end_sort=chrono::high_resolution_clock::now();
    double time_taken_sort=chrono::duration_cast<chrono::nanoseconds>(end_sort - start_sort).count();
    time_taken_sort *= 1e-9;


    ofstream outfile;

    string outfilename=argv[2];

    outfile.open(outfilename);


    int num_records=0;
    for (auto x : target_map)
    {
        outfile << x.first  << ": "<<endl;
        for ( auto y : x.second)
        {
            outfile << y.first << " ";
            outfile  << y.second << " ";
            num_records++;
        }
        outfile << endl;
    }
    outfile.close();

    cout << "Number of records: " << num_records << endl;
    cout << "Number of target sequences: " << target_map.size() << endl;

    cout << "Time taken by read is : " << fixed << time_taken_read << setprecision(9) << " sec" << endl;
    cout << "Time taken by sort is : " << fixed << time_taken_sort << setprecision(9) << " sec" << endl;

    return 0;
}
