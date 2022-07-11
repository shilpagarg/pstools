// Headers I definitely need
#include "centro_asm.hpp"
#include "bseq.h" // for reading condensed files
extern "C" {
    #include "ClassPro.h" 
	#include "FastK.h" 
	#include "libfastk.h" 
}
#include <fstream>
#include <iostream>
#include <vector>
#include <string.h>
#include <algorithm>
#include <set>

using namespace std;

int main_centro_asm(int argc, char* argv[]){
	// std::cout << "test\n";

	// // FastK
	// int argc_fk = 6;
	// char* argv_fk[argc_fk];
	// argv_fk[0] = "FastK";
	// argv_fk[2] = "-t1";
	// argv_fk[3] = "-p";
	// argv_fk[4] = "-v";
	// argv_fk[5] = "-k20";
	// argv_fk[1] = "/binf-isilon/sgarglab/people/graphSVTADs/analysis/test_data/subset.chr13.fastq";
	// main_fastk(argc_fk, argv_fk);
	// std::cout << "fastk done\n";

	// // ClassPro
	// int argc_cp = 2;
	// char* argv_cp[argc_cp];
	// argv_cp[0] = "ClassPro";
	// //argv_cp[1] = "-v";
	// argv_cp[1] = "/binf-isilon/sgarglab/people/graphSVTADs/centro/HG01891.chr5.fastq";
	// main_classpro(argc_cp, argv_cp);
	// std::cout << "classpro done\n";
	

	// // FastK Create kmertable of informative k-mers
	// int argc_fk = 5;
	// char* argv_fk[argc_fk];
	// argv_fk[0] = "FastK";
	// argv_fk[1] = "-t";
	// argv_fk[2] = "-v";
	// // argv_fk[3] = "-k20";
	// argv_fk[3] = "-T1";
	// argv_fk[4] = "/binf-isilon/sgarglab/people/graphSVTADs/centro/HG01891.chr5_inf_full.fastq";
	// main_fastk(argc_fk, argv_fk);
	// std::cout << "fastk informative table done\n";

	// int argc_fk = 5;
	// char* argv_fk[argc_fk];
	// argv_fk[0] = "FastK";
	// argv_fk[1] = "-p:/binf-isilon/sgarglab/people/graphSVTADs/centro/HG01891.chr5_inf_full";
	// argv_fk[2] = "-v";
	// argv_fk[3] = "-T1";
	// argv_fk[4] = "/binf-isilon/sgarglab/people/graphSVTADs/centro/HG01891.nano.chr5.fastq";
	// main_fastk(argc_fk, argv_fk);
	// std::cout << "fastk informative histogram done\n";

	// Display histogram
	int argc_pf = 3;
	// int counter = argc_pf;
	// char* argv_pf[argc_pf];
	// argv_pf[--counter] = "1-#";
	// argv_pf[--counter] = "/binf-isilon/sgarglab/people/graphSVTADs/centro/HG01891.nano.chr5";
	// argv_pf[--counter] = "Profex";

	const char *argv_pf[argc_pf] = { "Profex", "/binf-isilon/sgarglab/people/graphSVTADs/centro/HG01891.nano.chr5", "1-#"};
	std::cout << argv_pf[0] << argv_pf[1] << argv_pf[2] << "\n";
	main_prof(argc_pf, argv_pf);
	std::cout << "Displaying fastk informative histogram done\n";
	return 0;
}