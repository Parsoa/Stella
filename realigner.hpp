#ifndef REALIGNER_HPP
#define REALIGNER_HPP 

#include <omp.h>
#include <vector>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include <pthread.h>
#include <unordered_map>

#include "rle.h"
#include "rld0.h"
#include "mrope.h"

#include <zlib.h>

#include "bam.hpp"
#include "fastq.hpp"
#include "lprint.hpp"
#include "config.hpp"
#include "sfsutils.hpp"

using namespace std;

typedef std::pair<std::string, std::string> batch_atom_type ;

class Realigner {

public:
    void run();

private:

    int current_batch = 0 ;
    int last_dumped_batch = 0 ;

    samFile *bam_file ;
    bam_hdr_t *bam_header ;
    std::vector<std::vector<std::vector<bam1_t*>>> bam_entries ;
    bool load_batch_bam(int threads, int batch_size, int p) ;
    void load_input_sfs_batch() ;
    ofstream out_file ;
    ofstream tau_out_file ;

    void load_target_SFS_set() ;
    std::unordered_map<std::string, bool> target_sfs ;
    std::unordered_map<std::string, bool> target_sfs_hits ;
    
    std::unordered_map<std::string, std::vector<SFS>> SFSs ;
    std::vector<batch_atom_type> process_batch(int, int) ;
    std::vector<std::vector<std::vector<batch_atom_type>>> batches ;
    void output_batch(int) ;

    Configuration* config ;

    CIGAR rebuild_cigar(char*, char*, const std::vector<std::pair<int, int>> &alpairs) ;
};

#endif
