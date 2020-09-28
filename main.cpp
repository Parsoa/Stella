#include <omp.h>
#include <ctime>
#include <chrono>
#include <string>
#include <vector>
#include <thread>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include <pthread.h>
#include <unordered_map>

#include <zlib.h>

#include "rle.h"
#include "rld0.h"
#include "mrope.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

    using namespace std;

#ifdef DEBUG_MODE
#  define DEBUG(x) x
#  define NEBUG(x)
#else
#  define DEBUG(x)
#  define NEBUG(x) x
#endif

/** From fermi ***********/

#define fm6_comp(a) ((a) >= 1 && (a) <= 4? 5 - (a) : (a))

#define fm6_set_intv(e, c, ik) ((ik).x[0] = (e)->cnt[(int)(c)], (ik).x[2] = (e)->cnt[(int)(c)+1] - (e)->cnt[(int)(c)], (ik).x[1] = (e)->cnt[fm6_comp(c)], (ik).info = 0)

static unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
} ;

void seq_char2nt6(int l, unsigned char *s) {
    int i;
    for (i = 0; i < l; ++i)
        s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5;
}

/*************************/

static const vector<string> int2char ({"$", "A", "C", "G", "T", "N"});

struct fastq_entry_t {
    string head ;
    string seq ;
    string qual ;
    uint start ;
    uint len ;
    // constructor
    fastq_entry_t(const string &h, const string &s, const string &q, const uint st = 0, const uint l = 0) : head(h), seq(s), qual(q) {
        len = l ;
        start = st ;
    }
    bool operator==(const fastq_entry_t &o) const {
        return seq == o.seq ; 
    }
};

namespace std {
    template <> struct hash<fastq_entry_t> {
        std::size_t operator()(const fastq_entry_t& k) const {
            return std::hash<string>()(k.seq) ;
        }
    };
}

fastq_entry_t get_solution(fastq_entry_t fqe, int s, int l) {
    string S (fqe.seq, s, l) ;
    string Q (fqe.qual, s, l) ;
    return fastq_entry_t(fqe.head, S, Q, s, l) ;
}

string interval2str(rldintv_t sai) {
    return "[" + to_string(sai.x[0]) + "," + to_string(sai.x[1]) + "," + to_string(sai.x[2]) + "]";
}

bool backward_search(rld_t *index, const uint8_t *P, int p2) {
    rldintv_t sai ; // rldintv_t is the struct used to store a SA interval.
    fm6_set_intv(index, P[p2], sai) ;
    while(sai.x[2] != 0 && p2 > 0) {
        --p2;
        rldintv_t osai[6] ;
        rld_extend(index, &sai, osai, 1) ; //1: backward, 0: forward
        sai = osai[P[p2]] ;
    }
    return sai.x[2] != 0 ;
}
bool check_solution(rld_t* index, string S) {
    int l = S.length() ;
    uint8_t *P = (uint8_t*) S.c_str() ;
    seq_char2nt6(l, P);
    bool found_full = backward_search(index, P, l - 1) ;
    bool found_prefix = backward_search(index, P, l - 2) ;
    bool found_suffix = backward_search(index, P + 1, l - 2) ;
    //if (found_full) {
    //    cerr << "Found full." << endl ;
    //}
    //if (found_prefix) {
    //    cerr << "Found prefix." << endl ;
    //}
    //if (found_suffix) {
    //    cerr << "Found suffix." << endl ;
    //}
    return !found_full & (found_prefix || found_suffix) ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void ping_pong_search(rld_t *index, fastq_entry_t fqe, vector<fastq_entry_t>& solutions) {
    int l = fqe.seq.size() ;
    if (l == 0) {
        cerr << "Empty read." << endl ;
        return ;
    }
    DEBUG(cerr << "Read Length: " << l << endl ;)

    char *seq = new char[fqe.seq.size() + 1] ; // current sequence
    strcpy(seq, fqe.seq.c_str()) ; // seq
    uint8_t *P = (uint8_t*) seq ; 
    seq_char2nt6(l, P) ; // convert to integers
    rldintv_t sai ;

    int begin = l - 1 ;
    fm6_set_intv(index, P[begin], sai) ;
    while (begin >= 0) {
        // Backward search. Find a mismatching sequence. Stop at first mismatch.
        int bmatches = 0 ;
        DEBUG(cerr << "BS from " << int2char[P[begin]] << " (" << begin << "): " << interval2str(sai) << endl ;)
            bmatches = 0 ;
        while (sai.x[2] != 0 && begin > 0) {
            begin-- ;
            bmatches++ ;
            rldintv_t osai[6] ; // output SA intervals (one for each symbol between 0 and 5)
            rld_extend(index, &sai, osai, 1) ;
            sai = osai[P[begin]] ;
            DEBUG(cerr << "- BE with " << int2char[P[begin]] << " (" << begin << "): " << interval2str(sai) << endl ;)
        }
        //last sequence was a match
        if (begin == 0 && sai.x[2] != 0) {
            break ;
        }
        DEBUG(cerr << "Mismatch " << int2char[P[begin]] << " (" <<  begin << "). bmatches: " << to_string(bmatches) << endl ;)
            // Forward search: 
            int end = begin ;
        int fmatches = 0 ;
        fm6_set_intv(index, P[end], sai) ;
	rldintv_t sai_backup;
        DEBUG(cerr << "FS from " << int2char[P[end]] << " (" << end << "): " << interval2str(sai) << endl ;)
            while(sai.x[2] != 0) {
                sai_backup = sai;
                end++ ;
                fmatches++ ;
                rldintv_t osai[6] ;
                rld_extend(index, &sai, osai, 0) ;
                sai = osai[P[end] >= 1 && P[end] <= 4 ? 5 - P[end] : P[end]];
                DEBUG(cerr << "- FE with " << int2char[P[end]] << " (" <<  end << "): " << interval2str(sai) << endl ;)
            }
        DEBUG(cerr << "Mismatch " << int2char[P[end]] << " (" << end << "). fmatches: " << fmatches << endl ;)
        // add solution
        DEBUG(cerr << "Adding [" << begin << ", " << end << "]." << endl ;)
        //string S(fqe.seq, begin, end - begin + 1) ;
        //if (check_solution(index, S)) {
        //    DEBUG(cerr << "Found " << S << endl ;)
        //    exit(0) ;
        //}
        solutions.push_back(get_solution(fqe, begin, end - begin + 1)) ;
        //if (solutions.find(S) == solutions.end()) {
        //    solutions[S] = 0 ;
        //}
        //solutions[S] = solutions[S] + 1 ;
        //DEBUG(std::this_thread::sleep_for(std::chrono::seconds(1)) ;)
        // prepare for next round
        if (begin == 0) {
            break ;
        }
        // overlapping version:
	sai = sai_backup;
	
	// non-overlapping version:
	// --begin;
	// fm6_set_intv(index, P[begin], sai) ;
    }
    //DEBUG(std::this_thread::sleep_for(std::chrono::seconds(2)) ;)
    delete[] seq ;
}

gzFile fastq_file ;
kseq_t* fastq_iterator ;
vector<vector<vector<fastq_entry_t>>> fastq_entries ;
vector<unordered_map<fastq_entry_t, int>> search_solutions ;

bool load_batch_fastq(int threads, int batch_size, int p) {
    for (int i = 0; i < threads; i++) {
        fastq_entries[p][i].clear() ;
    }
    int l = 0 ;
    int n = 0 ;
    while ((l = kseq_read(fastq_iterator)) >= 0) {
        fastq_entries[p][n % threads].push_back(fastq_entry_t(fastq_iterator->name.s, fastq_iterator->seq.s, fastq_iterator->qual.s)) ;
        n += 1 ;
        if (n == batch_size) {
            return true ;
        }
    }
    if(fastq_entries[p].empty()) {
        return false ;
    }
    return true ;
}

vector<fastq_entry_t> process_batch_fastq(rld_t* index, vector<fastq_entry_t> fastq_entries) {
    vector<fastq_entry_t> solutions ;
    for (const auto fastq_entry : fastq_entries) {
        ping_pong_search(index, fastq_entry, solutions) ;
    }
    return solutions ;
}

int current_batch = 0 ;
struct OutputBatchArgs {
    int batch ;
} ;

void output_batch(void* args) {
    int batch = ((OutputBatchArgs*) args)->batch ;
    string path = "solution_batch_" + std::to_string(batch - 1) + ".fastq" ;
    std::ofstream o(path) ;
    for (const auto it : search_solutions[batch - 1]) {
        fastq_entry_t fastq_entry = it.first ;
        o << "@" << fastq_entry.head << ".css" << "_" << fastq_entry.start << ":"
            << fastq_entry.start + fastq_entry.len - 1 << ":" << it.second << endl
            << fastq_entry.seq << endl
            << "+" << endl
            << fastq_entry.qual << endl ;
    }
    search_solutions[batch - 1].clear() ;
    //pthread_exit(NULL) ;
}

int search_f3(int argc, char *argv[]) {
    // parse arguments
    char *index_path = argv[1] ;
    cout << "Restroing index.." << endl ;
    rld_t *index = rld_restore(index_path) ;
    cout << "Done." << endl ;
    char *sample_path = argv[2] ;
    fastq_file = gzopen(sample_path, "r") ;
    fastq_iterator = kseq_init(fastq_file) ;
    int threads = atoi(argv[3]) ;
    // load first batch
    unordered_map<fastq_entry_t, int> s ;
    search_solutions.push_back(s) ;
    vector<vector<vector<fastq_entry_t>>> batches ;
    for(int i = 0; i < 2; i++) {
        batches.push_back(vector<vector<fastq_entry_t>>(threads)) ; // previous and current output
        fastq_entries.push_back(vector<vector<fastq_entry_t>>(threads)) ; // current and next output
    }
    int p = 0 ;
    int batch_size = 10000 ;
    cerr << "Loading first batch" << endl ;
    load_batch_fastq(threads, batch_size, p) ;
    // main loop
    time_t t ;
    time(&t) ;
    int b = 0 ;
    uint64_t u = 0 ;
    while (true) {
        cerr << "Beginning batch " << b + 1 << endl ;
        uint64_t v = u ;
        for (int i = 0 ; i < threads ; i++) {
            u += fastq_entries[p][i].size() ;
        }
        if (v == u) {
            break ;
        }
        #pragma omp parallel for num_threads(threads + 2)
        for(int i = 0; i < threads + 2; i++) {
            if (i == 0) {
                // load next batch of entries
                load_batch_fastq(threads, batch_size, (p + 1) % 2) ;
                cerr << "Loaded." << endl ;
            } else if (i == 1) {
                // merge output of previous batch
                if (b >= 1) {
                    int y = 0 ;
                    for (const auto &batch : batches[(p + 1) % 2]) {
                        y += batch.size() ;
                        for (const auto fastq_entry : batch) {
                            if (search_solutions[current_batch].find(fastq_entry) == search_solutions[current_batch].end()) {
                                search_solutions[current_batch][fastq_entry] = 0 ;
                            }
                            search_solutions[current_batch][fastq_entry] += 1 ;
                        }
                    }
                    cerr << y << " total sequences." << endl ;
                }
                cerr << "Merged. " << search_solutions[current_batch].size() << " unique sequences." << endl ;
            } else {
                // process current batch
                batches[p][i - 2] = process_batch_fastq(index, fastq_entries[p][i - 2]) ;
            }
        }
        if (search_solutions[current_batch].size() >= 10000000) {
            cerr << "Memory limit reached, dumping output batch " << current_batch << ".." << endl ;
            current_batch += 1 ;
            unordered_map<fastq_entry_t, int> s ;
            search_solutions.push_back(s) ;
            OutputBatchArgs* b_args = new OutputBatchArgs() ;
            b_args->batch = current_batch ;
            output_batch((void*) b_args) ;
        }
        p += 1 ;
        p %= 2 ;
        b += 1 ;
        time_t s ;
        time(&s) ;
        if (s - t == 0) {
            s += 1 ;
        }
        cerr << "Processed batch " << std::left << std::setw(10) << b << ". Reads so far " << std::right << std::setw(12) << u << ". Reads per second: " <<  u / (s - t) << ". Time: " << std::setw(8) << std::fixed << s - t << "\n" ;
    }
    int y = 0 ;
    for (const auto &batch : batches[(p + 1) % 2]) {
        y += batch.size() ;
        for (const auto fastq_entry : batch) {
            if (search_solutions[current_batch].find(fastq_entry) == search_solutions[current_batch].end()) {
                search_solutions[current_batch][fastq_entry] = 0 ;
            }
            search_solutions[current_batch][fastq_entry] += 1 ;
        }
    }
    // last batch
    current_batch += 1 ;
    OutputBatchArgs* b_args = new OutputBatchArgs() ;
    b_args->batch = current_batch ;
    output_batch((void*) b_args) ;
    // cleanup
    kseq_destroy(fastq_iterator) ;
    gzclose(fastq_file) ;
    return u ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

int query(int argc, char *argv[]) {
    // parse arguments
    cout << "Single query mode.." << endl ;
    char *index_path = argv[1] ;
    rld_t *index = rld_restore(index_path) ;
    char *query_seq = argv[2] ;
    //
    vector<fastq_entry_t> solutions ;
    fastq_entry_t fastq_entry("Query", std::string(query_seq), std::string(query_seq)) ;
    ping_pong_search(index, fastq_entry, solutions) ;
    //
    for (const auto fastq_entry : solutions) {
        cout << "@" << fastq_entry.head << ".css" << "_" << fastq_entry.start << ":"
            << fastq_entry.start + fastq_entry.len - 1 << endl 
            << fastq_entry.seq << endl ;
    }
    return 0 ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

unordered_map<string, int> filter(rld_t* index, vector<string> lines) {
    unordered_map<string, int> solution ;
    for(const auto line: lines){
        int i = line.find_first_of(":") ;
        std::string seq(line, 0, i) ;
        int count = std::stoi(std::string(line, i + 2, line.length() - (i + 1)), nullptr, 10) ;
        if (!check_solution(index, seq)) {
            if (solution.find(seq) == solution.end()) {
                solution[seq] = 0 ;
            }
            solution[seq] += count ;
        }
    }
    return solution ;
}

vector<unordered_map<string, int>> pcheck_f3(const vector<vector<string>> &entries, rld_t *index, int threads) {
    DEBUG(threads = 1; )
    cout << "Running .. " << endl ;
    vector<unordered_map<string, int>> solutions (threads) ;
    #pragma omp parallel for num_threads(threads)
    for(int i = 0; i < threads; i++) {
        unordered_map<string, int> sol = filter(index, entries[i]) ;
        solutions[i] = sol ;
    }
    return solutions ;
}

int check_f3(int argc, char* argv[]) {
    char *index_path = argv[1] ;
    rld_t *index = rld_restore(index_path);
    std::ifstream ifs(argv[2]) ;
    int threads = atoi(argv[3]) ;

    int b = 0 ;
    int n = 0 ;
    int batch_size = 100000 ;
    std::string line ;
    vector<vector<string>> entries (threads) ;
    unordered_map<string, int> final_solution ;
    time_t t ;
    time(&t) ;
    while (std::getline(ifs, line)) {
        entries[n % threads].push_back(line);
        n += 1 ;
        if (n == batch_size) {
            n = 0 ;
            vector<unordered_map<string, int>> output = pcheck_f3(entries, index, threads);
            for (const auto &batch : output) {
                for (const auto sol : batch) {
                    if (final_solution.find(sol.first) == final_solution.end()) {
                        final_solution[sol.first] = 0 ;
                    }
                    //cout << sol.first << ": " << sol.second << endl ;
                    final_solution[sol.first] += sol.second ;
                }
            }
            b += 1 ;
            for (int i = 0 ; i < threads ; i++) {
                entries[i].clear() ;
            }
            time_t s ;
            time(&s) ;
            cerr << "Processed batch " << std::setw(10) << b << ". Reads so far " << std::setw(12) << b * batch_size << "Final solutions: " << setw(15) << final_solution.size() << ". Time: " << std::setw(8) << std::fixed << s - t << endl ;
        }
    }
    // last batch
    if(!entries.empty()) {
        vector<unordered_map<string, int>> output = pcheck_f3(entries, index, threads);
        for (const auto &batch : output) {
            for (const auto sol : batch) {
                if (final_solution.find(sol.first) == final_solution.end()) {
                    final_solution[sol.first] = 0 ;
                }
                cout << sol.first << ": " << sol.second << endl ;
                final_solution[sol.first] += sol.second ;
            }
        }
    }
    //while (std::getline(ifs, line)) {
    //    //cout << line << endl ;
    //    int i = line.find_first_of(":") ;
    //    std::string seq(line, 0, i) ;
    //    int count = std::stoi(std::string(line, i + 2, line.length() - (i + 1)), nullptr, 10) ;
    //    //cout << seq << ": " << count << endl ;
    //    if (!check_solution(index, seq)) {
    //        if (final_solution.find(seq) == final_solution.end()) {
    //            final_solution[seq] = 0 ;
    //        }
    //        final_solution[seq] += count ;
    //    }
    //    n += 1 ;
    //    if (n % 10000 == 0) {
    //        cout << n << " " << final_solution.size() << endl ;
    //    }
    //}
    return 0 ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

/** From ropebwt2 ********/

static inline int kputsn(const char *p, int l, kstring_t *s) {
    if (s->l + l + 1 >= s->m) {
        char *tmp;
        s->m = s->l + l + 2;
        kroundup32(s->m);
        if ((tmp = (char*)realloc(s->s, s->m))) s->s = tmp;
        else return EOF;
    }
    memcpy(s->s + s->l, p, l);
    s->l += l;
    s->s[s->l] = 0;
    return l;
}

/** Code adapted from ropebwt2 (main_ropebwt2 in main.c) **/
int main_index(int argc, char* argv[]) {
    // hardcoded parameters
    uint64_t m = (uint64_t)(.97 * 10 * 1024 * 1024 * 1024) + 1; // batch size for multi-string indexing
    int block_len = ROPE_DEF_BLOCK_LEN, max_nodes = ROPE_DEF_MAX_NODES, so = MR_SO_RCLO;
    int thr_min = 100; // switch to single thread when < 100 strings remain in a batch

    // the index
    mrope_t *mr = 0;

    // cmd line arguments
    bool binary_output = false;
    int c;
    while ((c = getopt (argc, argv, "a:bh")) != -1)
        switch (c) {
            case 'a':
                // append to existing index - we restore the index (it must be in binary FMR format)
                FILE *fp;
                if ((fp = fopen(optarg, "rb")) == 0) {
                    cerr << "fail to open file " << optarg << endl;
                    return 1;
                }
                mr = mr_restore(fp);
                fclose(fp);
                break;
            case 'b':
                binary_output = true;
                break;
            case 'h':
                cerr << "Usage: main index [-h] [-b] [-a index] <in.fq> > <index>" << endl;
                return 0;
            case '?':
                return 1;
            default:
                return 1;
        }
    if(optind == argc) {
        cerr << "Usage: main index [-h] [-b] [-a index] <in.fq> > <index>" << endl;
        return 1;
    }
    char *fpath = argv[optind]; // input sample

    // Initialize mr if not restored
    if (mr == 0) mr = mr_init(max_nodes, block_len, so);
    mr_thr_min(mr, thr_min);

    // Parsing the input sample
    gzFile fp = gzopen(fpath, "rb");
    kseq_t *ks = kseq_init(fp);
    kstring_t buf = { 0, 0, 0 }; // buffer, will contain the concatenation
    int l;
    uint8_t *s;
    int i;
    while ((l = kseq_read(ks)) >= 0) {
        s = (uint8_t*)ks->seq.s;

        // change encoding
        for (i = 0; i < l; ++i)
            s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5;

        // --- hard mask the sequence ---
        // if (ks->qual.l && min_q > 0)
        //   for (i = 0; i < l; ++i)
        // 	s[i] = ks->qual.s[i] - 33 >= min_q? s[i] : 5;

        // --- skip sequences containing ambiguous bases ---
        // if (flag & FLAG_NON) {
        //   for (i = 0; i < l; ++i)
        // 	if (s[i] == 5) break;
        //   if (i < l) continue;
        // }

        // --- cut at ambiguous bases and discard segment with length <INT
        // + cut one base if forward==reverse ---
        // if (flag & FLAG_CUTN) {
        //   int b, k;
        //   for (k = b = i = 0; i <= l; ++i) {
        // 	if (i == l || s[i] == 5) {
        // 	  int tmp_l = i - b;
        // 	  if (tmp_l >= min_cut_len) {
        // 	    if ((flag & FLAG_ODD) && is_rev_same(tmp_l, &s[k - tmp_l])) --k;
        // 	    s[k++] = 0;
        // 	  } else k -= tmp_l; // skip this segment
        // 	  b = i + 1;
        // 	} else s[k++] = s[i];
        //   }
        //   if (--k == 0) continue;
        //   ks->seq.l = l = k;
        // }
        // if ((flag & FLAG_ODD) && is_rev_same(l, s)) {
        //   ks->seq.s[--l] = 0;
        //   ks->seq.l = l;
        // }

        // Reverse the sequence
        for (i = 0; i < l>>1; ++i) {
            int tmp = s[l-1-i];
            s[l-1-i] = s[i];
            s[i] = tmp;
        }

        // Add forward to buffer
        kputsn((char*)ks->seq.s, ks->seq.l + 1, &buf);

        // Add reverse to buffer
        for (i = 0; i < l>>1; ++i) {
            int tmp = s[l-1-i];
            tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
            s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
            s[i] = tmp;
        }
        if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
        kputsn((char*)ks->seq.s, ks->seq.l + 1, &buf);

        if(buf.l >= m) {
            mr_insert_multi(mr, buf.l, (uint8_t*)buf.s, 1);
            buf.l = 0;
        }
    }

    if(buf.l) // last batch
        mr_insert_multi(mr, buf.l, (uint8_t*)buf.s, 1);

    free(buf.s);
    kseq_destroy(ks);
    gzclose(fp);

    // dump index to stdout
    if(binary_output)
        // binary FMR format
        mr_dump(mr, stdout);
    else {
        // FMD format
        mritr_t itr;
        const uint8_t *block;
        rld_t *e = 0;
        rlditr_t di;
        e = rld_init(6, 3);
        rld_itr_init(e, &di, 0);
        mr_itr_first(mr, &itr, 1);
        while ((block = mr_itr_next_block(&itr)) != 0) {
            const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
            while (q < end) {
                int c = 0;
                int64_t l;
                rle_dec1(q, c, l);
                rld_enc(e, &di, l, c);
            }
        }
        rld_enc_finish(e, &di);
        rld_dump(e, "-");
    }

    mr_destroy(mr);

    return 0;
}

/**********************************************************/

int main(int argc, char *argv[]) {
    string mode = argv[1] ;
    int retcode = 0 ;
    DEBUG(cerr << "DEBUG MODE" << endl ;)
        if(mode == "index") {
            retcode = main_index(argc - 1, argv + 1);
        } else if (mode == "sf3") {
            retcode = search_f3(argc - 1, argv + 1) ;
        } else if (mode == "cf3") {
            retcode = check_f3(argc - 1, argv + 1) ;
        } else if (mode == "query") {
            retcode = query(argc - 1, argv + 1) ;
        } else {
            retcode = 1 ;
        }
    return retcode ;
}
