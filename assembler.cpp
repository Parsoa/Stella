#include "assembler.hpp"

void Assembler::run() {
    auto c = Configuration::getInstance();
    int num_batches = c->aggregate_batches;
    int tau = c->cutoff;

    lprint({"Assembling high-abundance strings from", to_string(num_batches), "batches.."});
    #pragma omp parallel for num_threads(c->threads)
    for (int j = 0; j < num_batches; j++) {
        string s_j = std::to_string(j);
        string inpath = c->workdir + "/solution_batch_" + s_j + ".sfs";
        string outpath = c->workdir + "/solution_batch_" + s_j + ".assembled.sfs";
        ofstream outf(outpath);
        map<string, vector<SFS>> SFSs = parse_sfsfile(inpath, tau);
        for (map<string, vector<SFS>>::iterator it = SFSs.begin(); it != SFSs.end(); ++it) {
            string ridx = it->first;
            vector<SFS> sfs = it->second;
            vector<SFS> assembled_sfs = assemble(sfs);
            bool is_first = true;
            for (const SFS &sfs : assembled_sfs) {
                outf << (is_first ? ridx : "*") << "\t"
                    << "*"
                    << "\t" << sfs.s << "\t" << sfs.l << "\t" << sfs.c << endl;
                is_first = false;
            }
        }
        outf.close();
    }
}

vector<SFS> Assembler::assemble(vector<SFS> &sfs) {
    vector<SFS> assembled_sfs;
    sort(sfs.begin(), sfs.end());
    int i = 0;
    while (i < sfs.size()) {
        uint j;
        for (j = i + 1; j < sfs.size(); ++j) {
            if (sfs[j - 1].s + sfs[j - 1].l <= sfs[j].s) {
            // non-overlapping
                uint l = sfs[j - 1].s + sfs[j - 1].l - sfs[i].s;
                assembled_sfs.push_back(SFS(sfs[i].s, l, 1, sfs[i].isreversed));
                i = j;
                break;
            }
        }
        if (j == sfs.size()) {
            uint l = sfs[j - 1].s + sfs[j - 1].l - sfs[i].s;
            assembled_sfs.push_back(SFS(sfs[i].s, l, 1, sfs[i].isreversed));
            i = j;
        }
    }
    return assembled_sfs;
}
