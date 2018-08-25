#include <iostream>

using namespace std;

struct blocksinfo_t {
    int nblocks, nelectrons, ncsfstotal, norbitals, nvectotal, nvecsize;
};
struct block_t {
    int blockid, ncsfs, nevs, iatjp, iaspa;
    double eav;
    double * eigenstates, * eigenenergies;
};
extern "C" int mixread(const char *, blocksinfo_t*, block_t**);

#define COUT_MEMBER(s, f) cout << #s "." #f " = " << s.f << endl

int main(int argc, char * argv[]) {
    if(argc != 2) {
        cerr << "ERROR: Wrong number of arguments." << endl;
        return 1;
    }

    blocksinfo_t blocksinfo;
    block_t * blocks;
    int * ints;

    cout << "Reading state from: " << argv[1] << endl;

    int status = mixread(argv[1], &blocksinfo, &blocks);
    if(status != 0) {
        cout << "mixread returned with status = " << status << endl;
        return 1;
    }

    COUT_MEMBER(blocksinfo, nblocks);
    COUT_MEMBER(blocksinfo, nelectrons);
    COUT_MEMBER(blocksinfo, ncsfstotal);
    COUT_MEMBER(blocksinfo, norbitals);
    COUT_MEMBER(blocksinfo, nvectotal);
    COUT_MEMBER(blocksinfo, nvecsize);

    for(int ib=0; ib < blocksinfo.nblocks; ib++) {
        cout << "block: " << ib << endl;
        cout << "  "; COUT_MEMBER(blocks[ib], blockid);
        cout << "  "; COUT_MEMBER(blocks[ib], ncsfs);
        cout << "  "; COUT_MEMBER(blocks[ib], nevs);
        cout << "  "; COUT_MEMBER(blocks[ib], iatjp);
        cout << "  "; COUT_MEMBER(blocks[ib], iaspa);
        cout << "  "; COUT_MEMBER(blocks[ib], eav);

        cout << "  "; COUT_MEMBER(blocks[ib], eigenenergies);
        for(size_t i=0; i < blocks[ib].nevs; i++) {
            cout.width(20);
            cout << blocks[ib].eigenenergies[i] << " ";
        }
        cout << endl;

        cout << "  "; COUT_MEMBER(blocks[ib], eigenstates);
        for(size_t i=0; i < blocks[ib].ncsfs; i++) {
            for(size_t j=0; j < blocks[ib].nevs; j++) {
                //cout << blocks[i].eigenstates[j][i] << " ";
                cout.width(20);
                cout << blocks[ib].eigenstates[i + j*blocks[ib].ncsfs] << " ";
            }
            cout << endl;
        }
    }

    /*for(int i=0; i < 1000; i++) {
        cout << ints[i] << " ";
    }
    cout << endl;*/

    return 0;
}
