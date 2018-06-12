#include <iostream>
#include <iomanip>

using namespace std;

struct orbitals_t {
    // integer(c_int) :: npy
    // integer(c_int) :: naky
    // real(c_double) :: ey
    // integer(c_int) :: my
    // type(c_ptr) :: pa
    // type(c_ptr) :: qa
    // type(c_ptr) :: ra
    int npy, naky;
    double ey;
    int my;
    double *pa, *qa, *ra;
};
extern "C" void rwfnread(const char *, int*, orbitals_t**);

#define COUT_MEMBER(s, f) cout << #s "." #f " = " << s.f << endl

int main(int argc, char * argv[]) {
    if(argc != 2) {
        cerr << "ERROR: Wrong number of arguments." << endl;
        return 1;
    }

    int norbitals;
    orbitals_t * orbitals;

    cout << "Reading orbitals from: " << argv[1] << endl;
    cout << "Calling: rwfnread" << endl;
    rwfnread(argv[1], &norbitals, &orbitals);

    cout << "Number of orbitals: " << norbitals << endl;
    for(int i = 0; i < norbitals; i++) {
        cout << setw(2) << i << ":";
        cout << setw(3) << orbitals[i].npy << setw(3) << orbitals[i].naky;
        cout << setw(15) << orbitals[i].ey;
        cout << setw(6) << orbitals[i].my;
        cout << endl;
    }

    orbitals_t orb = orbitals[1];
    for(int i = 0; i < 5; i++) {
        cout << orb.ra[i] << " ";
        cout << orb.pa[i] << " ";
        cout << orb.qa[i] << " ";
        cout << endl;
    }

    return 0;
}
