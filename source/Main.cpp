#include "Defines.h"
#include "IO.h"
#include "XS.h"

Control Cntl;

int main(int argc, char* argv[]) {
	Cntl.input = argv[1];
	Readinput(Cntl.input);
	GC_t *GC = new GC_t[Cntl.num_GC_files];
	for (int iGC = 0; iGC < Cntl.num_GC_files; iGC++) {
		GC[iGC].filename = Cntl.GC_file[iGC];
		GC[iGC].asyname = Cntl.ASY_name[iGC];
		GC[iGC].ReadGCfile();
		GC[iGC].SetMatrix();
		cout << GC[iGC].asyname << " Matrix Done!" << endl;
		GC[iGC].WriteSVDfile();
		cout << GC[iGC].asyname << " SVD Done!" << endl;
	}
}