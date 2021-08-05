#include "IO.h"
#include "XS.h"
#include "Array.h"

GC_t GC;

#define Joule2Mev 62415064799.632
void GC_t::ReadGCfile() {
	hid_t file;
	Array<int> PinIdx, Bidx;

	
	Array<int> ibu_branch;
	ibu_branch.Create(3, 2, 4);
	for (int ik = 0, i = 0; ik < 4; ik++) {
		for (int ir = 0; ir < 3; ir++) {
			for (int ic = 0; ic < 2; ic++) {
				ibu_branch(ir, ic, ik) = i;
				i++;
			}
		}
	}


	for (int iGC = 0; iGC < 4; iGC++) {
		Array<int> ibu;
		ibu.Create(&ibu_branch(0, 0, iGC), 3, 2);
		for (int ir = 0, i = 10; ir < 3; ir++) {
			for (int ic = 0; ic < 2; ic++) {
				ibu(ir, ic) = i;
				i--;
			}
		}
	}

	file = H5Fopen(this->filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	this->num_groups = Cntl.num_groups;
	this->SVDRank = Cntl.SVDRank;

	// Read assembly type, 1 : oct-symmetry
	HDFread(file, "ASY_TYPE", H5T_NATIVE_INT, &this->asyType);
	HDFread(file, "NUM_PIN", H5T_NATIVE_INT, &this->num_pins);

	int narray = (int)sqrt((double)this->num_pins);
	MappingPinGC(narray, asyType, PinIdx);


	// Read # of burnup steps
	HDFread(file, "NUM_BURNUP_POINT", H5T_NATIVE_INT, &this->nburnup);

	// Branch point indicators
	lbranch.Create(this->nburnup);
	lTm.Create(this->nburnup);
	lrho.Create(this->nburnup);
	lTf.Create(this->nburnup);
	lppm.Create(this->nburnup);

	// Branch point data
	int nBranchPoint;
	HDFread(file, "NUM_BRANCH_POINT", H5T_NATIVE_INT, &nBranchPoint);
	this->nBranchPoint = nBranchPoint;

	Bidx.Create(nBranchPoint);
	HDFread(file, "BRANCH_POINT_IDX", H5T_NATIVE_INT, Bidx());
	this->Bidx = Bidx;
	
	Array<int> iso_namelist; iso_namelist.Create(this->num_nuclides);
	HDFread(file, "ISO_NAME_LIST", H5T_NATIVE_INT, iso_namelist());
	int num_nuclides = this->num_nuclides;
	vector<int> iso_idxlist;
	for (int iso = 0; iso < num_nuclides; iso++) {
		int isoidx = GetIsoIdx(iso_namelist[iso]);
		iso_idxlist.push_back(isoidx);
	}
	

	Array<int> nBranchType;
	nBranchType.Create(4, nBranchPoint);
	HDFread(file, "NUM_BRANCH_PER_TYPE", H5T_NATIVE_INT, nBranchType());

	this->nTmMax = 1;
	this->nrhoMax = 1;
	this->nTfMax = 1;
	this->nppmMax = 1;

	nstate = nburnup;
	nBranchState = nBranchPoint;
	for (int i = 0; i < nBranchPoint; i++) {
		this->nTmMax = max(nBranchType(0, i), this->nTmMax);
		this->nTfMax = max(nBranchType(1, i), this->nTfMax);
		this->nppmMax = max(nBranchType(2, i), this->nppmMax);
		this->nrhoMax = max(nBranchType(3, i), this->nrhoMax);

		nstate += nBranchType(0, i) + nBranchType(1, i) + nBranchType(2, i) + nBranchType(3, i);
		nBranchState += nBranchType(0, i) + nBranchType(1, i) + nBranchType(2, i) + nBranchType(3, i);
	}

	this->Num_Branches.Create(4);
	this->Num_Branches[0] = this->nTmMax;
	this->Num_Branches[1] = this->nTfMax;
	this->Num_Branches[2] = this->nppmMax;
	this->Num_Branches[3] = this->nrhoMax;

	// Allocate arrays
	lFuel.Create(this->nGcType);
	burnup.Create(this->nburnup, this->nGcType);
	pnum.Create(this->num_iso, this->nburnup, this->nGcType);
	pnum_avg.Create(this->num_iso, this->nGcType);
	flux_avg.Create(this->num_groups, this->nGcType);
	pnum_avg = 1E-30; flux_avg = 1E-30;

	Tm.Create(this->nTmMax + 1, this->nburnup, this->nGcType);
	Tf.Create(this->nTfMax + 1, this->nburnup, this->nGcType);
	ppm.Create(this->nppmMax + 1, this->nburnup, this->nGcType);
	rho.Create(this->nrhoMax + 1, this->nburnup, this->nGcType);

	Tm_list.Create(this->nTmMax + 1, this->nburnup, this->nGcType);
	Tf_list.Create(this->nTfMax + 1, this->nburnup, this->nGcType);
	ppm_list.Create(this->nppmMax + 1, this->nburnup, this->nGcType);
	rho_list.Create(this->nrhoMax + 1, this->nburnup, this->nGcType);

	nTm.Create(this->nburnup);
	nTf.Create(this->nburnup);
	nppm.Create(this->nburnup);
	nrho.Create(this->nburnup);

	b_Tm.Create(this->nburnup, this->nGcType);
	b_Tf.Create(this->nburnup, this->nGcType);
	b_ppm.Create(this->nburnup);
	b_rho.Create(this->nburnup);

	this->Branch_num_row = nGcType * nBranchState;
	this->nvar = nTmMax + nTfMax + nppmMax + nrhoMax;
	BaseMap.Create(this->nGcType, this->nburnup);
	BranchMap.Create(this->nGcType, this->nvar, this->nBranchPoint);

	// Base GCs and the derivatives
	Base.resize(this->nburnup);
	for (int ibu = 0; ibu < this->nburnup; ibu++) {
		AllocGC(Base[ibu], num_iso, num_groups, nGcType);
	}

	dTm.resize(nburnup);
	dTf.resize(nburnup);
	dppm.resize(nburnup);
	drho.resize(nburnup);

	for (int i = 0; i < nBranchPoint; i++) {
		int ibu = Bidx[i];
		int n_Tm, n_Tf, n_ppm, n_rho;
		n_Tm = nBranchType(0, i);
		n_Tf = nBranchType(1, i);
		n_ppm = nBranchType(2, i);
		n_rho = nBranchType(3, i);

		lbranch[ibu] = true;
		if (n_Tm > 0) {
			lTm[ibu] = true;
			nTm[ibu] = n_Tm;
			dTm[ibu].resize(n_Tm);
			for (int idev = 0; idev < n_Tm; idev++) {
				AllocGC(dTm[ibu][idev], num_iso, num_groups, nGcType);
			}
		}
		if (n_Tf > 0) {
			lTf[ibu] = true;
			nTf[ibu] = n_Tf;
			dTf[ibu].resize(n_Tf);
			for (int idev = 0; idev < n_Tf; idev++) {
				AllocGC(dTf[ibu][idev], num_iso, num_groups, nGcType);
			}
		}
		if (n_ppm > 0) {
			lppm[ibu] = true;
			nppm[ibu] = n_ppm;
			dppm[ibu].resize(n_ppm);
			for (int idev = 0; idev < n_ppm; idev++) {
				AllocGC(dppm[ibu][idev], num_iso, num_groups, nGcType);
			}
		}
		if (n_rho > 0) {
			lrho[ibu] = true;
			nrho[ibu] = n_rho;
			drho[ibu].resize(n_rho);
			for (int idev = 0; idev < n_rho; idev++) {
				AllocGC(drho[ibu][idev], num_iso, num_groups, nGcType);
			}
		}
	}

	stringstream s1, s2, s3;

	Array<float> dbuf;
	Array<float> d_Tm, d_Tf, d_ppm, d_rho;

	// Read burnup exposure at each step
	dbuf.Create(this->nburnup, this->num_pins);
	HDFread(file, "BUEXP", H5T_NATIVE_FLOAT, dbuf());

	for (int itype = 0; itype < nGcType; itype++) {
		int ipin = PinIdx[itype]; //// + 1 ?
		for (int ibu = 0; ibu < this->nburnup; ibu++) {
			burnup(ibu, itype) = dbuf(ibu, ipin);
		}
	}

	// Read and assign the base and branch conditions
	d_Tm.Create(nTmMax + 1, nburnup);
	HDFread(file, "COND_TMOD", H5T_NATIVE_FLOAT, d_Tm());

	d_Tf.Create(nTfMax + 1, nburnup);
	HDFread(file, "COND_TFUEL", H5T_NATIVE_FLOAT, d_Tf());

	d_ppm.Create(nppmMax + 1, nburnup);
	HDFread(file, "COND_BORON", H5T_NATIVE_FLOAT, d_ppm());

	d_rho.Create(nrhoMax + 1, nburnup);
	HDFread(file, "COND_RHO", H5T_NATIVE_FLOAT, d_rho());


	// Open the groups and containing the data set
	hid_t ISO_GC, ISO_NAME_LIST, PNUM_LIST, PIN_FLUX;
	ISO_GC = H5Gopen2(file, "ISO_GC", H5P_DEFAULT);
	PNUM_LIST = H5Gopen2(file, "ISO_ND_LIST", H5P_DEFAULT);
	PIN_FLUX = H5Gopen2(file, "PIN_FLUX", H5P_DEFAULT);

  for (int ibu = 0; ibu < nburnup; ibu++) {
  	s2.str("");
  	s2 << "D" << setw(3) << setfill('0') << to_string(ibu);
  	string burnup_name = s2.str();
  	Array<float> fluxbuf;
  	
  	int ncol = nTm(ibu) + nTf(ibu) + nppm(ibu) + nrho(ibu) + 1;


  	Array<float> bufGC, bufSS, bufHSS, bufD, bufpnum;
  
  	string tmp;
  
  	hid_t ISO_GC_D = H5Gopen2(ISO_GC, burnup_name.c_str(), H5P_DEFAULT);
  

  	bufGC.Create(6, num_groups, num_nuclides, num_pins);
  	HDFread(ISO_GC_D, "BASE_GC", H5T_NATIVE_FLOAT, bufGC());

  	bufpnum.Create(num_nuclides, num_pins);
  	HDFread(PNUM_LIST, burnup_name.c_str(), H5T_NATIVE_FLOAT, bufpnum());
  
		fluxbuf.Create(ncol, num_groups, num_pins);
  	HDFread(PIN_FLUX, burnup_name.c_str(), H5T_NATIVE_FLOAT, fluxbuf());
  
  	Deriv_t& myPt = Base[ibu];
		for (int itype = 0; itype < this->nGcType; itype++) {
			int ipin = PinIdx[itype];
			AssignGC(myPt, bufGC, iso_idxlist, itype, ipin);
		}

		for (int itype = 0; itype < this->nGcType; itype++) {
			int ipin = PinIdx[itype];
			for (int iso = 0; iso < num_nuclides; iso++) {
				int isoidx = iso_idxlist[iso];
				float _pnum = bufpnum(iso, ipin);
				pnum(isoidx, ibu, itype) = _pnum;
				pnum_avg(isoidx, itype) += _pnum;
			}
			for (int ig = 0; ig < num_groups; ig++) {
				flux_avg(ig, itype) += fluxbuf(0, ig, ipin);
			}
		}
  
  	int icol = 1;
  	if (lTm[ibu]) {
  		for (int ivar = 0; ivar < nTm(ibu); ivar++) {
  			s3.str("");
  			s3 << setw(2) << setfill('0') << to_string(ivar + 1);
  			string num_pad = s3.str();
  
				tmp = "VAR_GC_TMOD" + num_pad;
				HDFread(ISO_GC_D, tmp.c_str(), H5T_NATIVE_FLOAT, bufGC());
  
  			Deriv_t& myPt = dTm[ibu][ivar];
				for (int itype = 0; itype < nGcType; itype++) {
					int ipin = PinIdx[itype];
					AssignGC(myPt, bufGC, iso_idxlist, itype, ipin);

					for (int iso = 0; iso < num_nuclides; iso++) { ////
						int isoidx = iso_idxlist[iso];
						double _pnum = bufpnum(iso, ipin);
						pnum_avg(isoidx, itype) += _pnum;
					}
					for (int ig = 0; ig < num_groups; ig++) {
						flux_avg(ig, itype) += fluxbuf(icol, ig, ipin);
					}
				}
  			icol++;
  		}
  	}
  
  	if (lTf[ibu]) {
  		for (int ivar = 0; ivar < nTf(ibu); ivar++) {
  			s3.str("");
  			s3 << setw(2) << setfill('0') << to_string(ivar + 1);
  			string num_pad = s3.str();
  
				tmp = "VAR_GC_TFUEL" + num_pad;
				HDFread(ISO_GC_D, tmp.c_str(), H5T_NATIVE_FLOAT, bufGC());
  
  			Deriv_t& myPt = dTf[ibu][ivar];

				for (int itype = 0; itype < this->nGcType; itype++) {
					int ipin = PinIdx[itype];
					AssignGC(myPt, bufGC, iso_idxlist, itype, ipin);

					for (int iso = 0; iso < num_nuclides; iso++) {
						int isoidx = iso_idxlist[iso];
						double _pnum = bufpnum(iso, ipin);
						pnum_avg(isoidx, itype) += _pnum;
					}
					for (int ig = 0; ig < num_groups; ig++) {
						flux_avg(ig, itype) += fluxbuf(icol, ig, ipin);
					}
				}
  			icol++;
  		}
  	}
  
  	if (lppm[ibu]) {
  		for (int ivar = 0; ivar < nppm(ibu); ivar++) {
  			s3.str("");
  			s3 << setw(2) << setfill('0') << to_string(ivar + 1);
  			string num_pad = s3.str();
  
				tmp = "VAR_GC_BORON" + num_pad;
				HDFread(ISO_GC_D, tmp.c_str(), H5T_NATIVE_FLOAT, bufGC());
  
  			Deriv_t& myPt = dppm[ibu][ivar];

				for (int itype = 0; itype < this->nGcType; itype++) {
					int ipin = PinIdx[itype];
					AssignGC(myPt, bufGC, iso_idxlist, itype, ipin);
					for (int iso = 0; iso < num_nuclides; iso++) {
						int isoidx = iso_idxlist[iso];
						double _pnum = bufpnum(iso, ipin);
						pnum_avg(isoidx, itype) += _pnum;
					}
					for (int ig = 0; ig < num_groups; ig++) {
						flux_avg(ig, itype) += fluxbuf(icol, ig, ipin);
					}
				}
  			icol++;
  		}
  	}
  
  	if (lrho[ibu]) {
  		for (int ivar = 0; ivar < nrho(ibu); ivar++) {
  			s3.str("");
  			s3 << setw(2) << setfill('0') << to_string(ivar + 1);
  			string num_pad = s3.str();
  
				tmp = "VAR_GC_RHO" + num_pad;
				HDFread(ISO_GC_D, tmp.c_str(), H5T_NATIVE_FLOAT, bufGC());

  			Deriv_t& myPt = drho[ibu][ivar];
				for (int itype = 0; itype < nGcType; itype++) {
					int ipin = PinIdx[itype];
					AssignGC(myPt, bufGC, iso_idxlist, itype, ipin);

					for (int iso = 0; iso < num_nuclides; iso++) {
						int isoidx = iso_idxlist[iso];
						double _pnum = bufpnum(iso, ipin);
						pnum_avg(isoidx, itype) += _pnum;
					}
					for (int ig = 0; ig < num_groups; ig++) {
						flux_avg(ig, itype) += fluxbuf(icol, ig, ipin);
					}
				}
  			icol++;
  		}
  	}
  
  	H5Gclose(ISO_GC_D);
  
  }

	// Close the groups for the next iteration

	H5Gclose(ISO_GC);
	H5Gclose(PNUM_LIST);
	H5Gclose(PIN_FLUX);

	double inv_nstate = 1.0 / ((double)nstate);
	pnum_avg *= inv_nstate;
	flux_avg *= inv_nstate;

	// Close the HDF5 file
	H5Fclose(file);
}

void GC_t::MappingPinGC(const int narray, const int asyType, Array<int>& PinIdx)
{
	int narray2 = narray / 2 + narray % 2;
	int npin;
	Array<int> pinList, pinType;

	pinList.Create(narray, narray);

	for (int i = 0, idx = 0; i < narray; i++) {
		for (int j = 0; j < narray; j++) {
			pinList(j, i) = idx; idx++;
		}
	}

	pinType.Create(narray, narray);
	switch (asyType) {
	case 1:
		npin = narray2 * (narray2 + 1) / 2;
		for (int i = 0, pType = 0; i < narray2; i++)
			for (int j = i; j < narray2; j++) {
				int ib = narray - 1 - i, jb = narray - 1 - j;
				pinType(j, i) = pType; pinType(i, j) = pType;
				pinType(jb, i) = pType; pinType(i, jb) = pType;
				pinType(j, ib) = pType; pinType(ib, j) = pType;
				pinType(jb, ib) = pType; pinType(ib, jb) = pType;
				pType++;
			}
		break;
	case 2:
		npin = narray2 * narray;
		for (int i = 0, pType = 0; i < narray2; i++)
			for (int j = 0; j < narray; j++) {
				int ib = narray - 1 - i;
				pinType(j, i) = pType;
				pinType(j, ib) = pType;
				pType++;
			}
		break;
	case 3:
		npin = narray2 * narray;
		for (int i = 0, pType = 0; i < narray; i++)
			for (int j = 0; j < narray2; j++) {
				int jb = narray - 1 - j;
				pinType(j, i) = pType;
				pinType(jb, i) = pType;
				pType++;
			}
		break;
	case 4:
		npin = narray * narray;
		for (int i = 0, pType = 0; i < narray; i++)
			for (int j = 0; j < narray; j++) {
				pinType(j, i) = pType;
				pType++;
			}
		break;
	}

	this->nGcType = npin;

	PinIdx.Create(npin);

	bool *lread = new bool[npin];
	for (int i = 0; i < npin; i++) lread[i] = false;

	for (int i = 0; i < narray; i++)
		for (int j = 0; j < narray; j++) {
			int iGcType = pinType(j, i);
			if (!lread[iGcType]) {
				PinIdx[iGcType] = pinList(j, i);
				lread[iGcType] = true;
			}
		}

	delete[] lread;
}

void GC_t::AllocGC(Deriv_t& BP, int num_iso, int num_groups, int nGcType)
{
	BP.CHI.Create(num_groups, nGcType);
	BP.SPH.Create(num_groups, nGcType); BP.SPH = 1.0;

	BP.sigTr.Create(num_iso, num_groups, nGcType);
	BP.sigT.Create(num_iso, num_groups, nGcType);
	BP.sigA.Create(num_iso, num_groups, nGcType);
	BP.sigF.Create(num_iso, num_groups, nGcType);
	BP.sigNf.Create(num_iso, num_groups, nGcType);
	BP.sigKf.Create(num_iso, num_groups, nGcType);
	BP.chi.Create(num_iso, num_groups, nGcType);
	BP.sigS.Create(num_iso, num_groups, num_groups, nGcType);
	BP.n2n.Create(num_iso, num_groups, nGcType);

	BP.pnum135.Create(2, nGcType);

	int num_react = sigType::num_sigType - 4;
	BP.GC.Create(num_groups, num_iso, num_react, nGcType);
}

int GC_t::GetIsoIdx(int ZID) {
	for (int i = 0; i < num_nuclides; i++)
		if (nuclides[i] == ZID) return i;
};

void GC_t::AssignGC(Deriv_t& BP, Array<float>& GC, vector<int>& isolist, int iGcType, int ipin)
{
	int niso = isolist.size();

	for (int iso = 0; iso < niso; iso++) {
		int isoidx = isolist[iso];
		for (int ig = 0; ig < num_groups; ig++) {
			BP.GC(ig, isoidx, sigType::sigT,   iGcType)   = GC(0, ig, iso, ipin);
			BP.GC(ig, isoidx, sigType::sigA,   iGcType)   = GC(2, ig, iso, ipin);
			BP.GC(ig, isoidx, sigType::sigF,   iGcType)   = GC(3, ig, iso, ipin);
			BP.GC(ig, isoidx, sigType::sigNf,  iGcType)   = GC(4, ig, iso, ipin);
			BP.GC(ig, isoidx, sigType::sigCap, iGcType)   = GC(5, ig, iso, ipin);
		}
	}
}

void GC_t::SetMatrix()
{
	int num_react = sigType::num_sigType - 4; // micro form

	this->num_row = nGcType * nstate;
	this->num_col = num_iso * (num_react - 2) * num_groups; // except for sigTr, sigKf

	A.resize(num_row, num_col); 
	Aorg.resize(num_row, num_col);

	int ir = 0;

	for (int ibu = 0; ibu < nburnup; ibu++) {
		Deriv_t& myPt = Base[ibu];
		for (int iGcType = 0; iGcType < nGcType; iGcType++) {
			BaseMap(iGcType, ibu) = ir;
			int _react = 0;
			for (int react = 0; react < num_react; react++) {
				if (react == sigType::sigTr || react == sigType::sigKf) continue;
				for (int ig = 0; ig < num_groups; ig++) {
					double wt_flux = flux_avg(ig, iGcType);
					for (int iso = 0; iso < num_iso; iso++) {
						double wt_pnum = pnum_avg(iso, iGcType);
						int ic = _react * num_groups * num_iso + ig * num_iso + iso;
						Aorg(ir, ic) = myPt.GC(ig, iso, react, iGcType);
						A(ir, ic) = myPt.GC(ig, iso, react, iGcType) * wt_pnum * wt_flux;
					}
				}
				_react++;
			}
			ir++;
		}
		int bibu, ivar = 0;
		for (bibu = 0; bibu < nBranchPoint; bibu++) {
			if (this->Bidx(bibu) == ibu) break;
		}
		if (lTm[ibu]) {
			for (int idev = 0; idev < nTm[ibu]; idev++) {
				Deriv_t& myPt = dTm[ibu][idev];
				for (int iGcType = 0; iGcType < nGcType; iGcType++) {
					BranchMap(iGcType, ivar, bibu) = ir;
					int _react = 0;
					for (int react = 0; react < num_react; react++) {
						if (react == sigType::sigTr || react == sigType::sigKf) continue;
						for (int ig = 0; ig < num_groups; ig++) {
							double wt_flux = flux_avg(ig, iGcType);
							for (int iso = 0; iso < num_iso; iso++) {
								double wt_pnum = pnum_avg(iso, iGcType);
								int ic = _react * num_groups * num_iso + ig * num_iso + iso;
								Aorg(ir, ic) = myPt.GC(ig, iso, react, iGcType);
								A(ir, ic) = myPt.GC(ig, iso, react, iGcType) * wt_pnum * wt_flux;
							}
						}
						_react++;
					}
					ir++;
				}
				ivar++;
			}
		}
		if (lTf[ibu]) {
			for (int idev = 0; idev < nTf[ibu]; idev++) {
				Deriv_t& myPt = dTf[ibu][idev];
				for (int iGcType = 0; iGcType < nGcType; iGcType++) {
					BranchMap(iGcType, ivar, bibu) = ir;
					int _react = 0;
					for (int react = 0; react < num_react; react++) {
						if (react == sigType::sigTr || react == sigType::sigKf) continue;
						for (int ig = 0; ig < num_groups; ig++) {
							double wt_flux = flux_avg(ig, iGcType);
							for (int iso = 0; iso < num_iso; iso++) {
								double wt_pnum = pnum_avg(iso, iGcType);
								int ic = _react * num_groups * num_iso + ig * num_iso + iso;
								Aorg(ir, ic) = myPt.GC(ig, iso, react, iGcType);
								A(ir, ic) = myPt.GC(ig, iso, react, iGcType) * wt_pnum * wt_flux;
							}
						}
						_react++;
					}
					ir++;
				}
				ivar++;
			}
		}
		if (lppm[ibu]) {
			for (int idev = 0; idev < nppm[ibu]; idev++) {
				Deriv_t& myPt = dppm[ibu][idev];
				for (int iGcType = 0; iGcType < nGcType; iGcType++) {
					BranchMap(iGcType, ivar, bibu) = ir;
					int _react = 0;
					for (int react = 0; react < num_react; react++) {
						if (react == sigType::sigTr || react == sigType::sigKf) continue;
						for (int ig = 0; ig < num_groups; ig++) {
							double wt_flux = flux_avg(ig, iGcType);
							for (int iso = 0; iso < num_iso; iso++) {
								double wt_pnum = pnum_avg(iso, iGcType);
								int ic = _react * num_groups * num_iso + ig * num_iso + iso;
								Aorg(ir, ic) = myPt.GC(ig, iso, react, iGcType);
								A(ir, ic) = myPt.GC(ig, iso, react, iGcType) * wt_pnum * wt_flux;
							}
						}
						_react++;
					}
					ir++;
				}
				ivar++;
			}
		}
		if (lrho[ibu]) {
			for (int idev = 0; idev < nrho[ibu]; idev++) {
				Deriv_t& myPt = drho[ibu][idev];
				for (int iGcType = 0; iGcType < nGcType; iGcType++) {
					BranchMap(iGcType, ivar, bibu) = ir;
					int _react = 0;
					for (int react = 0; react < num_react; react++) {
						if (react == sigType::sigTr || react == sigType::sigKf) continue;
						for (int ig = 0; ig < num_groups; ig++) {
							double wt_flux = flux_avg(ig, iGcType);
							for (int iso = 0; iso < num_iso; iso++) {
								double wt_pnum = pnum_avg(iso, iGcType);
								int ic = _react * num_groups * num_iso + ig * num_iso + iso;
								Aorg(ir, ic) = myPt.GC(ig, iso, react, iGcType);
								A(ir, ic) = myPt.GC(ig, iso, react, iGcType) * wt_pnum * wt_flux;
								if (react == sigType::sigKf) A(ir, ic) *= Joule2Mev;
							}
						}
						_react++;
					}
					ir++;
				}
				ivar++;
			}
		}
	}
}

//void GC_t::SetBaseMatrix()
//{
//	int num_react = sigType::num_sigType - 3;
//
//	Base_num_row = nGcType * nburnup;
//	num_col = num_iso * num_react * num_groups;
//
//	A_base.resize(Base_num_row, num_col);
//
//	for (int ibu = 0; ibu < nburnup; ibu++) {
//		Deriv_t& myPt = Base[ibu];
//		for (int iGcType = 0; iGcType < nGcType; iGcType++) {
//			int ir = ibu * nGcType + iGcType;
//			for (int react = 0; react < num_react; react++) {
//				for (int iso = 0; iso < num_iso; iso++) {
//					double wt_pnum = pnum_avg(iso, iGcType);
//					for (int ig = 0; ig < num_groups; ig++) {
//						int ic = react * num_iso * num_groups + iso * num_groups + ig;
//						A_base(ir, ic) = myPt.GC(ig, iso, react, iGcType) * wt_pnum;
//					}
//				}
//			}
//		}
//	}
//}
//
//void GC_t::SetBranchMatrix()
//{
//	int num_react = sigType::num_sigType - 3;
//
//	Branch_num_row = nGcType * nBranchState;
//
//	A_branch.resize(Branch_num_row, num_col);
//
//	int ir = 0;
//
//	for (int bibu = 0; bibu < nBranchPoint; bibu++) {
//		int ibu = Bidx(bibu);
//		Deriv_t& myPt = Base[ibu];
//		for (int iGcType = 0; iGcType < nGcType; iGcType++) {
//			for (int react = 0; react < num_react; react++) {
//				for (int iso = 0; iso < num_iso; iso++) {
//					double wt_pnum = pnum_avg(iso, iGcType);
//					for (int ig = 0; ig < num_groups; ig++) {
//						int ic = react * num_iso * num_groups + iso * num_groups + ig; 
//						A_branch(ir, ic) = myPt.GC(ig, iso, react, iGcType) * wt_pnum;
//					}
//				}
//			}
//			ir++;
//		}
//
//		for (int idev = 0; idev < nTm[ibu]; idev++) {
//			Deriv_t& myPt = dTm[ibu][idev];
//			for (int iGcType = 0; iGcType < nGcType; iGcType++) {
//				for (int react = 0; react < num_react; react++) {
//					for (int iso = 0; iso < num_iso; iso++) {
//						double wt_pnum = pnum_avg(iso, iGcType);
//						for (int ig = 0; ig < num_groups; ig++) {
//							int ic = react * num_iso * num_groups + iso * num_groups + ig;
//							A_branch(ir, ic) = myPt.GC(ig, iso, react, iGcType) * wt_pnum;
//						}
//					}
//				}
//				ir++;
//			}
//		}
//		for (int idev = 0; idev < nTf[ibu]; idev++) {
//			Deriv_t& myPt = dTf[ibu][idev];
//			for (int iGcType = 0; iGcType < nGcType; iGcType++) {
//				for (int react = 0; react < num_react; react++) {
//					for (int iso = 0; iso < num_iso; iso++) {
//						double wt_pnum = pnum_avg(iso, iGcType);
//						for (int ig = 0; ig < num_groups; ig++) {
//							int ic = react * num_iso * num_groups + iso * num_groups + ig;
//							A_branch(ir, ic) = myPt.GC(ig, iso, react, iGcType) * wt_pnum;
//						}
//					}
//				}
//				ir++;
//			}
//		}
//		for (int idev = 0; idev < nppm[ibu]; idev++) {
//			Deriv_t& myPt = dppm[ibu][idev];
//			for (int iGcType = 0; iGcType < nGcType; iGcType++) {
//				for (int react = 0; react < num_react; react++) {
//					for (int iso = 0; iso < num_iso; iso++) {
//						double wt_pnum = pnum_avg(iso, iGcType);
//						for (int ig = 0; ig < num_groups; ig++) {
//							int ic = react * num_iso * num_groups + iso * num_groups + ig;
//							A_branch(ir, ic) = myPt.GC(ig, iso, react, iGcType) * wt_pnum;
//						}
//					}
//				}
//				ir++;
//			}
//		}
//		for (int idev = 0; idev < nrho[ibu]; idev++) {
//			Deriv_t& myPt = drho[ibu][idev];
//			for (int iGcType = 0; iGcType < nGcType; iGcType++) {
//				for (int react = 0; react < num_react; react++) {
//					for (int iso = 0; iso < num_iso; iso++) {
//						double wt_pnum = pnum_avg(iso, iGcType);
//						for (int ig = 0; ig < num_groups; ig++) {
//							int ic = react * num_iso * num_groups + iso * num_groups + ig;
//							A_branch(ir, ic) = myPt.GC(ig, iso, react, iGcType) * wt_pnum;
//						}
//					}
//				}
//				ir++;
//			}
//		}
//	}
//}

void GC_t::WriteSVDfile()
{
	string filename;
	filename = this->asyname + "_SVD.h5";

	svdfile = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 

	hsize_t dims[3];

	dims[0] = 1;
	HDFwrite(svdfile, "NUM_STATES", 1, dims, H5T_NATIVE_INT, &Branch_num_row);
	HDFwrite(svdfile, "SVD_RANK", 1, dims, H5T_NATIVE_INT, &SVDRank);
	HDFwrite(svdfile, "NUM_CROSS_SECTIONS", 1, dims, H5T_NATIVE_INT, &num_col);
	HDFwrite(svdfile, "NUM_BURNUP", 1, dims, H5T_NATIVE_INT, &nburnup);
	HDFwrite(svdfile, "NUM_BRANCH_BURNUP", 1, dims, H5T_NATIVE_INT, &nBranchPoint);
	HDFwrite(svdfile, "NUM_PINTYPE", 1, dims, H5T_NATIVE_INT, &nGcType);
	HDFwrite(svdfile, "NUM_NUCLIDES", 1, dims, H5T_NATIVE_INT, &num_iso);

	dims[0] = 4;
	HDFwrite(svdfile, "NUM_BRANCHES", 1, dims, H5T_NATIVE_INT, Num_Branches());

	dims[0] = nGcType;
	dims[1] = num_iso;
	HDFwrite(svdfile, "WT_PNUM", 2, dims, H5T_NATIVE_DOUBLE, pnum_avg());

	dims[1] = num_groups;
	HDFwrite(svdfile, "WT_FLUX", 2, dims, H5T_NATIVE_DOUBLE, flux_avg());

	GetSVD();

	H5Fclose(svdfile);
}

void GC_t::GetSVD()
{

//	if (GCType == GC_TYPE::BASE) {
//		num_row = this->Base_num_row;
//		num_col = this->num_col;
//		Type = H5Gcreate2(svdfile, "Base", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//	}
//	else {
//		num_row = this->Branch_num_row;
//		num_col = this->num_col;
//		Type = H5Gcreate2(svdfile, "Branch", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//	}

	Eigen::MatrixXd U, V, F, v;
	Eigen::VectorXd S, s;
	F.resize(num_row, SVDRank);
	s.resize(SVDRank);

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::DecompositionOptions::ComputeThinU | 
		Eigen::DecompositionOptions::ComputeThinV);

	S = svd.singularValues();
	U = svd.matrixU();
	V = svd.matrixV();


	for (int ic = 0; ic < SVDRank; ic++) {
		for (int ir = 0; ir < num_row; ir++) {
			F(ir, ic) = U(ir, ic) * S(ic);
		}
	}
	
////	v.resize(SVDRank, num_col);
	v.resize(num_col, SVDRank);
	for (int ic = 0; ic < SVDRank; ic++) {
	////	v.row(ic) = V.col(ic);
		v.col(ic) = V.col(ic);
		s(ic) = S(ic);
	}

	Eigen::MatrixXd F_base, F_branch;
	F_base.resize(SVDRank, nburnup * nGcType);
	F_branch.resize(SVDRank, Branch_num_row);

	int ir_base = 0, ir_branch = 0;
	for (int ibu = 0; ibu < nburnup; ibu++) {
		for (int iGcType = 0; iGcType < nGcType; iGcType++) {
			int ir = BaseMap(iGcType, ibu);
			F_base.col(ir_base) = F.row(ir);
			ir_base++;
		}
	}
	for (int bibu = 0; bibu < nBranchPoint; bibu++) {
		int ibu = this->Bidx(bibu);
		for (int iGcType = 0; iGcType < nGcType; iGcType++) {
			int ir = BaseMap(iGcType, ibu);
			F_branch.col(ir_branch) = F.row(ir);
			ir_branch++;
		}
		for (int ivar = 0; ivar < nvar; ivar++) {
			for (int iGcType = 0; iGcType < nGcType; iGcType++) {
				int ir = BranchMap(iGcType, ivar, bibu);
				F_branch.col(ir_branch) = F.row(ir);
				ir_branch++;
			}
		}
	}

	hsize_t dims[2];

////	dims[0] = num_col;
////	dims[1] = num_row;
////	HDFwrite(svdfile, "A", 2, dims, H5T_NATIVE_DOUBLE, A.data());
////	HDFwrite(svdfile, "Aorg", 2, dims, H5T_NATIVE_DOUBLE, Aorg.data()); ////
////
////	dims[0] = min(num_row, num_col);
////	dims[1] = num_row;
////	HDFwrite(svdfile, "U", 2, dims, H5T_NATIVE_DOUBLE, U.data());
////	dims[0] = SVDRank;
////	dims[1] = num_row;
////	HDFwrite(svdfile, "F", 2, dims, H5T_NATIVE_DOUBLE, F.data());
////	
////	dims[1] = num_col;
////	dims[0] = SVDRank;
////	HDFwrite(svdfile, "V", 2, dims, H5T_NATIVE_DOUBLE, V.data());

////	dims[1] = SVDRank;
////	dims[0] = num_col;
	dims[0] = SVDRank;
	dims[1] = num_col;

	HDFwrite(svdfile, "v", 2, dims, H5T_NATIVE_DOUBLE, v.data()); //////// check !!!!!!!!!!!!!

	dims[0] = min(num_row, num_col);
////	HDFwrite(svdfile, "S", 1, dims, H5T_NATIVE_DOUBLE, S.data());
////	dims[0] = SVDRank;
////	HDFwrite(svdfile, "s", 1, dims, H5T_NATIVE_DOUBLE, s.data());

	dims[1] = SVDRank;
	dims[0] = nburnup * nGcType;
	HDFwrite(svdfile, "F_base", 2, dims, H5T_NATIVE_DOUBLE, F_base.data());

	dims[0] = Branch_num_row;
	HDFwrite(svdfile, "F_branch", 2, dims, H5T_NATIVE_DOUBLE, F_branch.data());


//	hid_t SingularVector;
	
//	SingularVector = H5Gcreate2(Type, "v", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//	dims[0] = num_col;
//	for (int i = 0; i < min(num_row, num_col); i++) {
//		stringstream ivec;
//		ivec << setfill('0') << setw(5) << i;
//		HDFwrite(SingularVector, (ivec.str()).c_str(), 1, dims, H5T_NATIVE_DOUBLE, v[i].data());
//	}

}




//A.resize(4, 2);
//int i = 0;
//for (int ir = 0; ir < 4; ir++) {
//	for (int ic = 0; ic < 2; ic++) {
//		A(ir, ic) = i;
//		i++;
//	}
//}
//Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::DecompositionOptions::ComputeThinU | //Eigen::DecompositionOptions::ComputeThinV);
//
//Eigen::MatrixXd U, V;
//Eigen::VectorXd S;
//
//S = svd.singularValues();
//U = svd.matrixU();
//V = svd.matrixV();
//
//ofstream fout("SVD.txt");
//fout.setf(ios::scientific, ios::floatfield);
//fout << setprecision(12) << S << endl;
//fout << setprecision(12) << U << endl;
//fout << setprecision(12) << V << endl;
//
//hid_t svdfile;
//string filename;
//
//filename = "SVD.h5";
//
//svdfile = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//hsize_t dims[2];
//dims[0] = 4; dims[1] = 2;
//HDFwrite(svdfile, "U", 2, dims, H5T_NATIVE_DOUBLE, U.data());
//dims[0] = 2; dims[1] = 2;
//HDFwrite(svdfile, "V", 2, dims, H5T_NATIVE_DOUBLE, V.data());
//HDFwrite(svdfile, "S", 1, dims, H5T_NATIVE_DOUBLE, S.data());
//
//H5Fclose(svdfile);