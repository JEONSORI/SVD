#pragma once
#include "Defines.h"
#include "IO.h"
#include "Array.h"
#include "Eigen/Dense"

extern GC_t GC;

enum sigType {
	sigTr,
	sigA,
	sigF,
	sigNf,
	sigKf,
	sigT,
	Chi,
	n2n,
	sigS,
	CHI,
	SPH,
	num_sigType
};

struct Deriv_t { // chi: isotope-wise chi / CHI: pin-averaged chi
	Array<double> D, sigTr, sigA, sigR, sigF, sigNf, sigKf, sigT, chi, sigS, n2n;
	Array<double> Flux, CHI, SPH, pnum135, surf_flux, surf_jout;
	Array<double> GC;
};

struct GC_t
{
	hid_t svdfile;
	string filename, asyname;
	int asyType;
	int num_pins, num_groups, num_iso = 30, nburnup, nstate = 0, nBranchState = 0;
	int nGcType;
	int nTmMax, nTfMax, nppmMax, nrhoMax;
	int num_nuclides = 30;
	int nuclides[30] = {
			0,   5000,   5010,  53135,  54135,  60147,  61147,  61148,
		61149,  62149,  64000,  92234,  92235,  92236,  92237,  92238,
		93237,  93238,  93239,  94238,  94239,  94240,  94241,  94242,
		95241,  95242,  95243,  96242,  96243,  96244
	};
	int nBranchPoint, nvar, SVDRank;
	int Base_num_row, Branch_num_row, num_row, num_col;
	Array<bool> lFuel, lbranch, lTm, lrho, lTf, lppm;
	Array<double> burnup, pnum, Tm, Tf, rho, ppm, pnum_avg, flux_avg;
	Array<double> b_Tm, b_Tf, b_rho, b_ppm;
	vector<Deriv_t> Base;
	vector<vector<Deriv_t>> dTm, dTf, drho, dppm;
	Array<int> nTm, nTf, nppm, nrho;
	Array<int> Num_Branches, BaseMap, BranchMap;
	Array<int> Tm_list, Tf_list, ppm_list, rho_list;
	Array<int> Bidx;
	Array<int> xs_count;
	Array<double> xs_avg;
	Eigen::MatrixXd A, Aorg;
	void ReadGCfile();
	void MappingPinGC(const int narray, const int asyType, Array<int>& PinIdx);
	void AllocGC(Deriv_t& BP, int num_iso, int num_groups, int nGcType);
	int GetIsoIdx(int ZID);
	void AssignGC(Deriv_t& BP, Array<double>& GC, vector<int>& isolist, int iGcType);
	void SetMatrix();
	void SetBaseMatrix();
	void SetBranchMatrix();
	void WriteSVDfile();
	void GetSVD();
	void SetCoeffMatrix();
};
