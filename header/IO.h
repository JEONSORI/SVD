#pragma once

#include "Defines.h"
#include <string.h>

extern Control Cntl;

struct Control {
	string input;
	int num_groups, SVDRank;
	int num_GC_files = 0;
	vector<string> GC_file;
	vector<string> ASY_name;
};

const char BLANK = ' ';

const string CardName[3] 
{
	"GROUP",
	"SVD_RANK",
	"GC_FILE"
};

enum class XSEC
{
	GROUP,
	SVD_RANK,
	GC_FILE
};

enum GC_TYPE
{
	BASE,
	BRANCH
};

void Readinput(string filename);
int SplitFields(string oneline, string*& fields, const char *delimiter, bool repeat);
template <typename T>
inline void HDFread(hid_t &file, const char *name, hid_t type, T *data);