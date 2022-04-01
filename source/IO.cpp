#include "Defines.h"
#include "IO.h"


void Uppercase(string& oneline)
{
	transform(oneline.begin(), oneline.end(), oneline.begin(), ::toupper);
}

int GetCardID(string oneline)
{
	int num_cards = 3;

	int pos_end = oneline.find(BLANK, 0);
	string block = oneline.substr(0, pos_end);

	Uppercase(block);
	for (int i = 0; i < num_cards; i++)
		if (!block.compare(CardName[i])) return i;
}



int Repeat(string& field)
{
	if (field.find('*') == string::npos) return 1;
	else {
		int pos = field.find('*');
		int rpt = atoi(field.substr(0, pos).c_str());
		field = field.substr(pos + 1, field.size() - pos - 1);
		return rpt;
	}
}

int CountFields(string oneline, const char *delimiter, bool repeat = true)
{
	int num_fields = 0;
	int rpt;
	char *pch;

	pch = strtok((char*)oneline.c_str(), delimiter);
	if (pch != NULL) {
		string field = pch;
		rpt = 1;
		if (repeat) rpt = Repeat(field);
		num_fields += rpt;
	}

	while (pch != NULL) {
		pch = strtok(NULL, delimiter);
		if (pch != NULL) {
			string field = pch;
			rpt = 1;
			if (repeat) rpt = Repeat(field);
			num_fields += rpt;
		}
	}

	return num_fields;
}

int SplitFields(string oneline, string*& fields, const char *delimiter, bool repeat)
{
	int num_fields = 0;
	int rpt, idx = 0;
	char *pch;

	num_fields = CountFields(oneline, delimiter, repeat);

	if (!num_fields) return 0;

	if (fields != NULL) delete[] fields;
	fields = new string[num_fields];

	num_fields = 0;

	pch = strtok((char*)oneline.c_str(), delimiter);
	if (pch != NULL) {
		string field = pch;
		rpt = 1;
		if (repeat) rpt = Repeat(field);
		num_fields += rpt;
		for (int i = 0; i < rpt; i++) {
			fields[idx] = field;
			idx++;
		}
	}

	while (pch != NULL) {
		pch = strtok(NULL, delimiter);
		if (pch != NULL) {
			string field = pch;
			rpt = 1;
			if (repeat) rpt = Repeat(field);
			num_fields += rpt;
			for (int i = 0; i < rpt; i++) {
				fields[idx] = field;
				idx++;
			}
		}
	}

	return num_fields;
}



void Readinput(string filename) {
	ifstream fin;
	string oneline;
	string *fields = NULL;
	int card;
	int num_fields;

	fin.open(filename);

	while (!fin.eof()) {
		getline(fin, oneline);
		if (oneline.empty()) continue;
		card = GetCardID(oneline);
		num_fields = SplitFields(oneline, fields, " ", true);
		switch (card) {
		case (int)XSEC::GROUP:
			Cntl.num_groups = stoi(fields[1]);
			break;
		case (int)XSEC::SVD_RANK:
			Cntl.SVDRank = stoi(fields[1]);
			break;
		case (int)XSEC::GC_FILE:
			Cntl.num_GC_files++;
			Cntl.GC_file.push_back(fields[2]);
			Cntl.ASY_name.push_back(fields[1]);
			break;
		}
	}

	fin.close();

}

template <typename T>
inline void HDFread(hid_t &file, const char *name, hid_t type, T *data)
{
	hid_t dataset;

	dataset = H5Dopen2(file, name, H5P_DEFAULT);
	H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	H5Dclose(dataset);
}