#pragma once

#include "Defines.h"
// Termination

enum errorCode {
	FILE_NOT_FOUND,
	LABEL_NOT_FOUND,
	INVALID_CARD,
	INVALID_BLOCK,
	INVALID_INTEGER,
	INVALID_FLOATING_POINT,
	INVALID_LOGICAL,
	INVALID_ARGUMENT,
	INVALID_ARGUMENT_COUNT,
	INCONSISTENT_CGROUP_NUMBER,
	INCONSISTENT_FGROUP_NUMBER,
	UNSPECIFIED_GC_FILE,
	UNKNOWN_GC_FORMAT,
	UNSPECIFIED_ERROR = -1
};

static void Error(string message)
{
	cout << "                                                                                      " << endl;
	cout << "=============================== VANGARD Internal Error ===============================" << endl;
	cout << "                                                                                      " << endl;
	cout << message << endl;
	cout << "                                                                                      " << endl;
	cout << "======================================================================================" << endl;
	cout << "                                                                                      " << endl;

	exit(EXIT_FAILURE);
}

static void Error(int code, string info = "", string message = "")
{
	cout << "                                                                                      " << endl;
	cout << "=============================== VANGARD Internal Error ===============================" << endl;
	cout << "                                                                                      " << endl;

	switch (code)
	{
	case errorCode::FILE_NOT_FOUND:
		cout << info;
		cout << ": File not found!" << endl;
		break;
	case errorCode::LABEL_NOT_FOUND:
		cout << info;
		cout << ": Label not found in the input!" << endl;
		break;
	case errorCode::INVALID_CARD:
		cout << info;
		cout << ": Invalid card!" << endl;
		break;
	case errorCode::INVALID_BLOCK:
		cout << info;
		cout << ": Invalid block!" << endl;
		break;
	case errorCode::INVALID_INTEGER:
		cout << info;
		cout << ": Invalid integer expression!" << endl;
		break;
	case errorCode::INVALID_FLOATING_POINT:
		cout << info;
		cout << ": Invalid floating point expression!" << endl;
		break;
	case errorCode::INVALID_LOGICAL:
		cout << info;
		cout << ": Invalid logical expression!" << endl;
		break;
	case errorCode::INVALID_ARGUMENT:
		cout << info;
		cout << ": Invalid argument!" << endl;
		break;
	case errorCode::INVALID_ARGUMENT_COUNT:
		cout << info;
		cout << ": Invalid number of arguments!" << endl;
		break;
	case errorCode::INCONSISTENT_CGROUP_NUMBER:
		cout << info;
		cout << ": Number of condensed group is inconsistent!" << endl;
		break;
	case errorCode::INCONSISTENT_FGROUP_NUMBER:
		cout << info;
		cout << ": Number of fine group is inconsistent!" << endl;
		break;
	case errorCode::UNSPECIFIED_GC_FILE:
		cout << info;
		cout << ": Group constant file is not specified!" << endl;
		break;
	case errorCode::UNKNOWN_GC_FORMAT:
		cout << info;
		cout << ": Group constant file format is unknown!" << endl;
		break;
	case errorCode::UNSPECIFIED_ERROR:
		cout << info;
		cout << ": " << message << endl;
	}

	cout << "                                                                                      " << endl;
	cout << "======================================================================================" << endl;
	cout << "                                                                                      " << endl;

	exit(EXIT_FAILURE);
}

// Filesystem

template <typename T>
inline void HDFwrite(hid_t& loc, const char *name, int dim, hsize_t* dims, hid_t type, T *data, bool overwrite = false)
{
	hid_t dataset;

	if (overwrite) dataset = H5Dopen2(loc, name, H5P_DEFAULT);
	else {
		hid_t dataspace;

		dataspace = H5Screate_simple(dim, dims, NULL);
		dataset = H5Dcreate2(loc, name, type, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		H5Sclose(dataspace);
	}

	H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	H5Dclose(dataset);
}

template <typename T>
inline void HDFread(hid_t &file, const char *name, hid_t type, T *data)
{
	hid_t dataset;

	dataset = H5Dopen2(file, name, H5P_DEFAULT);
	H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	H5Dclose(dataset);
}
