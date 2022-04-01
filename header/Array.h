#pragma once

#include "Defines.h"
#include "Util.h"
#include <string.h>
namespace ArrayOp {

	template <typename T>
	T dot(Array<T>& a, Array<T>& b) {
		return inner_product(a(), a() + a.size(), b(), 0);
	}

	template <typename T>
	T norm(Array<T>& a) {
		return sqrt(dot(a, a));
	}

};

template <typename T>
class Array
{
protected:
	T *Entry, *d_Entry;
	int n = 0, nx = 0, ny = 0, nz = 0, nw = 0, nv = 0, nxy = 0, nxyz = 0, nxyzw = 0;
	int dx = 0, dy = 0, dz = 0, dw = 0, dv = 0;
	bool lAlloc = false;
	bool lAlias = false;
	bool lDevice = false;
	bool lStatic = false; /// Whether to invoke destructor
	Array<T> *dev_Handle;

public:
	~Array();
	Array() {};
	Array(const Array<T>& arr);
	explicit Array(int nx, int ny = 0, int nz = 0, int nw = 0, int nv = 0) { Create(nx, ny, nz, nw, nv); }
	void Offset(int x, int y = 0, int z = 0, int w = 0, int v = 0) { dx = x; dy = y; dz = z; dw = w; dv = v; }
	void Create(int nx, int ny = 0, int nz = 0, int nw = 0, int nv = 0);
	void Create(T *ptr, int nx, int ny = 0, int nz = 0, int nw = 0, int nv = 0);
	void CreateStatic(int nx, int ny = 0, int nz = 0, int nw = 0, int nv = 0);
	void CreateDeviceStatic(int nx, int ny = 0, int nz = 0, int nw = 0, int nv = 0);
	void Alloc(int nx, int ny = 0, int nz = 0, int nw = 0, int nv = 0); /// Do not initialize
	void Alias(Array<T>& arr);
	void Destroy();
	void Clear();
	int size() { return n; }
	T* begin() { return Entry; }
	T* end() { return Entry + n; }
	T* host_ptr() { return Entry; }
	T* device_ptr() { return d_Entry; }
	Array<T>* Handle() { return dev_Handle; }
	 T& operator[] (int i) {
#ifdef __CUDA_ARCH__
		return d_Entry[i];
#else
		return Entry[i];
#endif
	}
	T*& operator() () {
#ifdef __CUDA_ARCH__
		return d_Entry;
#else
		return Entry;
#endif
	}
	 T& operator() (int ix) {
#ifdef __CUDA_ARCH__
		return d_Entry[ix - dx];
#else
		return Entry[ix - dx];
#endif
	}
	 T& operator() (int ix, int iy) {
#ifdef __CUDA_ARCH__
		return d_Entry[ix - dx + nx * (iy - dy)];
#else
		return Entry[ix - dx + nx * (iy - dy)];
#endif
	}
	T& operator() (int ix, int iy, int iz) {
#ifdef __CUDA_ARCH__
		return d_Entry[ix - dx + nx * (iy - dy) + nxy * (iz - dz)];
#else
		return Entry[ix - dx + nx * (iy - dy) + nxy * (iz - dz)];
#endif
	}
	 T& operator() (int ix, int iy, int iz, int iw) {
#ifdef __CUDA_ARCH__
		return d_Entry[ix - dx + nx * (iy - dy) + nxy * (iz - dz) + nxyz * (iw - dw)];
#else
		return Entry[ix - dx + nx * (iy - dy) + nxy * (iz - dz) + nxyz * (iw - dw)];
#endif
	}
	T& operator() (int ix, int iy, int iz, int iw, int iv) {
#ifdef __CUDA_ARCH__
		return d_Entry[ix - dx + nx * (iy - dy) + nxy * (iz - dz) + nxyz * (iw - dw) + nxyzw * (iv - dv)];
#else
		return Entry[ix - dx + nx * (iy - dy) + nxy * (iz - dz) + nxyz * (iw - dw) + nxyzw * (iv - dv)];
#endif
	}
	Array<T>& operator= (const T val);
	Array<T>& operator= (const T* arr);
	Array<T>& operator= (const Array<T>& arr);
	Array<T>  operator+ (const Array<T>& arr);
	Array<T>  operator- (const Array<T>& arr);
	Array<T>  operator* (const Array<T>& arr);
	Array<T>  operator/ (const Array<T>& arr);
	Array<T>& operator+=(const T val);
	Array<T>& operator+=(const Array<T>& arr);
	Array<T>& operator-=(const T val);
	Array<T>& operator-=(const Array<T>& arr);
	Array<T>& operator*=(const T val);
	Array<T>& operator*=(const Array<T>& arr);
	Array<T>& operator/=(const T val);
	Array<T>& operator/=(const Array<T>& arr);
};

template <typename T>
inline Array<T>::~Array()
{
	if (!lStatic && lAlloc) Destroy();
}

template <typename T>
inline Array<T>::Array(const Array<T>& arr)
{
	Create(arr.nx, arr.ny, arr.nz, arr.nw, arr.nv);

	memcpy(Entry, arr.Entry, n * sizeof(T));
}

template <typename T>
inline void Array<T>::Create(int nx, int ny, int nz, int nw, int nv)
{
	if (lAlloc) delete[] Entry;

	this->nx = nx;
	this->ny = ny;
	this->nz = nz;
	this->nw = nw;
	this->nv = nv;
	this->nxy = nx * ny;
	this->nxyz = nx * ny * nz;
	this->nxyzw = nx * ny * nz * nw;

	this->n = nx * max(ny, 1) * max(nz, 1) * max(nw, 1) * max(nv, 1);

	Entry = new T[n];
	Clear();

	lAlloc = true;
}

template <typename T>
inline void Array<T>::Create(T *ptr, int nx, int ny, int nz, int nw, int nv)
{
	if (lAlloc) Destroy();

	this->nx = nx;
	this->ny = ny;
	this->nz = nz;
	this->nw = nw;
	this->nv = nv;
	this->nxy = nx * ny;
	this->nxyz = nx * ny * nz;
	this->nxyzw = nx * ny * nz * nw;

	this->n = nx * max(ny, 1) * max(nz, 1) * max(nw, 1) * max(nv, 1);

	Entry = ptr;

	lAlias = true;
}


template <typename T>
inline void Array<T>::CreateStatic(int nx, int ny, int nz, int nw, int nv)
{
	Create(nx, ny, nz, nw, nv);

	lStatic = true;
}

template <typename T>
inline void Array<T>::Alloc(int nx, int ny, int nz, int nw, int nv)
{
	if (lAlloc) delete[] Entry;

	this->nx = nx;
	this->ny = ny;
	this->nz = nz;
	this->nw = nw;
	this->nv = nv;
	this->nxy = nx * ny;
	this->nxyz = nx * ny * nz;
	this->nxyzw = nx * ny * nz * nw;

	this->n = nx * max(ny, 1) * max(nz, 1) * max(nw, 1) * max(nv, 1);

	Entry = new T[n];

	lAlloc = true;
}

template <typename T>
inline void Array<T>::Alias(Array<T>& arr)
{
	if (lAlloc) Destroy();

	this->nx = arr.nx;
	this->ny = arr.ny;
	this->nz = arr.nz;
	this->nw = arr.nw;
	this->nv = arr.nv;

	this->nxy = arr.nxy;
	this->nxyz = arr.nxyz;
	this->nxyzw = arr.nxyzw;

	this->n = arr.n;

	this->dx = arr.dx;
	this->dy = arr.dy;
	this->dz = arr.dz;
	this->dw = arr.dw;
	this->dv = arr.dv;

	Entry = arr.Entry;
	d_Entry = arr.d_Entry;
	dev_Handle = arr.dev_Handle;

	lAlias = true;
	lDevice = arr.lDevice;
}

template <typename T>
inline void Array<T>::Destroy()
{
	if (!lAlias) delete[] Entry;
	else		 Entry = NULL;

	if (!lAlias && lDevice) {
	}
	else {
		d_Entry = NULL;
		dev_Handle = NULL;
	}

	lAlloc = false;
	lDevice = false;
	lStatic = false;
}

template <typename T>
inline void Array<T>::Clear()
{
	memset(Entry, 0, n * sizeof(T));
}

template <typename T>
Array<T>& Array<T>::operator=(const Array<T>& arr)
{
	if (!lAlloc && !lAlias) Create(arr.nx, arr.ny, arr.nz, arr.nw, arr.nv);

	memcpy(Entry, arr.Entry, n * sizeof(T));
	return *this;
}