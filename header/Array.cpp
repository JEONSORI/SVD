#include "Array.h"

template <typename T>
Array<T>& Array<T>::operator=(const T val)
{
	fill_n(Entry, n, val);
	return *this;
}

template <typename T>
Array<T>& Array<T>::operator=(const T* arr)
{
	memcpy(Entry, arr, n * sizeof(T));
	return *this;
}

template <typename T>
Array<T> Array<T>::operator+(const Array<T>& arr)
{
	Array<T> result = *this;
	result += arr;

	return result;
}

template <typename T>
Array<T> Array<T>::operator-(const Array<T>& arr)
{
	Array<T> result = *this;
	result -= arr;

	return result;
}

template <typename T>
Array<T> Array<T>::operator*(const Array<T>& arr)
{
	Array<T> result = *this;
	result *= arr;

	return result;
}

template <typename T>
Array<T> Array<T>::operator/(const Array<T>& arr)
{
	Array<T> result = *this;
	result /= arr;

	return result;
}

template <typename T>
Array<T>& Array<T>::operator+=(const T val)
{
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n; i++)
		Entry[i] = Entry[i] + val;

	return *this;
}

template <typename T>
Array<T>& Array<T>::operator+=(const Array<T>& arr)
{
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n; i++)
		Entry[i] = Entry[i] + arr.Entry[i];

	return *this;
}

template <typename T>
Array<T>& Array<T>::operator-=(const T val)
{
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n; i++)
		Entry[i] = Entry[i] - val;

	return *this;
}

template <typename T>
Array<T>& Array<T>::operator-=(const Array<T>& arr)
{
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n; i++)
		Entry[i] = Entry[i] - arr.Entry[i];

	return *this;
}

template <typename T>
Array<T>& Array<T>::operator*=(const T val)
{
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n; i++)
		Entry[i] = Entry[i] * val;

	return *this;
}

template <typename T>
Array<T>& Array<T>::operator*=(const Array<T>& arr)
{
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n; i++)
		Entry[i] = Entry[i] * arr.Entry[i];

	return *this;
}

template <typename T>
Array<T>& Array<T>::operator/=(const T val)
{
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n; i++)
		Entry[i] = Entry[i] / val;

	return *this;
}

template <typename T>
Array<T>& Array<T>::operator/=(const Array<T>& arr)
{
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < n; i++)
		Entry[i] = Entry[i] / arr.Entry[i];

	return *this;
}

template class Array<bool>;
template class Array<int>;
template class Array<float>;
template class Array<double>;