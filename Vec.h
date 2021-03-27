// W.H.Press et al "Numerical Recipes in C++"
//  http://numerical.recipes

#ifndef __Vec_h__
#define __Vec_h__

template <class T>
class Vec {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	Vec();
	explicit Vec(int n);		// Zero-based array
	Vec(const T &a, int n);	//initialize to constant value
	Vec(const T *a, int n);	// Initialize to array
	Vec(const Vec &rhs);	// Copy constructor
	Vec & operator=(const Vec &rhs);	//assignment
	Vec & operator=(const T &a);	//assign a to every element
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	~Vec();
};

template <class T>
Vec<T>::Vec() : nn(0), v(0) {}

template <class T>
Vec<T>::Vec(int n) : nn(n), v(new T[n]) {}

template <class T>
Vec<T>::Vec(const T& a, int n) : nn(n), v(new T[n])
{
	for(int i=0; i<n; i++)
		v[i] = a;
}

template <class T>
Vec<T>::Vec(const T *a, int n) : nn(n), v(new T[n])
{
	for(int i=0; i<n; i++)
		v[i] = *a++;
}

template <class T>
Vec<T>::Vec(const Vec<T> &rhs) : nn(rhs.nn), v(new T[nn])
{
	for(int i=0; i<nn; i++)
		v[i] = rhs[i];
}

template <class T>
Vec<T> & Vec<T>::operator=(const Vec<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			if (v != 0) delete [] (v);
			nn=rhs.nn;
			v= new T[nn];
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

template <class T>
Vec<T> & Vec<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i<nn; i++)
		v[i]=a;
	return *this;
}

template <class T>
inline T & Vec<T>::operator[](const int i)	//subscripting
{
	return v[i];
}

template <class T>
inline const T & Vec<T>::operator[](const int i) const	//subscripting
{
	return v[i];
}

template <class T>
inline int Vec<T>::size() const
{
	return nn;
}

template <class T>
Vec<T>::~Vec()
{
	if (v != 0)
		delete[] (v);
}

template <class T>
class Mat {
private:
	int nn;
	int mm;
	T **v;
public:
	Mat();
	Mat(int n, int m);			// Zero-based array
	Mat(const T &a, int n, int m);	//Initialize to constant
	Mat(const T *a, int n, int m);	// Initialize to array
	Mat(const Mat &rhs);		// Copy constructor
	Mat & operator=(const Mat &rhs);	//assignment
	Mat & operator=(const T &a);		//assign a to every element
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline int nrows() const;
	inline int ncols() const;
	~Mat();
};

template <class T>
Mat<T>::Mat() : nn(0), mm(0), v(0) {}

template <class T>
Mat<T>::Mat(int n, int m) : nn(n), mm(m), v(new T*[n])
{
	v[0] = new T[m*n];
	for (int i=1; i< n; i++)
		v[i] = v[i-1] + m;
}

template <class T>
Mat<T>::Mat(const T &a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
	int i,j;
	v[0] = new T[m*n];
	for (i=1; i< n; i++)
		v[i] = v[i-1] + m;
	for (i=0; i< n; i++)
		for (j=0; j<m; j++)
			v[i][j] = a;
}

template <class T>
Mat<T>::Mat(const T *a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
	int i,j;
	v[0] = new T[m*n];
	for (i=1; i< n; i++)
		v[i] = v[i-1] + m;
	for (i=0; i< n; i++)
		for (j=0; j<m; j++)
			v[i][j] = *a++;
}

template <class T>
Mat<T>::Mat(const Mat &rhs) : nn(rhs.nn), mm(rhs.mm), v(new T*[nn])
{
	int i,j;
	v[0] = new T[mm*nn];
	for (i=1; i< nn; i++)
		v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++)
		for (j=0; j<mm; j++)
			v[i][j] = rhs[i][j];
}

template <class T>
Mat<T> & Mat<T>::operator=(const Mat<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i,j;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != 0) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			v = new T*[nn];
			v[0] = new T[mm*nn];
			for (i=1; i< nn; i++)
				v[i] = v[i-1] + mm;
		}
		for (i=0; i< nn; i++)
			for (j=0; j<mm; j++)
				v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
Mat<T> & Mat<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i< nn; i++)
		for (int j=0; j<mm; j++)
			v[i][j] = a;
	return *this;
}

template <class T>
inline T* Mat<T>::operator[](const int i)	//subscripting: pointer to row i
{
	return v[i];
}

template <class T>
inline const T* Mat<T>::operator[](const int i) const
{
	return v[i];
}

template <class T>
inline int Mat<T>::nrows() const
{
	return nn;
}

template <class T>
inline int Mat<T>::ncols() const
{
	return mm;
}

template <class T>
Mat<T>::~Mat()
{
	if (v != 0) {
		delete[] (v[0]);
		delete[] (v);
	}
}

typedef Vec<double> Vec_DP;
typedef Mat<double> Mat_DP;

#endif // __Vec_h__