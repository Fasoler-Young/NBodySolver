#pragma once


template<class T>
class Point3
{
public:
	typedef size_t size_type;
	typedef T value_type;
	T x, y, z;
	Point3(void) { x = 0; y = 0; z = 0; }
	Point3(T value) { x = y = z = value; }
	// ћой старый неработающий конструктор
	//Point3(Point3<T>& V) { x = V.x; y = V.y; z = V.z; }
	template<class V>
	Point3(const V& copy) : x((T)copy.x), y((T)copy.y), z((T)copy.z) {}
	template<class V>
	Point3& operator = (const V& copy) { x = (T)copy.x; y = (T)copy.y; z = (T)copy.z; return *this; }
	Point3(const T& x_, const  T& y_, const  T& z_) : x(x_), y(y_), z(z_) {}
	// ¬ектор со случайными значени€ми из диапазона
	Point3(T min, T max) { x = rand_in_range(min, max); y = rand_in_range(min, max); z = 0;/* rand_in_range(min, max);*/ }


	T* data(void) { return &x; }
	const T* data(void) const { return &x; }
	const T& operator [] (int n) const
	{
		return data()[n];
	}
	T& operator [] (int n)
	{
		return data()[n];
	}
	const T& operator [] (size_t n) const
	{
		return data()[n];
	}
	T& operator [] (size_t n)
	{
		return data()[n];
	}



	bool operator == (const Point3 & V) { return x == V.x && y = V.y && z = V.z; }
	bool operator != (const Point3& V) { return !this == V; }

	Point3<T> operator + (const Point3& V)const{
		return Point3(x + V.x, y + V.y, z + V.z);
	}

	Point3<T>& operator += (const Point3& V) {
		x += V.x;
		y += V.y;
		z += V.z;
		return *this;
	}

	Point3<T> operator - (const Point3& V)const {
		return Point3<T>(x - V.x, y - V.y, z - V.z);
	}

	Point3<T>& operator -= (const Point3& V) {
		x -= V.x;
		y -= V.y;
		z -= V.y;
		return *this;
	}

	Point3<T> operator / (const Point3& V) const {
		return Point3(x / V.x, y / V.y, z / V.z);
	}

	Point3<T> operator / (const T V)const {
		return Point3<T>(x / V, y / V, z / V);
	}

	Point3<T>& operator /= (const T V) {
		x /= V;
		y /= V;
		z /= V;
		return *this;
	}

	Point3<T>& operator *= (const T V) {
		x *= V;
		y *= V;
		z *= V;
		return *this;
	}

	Point3<T> operator * (const T v) const {
		return Point3(x * v, y * v, z * v);
	}

	Point3<T> operator * (const Point3& V)const {
		return Point3(x * V.x, y * V.y, z * V.z);
	}

	Point3<T> operator ^ (const Point3<T>& V) const
	{
		return Point3<T>(	y * V.z - z * V.y,
							z * V.x - x * V.z,
							x * V.y - y * V.x);
	}
	//T operator * (const Point3& V)const {
	//	return (x * V.x + y * V.y + z * V.z);
	//}

	T distance(const Point3& V) {
		Point3<T> difference = (*this) - V;
		return difference.length();
	}

	T length(void) {
		return sqrt(this->norm());
	}

	T norm() {
		return x*x + y*y + z*z;
	}

private:
	T rand_in_range(T min, T max) {
		return min + (max - min) * (T) rand() / (T) RAND_MAX;
	}

};
















//
//template<class T>
//struct Point3
//{
//	typedef size_t		size_type;
//	typedef T			value_type;
//	T x, y, z;
//	Point3(void) { x = y = z = 0; }
//	Point3(T value) { x = y = z = value; }
//	template<class V>
//	Point3(const V& copy) : x((T)copy.x), y((T)copy.y), z((T)copy.z) {}
//	template<class V>
//	Point3& operator = (const V& copy) { x = (T)copy.x; y = (T)copy.y; z = (T)copy.z; return *this; }
//	Point3(const T& x_, const  T& y_, const  T& z_) : x(x_), y(y_), z(z_) {}
//
//	T* data(void) { return &x; }
//	const T* data(void) const { return &x; }
//	const T& operator [] (int n) const
//	{
//		return data()[n];
//	}
//	T& operator [] (int n)
//	{
//		return data()[n];
//	}
//	const T& operator [] (size_t n) const
//	{
//		return data()[n];
//	}
//	T& operator [] (size_t n)
//	{
//		return data()[n];
//	}
//	bool operator == (const Point3<T>& V) const
//	{
//		return (x == V.x && y == V.y && z == V.z);
//	}
//	bool operator != (const Point3<T>& V) const
//	{
//		return (x != V.x || y != V.y || z != V.z);
//	}
//
//	//************************************************************************************************
//	// Vector Sum
//	//************************************************************************************************
//
//	Point3<T> operator + (const Point3<T>& V) const
//	{
//		return Point3<T>(x + V.x, y + V.y, z + V.z);
//	}
//	Point3<T>& operator += (const Point3<T>& V)
//	{
//		x += V.x;
//		y += V.y;
//		z += V.z;
//		return *this;
//	}
//	//************************************************************************************************
//	//Vector Sub
//	//************************************************************************************************
//	Point3<T> operator - (const Point3<T>& V) const
//	{
//		return Point3<T>(x - V.x, y - V.y, z - V.z);
//	}
//	Point3<T>& operator -= (const Point3<T>& V)
//	{
//		x -= V.x;
//		y -= V.y;
//		z -= V.z;
//		return *this;
//	}
//	//************************************************************************************************
//	//Mul
//	//************************************************************************************************
//	Point3<T> operator * (const T a) const
//	{
//		return Point3<T>(x * a, y * a, z * a);
//	}
//	Point3<T>& operator *= (const T a)
//	{
//		x *= a;
//		y *= a;
//		z *= a;
//		return *this;
//	}
//	//************************************************************************************************
//	//Div
//	//************************************************************************************************
//	Point3<T> operator / (const T a) const
//
//	{
//		return Point3<T>(x / a, y / a, z / a);
//	}
//	Point3<T>& operator /=(const T a)
//	{
//		x /= a;
//		y /= a;
//		z /= a;
//		return *this;
//	}
//	Point3<T>& operator /= (const Point3<T>& V)
//	{
//		x /= V.x;
//		y /= V.y;
//		z /= V.z;
//		return *this;
//	}
//	//************************************************************************************************
//	//Unary minus
//	//************************************************************************************************
//	Point3<T> operator - (void) const
//	{
//		return Point3<T>(-x, -y, -z);
//	}
//	//************************************************************************************************
//	//Scalar mul
//	//************************************************************************************************
//	T operator * (const Point3<T>& V) const
//	{
//		return (x * V.x + y * V.y + z * V.z);
//	}
//	//************************************************************************************************
//	//Vector mul
//	//************************************************************************************************
//	Point3<T> operator ^ (const Point3<T>& V) const
//	{
//		return Point3<T>(y * V.z - z * V.y,
//			z * V.x - x * V.z,
//			x * V.y - y * V.x);
//	}
//	//************************************************************************************************
//	//	Helpers
//	//************************************************************************************************
//	void normalize(void)
//	{
//		(*this) = (*this) / length();
//	}
//	T distance(const Point3<T>& Vd) const
//	{
//		Point3<T> V = (*this) - Vd;
//		return V.length();
//	}
//	T norm(void) const
//	{
//		return (x * x + y * y + z * z);
//	}
//	T length(void) const
//	{
//		return sqrt((*this) * (*this));
//	}
//	Point3<T> mirror(const Point3<T>& V) const
//	{
//		return (V + ((*this) - V) * ((T)2));
//	}
//};

