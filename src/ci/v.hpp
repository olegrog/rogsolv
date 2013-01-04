#ifndef _V_H_
#define _V_H_

#include <iostream>
#include <vector>
#include <cmath>

namespace dod_vector {

template <typename T> class V2;

template <typename T> class V3 {
	public:
		V3(T y = 0) {
			x[0] = y; x[1] = y; x[2] = y;
		}
		V3(T y1, T y2, T y3) {
			x[0] = y1; x[1] = y2; x[2] = y3;
		}
		V3(V2<T> y, T y2 = 0) {
			x[0] = y[0]; x[1] = y[1]; x[2] = y2;
		}
		template <typename T1> V3(const V3<T1>& y) {
			x[0] = y[0]; x[1] = y[1]; x[2] = y[2];
		}

		template <typename T2> 
		V3<T>& operator=(const V3<T2>& u) {
			if (static_cast<void*>(this) != static_cast<void*>(&u)) {
				x[0] = u[0]; x[1] = u[1]; x[2] = u[2];
			}
			return *this;
		}

		bool operator==(const V3<T>& u) {
			return (x[0] == u[0]) && (x[1] == u[1]) && (x[2] == u[2]);
		}
		bool operator!=(const V3<T>& u) {
			return !(*this == u);
		}

		const V3<T> operator-() const {
			return *this;
		}

		V3<T>& operator+=(const V3<T>& u) {
			x[0] += u[0]; x[1] += u[1]; x[2] += u[2];
			return *this;
		}
		V3<T>& operator-=(const V3<T>& u) {
			x[0] -= u[0]; x[1] -= u[1]; x[2] -= u[2];
			return *this;
		}
		V3<T>& operator*=(const V3<T>& u) {
			x[0] *= u[0]; x[1] *= u[1]; x[2] *= u[2];
			return *this;
		}
		V3<T>& operator/=(const V3<T>& u) {
			x[0] /= u[0]; x[1] /= u[1]; x[2] /= u[2];
			return *this;
		}

		T& operator[](size_t i) { return x[i]; }
		const T operator[](size_t i) const { return x[i]; }

		friend const V3 operator+(const V3& u, const V3& v) {
			return V3(u[0]+v[0], u[1]+v[1], u[2]+v[2]);
		}
		friend const V3 operator-(const V3& u, const V3& v) {
			return V3(u[0]-v[0], u[1]-v[1], u[2]-v[2]);
		}
		friend const V3 operator*(const V3& u, const V3& v) {
			return V3(u[0]*v[0], u[1]*v[1], u[2]*v[2]);
		}
		friend const V3 operator/(const V3& u, const V3& v) {
			return V3(u[0]/v[0], u[1]/v[1], u[2]/v[2]);
		}

		friend bool operator>(const V3& u, const V3& v) {
			return (u[0] > v[0]) && (u[1] > v[1]) && (u[2] > v[2]);
		}
		friend bool operator>=(const V3& u, const V3& v) {
			return (u[0] >= v[0]) && (u[1] >= v[1]) && (u[2] >= v[2]);
		}
		friend bool operator<(const V3& u, const V3& v) {
			return (u[0] < v[0]) && (u[1] < v[1]) && (u[2] < v[2]);
		}
		friend bool operator<=(const V3& u, const V3& v) {
			return (u[0] <= v[0]) && (u[1] <= v[1]) && (u[2] <= v[2]);
		}

	private:
		T x[3];
};

template <typename T> inline T dot(const V3<T>& u, const V3<T>& v) {
	return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
}

template <typename T> inline T sqr(const V3<T>& u) {
	return dot(u, u);
}

template <typename T> inline T norm(const V3<T>& u) {
	return std::sqrt(sqr(u));
}

template <typename T> inline V3<T> cross(const V3<T>& u, const V3<T>& v) {
	return V3<T>(	u[1]*v[2]-u[2]*v[1],
					u[2]*v[0]-u[0]*v[2],
					u[0]*v[1]-u[1]*v[0]);
}

template <typename T> inline V3<T> rotate(const V3<T>& v, const V3<T>& n, double theta) {
	double ct = std::cos(theta);
	double st = std::sin(theta);
	return 	v * ct + n*dot(n, v)*(1 - ct) + cross(n, v)*st;
}

template <typename T> inline double angle(const T& u, const T& v) {
	double r2 = sqr(u)*sqr(v);
	if (r2 != 0.)
		return std::acos(dot(u, v)/std::sqrt(r2));
	else
		return 0.;
}

typedef V3<double> V3d;
typedef V3<int> V3i;

template <class T> inline
std::istream& operator>>(std::istream& from, V3<T>& v) {
	V3<T> u;
	char c = 0;
	from >> c;
	if (c == '(') 
		from >> u[0] >> u[1] >> u[2] >> c;
	else { /* TODO excepton	*/ }
	v = u;
	return from;
}

template <class T> inline
std::ostream& operator<<(std::ostream& to, const V3<T>& v) {
	to << '(' << v[0] << ' ' << v[1] << ' ' << v[2] << ')';
	return to;
}

template <typename T> class V2 {
	public:
		V2(T y = 0) {
			x[0] = y; x[1] = y;
		}
		template <typename T1> V2(T1 y1, T1 y2) {
			x[0] = y1; x[1] = y2;
		}
		template <typename T1> V2(const V2<T1>& y) {
			x[0] = y[0]; x[1] = y[1];
		}

		template <typename T2> 
		V2& operator=(const V2<T2>& u) {
			if (static_cast<void*>(this) != static_cast<void*>(&u)) {
				x[0] = u[0]; x[1] = u[1];
			}
			return *this;
		}

		bool operator==(const V2<T>& u) {
			return (x[0] == u[0]) && (x[1] == u[1]);
		}
		bool operator!=(const V2<T>& u) {
			return !(*this == u);
		}

		V2<T>& operator+=(const V2<T>& u) {
			x[0] += u[0]; x[1] += u[1];
			return *this;
		}
		V2<T>& operator-=(const V2<T>& u) {
			x[0] -= u[0]; x[1] -= u[1];
			return *this;
		}
		V2<T>& operator*=(const V2<T>& u) {
			x[0] *= u[0]; x[1] *= u[1];
			return *this;
		}
		V2<T>& operator/=(const V2<T>& u) {
			x[0] /= u[0]; x[1] /= u[1];
			return *this;
		}

		T& operator[](size_t i) { return x[i]; }
		const T operator[](size_t i) const { return x[i]; }

		friend const V2 operator+(const V2& u, const V2& v) {
			return V2(u[0]+v[0], u[1]+v[1]);
		}
		friend const V2 operator-(const V2& u, const V2& v) {
			return V2(u[0]-v[0], u[1]-v[1]);
		}
		friend const V2 operator*(const V2& u, const V2& v) {
			return V2(u[0]*v[0], u[1]*v[1]);
		}
		friend const V2 operator/(const V2& u, const V2& v) {
			return V2(u[0]/v[0], u[1]/v[1]);
		}

		friend bool operator>(const V2& u, const V2& v) {
			return (u[0] > v[0]) && (u[1] > v[1]);
		}
		friend bool operator>=(const V2& u, const V2& v) {
			return (u[0] >= v[0]) && (u[1] >= v[1]);
		}
		friend bool operator<(const V2& u, const V2& v) {
			return (u[0] < v[0]) && (u[1] < v[1]);
		}
		friend bool operator<=(const V2& u, const V2& v) {
			return (u[0] <= v[0]) && (u[1] <= v[1]);
		}

	private:
		T x[3];
};

template <typename T1, typename T2> inline const T1 dot(const V2<T1>& u, const V2<T2>& v) {
	return u[0]*v[0]+u[1]*v[1];
}

template <typename T> inline const T sqr(const V2<T>& u) {
	return dot(u, u);
}

template <typename T> inline const T norm(const V2<T>& u) {
	return std::sqrt(sqr(u));
}

template <typename T> inline const T arg(const V2<T>& v) {
	double n = norm(v);
	if (v[1] > 0.)
		return std::acos(v[0]/n);
	else
		return -std::acos(v[0]/n);
}

template <typename T> inline const V2<T> rotate(const V2<T>& v, double theta) {
	double ct = std::cos(theta);
	double st = std::sin(theta);
	return V2<T>(	v[0]*ct - v[1]*st	,
					v[0]*st + v[1]*ct	);
}



typedef V2<double> V2d;
typedef V2<int> V2i;

template <class T> inline
std::istream& operator>>(std::istream& from, V2<T>& v) {
	V2<T> u;
	char c = 0;
	from >> c;
	if (c == '(') 
		from >> u[0] >> u[1] >> c;
	else { /* TODO excepton	*/ }
	v = u;
	return from;
}

template <class T> inline
std::ostream& operator<<(std::ostream& to, const V2<T>& v) {
	to << '(' << v[0] << ' ' << v[1] << ')';
	return to;
}

}

#endif
