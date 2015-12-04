#ifndef UTIL_H
#define UTIL_H

#include <cstddef>
#include <string>

namespace util {
	struct vec_t
	{
		double x, y, z;
	};

	vec_t operator+(const util::vec_t& v1, const util::vec_t& v2);
	double operator*(const util::vec_t& v1, const util::vec_t& v2);
	std::ostream& operator<<(std::ostream &str, util::vec_t &arg);

	std::string pretty_size(size_t size);
};

#endif // UTIL_H
