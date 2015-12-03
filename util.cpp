#include <sstream>
#include "util.h"

namespace util {
	std::string pretty_size(size_t size)
	{
		std::stringstream pretty_string_stream;

		if(size > 1024*1024*1024*1024UL)
			pretty_string_stream << size / (1024*1024*1024*1024.) << " TiB";
		else if(size > 1024*1024*1024UL)
			pretty_string_stream << size / (1024*1024*1024.) << " GiB";
		else if(size > 1024*1024UL)
			pretty_string_stream << size / (1024*1024.) << " MiB";
		else if(size > 1024UL)
			pretty_string_stream << size / (1024.) << " kiB";
		else
			pretty_string_stream << size << " B";

		return pretty_string_stream.str();
	}

	vec_t operator+(const vec_t& v1, const vec_t& v2)
	{
		vec_t res;
		res.x = v1.x + v2.x;
		res.y = v1.y + v2.y;
		res.z = v1.z + v2.z;

		return res;
	}

	double operator*(const vec_t& v1, const vec_t& v2)
	{
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	}

};
