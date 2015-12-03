#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "solver.h"
#include "util.h"

size_t availableMemory()
{
	std::ifstream meminfo("/proc/meminfo");
	if(meminfo.good())
	{
		char buffer[128];
		while(!meminfo.eof())
		{
			meminfo.getline(buffer, 128);
			std::string line(buffer, 128);
			if(line.substr(0, 8) == "MemFree:")
			{
				std::stringstream strstream(line.substr(9));
				size_t available_memory;
				strstream >> available_memory;

				return available_memory * 1024;
			}
		}
	}

	return 0;
}

int main(int argc, char **argv)
{
	bool override_memory_error = false;
	if(argc > 1)
	{
		if(std::string(argv[1]) == "--disable-memory-check")
			override_memory_error = true;
	}

	Solver solver(100, .1, 100, .1, 100, .1, 500, .05);
//	Solver solver(10, .1, 10, .1, 10, .1, 2, .08);
	
	std::cout << "Memory required : " << util::pretty_size(solver.memoryRequired()) << "." << std::endl;
	std::cout << "Available memory: " << util::pretty_size(availableMemory()) << "." << std::endl;
	if(!override_memory_error && availableMemory() < solver.memoryRequired())
	{
		std::cerr << "Not enough memory available." << std::endl;
		return EXIT_FAILURE;
	}

	system("rm -f *.dat");
	
	solver.setBoundaryConditions();
	solver.allocate();
	solver.fillMedium(0);
	solver.solve();

	return EXIT_SUCCESS;
}
