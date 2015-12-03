#include <cmath>
#include <stdexcept>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <cassert>
#include "solver.h"
#include "medium.h"

// For readability in the formulas.
#define PLUS_HALF	+1
#define MINUS_HALF

/*
 * Yee cell:
 * 
 *     i integer          i half-integer
 * 
 *    o   Hy   o          Hx  Ez  Hx
 *    Hz  Ex  Hz          Ey   o  Ey           etc.
 *    o   Hy   o          Hx  Ez  Hx
 */

Solver::Solver(int xsize, double dx, int ysize, double dy, int zsize, double dz, int num_time_steps, double dt)
 : mXsize(xsize), mYsize(ysize), mZsize(zsize), mNumTimeSteps(num_time_steps), mDx(dx), mDy(dy), mDz(dz), mDt(dt), mCurrentTimeStep(0)
{
	mEx_current = NULL;
	mEy_current = NULL;
	mEz_current = NULL;
	mHx_current = NULL;
	mHy_current = NULL;
	mHz_current = NULL;
	mEx_future = NULL;
	mEy_future = NULL;
	mEz_future = NULL;
	mHx_future = NULL;
	mHy_future = NULL;
	mHz_future = NULL;
	mPermittivity = NULL;

	mBoundaryLowX = unknown;
	mBoundaryHighX = unknown;
	mBoundaryLowY = unknown;
	mBoundaryHighY = unknown;
	mBoundaryLowZ = unknown;
	mBoundaryHighZ = unknown;

	mFullXsize = 0;
	mFullYsize = 0;
	mFullZsize = 0;
}

Solver::~Solver()
{
	delete[] mEx_current;
	delete[] mEy_current;
	delete[] mEz_current;
	delete[] mHx_current;
	delete[] mHy_current;
	delete[] mHz_current;
	delete[] mEx_future;
	delete[] mEy_future;
	delete[] mEz_future;
	delete[] mHx_future;
	delete[] mHy_future;
	delete[] mHz_future;
	delete[] mPermittivity;
}

size_t Solver::memoryRequired()
{
	return
	(
		2 * ((size_t)
			(mXsize+2) * (mYsize+2) * (mZsize+2) +
			(mXsize+2) * (mYsize+2) * (mZsize+2) +
			(mXsize+2) * (mYsize+2) * (mZsize+2) +
			(mXsize+2) * (mYsize+2) * (mZsize+2) +
			(mXsize+2) * (mYsize+2) * (mZsize+2) +
			(mXsize+2) * (mYsize+2) * (mZsize+2) 
			) +
		(mXsize+2) * (mYsize+2) * (mZsize+2)
	) * sizeof(double);
}

void Solver::allocate()
{
	if(mFullXsize == 0)
		throw std::length_error("mFullXsize == 0");
	if(mFullYsize == 0)
		throw std::length_error("mFullYsize == 0");
	if(mFullZsize == 0)
		throw std::length_error("mFullZsize == 0");

	size_t num_grid_cells = (mXsize+2) * (mYsize+2) * (mZsize+2);
	mEx_current   = new double[num_grid_cells];
	mEy_current   = new double[num_grid_cells];
	mEz_current   = new double[num_grid_cells];
	mHx_current   = new double[num_grid_cells];
	mHy_current   = new double[num_grid_cells];
	mHz_current   = new double[num_grid_cells];

	mEx_future    = new double[num_grid_cells];
	mEy_future    = new double[num_grid_cells];
	mEz_future    = new double[num_grid_cells];
	mHx_future    = new double[num_grid_cells];
	mHy_future    = new double[num_grid_cells];
	mHz_future    = new double[num_grid_cells];

	mPermittivity = new double[num_grid_cells];

	std::memset(mEx_current, 0, sizeof(double) * num_grid_cells);
	std::memset(mEy_current, 0, sizeof(double) * num_grid_cells);
	std::memset(mEz_current, 0, sizeof(double) * num_grid_cells);
	std::memset(mHx_current, 0, sizeof(double) * num_grid_cells);
	std::memset(mHy_current, 0, sizeof(double) * num_grid_cells);
	std::memset(mHz_current, 0, sizeof(double) * num_grid_cells);

	std::memset(mEx_future, 0, sizeof(double) * num_grid_cells);
	std::memset(mEy_future, 0, sizeof(double) * num_grid_cells);
	std::memset(mEz_future, 0, sizeof(double) * num_grid_cells);
	std::memset(mHx_future, 0, sizeof(double) * num_grid_cells);
	std::memset(mHy_future, 0, sizeof(double) * num_grid_cells);
	std::memset(mHz_future, 0, sizeof(double) * num_grid_cells);

	const double wavenumber = 5 * (2 * M_PI) / 10;
	for(int i = 0; i < mXsize+2; i++)
	{
		for(int j = 0; j < mYsize+2; j++)
		{
			for(int k = 0; k < mZsize+2; k++)
			{
				Ex(i, j, k, 0) = std::cos(wavenumber * k*mDz);
			}
		}
	}
	for(int i = 0; i < mXsize+2; i++)
	{
		for(int j = 0; j < mYsize+2; j++)
		{
			for(int k = 0; k < mZsize+2; k++)
			{
				Hy(i, j, k, 0) = std::cos(wavenumber * k*mDz);
			}
		}
	}
}

void Solver::fillMedium(double chi1)
{
	for(int i = 0; i < mXsize+2; ++i)
	{
		for(int j = 0; j < mYsize+2; ++j)
		{
			for(int k = 0; k < mZsize+2; ++k)
			{
				epsilon(i, j, k) = 1 + chi1;
			}
		}
	}

	SphereMedium sphere({0, 0, 0}, 10, 2);

	// Sphere position in grid points with origin in lowest corner.
	util::vec_t spherepos = {
		sphere.pos().x / mDx + mXsize / 2,
		sphere.pos().y / mDy + mYsize / 2,
		sphere.pos().z / mDz + mZsize / 2
	};

	for(int i = spherepos.x - sphere.r() - 1; i <= spherepos.x + sphere.r() + 1; ++i)
	{
		for(int j = spherepos.y - sphere.r() - 1; j <= spherepos.y + sphere.r() + 1; ++j)
		{
			for(int k = spherepos.z - sphere.r() - 1; k <= spherepos.z + sphere.r() + 1; ++k)
			{
				const int num_subpixels = 3;
				int subpixels_inside_sphere = 0;

				util::vec_t pos = {
					(i - mXsize / 2) * mDx - spherepos.x,
					(i - mYsize / 2) * mDy - spherepos.y,
					(i - mZsize / 2) * mDz - spherepos.z
				};
				// Use subpixels to avoid hard boundaries.
				util::vec_t subpixel_offset;
				for(int i1 = 0; i1 < 3; ++i1)
				{
					subpixel_offset.x = i1 * mDx / num_subpixels;
					for(int i2 = 0; i2 < 3; ++i2)
					{
						subpixel_offset.y = i2 * mDy / num_subpixels;
						for(int i3 = 0; i3 < 3; ++i3)
						{
							subpixel_offset.z = i3 * mDz / num_subpixels;

							if((pos + subpixel_offset) * (pos + subpixel_offset) <= sphere.r())
								subpixels_inside_sphere++;
						}
					}
				}
				
				epsilon(i, j, k) = 1 + sphere.chi1();// * subpixels_inside_sphere / (num_subpixels*num_subpixels*num_subpixels);
			}
		}
	}
	/*
	for(int i = 3*mXsize/8+1; i < 5*mXsize/8+1; ++i)
		for(int j = 3*mYsize/8+1; j < 5*mYsize/8+1; ++j)
			for(int k = 3*mZsize/8+1; k < 5*mZsize/8+1; ++k)
			{
				if((i-mXsize/2)*(i-mXsize/2) + (j-mYsize/2)*(j-mYsize/2) + (k-mZsize/2)*(k-mZsize/2) <= (mXsize/8)*(mXsize/8))
					epsilon(i, j, k) = 1 + 2;
			}
*/
	std::ofstream eps("epsilon.dat");
	for(int i = 0; i < (mXsize+2)*(mYsize+2)*(mZsize+2); ++i)
		eps << mPermittivity[i] << '\t';
	throw 2;
}

void Solver::setBoundaryConditions()
{
	mBoundaryLowX = periodic;
	mBoundaryHighX = periodic;
	mBoundaryLowY = periodic;
	mBoundaryHighY = periodic;
	mBoundaryLowZ = periodic;
	mBoundaryHighZ = periodic;

	if(mBoundaryLowX == periodic || mBoundaryHighX == periodic)
	{
		if(mBoundaryLowX != periodic || mBoundaryHighX != periodic)
			throw std::range_error("x boundaries conflict.");

		mFullXsize = mXsize + 1;
	}
	if(mBoundaryLowY == periodic || mBoundaryHighY == periodic)
	{
		if(mBoundaryLowY != periodic || mBoundaryHighY != periodic)
			throw std::range_error("y boundaries conflict.");

		mFullYsize = mYsize + 1;
	}
	if(mBoundaryLowZ == periodic || mBoundaryHighZ == periodic)
	{
		if(mBoundaryLowZ != periodic || mBoundaryHighZ != periodic)
			throw std::range_error("z boundaries conflict.");

		mFullZsize = mZsize + 1;
	}
}

void Solver::solve()
{
	if(mBoundaryLowX == unknown)
		throw std::range_error("Low x boundary is unknown.");
	if(mBoundaryHighX == unknown)
		throw std::range_error("High x boundary is unknown.");
	if(mBoundaryLowY == unknown)
		throw std::range_error("Low y boundary is unknown.");
	if(mBoundaryHighY == unknown)
		throw std::range_error("High y boundary is unknown.");
	if(mBoundaryLowZ == unknown)
		throw std::range_error("Low z boundary is unknown.");
	if(mBoundaryHighZ == unknown)
		throw std::range_error("High z boundary is unknown.");

	for(int n = 0; n < mNumTimeSteps; ++n)
	{
		mCurrentTimeStep = n;

		std::cout << "t = " << (mCurrentTimeStep * mDt) << std::endl;

		applyBoundaryConditionsForE();

		// Evaluate the Hx-fields.
#pragma omp parallel for
		for(int i = 1; i < mXsize+1; ++i)
		{
			for(int j = 1; j < mYsize+1; ++j)
			{
				for(int k = 1; k < mZsize+1; ++k)
				{
					Hx(i, j, k, n + 1) = Hx(i, j, k, n) + (mDt / mu(i, j, k)) * (
						(Ey(i, j, k, n) - Ey(i, j, k - 1, n)) / mDz -
						(Ez(i, j, k, n) - Ez(i, j - 1, k, n)) / mDy
						);
				}
			}
		}

		// Evaluate the Hy-fields.
#pragma omp parallel for
		for(int i = 1; i < mXsize+1; ++i)
		{
			for(int j = 1; j < mYsize+1; ++j)
			{
				for(int k = 1; k < mZsize+1; ++k)
				{
					Hy(i, j, k, n + 1) = Hy(i, j, k, n) + (mDt / mu(i, j, k)) * (
						(Ez(i, j, k, n) - Ez(i - 1, j, k, n)) / mDx -
						(Ex(i, j, k, n) - Ex(i, j, k - 1, n)) / mDz
						);
				}
			}
		}

		// Evaluate the Hz-fields.
#pragma omp parallel for
		for(int i = 1; i < mXsize+1; ++i)
		{
			for(int j = 1; j < mYsize+1; ++j)
			{
				for(int k = 1; k < mZsize+1; ++k)
				{
					Hz(i, j, k, n + 1) = Hz(i, j, k, n) + (mDt / mu(i, j, k)) * (
						(Ex(i, j, k, n) - Ex(i, j - 1, k, n)) / mDy -
						(Ey(i, j, k, n) - Ey(i - 1, j, k, n)) / mDx
						);
				}
			}
		}

		applyBoundaryConditionsForH();

		// Evaluate the Ex-fields.
#pragma omp parallel for
		for(int i = 1; i < mXsize+1; ++i)
		{
			for(int j = 1; j < mYsize+1; ++j)
			{
				for(int k = 1; k < mZsize+1; ++k)
				{
					Ex(i, j, k, n+1) = Ex(i, j, k, n) + (mDt / epsilon(i, j, k)) * (
						(Hz(i, j + 1, k, n + 1) - Hz(i, j, k, n + 1)) / mDy -
						(Hy(i, j, k + 1, n + 1) - Hy(i, j, k, n + 1)) / mDz
						);

					if(std::isnan(Ex(i, j, k, n+1)))
					{
						std::cerr << "(" << i << ", " << j << ", " << k << ")" << std::endl;
					}
				}
			}
		}

		// Evaluate the Ey-fields.
#pragma omp parallel for
		for(int i = 1; i < mXsize+1; ++i)
		{
			for(int j = 1; j < mYsize+1; ++j)
			{
				for(int k = 1; k < mZsize; ++k)
				{
					Ey(i, j, k, n+1) = Ey(i, j, k, n) + (mDt / epsilon(i, j, k)) * (
						(Hx(i, j, k + 1, n + 1) - Hx(i, j, k, n + 1)) / mDz -
						(Hz(i + 1, j, k, n + 1) - Hz(i, j, k, n + 1)) / mDx
						);
				}
			}
		}

		// Evaluate the Ez-fields.
#pragma omp parallel for
		for(int i = 1; i < mXsize+1; ++i)
		{
			for(int j = 1; j < mYsize+1; ++j)
			{
				for(int k = 1; k < mZsize+1; ++k)
				{
					Ez(i, j, k, n+1) = Ez(i, j, k, n) + (mDt / epsilon(i, j, k)) * (
						(Hy(i + 1, j, k, n + 1) - Hy(i, j, k, n + 1)) / mDx -
						(Hx(i, j + 1, k, n + 1) - Hx(i, j, k, n + 1)) / mDy
						);
				}
			}
		}

		processData();
		swapBuffers();
	}
}

void Solver::processData()
{
	if(mCurrentTimeStep % 10 != 0)
		return;

	std::stringstream strstream;
	strstream << "Ex_";
	strstream << std::setw(5) << std::setfill('0') << (mCurrentTimeStep+1);
	strstream << ".dat";
	std::ofstream outEx(strstream.str().c_str());

	strstream.str("");
	strstream.clear();
	strstream << "Ey_";
	strstream << std::setw(5) << std::setfill('0') << (mCurrentTimeStep+1);
	strstream << ".dat";
	std::ofstream outEy(strstream.str().c_str());

	strstream.str("");
	strstream.clear();
	strstream << "Ez_";
	strstream << std::setw(5) << std::setfill('0') << (mCurrentTimeStep+1);
	strstream << ".dat";
	std::ofstream outEz(strstream.str().c_str());

	strstream.str("");
	strstream.clear();
	strstream << "Hx_";
	strstream << std::setw(5) << std::setfill('0') << (mCurrentTimeStep+1);
	strstream << ".dat";
	std::ofstream outHx(strstream.str().c_str());

	strstream.str("");
	strstream.clear();
	strstream << "Hy_";
	strstream << std::setw(5) << std::setfill('0') << (mCurrentTimeStep+1);
	strstream << ".dat";
	std::ofstream outHy(strstream.str().c_str());

	strstream.str("");
	strstream.clear();
	strstream << "Hz_";
	strstream << std::setw(5) << std::setfill('0') << (mCurrentTimeStep+1);
	strstream << ".dat";
	std::ofstream outHz(strstream.str().c_str());

	if(outEx.good() && outEy.good() && outEz.good() && outHx.good() && outHy.good() && outHz.good())
	{
		outEx.precision(16);
		outHz.precision(16);
		for(int i = 0; i < mXsize+2; ++i)
		{
			for(int j = 0; j < mYsize+2; ++j)
			{
				for(int k = 0; k < mZsize+2; ++k)
				{
					outEx << Ex(i, j, k, mCurrentTimeStep+1) << '\n';
					outEy << Ey(i, j, k, mCurrentTimeStep+1) << '\n';
					outEz << Ez(i, j, k, mCurrentTimeStep+1) << '\n';
					outHx << Hx(i, j, k, mCurrentTimeStep+1) << '\n';
					outHy << Hy(i, j, k, mCurrentTimeStep+1) << '\n';
					outHz << Hz(i, j, k, mCurrentTimeStep+1) << '\n';
				}
			}
		}
	}
	else
		std::cerr << "Couldn't write to files" << std::endl;
}

void Solver::swapBuffers()
{
	double *tmp;

	tmp = mEx_current;
	mEx_current = mEx_future;
	mEx_future = tmp;

	tmp = mEy_current;
	mEy_current = mEy_future;
	mEy_future = tmp;

	tmp = mEz_current;
	mEz_current = mEz_future;
	mEz_future = tmp;

	tmp = mHx_current;
	mHx_current = mHx_future;
	mHx_future = tmp;

	tmp = mHy_current;
	mHy_current = mHy_future;
	mHy_future = tmp;

	tmp = mHz_current;
	mHz_current = mHz_future;
	mHz_future = tmp;
}




void Solver::applyBoundaryConditionsForE()
{
	if(mBoundaryHighX == periodic)
	{
		for(int j = 0; j < mYsize+2; ++j)
		{
			for(int k = 0; k < mZsize+2; ++k)
			{
				Ex(mXsize + 1, j, k, mCurrentTimeStep) = Ex(1     , j, k, mCurrentTimeStep);
				Ex(0         , j, k, mCurrentTimeStep) = Ex(mXsize, j, k, mCurrentTimeStep);
				Ey(mXsize + 1, j, k, mCurrentTimeStep) = Ey(1     , j, k, mCurrentTimeStep);
				Ey(0         , j, k, mCurrentTimeStep) = Ey(mXsize, j, k, mCurrentTimeStep);
				Ez(mXsize + 1, j, k, mCurrentTimeStep) = Ez(1     , j, k, mCurrentTimeStep);
				Ez(0         , j, k, mCurrentTimeStep) = Ez(mXsize, j, k, mCurrentTimeStep);
			}
		}
	}

	if(mBoundaryHighY == periodic)
	{
		for(int i = 0; i < mXsize+2; ++i)
		{
			for(int k = 0; k < mZsize+2; ++k)
			{
				Ex(i, mYsize + 1, k, mCurrentTimeStep) = Ex(i, 1     , k, mCurrentTimeStep);
				Ex(i, 0         , k, mCurrentTimeStep) = Ex(i, mYsize, k, mCurrentTimeStep);
				Ey(i, mYsize + 1, k, mCurrentTimeStep) = Ey(i, 1     , k, mCurrentTimeStep);
				Ey(i, 0         , k, mCurrentTimeStep) = Ey(i, mYsize, k, mCurrentTimeStep);
				Ez(i, mYsize + 1, k, mCurrentTimeStep) = Ez(i, 1     , k, mCurrentTimeStep);
				Ez(i, 0         , k, mCurrentTimeStep) = Ez(i, mYsize, k, mCurrentTimeStep);
			}
		}
	}

	if(mBoundaryHighZ == periodic)
	{
		for(int i = 0; i < mXsize+2; ++i)
		{
			for(int j = 0; j < mYsize+2; ++j)
			{
				Ex(i, j, mZsize + 1, mCurrentTimeStep) = Ex(i, j, 1     , mCurrentTimeStep);
				Ex(i, j, 0         , mCurrentTimeStep) = Ex(i, j, mZsize, mCurrentTimeStep);
				Ey(i, j, mZsize + 1, mCurrentTimeStep) = Ey(i, j, 1     , mCurrentTimeStep);
				Ey(i, j, 0         , mCurrentTimeStep) = Ey(i, j, mZsize, mCurrentTimeStep);
				Ez(i, j, mZsize + 1, mCurrentTimeStep) = Ez(i, j, 1     , mCurrentTimeStep);
				Ez(i, j, 0         , mCurrentTimeStep) = Ez(i, j, mZsize, mCurrentTimeStep);
			}
		}
	}


}

void Solver::applyBoundaryConditionsForH()
{
	if(mBoundaryHighX == periodic)
	{
		for(int j = 0; j < mYsize+2; ++j)
		{
			for(int k = 0; k < mZsize+2; ++k)
			{
				Hx(mXsize + 1, j, k, mCurrentTimeStep+1) = Hx(1     , j, k, mCurrentTimeStep+1);
				Hx(0         , j, k, mCurrentTimeStep+1) = Hx(mXsize, j, k, mCurrentTimeStep+1);
				Hy(mXsize + 1, j, k, mCurrentTimeStep+1) = Hy(1     , j, k, mCurrentTimeStep+1);
				Hy(0         , j, k, mCurrentTimeStep+1) = Hy(mXsize, j, k, mCurrentTimeStep+1);
				Hz(mXsize + 1, j, k, mCurrentTimeStep+1) = Hz(1     , j, k, mCurrentTimeStep+1);
				Hz(0         , j, k, mCurrentTimeStep+1) = Hz(mXsize, j, k, mCurrentTimeStep+1);
			}
		}
	}

	if(mBoundaryHighY == periodic)
	{
		for(int i = 0; i < mXsize+2; ++i)
		{
			for(int k = 0; k < mZsize+2; ++k)
			{
				Hx(i, mYsize + 1, k, mCurrentTimeStep+1) = Hx(i, 1     , k, mCurrentTimeStep+1);
				Hx(i, 0         , k, mCurrentTimeStep+1) = Hx(i, mYsize, k, mCurrentTimeStep+1);
				Hy(i, mYsize + 1, k, mCurrentTimeStep+1) = Hy(i, 1     , k, mCurrentTimeStep+1);
				Hy(i, 0         , k, mCurrentTimeStep+1) = Hy(i, mYsize, k, mCurrentTimeStep+1);
				Hz(i, mYsize + 1, k, mCurrentTimeStep+1) = Hz(i, 1     , k, mCurrentTimeStep+1);
				Hz(i, 0         , k, mCurrentTimeStep+1) = Hz(i, mYsize, k, mCurrentTimeStep+1);
			}
		}
	}

	if(mBoundaryHighZ == periodic)
	{
		for(int i = 0; i < mXsize+2; ++i)
		{
			for(int j = 0; j < mYsize+2; ++j)
			{
				Hx(i, j, mZsize + 1, mCurrentTimeStep+1) = Hx(i, j, 1     , mCurrentTimeStep+1);
				Hx(i, j, 0         , mCurrentTimeStep+1) = Hx(i, j, mZsize, mCurrentTimeStep+1);
				Hy(i, j, mZsize + 1, mCurrentTimeStep+1) = Hy(i, j, 1     , mCurrentTimeStep+1);
				Hy(i, j, 0         , mCurrentTimeStep+1) = Hy(i, j, mZsize, mCurrentTimeStep+1);
				Hz(i, j, mZsize + 1, mCurrentTimeStep+1) = Hz(i, j, 1     , mCurrentTimeStep+1);
				Hz(i, j, 0         , mCurrentTimeStep+1) = Hz(i, j, mZsize, mCurrentTimeStep+1);
			}
		}
	}
}

void Solver::checkBounds(int i, int j, int k, const std::string &field)
{
	if(i < 0)
		throw std::out_of_range(field + " space (i < 0)");
	if(i >= mXsize+2)
		throw std::out_of_range(field + " space (i >= mXsize+2)");
	if(j < 0)
		throw std::out_of_range(field + " space (j < 0)");
	if(j >= mYsize+2)
		throw std::out_of_range(field + " space (j >= mYsize+2)");
	if(k < 0)
		throw std::out_of_range(field + " space (k < 0)");
	if(k >= mZsize+2)
		throw std::out_of_range(field + " space (k >= mZsize+2)");
}

double &Solver::Ex(int i, int j, int k, int n)
{
	checkBounds(i, j, k, "Ex");

	if(n == mCurrentTimeStep)
		return mEx_current[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
	else if(n == mCurrentTimeStep + 1)
		return mEx_future[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
	else
		throw std::out_of_range("Ex time");
}

double &Solver::Ey(int i, int j, int k, int n)
{
	if(i >= 0 && i < mXsize+2 && j >= 0 && j < mYsize+2 && k >= 0 && k < mZsize+2)
	{
		if(n == mCurrentTimeStep)
			return mEy_current[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
		else if(n == mCurrentTimeStep + 1)
			return mEy_future[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
	}

	static double dummy;
	dummy = 0;
	return dummy;
}

double &Solver::Ez(int i, int j, int k, int n)
{
	if(i >= 0 && i < mXsize+2 && j >= 0 && j < mYsize+2 && k >= 0 && k < mZsize+2)
	{
		if(n == mCurrentTimeStep)
			return mEz_current[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
		else if(n == mCurrentTimeStep + 1)
			return mEz_future[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
	}

	static double dummy;
	dummy = 0;
	return dummy;
}

double &Solver::Hx(int i, int j, int k, int n)
{
	if(i >= 0 && i < mXsize+2 && j >= 0 && j < mYsize+2 && k >= 0 && k < mZsize+2)
	{
		if(n == mCurrentTimeStep)
			return mHx_current[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
		else if(n == mCurrentTimeStep + 1)
			return mHx_future[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
	}

	static double dummy;
	dummy = 0;
	return dummy;
}

double &Solver::Hy(int i, int j, int k, int n)
{
	if(i >= 0 && i < mXsize+2 && j >= 0 && j < mYsize+2 && k >= 0 && k < mZsize+2)
	{
		if(n == mCurrentTimeStep)
			return mHy_current[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
		else if(n == mCurrentTimeStep + 1)
			return mHy_future[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
	}

	static double dummy;
	dummy = 0;
	return dummy;
}

double &Solver::Hz(int i, int j, int k, int n)
{
	if(i >= 0 && i < mXsize+2 && j >= 0 && j < mYsize+2 && k >= 0 && k < mZsize+2)
	{
		if(n == mCurrentTimeStep)
			return mHz_current[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
		else if(n == mCurrentTimeStep + 1)
			return mHz_future[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
	}

	static double dummy;
	dummy = 0;
	return dummy;
}

double &Solver::epsilon(int i, int j, int k)
{
	checkBounds(i, j, k, "epsilon");

	return mPermittivity[k + (mZsize+2) * j + (mZsize+2) *(mYsize+2) * i];
}
