#ifndef SOLVER_H
#define SOLVER_H

#include <cstddef>
#include <string>

class Solver
{
public:
	enum EBoundaryCondition
	{
		unknown,
		absorbing,
		reflecting,
		periodic
	};

public:
	Solver(int xsize, double dx, int ysize, double dy, int zsize, double dz, int num_time_steps, double dt);
	~Solver();

	void allocate();
	size_t memoryRequired();
	void fillMedium(double chi1);
	void setBoundaryConditions();
	void solve();

private:
	// The fields at time-step n
	double *mEx_current, *mEy_current, *mEz_current;
	double *mHx_current, *mHy_current, *mHz_current;

	// The fields at time-step n+1
	double *mEx_future, *mEy_future, *mEz_future;
	double *mHx_future, *mHy_future, *mHz_future;

	// The linear medium
	double *mPermittivity;

	// Dimensions of the simulation domain excluding boundary conditions
	int mXsize, mYsize, mZsize, mNumTimeSteps;

	// Dimensions of the simulation domain including boundary conditions
	int mFullXsize, mFullYsize, mFullZsize;

	// Step sizes
	double mDx, mDy, mDz, mDt;
	
	// Boundary conditions
	EBoundaryCondition mBoundaryLowX, mBoundaryHighX, mBoundaryLowY, mBoundaryHighY, mBoundaryLowZ, mBoundaryHighZ;

	int mCurrentTimeStep;

	inline double& Ex(int i, int j, int k, int n);
	inline double& Ey(int i, int j, int k, int n);
	inline double& Ez(int i, int j, int k, int n);
	inline double& Hx(int i, int j, int k, int n);
	inline double& Hy(int i, int j, int k, int n);
	inline double& Hz(int i, int j, int k, int n);
	inline double& epsilon(int i, int j, int k);
	inline double& mu(int, int, int) {static double dummy; dummy = 1; return dummy;}
	void checkBounds(int i, int j, int k, const std::string& field);
    void processData();
    void swapBuffers();
    void applyBoundaryConditionsForE();
    void applyBoundaryConditionsForH();
};

#endif // SOLVER_H
