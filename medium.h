#ifndef MEDIUM_H
#define MEDIUM_H

#include "util.h"

class Medium
{
public:
	Medium();
	virtual ~Medium();

	double chi1() const {return mChi1;}

protected:
	double mChi1;
};


class SphereMedium : public Medium
{
public:
	SphereMedium(util::vec_t pos, double r, double chi1);
	~SphereMedium();

	double r() const {return mR;}
	util::vec_t pos() const {return mPos;}

private:
	util::vec_t mPos;
	double mR;
};


#endif // MEDIUM_H
