#include "medium.h"

Medium::Medium()
{

}

Medium::~Medium()
{

}


SphereMedium::SphereMedium(util::vec_t pos, double r, double chi1)
 : mPos(pos), mR(r)
{
	mChi1 = chi1;
}

SphereMedium::~SphereMedium()
{

}
