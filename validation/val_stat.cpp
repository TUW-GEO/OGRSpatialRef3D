#include <iostream>
#include <iomanip>

#include "validate.h"

#include "cpl_conv.h"

void
	SummStat::add(double value)
{
	n += 1;
	min = MIN(min, value);
	max = MAX(max, value);
	sum += value;
	ssq += value*value;
}

void
	SummStat::printout(int shortprec, int longprec)
{
	std::cout << std::setprecision(shortprec);
	std::cout << "\tmin : " << min << std::endl; 
	std::cout << "\tmax : " << max << std::endl; 
	std::cout << "\tavg : " << sum / (double)n << std::endl; 
	if(n>1){
		std::cout << std::setprecision(longprec);
		std::cout << "\tvar : " << ((double)n*ssq - sum*sum) / (double)(n*(n-1)) << std::endl; 
	}
}