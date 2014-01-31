/******************************************************************************
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  library for doing interpolation
 * Authors:  Gottfried Mandlburger, Johannes Otepka, Peb Ruswono Aryan
 *
 ******************************************************************************
 * Copyright (c) 2012-2013,  I.P.F., TU Vienna.
  *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ****************************************************************************/
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