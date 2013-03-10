/* 
 *  Copyright (c) 2008-2010  Noah Snavely (snavely (at) cs.cornell.edu)
 *    and the University of Washington
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 */

/* ParameterBound.h */

#ifndef __parameter_bound_h__
#define __parameter_bound_h__

class ParameterBound 
{
public:
    ParameterBound() { }
    ParameterBound(float min_x, float min_y, float max_x, float max_y) :
	m_min_x(min_x), m_min_y(min_y), m_max_x(max_x), m_max_y(max_y)
    { }
    
    float m_min_x, m_min_y;
    float m_max_x, m_max_y;
};

#endif /* __parameter_bound_h__ */
