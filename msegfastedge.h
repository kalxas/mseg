/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * msegstd.h: MSEG Fast mode class with edge compensation                 *
 * Version: 0.9.x                                                         *
 * Last revised: 31/03/2010                                               *
 *                                                                        *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * MSEG algorithm by Angelos Tzotsos (tzotsos@gmail.com)		  *
 * Remote Sensing Lab NTUA - GCpp                         March 2010      *
 *									  *
 *   Copyright (C) Angelos Tzotsos <tzotsos@gmail.com>	  		  *
 *									  *
 * This program is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU General Public License as published by   *
 * the Free Software Foundation, either version 2 of the License, or      *
 * (at your option) any later version.					  *
 *									  *
 * This program is distributed in the hope that it will be useful,	  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 	  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  *
 * GNU General Public License for more details.				  *
 *									  *
 * You should have received a copy of the GNU General Public License	  *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *									  *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * */


#ifndef MSEGFASTEDGE_H
#define MSEGFASTEDGE_H

#include <cmath>
#include "msegcore.h"

	class FastEdge_Mode{
	    public:

     	int firstpass(Level& pass1, ERS_Image& Img, ERS_Image& ImgEdge, MSEG_Init& mseg, float level_ID, MSEG_Param& Param);
     	int firstpass(Level& pass1, ERS_Image& Img, ERS_Image& ImgEdge, MSEG_Init& mseg, float level_ID, MSEG_Param& Param, Level& super_level);
     	int secondpass(Level& pass2, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass);
     	int secondpass(Level& pass2, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass, Level& super_level);
     	int nthpass(Level& passn, ERS_Image& Img, float level_ID, MSEG_Init& mseg, MSEG_Param& Param, Level& sub_level);
     	int nthpass(Level& passn, ERS_Image& Img, float level_ID, MSEG_Init& mseg, MSEG_Param& Param, Level& sub_level, Level& super_level);
     	int nthpass(Level& passn, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass);
     	int nthpass(Level& passn, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass, Level& super_level);
     	int lastpass(Level& passn, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass);
     	int lastpass(Level& passn, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass, Level& super_level);
	};

#endif /* MSEGFASTEDGE_H */
