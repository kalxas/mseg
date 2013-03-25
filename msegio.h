/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * msegio.h: Image base classes		                                  *
 * Version: 0.9.x                                                         *
 * Last revised: 05/08/2009                                               *
 *                                                                        *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * MSEG algorithm by Angelos Tzotsos (tzotsos@gmail.com)		  *
 * Remote Sensing Lab NTUA - GCpp                         August 2009     *
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


#ifndef MSEGIO_H
#define MSEGIO_H

#define TRUE 1
#define FALSE 0

#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <math.h>
#include <map>
#include <queue>
#include "FreeImage.h"
#include <stdio.h>
#include <string.h>


using namespace std;
//typedef unsigned char CELL_TYPE;
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef signed short sshort;

string find_cell_type(string HeaderName);
int D1(int imglines, int imgcolumns, int line, int column);

    struct Coordinate {
        string Datum;
        string Projection;
        string CoordinateType;
        string Units;
        string Rotation;
    };

    struct Macroblock {
        int size;
        int LineStart;
        int LineEnd;
        int ColumnStart;
        int ColumnEnd;
    };


    struct pixel {
    	int line;
    	int column;
    };

    struct point {
    double x;
    double y;
    double z;
    };

    struct CoPair {
        int first;
        int second;
    };

    typedef vector<CoPair> CoDistribution;

pixel D2(int d, int imgcolumns);

    class Image {											//Virtual class
        public:

        //Image basic components
        int Bands;											//Number of Bands
        int Lines;											//Number of Lines
        int Columns;                                        //Number of Columns
        int Bytes;
        //size_t CELL_TYPE;									//Spectral depth of image data


        string FileName;
        string CellType;

        //Raw data file
        ifstream file;											//or FILE* if necesary

        //Memory Buffer
        vector<uchar> Buffer8;								//in case of 1D array
        vector<ushort> Buffer16;
        vector<sshort> Buffer16s;
        vector<float> Buffer32;
        //CELL_TYPE ***Buffer;									//in case of 3D array
        int BuffLineStart;
        int BuffLineEnd;

        //Band Weights for heterogeneity calculation
        vector<float> BandWeight;

        //Macroblock list
        int MacroSize;
        int MacroLines;
        int MacroColumns;

	vector<struct Macroblock> Macroblocks;

        Image();											//Constructor
        ~Image(); 											//Destructor
    };//Image End

    class ERS_Image : public Image {
        public:

        int NullValue;
        //ERS Header strings
        string HeaderName;
        string ERSVersion;
        string DataType;
        string ByteOrder;
        string LastUpdated;


        //ERS Band names
        vector<string> BandID;

        //ERS Geocoding Information
        Coordinate CoordinateSpace;
        point RegistrationPoint;
	point RegistrationCell;
        double GroundSizeX;
        double GroundSizeY;

        //Data Input Functions
        float Buffer(int line, int column, int band);//Returns one pixel value in 3D representation
        float Buffer(int indx);//Returns one pixel value in 1D representation
        //void FillBuffer(int LineStart, int LineEnd);	//Fills the Buffer with several lines
        //void FillBuffer(Macroblock* macro);				//Fills the buffer with one macroblock
        void FullBuffer();								//Loads the whole image on memory

        void SaveAs(string);							//Saves the image in ERS format
        void SaveAs(string, int);						//Saves the image with null value
        int D1(int line, int column, int band);

        //Constructors and Destructors
        ERS_Image(string, int);		//Constructor accepting the Header Filename as an argument

        ERS_Image(string NAME, int BANDS, int LINES, int COLUMNS, string cell_type);
        ERS_Image(string NAME, int BANDS, int LINES, int COLUMNS, int BYTES);

        	//Generic constructor for output images
        ~ERS_Image();

    };//ERS_Image End

    class Texture_Image : public Image {
        public:

        //Buffer stuff
        vector<int> YBuffer;
        int BuffMin;
        int BuffMax;
        int bitdepth;
        int quant;
        int distance;
        ERS_Image *Img;

        void FullYBuffer();
        int Transform(int p);
        vector<unsigned long> hist;
        //unsigned long hist[256];
        void BuildHistogram();

        //LUT stuff

        //vector<CELL_TYPE> LUT;
        uchar LUT[256];

        void CalculateLUT();

        //co-occurence cube stuff
        vector<CoDistribution> CoCube;

        void BuildCoCube();
        void CheckCoCube();
        bool CoCubeBuilt;

        //Constructors and Destructors
        Texture_Image();            //Constructor
        Texture_Image(ERS_Image& Img, int quantizer, int distance);
        ~Texture_Image();           //Destructor

    };
#endif /* MSEGIO_H */

