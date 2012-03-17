/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * libmseg.h: Multiscale SEGmentation library                             *
 * Version: 0.5.0 Final                                                   *
 * Last revised: 15/02/2005                                               *
 *                                                                        *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * MSEG algorithm by Angelos Tzotsos (tzotsos@gmail.com)		  *
 * Remote Sensing Lab NTUA - GCpp                         February 2005   *
 *									  *
 *   Copyright (C) February 2005 Angelos Tzotsos <tzotsos@gmail.com>	  *
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

#ifndef MSEG_H
#define MSEG_H

#define TRUE 1
#define FALSE 0

#include <stdio.h>
#include <string.h>
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


using namespace std;
typedef unsigned char CELL_TYPE;

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

pixel D2(int d, int imgcolumns);
    
    class Image {											//Virtual class
        public:
        
        //Image basic components
        int Bands;											//Number of Bands
        int Lines;											//Number of Lines
        int Columns;										//Number of Columns
        //size_t CELL_TYPE;									//Spectral depth of image data

        
        string FileName;
        string CellType;
        
        //Raw data file
        ifstream file;											//or FILE* if necesary
        
        //Memory Buffer
        vector<CELL_TYPE> Buffer;								//in case of 1D array
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
        
        //Data Input Functions
        CELL_TYPE GetData (int line, int column, int band);//Returns one pixel value
        void FillBuffer(int LineStart, int LineEnd);	//Fills the Buffer with several lines
        //void FillBuffer(Macroblock* macro);				//Fills the buffer with one macroblock
        void FullBuffer();								//Loads the whole image on memory
        
        void SaveAs(string);							//Saves the image in ERS format
        void SaveAs(string, int);						//Saves the image with null value
        int D1(int line, int column, int band);
                
        //Constructors and Destructors
        ERS_Image(string, int);		//Constructor accepting the Header Filename as an argument
        
        ERS_Image(string NAME, int BANDS, int LINES, int COLUMNS, string cell_type);			
        	//Generic constructor for output images
        ~ERS_Image();
        
    };//ERS_Image End

    struct BoundaryLine {
        int Line;
        int ColumnStart;
        int ColumnEnd;
        
        //bool operator< (const BoundaryLine& Second);
        //bool operator> (const BoundaryLine& Second);
    };
    
    bool operator < (const BoundaryLine& First, const BoundaryLine& Second);
    bool operator > (const BoundaryLine& First, const BoundaryLine& Second);
    
    struct Neighbor {
        int id;
        float merging_heterogeneity;
        //bool operator == (const Neighbor& Second);
        //bool operator < (const Neighbor& Second);
        //bool operator > (const Neighbor& Second);
    };
    
    bool operator == (const Neighbor& First, const Neighbor& Second);
    bool operator < (const Neighbor& First, const Neighbor& Second);
	bool operator > (const Neighbor& First, const Neighbor& Second);

    class Object {
    	public:
	    int id;		//perhaps double???
	    int area;
	    float perimeter;
	    //int cycle;			//number of segmentation cycle that created this object
     	bool merged;		/*Shows if the object was merged during current cycle 
      						of segmentation*/
      	bool virtual_merged;//Shows if the virtual calculations have been performed over a cyrcle
      	bool priority_loaded;/*Shows if the object has been loaded to the priority list
      						in order to avoid duplicate loadings*/
       	vector<BoundaryLine> Boundary;		//list??
	    
	    vector<float> MeanValue;//Holds the mean value in order to create the output object image
    	
    	vector<Neighbor> Neighbors;
        /*List containing the total neighbors of the object, plus the heterogeneity 
         "cost" of a posible merge*/
        bool check_neighbor(int y);//Check if an id is present into the Neighbor list
        void merge_neighbors(const Object& First, const Object& Second);
        void fix_neighbors();
        
        Object(int x);
        Object();
        Object(const Object& O);
        ~Object();
        void clear();
        
        //friend Object operator+ (const Object& Left, const Object& Right);
        //will be used with the copy constructor. SOS! 1 object more to be created
        
        void merge(const Object& First, const Object& Second);
        //will not use temporary object
        
        void Calculate_MeanValue(ERS_Image& Img);
        float Color_Heterogeneity(ERS_Image& Img);
        float Compactness();
        float Smoothness();
        //float Edge_Compensation(ERS_Image& img);
        
    };//Object End

    struct PropertyRow {
        float mean_value;
        float stdev;
    };   

	class Level {
    	public:
    	map<int, Object> Objects;
    	//map<long int, PropertyTable> PropertyMap;	//TODO
     	
      	float HierarchyID;
    	vector<int> raster;
    	vector<bool> BoundaryMap;//Mark boudary pixels as 1, rest as 0
     
        int Lines;
        int Columns;
        
        float ScaleParameter;
        float Color;
        float Compactness;
        float Edge;
        
        int cycles;						//Number of segmentation passes
        
        //Multiresolution parameters
        bool ExistenceOfSuperlevel;
        bool ExistenceOfSublevel;
        Level* Superlevel;
        Level* Sublevel;
        
        Level();
        Level(const Level& L);        //Copy constructor    
        ~Level();
        
        void UpdateRaster(int z);
        void CalculatePerimeter(int id);
        int CalculatePerimeter(int id1, int id2);
        
        //void Save(string LevelFile);//TODO
        //void Load(string LevelFile);//TODO
        void SaveRaster(string RasterFile);
        //void CreateERS(string ImageFileName, ERS_Image& Img, int BandNumber);//TODO
        void SaveAsERS(string ImageFileName, ERS_Image& Img);
        
        void CreateBoundaryMap();
        void SaveBoundaryMapERS(string ImageFileName);
        void clear();
        
    };//Level End

    class PropertyTable {//It's a class because some functions will be implemented later
    	public:
     	int id;
     	
     	vector<PropertyRow> BandProperties; 
        //Property vectors. Each vector holds properties for one band
      	
       /*The property table has limited records and is not dynamic. A vector, 
       however is necesary due to the unknown number of bands. It is expected to 
       be replaced with Vector<float[]>. Perhaps a new class caled PropertyRow so 
       to get vector<PropertyRow>.*/
      	
      	PropertyTable(int ID, int BANDS); 	//Constructor
        ~PropertyTable();						//Destructor
    };//PropertyTable End    

    struct TopologyObject {
        int id;
        bool merged;
        int new_id;
    };

    class PriorityObject {
    	public:
     	int primary;
    	int secondary;
    	int id;
    	//friend bool operator < (const PriorityObject& x, const PriorityObject& y);//TEST ME: Eckel's proposition
    	//bool operator < (const PriorityObject& y);
    	PriorityObject();
        PriorityObject(const PriorityObject& PO); 	
        PriorityObject(int pri, int sec, int x);
        ~PriorityObject();    
    };

    class MSEG_Param {
        public:
        
        float Scale;// Heterogeneity Threshold
        float Color;
        float Compact;
        float Edge;
        
        bool EC;//Edge Compensation
        short GHH; //Global Heterogeneity Heuristic: 0=off, 1=mode local, 2=mode global
        bool MA;//Multiresolution Algorithm
        
        MSEG_Param();
        MSEG_Param(float Sc, float Clr, float Cmp, float Edg, bool ec, bool ma, short ghh);
        MSEG_Param(float Sc, float Clr, float Cmp);
        ~MSEG_Param();
        
    };//MSEG_Param End
	
 	struct LevelQueueObj{
	    int number;
		MSEG_Param param;
		int sub_level_number;
    };
    
    struct TopologyPair{
        int from;
        int to;
    };    

    class MSEG_Init {			//Initializes the basic structures for the MSEG algorithm
        public:
        
        map<int, TopologyObject> Topology;
        /*Holds all the object ID's created and if it is merged, points to the new object*/
        
        vector<TopologyPair> Temp_Topology;
        /*Temporarily holds topology until a level is formed and then copies objects to the
        main Topology map*/
        
        priority_queue<PriorityObject> PriorityList;
        /*Holds the queue, in which the objects will be evaluated (virtualy merged) 
        by the algorithm*/
        
        map<float, Level*> LevelHierarchy;
                //keeps the real hierarchy of the levels after they are created
                //float keeps the id of the level
                //Level* is the pointer for the memory where the level is deployed
                
        vector<LevelQueueObj> LevelQueue;
                //keeps the flow from the input txt file, before the levels are created
                //int is the Level Queue number, the priority to be created
                //MSEG_Param stores the parameters from the txt file
                //int is the level int from which is sub level
        
        map<int, float> Queue_to_Hierarchy;
                //map between LevelHierarchy and LevelQueue
                //computed after the Queue is completed and during creation of levels
        
        map<int, short> Multi_Pushed;
        		//This map holds id's that occur many times in Priority list, in order
        		//to avoid infinite loops of the list
        
        vector<pixel> StartingPoints;
        /*The results of the SPE algorithm*/
        void SaveSPE(string ImageFileName, ERS_Image& Img);
        
        int last_id;
        
        MSEG_Init();
        ~MSEG_Init();
        
    };//MSEG_Init End
   
	void update_neighbors(Object& O, MSEG_Init& mseg);
	bool operator < (const PriorityObject& x, const PriorityObject& y);
	bool operator > (const PriorityObject& x, const PriorityObject& y);

int firstpass(Level& pass1, ERS_Image& Img, MSEG_Init& mseg, float level_ID, MSEG_Param& Param);

//void firstpass(Level& pass1, ERS_Image& Img, MSEG_Init& mseg, float level_ID, MSEG_Param& Param, Level& super_level);

int secondpass(Level& pass2, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass);

//void secondpass(Level& pass2, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass, Level& super_level);

int nthpass(Level& passn, ERS_Image& Img, float level_ID, MSEG_Init& mseg, MSEG_Param& Param, Level& sub_level);

//void nthpass(Level& passn, ERS_Image& Img, float level_ID, MSEG_Init& mseg, MSEG_Param& Param, Level& sub_level, Level& super_level);

int nthpass(Level& passn, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass);

//void nthpass(Level& passn, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass, Level& super_level);

void SPE (ERS_Image& Img, MSEG_Init& mseg, int SPE_mode);


#endif /* MSEG_H */
