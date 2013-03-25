/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * msegcore.h: MSEG base modules	                                  *
 * Version: 0.9.x                                                         *
 * Last revised: 31/03/2010                                               *
 *                                                                        *
 * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * MSEG algorithm by Angelos Tzotsos (tzotsos@gmail.com)		  *
 * Remote Sensing Lab NTUA - GCpp                         March 2010      *
 *									  *
 *   Copyright (C) Angelos Tzotsos <tzotsos@gmail.com>          	  *
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


#ifndef MSEGCORE_H
#define MSEGCORE_H

#include "msegio.h"
#include <sstream>
#include "tinyxml.h"

	string ftoa(float);
	string itoa(int);
    double sqr(double x);
	typedef vector<int> CoAngle;

	struct BoundaryLine {
        int Line;
        int ColumnStart;
        int ColumnEnd;
    };

    bool operator < (const BoundaryLine& First, const BoundaryLine& Second);
    bool operator > (const BoundaryLine& First, const BoundaryLine& Second);

    struct Neighbor {
        int id;
        float merging_heterogeneity; // Why do we need that?? TODO: Remove
    };

    bool operator == (const Neighbor& First, const Neighbor& Second);
    bool operator < (const Neighbor& First, const Neighbor& Second);
	bool operator > (const Neighbor& First, const Neighbor& Second);

	struct ClassAttribute {
	    int id;
	    int R;
	    int G;
	    int B;
	    string Name;
	};

    class Object {
    	public:
	    int id;		//32 bit in x86, 64 bit in x64
	    int area;
	    int perimeter;
	    //int cycle;			//number of segmentation cycle that created this object
     	bool merged;		/*Shows if the object was merged during current cycle
      						of segmentation*/
      	bool virtual_merged;//Shows if the virtual calculations have been performed over a cyrcle
      	bool priority_loaded;/*Shows if the object has been loaded to the priority list
      						in order to avoid duplicate loadings*/
        bool edge;
        bool properties_exported;//Check if properties have been calculated for this object
       	vector<BoundaryLine> Boundary;		//list??

	    vector<float> Sum;//Holds the sum of all pixels in order to calculate statistics of the Object
    	vector<float> SumSq;//Holds the sum of squares of all the pixels to calculate statistics of the Object

    	vector<Neighbor> Neighbors;
        /*List containing the total neighbors of the object, plus the heterogeneity
         "cost" of a posible merge*/
        bool check_neighbor(int y);//Check if an id is present into the Neighbor list
        void merge_neighbors(const Object& First, const Object& Second);// Join the Neighbors table after a merge
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

        //TODO
        //void Calculate_Sums();
        //void Calculate_MeanValue();// we don't need the image anymore
        //void Calculate_Std();// we don't need the image anymore

        float Color_Heterogeneity(ERS_Image& Img);// we keep the image for band weights only
        float Compactness();
        float Smoothness();

        //TODO
        //float Edge_Compensation(ERS_Image& img);
        //float Width_to_Length();
        //void BoundingBox();

        //Texture stuff
    	vector<CoAngle> SDM;//Co-Occurence Matrix
    	vector<double> ASM;//Angular Second Moment
    	vector<long int> R;
    	//void InitSDM(Texture_Image& TexImg);//not to be used...
    	void BuildSDM(Texture_Image& TexImg);
    	void ClearSDM();
    	void CalculateASM(Texture_Image& TexImg);

    };//Object End

    class PropertyObject {//It's a class because some functions will be implemented later
    	public:
     	//int id;
        int classid;
     	vector<float> MeanValue;
		vector<float> StdDev;
        //Property vectors. Each vector holds properties for one band

        //Second Order Texture Properties
        vector<double> ASM;
        vector<double> Contrast;
        vector<double> Correlation;
        vector<double> Variance;


       /*The property table has limited records and is not dynamic. A vector,
       however is necesary due to the unknown number of bands. It is expected to
       be replaced with Vector<float[]>. Perhaps a new class caled PropertyRow so
       to get vector<PropertyRow>.*/
      	void clear();
      	PropertyObject(); 					//Constructor
        ~PropertyObject();					//Destructor
    };//PropertyObject End

/*    struct PropertyBand {
        float mean;
        float stdev;
    };
*/
	class Level {
    	public:
    	map<int, Object> Objects;
    	map<int, PropertyObject> Properties;

      	float HierarchyID;
    	vector<int> raster;
    	vector<bool> BoundaryMap;//Mark boudary pixels as 1, rest as 0
    	vector<int> SampleMap;
    	vector<int> ClassificationMap;
    	map<int, ClassAttribute> SampleAttributes;

        int Lines;
        int Columns;

        float ScaleParameter;
        float Color;
        float Compactness;
        float Edge;
        float Texture;

        int cycles;						//Number of segmentation passes

        //Multiresolution parameters
        bool ExistenceOfSuperlevel;
        bool ExistenceOfSublevel;
        Level* Superlevel;
        Level* Sublevel;

        Level();
        Level(const Level& L);        //Copy constructor
        ~Level();

        void UpdateRaster(int z);   //Updates the raster values according to the new merged object with id z. Uses the boundary lines of the object.
        void UpdateRaster(int z, int nz);
        void CalculatePerimeter(int id);
        int CalculatePerimeter(int id1, int id2);

        //void SaveLevel(string XMLFile, string RasterFile);//Combine XML and Raster Export for full Export Ability
	//void LoadLevel(string XMLFile, string RasterFile);//Combine XML and Raster Import for full Import Ability

        void SaveXML(string XMLFile, ERS_Image& Img);//Will save an XML document for the Level class
        //void LoadXML(string XMLFile, int bands);//TODO: Will load the above XML document
        void SaveMiniXML(string XMLFile);
        void LoadMiniXML(string XMLFile);
        void CalculateProperties(ERS_Image& Img);// may need image in the future...
        //TODO
        void CalculateNeighbors(); //Calculates iterativelly all the neighbors for each object - Load Mode
        void DeleteProperties();//Unloads the Properties from memory
        //void CalculateTextureProperties(Texture_Image& TexImg);
	void SaveProperties(string CSVFile);
	void SaveSampleProperties(string CSVFile);
	void SaveSVMTraining(string CSVFile);
	void SaveSVMTrainingVoting(string CSVFile, double percent);//Exports samples only if percentage of sample pixels is found in the region
	void SaveSVMTesting(string CSVFile);
	void LoadTTA(string TTAName);//Loads Sample Map from eCognition's TTAMask file
	//TODO
	//void LoadERSSamples();//Loads Sample Map from an ERS Image
	void LoadAttributes(string TTAName);//Loads sample classes from eCognition's TTAMask attribute file
	void SaveRaster(string RasterFile, ERS_Image& Img);//Saves the id raster vector.
        void LoadRaster(string RasterFile);//Loads the id raster vector.
        void LoadObjects();//Loads the Objects to memory from the id raster vector - Load Mode
	//TODO
	void LoadObjectSums(ERS_Image& Img); //Restores the primitive spectral statistics for the Objects after they are loaded from a raster - Load Mode
	void DeleteObjects();//Unloads the Objects from memory
	//void CreateERS(string ImageFileName, ERS_Image& Img, int BandNumber);//TODO -> Discontinued for api 1.0
        void SaveMeanAsERS(string ImageFileName, ERS_Image& Img);//Saves an image with object means

        void SVM2ClassificationMap(string TestFile, string PredictionFile);
        void SaveClassificationMap(string ClassFile);
        void SaveClassificationAsERS(string ClassFile);

        void CreateBoundaryMap();//Obsolete since gdal_polygonize routine
        void SaveBoundaryMapERS(string ImageFileName);//Obsolete since gdal_polygonize routine
        void clear();

    };//Level End

/*    struct TopologyObject {	//used in old api
        int id;
        bool merged;
        int new_id;
    };
*/

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
        float Texture;
        int quantizer;
        int distance;

        bool EC;//Edge Compensation
        short GHH; //Global Heterogeneity Heuristic: 0=off, 1=mode local, 2=mode global
        bool MA;//Multiresolution Algorithm
        bool TH;//2nd order Texture Heterogeneity

        MSEG_Param();
        MSEG_Param(float Sc, float Clr, float Cmp, float Edg,float Tex, int quant, int dist,bool ec, bool ma, short ghh, bool tx);
        MSEG_Param(float Sc, float Clr, float Cmp);
        ~MSEG_Param();

    };//MSEG_Param End

 	struct LevelQueueObj{
	    int number;
		MSEG_Param param;
		int sub_level_number;
    };

/*    struct TopologyPair{	//used in old api
        int from;
        int to;
    };
*/

    class MSEG_Init {			//Initializes the basic structures for the MSEG algorithm
        public:

        //map<int, TopologyObject> Topology;
        /*Holds all the object ID's created and if it is merged, points to the new object*/

        //vector<TopologyPair> Temp_Topology;
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

	//void update_neighbors(Object& O, MSEG_Init& mseg);
	bool operator < (const PriorityObject& x, const PriorityObject& y);
	bool operator > (const PriorityObject& x, const PriorityObject& y);

	class Starting_Points_Estimation{
	    public:

     	int SPE_mode;
     	void HSI(ERS_Image& Img, MSEG_Init& mseg);
     	void YUV(ERS_Image& Img, MSEG_Init& mseg);
     	void PCA(ERS_Image& Img, MSEG_Init& mseg);
     	void Dithering(ERS_Image& Img, MSEG_Init& mseg, string Method);

     	Starting_Points_Estimation(int mode);
     	~Starting_Points_Estimation();
	};//Starting_Points_Estimation End

#endif /* MSEGCORE_H */
