/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * libmseg.cpp: Multiscale SEGmentation library                           *
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

#include "libmseg.h"

string find_cell_type(string HeaderName){
    
    cout << "Finding bit depth from image header..." << endl;
    
    ifstream infile;
    infile.open(HeaderName.c_str());
    
    if (!infile){								//Error handling 
	    cout << "Error opening file " << HeaderName << endl;
	    //throw FileError;
     	
      	cout << "Press any key to quit..." << endl;
	    cin.get();
	    exit(1);
	}
	
	string CellType;
	string word;
	
	//Header file parser
 	while (!infile.eof()){
 	    infile >> word;
 	    
 	    if (strncmp(word.c_str(), "CellType", 8) == 0){
           infile >> word;		/*Skips the "="*/
           infile >> CellType;
        }
    }
    
    infile.close();
    
    cout << "Found ERS Cell Type : " << CellType << endl;
    
    return CellType;

}

pixel D2(int d, int imgcolumns){
    pixel x;
    int k = d % imgcolumns;
    x.column = (k) ? k : imgcolumns;
    x.line = (k) ? ((int)(d/imgcolumns) + 1) : ((int)(d/imgcolumns));
    return x;
}

int D1(int imglines, int imgcolumns, int line, int column){
    return (((line-1)*imgcolumns) + column -1);
}      

void ERS_Image::FillBuffer (int LineStart, int LineEnd){
    
    if(BuffLineStart != BuffLineEnd){
        cout << "Buffer is not empty" << endl;
        cout << "Deleting buffer..." << endl;
        Buffer.clear();
    }
    
    BuffLineStart = LineStart;
    BuffLineEnd = LineEnd;
    
    CELL_TYPE tmp;
    file.seekg(0, ios::end);
    streampos file_end = file.tellg(); // file_end keeps the file size in bytes
    cout << "File size is: " << file_end << endl;
    file.seekg( ((BuffLineStart - 1)*Bands*Columns*sizeof(CELL_TYPE)), ios::beg);
    streampos buff_start = file.tellg();
    
    for(int i=0;i<((BuffLineEnd - BuffLineStart + 1)*Bands*Columns);i++){
        file >> tmp;
        Buffer.push_back(tmp);
    }
    streampos buff_end = file.tellg();
    file.seekg(0, ios::beg);
    
    cout << "Buffer write succesfull" << endl;
    
}

/*void ERS_Image::FillBuffer (Macroblock* macro){
    
}*/

void ERS_Image::FullBuffer (){
    
    if(BuffLineStart != BuffLineEnd){
        cout << "Buffer is not empty" << endl;
        cout << "Deleting buffer..." << endl;
        Buffer.clear();
    }    
    
    CELL_TYPE tmp;
    
    file.seekg(0, ios::beg);
    /*
    while (file){
        file >> tmp;
        if(!file) break;
        Buffer.push_back(tmp);
    }
    */
    
	for(int i=0; i<(Lines*Columns*Bands); i++){
	    file.read(reinterpret_cast<char *>(&tmp), sizeof(CELL_TYPE));
	    Buffer.push_back(tmp);
	}
	
	BuffLineStart = 1;
    BuffLineEnd = Lines;
    file.seekg(0, ios::beg);     
    cout << "Buffer Full!" << endl;
    
}

CELL_TYPE ERS_Image::GetData (int line, int column, int band){
   //TODO: Error handling if line, column, band > image's 
   file.seekg( ((((line-1)*Bands*Columns) + ((band-1)*Columns) + column)*sizeof(CELL_TYPE)), ios::beg);
   CELL_TYPE tmp;
   file >> tmp;
   return tmp;
   
}

int ERS_Image::D1(int line, int column, int band){
    return (((line-1)*Bands*Columns) + ((band-1)*Columns) + column -1);
}

void ERS_Image::SaveAs(string OutputName){
    ofstream head_out, file_out;
    string tmp;
    tmp = OutputName + ".ers";
    
    head_out.open(tmp.c_str());
    file_out.open(OutputName.c_str(), ios::binary);
    
    //cout << "Buffer size is: " << Buffer.size() << endl;
    
	//CELL_TYPE t;
	//save binary data
    for(int i=0; i<Buffer.size(); i++){
        //t = Buffer[i];
		file_out.write(reinterpret_cast<char *> (&Buffer[i]), sizeof(CELL_TYPE));
    }
    //Coordinate CoordinateSpace;
    //save header
    head_out << "DatasetHeader Begin" << endl;
    head_out << "\tVersion \t= \"" << ERSVersion << "\"" << endl;
    //head_out << "\tName\t\t= \"" << HeaderName << "\"" << endl;
    head_out << "\tLastUpdated\t= " << LastUpdated << endl;
    head_out << "\tDataSetType\t= ERStorage" << endl;
    head_out << "\tDataType\t= " << DataType << endl;
    head_out << "\tByteOrder\t= " << ByteOrder << endl;
    head_out << "\tCoordinateSpace Begin" << endl;
    head_out << "\t\tDatum\t\t= \"" << CoordinateSpace.Datum << "\"" << endl;
    head_out << "\t\tProjection\t= \"" << CoordinateSpace.Projection << "\"" << endl;
    head_out << "\t\tCoordinateType\t= " << CoordinateSpace.CoordinateType << endl;
    head_out << "\t\tUnits\t\t= \"" << CoordinateSpace.Units << "\"" << endl;
    head_out << "\t\tRotation\t= " << CoordinateSpace.Rotation << endl;
    head_out << "\tCoordinateSpace End" << endl;
    head_out << "\tRasterInfo Begin" << endl;
    head_out << "\t\tCellType\t= " << CellType << endl;
    head_out << "\t\tNrOfLines\t= " << Lines << endl;
    head_out << "\t\tNrOfCellsPerLine\t= " << Columns << endl;
    head_out << "\t\tNrOfBands\t= " << Bands << endl;
    for(int i=0;i<BandID.size();i++){
        head_out << "\t\tBandId Begin" << endl;
        head_out << "\t\t\tValue\t\t= " << BandID[i] << endl;
        head_out << "\t\tBandId End" << endl;
    }
    head_out << "\t\tRegionInfo Begin" << endl;
    head_out << "\t\t\tType\t\t= Polygon" << endl;
    head_out << "\t\t\tRegionName\t= \"All\"" << endl;
    //head_out << "\t\t\tSourceDataset\t= \"" << FileName << "\"" << endl;
    head_out << "\t\t\tRGBcolour Begin" << endl;
    head_out << "\t\t\t\tRed\t\t= 65535" << endl;
    head_out << "\t\t\t\tGreen\t\t= 65535" << endl;
    head_out << "\t\t\t\tBlue\t\t= 65535" << endl;
    head_out << "\t\t\tRGBcolour End" << endl;
    head_out << "\t\t\tSubRegion\t= {" << endl;
    head_out << "\t\t\t\t0\t0" << endl;
    head_out << "\t\t\t\t0\t" << Lines << endl;
    head_out << "\t\t\t\t" << Columns << "\t" << Lines << endl;
    head_out << "\t\t\t\t" << Columns << "\t0" << endl;
    head_out << "\t\t\t}" << endl;
    head_out << "\t\tRegionInfo End" << endl;
    head_out << "\tRasterInfo End" << endl;
    head_out << "DatasetHeader End" << endl;
    
    //Close files
    file_out.close();
    head_out.close();
}

void ERS_Image::SaveAs(string OutputName, int nullval){
    NullValue = nullval;
    ofstream head_out, file_out;
    string tmp;
    tmp = OutputName + ".ers";
    
    head_out.open(tmp.c_str());
    file_out.open(OutputName.c_str(), ios::binary);
    
    //save binary data
    for(int i=0; i<Buffer.size(); i++){
        file_out.write(reinterpret_cast<char *> (&Buffer[i]), sizeof(CELL_TYPE));
    }
    //Coordinate CoordinateSpace;
    //save header
    head_out << "DatasetHeader Begin" << endl;
    head_out << "\tVersion \t= \"" << ERSVersion << "\"" << endl;
    //head_out << "\tName\t\t= \"" << HeaderName << "\"" << endl;
    head_out << "\tLastUpdated\t= " << LastUpdated << endl;
    head_out << "\tDataSetType\t= ERStorage" << endl;
    head_out << "\tDataType\t= " << DataType << endl;
    head_out << "\tByteOrder\t= " << ByteOrder << endl;
    head_out << "\tCoordinateSpace Begin" << endl;
    head_out << "\t\tDatum\t\t= \"" << CoordinateSpace.Datum << "\"" << endl;
    head_out << "\t\tProjection\t= \"" << CoordinateSpace.Projection << "\"" << endl;
    head_out << "\t\tCoordinateType\t= " << CoordinateSpace.CoordinateType << endl;
    head_out << "\t\tUnits\t\t= \"" << CoordinateSpace.Units << "\"" << endl;
    head_out << "\t\tRotation\t= " << CoordinateSpace.Rotation << endl;
    head_out << "\tCoordinateSpace End" << endl;
    head_out << "\tRasterInfo Begin" << endl;
    head_out << "\t\tCellType\t= " << CellType << endl;
    head_out << "\t\tNullCellValue\t= " << NullValue << endl;
    head_out << "\t\tNrOfLines\t= " << Lines << endl;
    head_out << "\t\tNrOfCellsPerLine\t= " << Columns << endl;
    head_out << "\t\tNrOfBands\t= " << Bands << endl;
    for(int i=0;i<BandID.size();i++){
        head_out << "\t\tBandId Begin" << endl;
        head_out << "\t\t\tValue\t\t= " << BandID[i] << endl;
        head_out << "\t\tBandId End" << endl;
    }
    head_out << "\t\tRegionInfo Begin" << endl;
    head_out << "\t\t\tType\t\t= Polygon" << endl;
    head_out << "\t\t\tRegionName\t= \"All\"" << endl;
    //head_out << "\t\t\tSourceDataset\t= \"" << FileName << "\"" << endl;
    head_out << "\t\t\tRGBcolour Begin" << endl;
    head_out << "\t\t\t\tRed\t\t= 65535" << endl;
    head_out << "\t\t\t\tGreen\t\t= 65535" << endl;
    head_out << "\t\t\t\tBlue\t\t= 65535" << endl;
    head_out << "\t\t\tRGBcolour End" << endl;
    head_out << "\t\t\tSubRegion\t= {" << endl;
    head_out << "\t\t\t\t0\t0" << endl;
    head_out << "\t\t\t\t0\t" << Lines << endl;
    head_out << "\t\t\t\t" << Columns << "\t" << Lines << endl;
    head_out << "\t\t\t\t" << Columns << "\t0" << endl;
    head_out << "\t\t\t}" << endl;
    head_out << "\t\tRegionInfo End" << endl;
    head_out << "\tRasterInfo End" << endl;
    head_out << "DatasetHeader End" << endl;
    
    //Close files
    file_out.close();
    head_out.close();
}    

Image::Image(){
    
}

Image::~Image(){
    
}    

ERS_Image::ERS_Image(string InputHeaderName, int msize){

	cout << "Loading image " << InputHeaderName << " ..." << endl;
	//system("Pause");
	
 	if(msize){
     	MacroSize = msize;
	}
 	else{
      	MacroSize = 16;
 	}
         	
 	ifstream header;
	header.open(InputHeaderName.c_str());		//It's an input fstream, explicit call.
	
 	if (!header){								//Error handling 
	    cout << "Error opening file " << InputHeaderName << endl;
	    //throw FileError;
     	
      	cout << "Press any key to quit..." << endl;
	    cin.get();
	    exit(1);
	}
 	
  	HeaderName = InputHeaderName;/*Don't want to aquire the HeaderName from within,
								 it's better to trust the filesystem's name for the .ers file*/
 	FileName = HeaderName.substr(0, HeaderName.size()-4);	/*Removes ".ers" */
 	cout << "Raw Data File Name : " << FileName << endl;
 	
   	string word;				//temporary string for input
 	
 	//Header file parser
 	
 	while (header){
 	    header >> word;
 	    
 	    if (strncmp(word.c_str(), "Version", 6) == 0){
 	        header >> word;		/*Skips the "="*/
 	        header >> word;
 	        ERSVersion = word.substr(1,word.size()-2);		/*Removes the "___"*/
 	        //string ERSVersion(word, 2, length-1);				/*Alternative*/ 
          	cout << "Found ERS Version : " << ERSVersion << endl;
       }
       
       if (strncmp(word.c_str(), "LastUpdated", 11) == 0){
           header >> word;		/*Skips the "="*/
           for (int i=0; i<5; i++){
               header >> word;
               LastUpdated = LastUpdated + word + " ";/*Reads 5 words*/
           }
           header >> word;
           LastUpdated = LastUpdated + word; /*Does not add the space after the sixth word*/
           cout << "Last ERS Update : " << LastUpdated << endl;
       }
       
       if (strncmp(word.c_str(), "DataType", 8) == 0){
           header >> word;		/*Skips the "="*/
           header >> DataType;
           cout << "ERS Data Type : " << DataType << endl;
       }
       
       if (strncmp(word.c_str(), "ByteOrder", 9) == 0){
           header >> word;		/*Skips the "="*/
           header >> ByteOrder;
           cout << "ERS Byte Order : " << ByteOrder << endl;
       }
       
       if (strncmp(word.c_str(), "Datum", 5) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           CoordinateSpace.Datum = word.substr(1,word.size()-2);		/*Removes the "___"*/
           cout << "ERS Datum : " << CoordinateSpace.Datum << endl;
       }
       
       if (strncmp(word.c_str(), "Projection", 10) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           CoordinateSpace.Projection = word.substr(1,word.size()-2);	/*Removes the "___"*/
           cout << "ERS Projection : " << CoordinateSpace.Projection << endl;
       }
       
       if (strncmp(word.c_str(), "CoordinateType", 14) == 0){
           header >> word;		/*Skips the "="*/
           header >> CoordinateSpace.CoordinateType;
           cout << "ERS Coordinate Type : " << CoordinateSpace.CoordinateType << endl;
       }
       
       if (strncmp(word.c_str(), "Units", 5) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           CoordinateSpace.Units = word.substr(1,word.size()-2);	/*Removes the "___"*/
           cout << "ERS Units : " << CoordinateSpace.Units << endl;
       }
       
       if (strncmp(word.c_str(), "Rotation", 8) == 0){
           header >> word;		/*Skips the "="*/
           header >> CoordinateSpace.Rotation;
           cout << "ERS Rotation : " << CoordinateSpace.Rotation << endl;
       }
       
       if (strncmp(word.c_str(), "CellType", 8) == 0){
           header >> word;		/*Skips the "="*/
           header >> CellType;
           cout << "ERS Cell Type : " << CellType << endl;
       }
       
       if (strncmp(word.c_str(), "NullCellValue", 13) == 0){
           header >> word;		/*Skips the "="*/
           header >> NullValue;
           cout << "ERS Null Value : " << NullValue << endl;
       }
       
       if (strncmp(word.c_str(), "NrOfLines", 9) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           Lines = atoi(word.c_str());
           if (Lines < 1){
               cout << "Wrong number of lines!" << endl;
               cout << "Press any key to quit..." << endl;
               cin.get();
               exit(1);
           }   
           cout << "ERS Number Of Lines : " << Lines << endl;
       }
       
       if (strncmp(word.c_str(), "NrOfCellsPerLine", 16) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           Columns = atoi(word.c_str());
           if (Columns < 1){
               cout << "Wrong number of pixels per line!" << endl;
               cout << "Press any key to quit..." << endl;
               cin.get();
               exit(1);
           }   
           cout << "ERS Number Of Pixels per Line : " << Columns << endl;
       }
       
       if (strncmp(word.c_str(), "NrOfBands", 9) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           Bands = atoi(word.c_str());
           if (Bands < 1){
               cout << "Wrong number of bands!" << endl;
               cout << "Press any key to quit..." << endl;
               cin.get();
               exit(1);
           }    
           cout << "ERS Number Of Bands : " << Bands << endl;
           //vector<string> BandID;
           for (int i=0; i<5; i++){
               header >> word;
           }
           BandID.push_back(word);
           cout << "ERS Band Names:" << endl << word << endl;
           for (int k=0;k < Bands-1; k++){
               for (int j=0; j<7; j++){
               		header >> word;
     		   }
          	   BandID.push_back(word);
          	   cout << word << endl;
           }
       }
   }//Header file parser END
   header.close();
   cout << "Header closed" << endl;   
   
   //Opening binary data file for input
   file.open(FileName.c_str(), ios::binary);
   if (!file){								//Error handling 
	    cout << "Error opening file " << FileName << endl;
	    //throw FileError;
     	
      	cout << "Press any key to quit..." << endl;
	    cin.get();
	    exit(1);
   }
   cout << "Image File Opened... " << endl;
   
   //Buffer Initialization
   cout << "Buffer initialized..." << endl;
   BuffLineStart = 0;
   BuffLineEnd = 0;   // if BuffLineStart == BuffLineEnd ---> Buffer empty
   //system("Pause");
   
   /*Default Band Weight values can be modified through the main routine from
   an outside *.txt file */
   BandWeight.push_back(0.0);
   for (int i=1; i<=Bands; i++){
       BandWeight.push_back(1.0);
   }
   cout << "Band weights initialized..." << endl;
   
   //Macroblock calculation and initialization
   
   int Line_Remainder;
   int Row_Remainder;
   
   /*TODO: Check here if MacroSize is divided by 2*/
   
   Line_Remainder = Lines % MacroSize;
   Row_Remainder = Columns % MacroSize;
   
   /*TODO: Check the rounding procedure int(...)*/
   
   if(Line_Remainder == 0){
       MacroLines = Lines / MacroSize;
   } 
   else{
       MacroLines = int(Lines/MacroSize) + 1;
   }
   
   if (Row_Remainder == 0){
       MacroColumns = Columns / MacroSize;
   } 
   else{
       MacroColumns = int(Columns/MacroSize) + 1;
   }
   
   cout << "Macrolines = " << MacroLines << "    Macrocolumns = " << MacroColumns << endl;
   
   Macroblock TempMb;
   TempMb.size=MacroSize;
   //int NumOfMb = MacroRows * MacroLines;
   
   //Adds a null macroblock in place Macroblocks[0]
   /*TempMB.LineEnd = 0;
   TempMB.LineStart = 0;
   TempMb.ColumnEnd = 0;
   TempMb.ColumnStart = 0;
   Macroblocks.push_back(TempMb);
   */
   
   for (int i=1;i<=MacroLines;i++){
       for (int j=1;j<=MacroColumns;j++){
           if (i*MacroSize >= Lines){
               TempMb.LineEnd = Lines;
           }
           else {
               TempMb.LineEnd = i*MacroSize;
           }
           TempMb.LineStart = ((i-1)* MacroSize) + 1;
           
           if (j*MacroSize >= Columns){
               TempMb.ColumnEnd = Columns;
           }
           else {
               TempMb.ColumnEnd = j*MacroSize;
           }
           TempMb.ColumnStart = ((j-1)* MacroSize) + 1;
           
           Macroblocks.push_back(TempMb);
       }
   }              
   cout << "Macroblocks loaded = " << Macroblocks.size() << endl;
   //Macroblock calculation End
   cout << "ERS_Image constructor succesfully finished..." << endl;
   //system("Pause");
}

ERS_Image::ERS_Image(string Name, int BANDS, int LINES, int COLUMNS, string cell_type){
    
    Bands = BANDS;
    Lines = LINES;
    Columns = COLUMNS;
    
    FileName = Name;
    CellType = cell_type;
    HeaderName = FileName + ".ers";
    
    //struct Coordinate CoordinateSpace;
    
    ERSVersion = "6.3";
    LastUpdated = "Tue Mar 09 17:41:16 GMT 2004";
    DataType = "Raster";
    ByteOrder = "LSBFirst";
    CoordinateSpace.Datum = "RAW";
    CoordinateSpace.Projection = "RAW";
    CoordinateSpace.CoordinateType = "RAW";
    CoordinateSpace.Units = "natural";
    CoordinateSpace.Rotation = "0:0:0.0";
    
    //vector<string> BandID;
    //BandID.push_back("Level 1");
    
    //vector<CELL_TYPE> Buffer;
    //Buffer.reserve(Lines*Columns*Bands);
    BuffLineStart = 0;
    BuffLineEnd = 0;
    
    //vector<float> BandWeight;
    //BandWeight[0] = 1.0;
    //vector<Macroblock> Macroblocks;
    MacroSize = 0;
    MacroColumns = 0;
    MacroLines = 0;
    
}

ERS_Image::~ERS_Image(){

}

Object::Object(int x){
    id = x;
    merged = FALSE;
    priority_loaded = FALSE;
    virtual_merged = FALSE;
    area = 0;
    perimeter = 0.0;
}

Object::Object(){
    id = 0;
    merged = FALSE;
    priority_loaded = FALSE;
    virtual_merged = FALSE;
    area = 0;
    perimeter = 0.0;
}

Object::Object(const Object& O){
    id = O.id;
    merged = O.merged;
    priority_loaded = O.priority_loaded;
    virtual_merged = O.virtual_merged;
    area = O.area;
    perimeter = O.perimeter;
    
    for(int i=0; i<O.Boundary.size(); i++){
        Boundary.push_back(O.Boundary[i]);
    }
    for(int i=0; i<O.MeanValue.size(); i++){
        MeanValue.push_back(O.MeanValue[i]);
    }
    for(int i=0; i<O.Neighbors.size(); i++){
        Neighbors.push_back(O.Neighbors[i]);
    }            
}    

Object::~Object(){
}

/*Object operator+ (const Object& Left, const Object& Right){
    
}*/

void Object::clear(){
    id = 0;
    merged = FALSE;
    priority_loaded = FALSE;
    virtual_merged = FALSE;
    area = 0; //area will be 1 or 2 in 1st pass where this is useful
    perimeter = 0.0;
    //cycle = 0;
    Boundary.clear();
    MeanValue.clear();
    Neighbors.clear();
}

bool Object::check_neighbor(int y){
    for(int i=0; i<Neighbors.size(); i++){
        if(Neighbors[i].id == y)
            return TRUE;
    }
    return FALSE;  
}    

void Object::merge(const Object& First, const Object& Second){
	 //First do the easy stuff :)
	 merged = FALSE;
	 priority_loaded = FALSE;
	 virtual_merged = FALSE;
	 perimeter = 0;
	 area = First.area + Second.area;
     //Boundary.clear();
     
     //TODO: Do the Neigbors stuff here
     
     //The Boundary stuff begin here
     BoundaryLine tmp;
     tmp.Line = -1;
     tmp.ColumnStart = -1;
     tmp.ColumnEnd = -1;
     //bool in =FALSE;
     //Create a temporary vector to merge all the BoundaryLines
	 vector<BoundaryLine> Boundtemp;
	 //Push everything to temp vector
	 for(int i=0; i<First.Boundary.size(); i++){
	     Boundtemp.push_back(First.Boundary[i]);
	 }
	 for(int i=0; i<Second.Boundary.size(); i++){
	     Boundtemp.push_back(Second.Boundary[i]);
	 }
	 sort(Boundtemp.begin(), Boundtemp.end());
	 Boundtemp.push_back(tmp);//Ads an "EOF" object at the end of list in order to break correctly
	 for(int i=0; i<Boundtemp.size(); i++){
	     tmp.Line = Boundtemp[i].Line;
         tmp.ColumnStart = Boundtemp[i].ColumnStart;
         tmp.ColumnEnd = Boundtemp[i].ColumnEnd;
         while((Boundtemp[i].Line == Boundtemp[i+1].Line)&&((Boundtemp[i].ColumnEnd + 1) == Boundtemp[i+1].ColumnStart)){
             tmp.ColumnEnd = Boundtemp[i+1].ColumnEnd;
             i++;
             //in = TRUE;
         }
         if(Boundtemp[i].Line == -1){
             break;
         }
         Boundary.push_back(tmp);
         /*if(in){
            i--;
            in = FALSE;
         }*/
     }            
}        

float Object::Color_Heterogeneity(ERS_Image& Img){
    if(area == 1) return 0.0;
    //cout << "area is " << area << endl;
    float std[Img.Bands];
    float s[Img.Bands];
    float ss[Img.Bands];
    for(int b=0; b<Img.Bands; b++){
        s[b] = 0.0;
        ss[b] = 0.0;
    }    
    int l, cur_pos;
    float h = 0.0;
    MeanValue.clear();
    
    for(int i=0; i<Boundary.size(); i++){
        l = Boundary[i].Line;
        for(int b=0; b<Img.Bands; b++){
            for(int c = Boundary[i].ColumnStart; c<=Boundary[i].ColumnEnd; c++){
                cur_pos = Img.D1(l, c, b+1);
                s[b]+=Img.Buffer[cur_pos];
                ss[b]+=(Img.Buffer[cur_pos]*Img.Buffer[cur_pos]);
            }
        }
    }
    for(int b=0; b<Img.Bands; b++){
        MeanValue.push_back((s[b]/area));
        std[b] = (sqrt(ss[b]/area-s[b]*s[b]/area/area));
        h+=(Img.BandWeight[b+1]*area*std[b]);
    }    
    return h;
}

float Object::Compactness(){
    //SOS:must have calculated the perimeter...
    return (area * perimeter)/(sqrt((float)area));
}

/*float Object::Edge_Compensation(){
    return 0.0;
}*/

void Object::merge_neighbors(const Object& First, const Object& Second){
	 for(int i=0; i<First.Neighbors.size(); i++){
	     if(First.Neighbors[i].id == Second.id) continue;//Skip because of merge 
		 Neighbors.push_back(First.Neighbors[i]);
	 }
	 for(int i=0; i<Second.Neighbors.size(); i++){
	     if(Second.Neighbors[i].id == First.id) continue;//Skip because of merge
		 Neighbors.push_back(Second.Neighbors[i]);
	 }
	 sort(Neighbors.begin(), Neighbors.end());
	 vector<Neighbor>::iterator p = unique(Neighbors.begin(), Neighbors.end());
	 Neighbors.erase(p, Neighbors.end());
}

void Object::fix_neighbors(){
	 sort(Neighbors.begin(), Neighbors.end());
	 vector<Neighbor>::iterator p = unique(Neighbors.begin(), Neighbors.end());
	 Neighbors.erase(p, Neighbors.end());
}

void Object::Calculate_MeanValue(ERS_Image& Img){
    float s[Img.Bands];
    int l;
    MeanValue.clear();
    for(int i=0; i<Boundary.size(); i++){
        l = Boundary[i].Line;
        for(int b=0; b<Img.Bands; b++){
            for(int c = Boundary[i].ColumnStart; c<=Boundary[i].ColumnEnd; c++){
                s[b]+=Img.Buffer[Img.D1(l, c, b+1)];
            }
        }
    }
    for(int b=0; b<Img.Bands; b++){
        MeanValue.push_back(s[b]/area);
    }
}    

float Object::Smoothness(){
    int mincol, maxcol, minline, maxline;
    minline = Boundary[0].Line;
    maxline = Boundary[Boundary.size()-1].Line;
    mincol = Boundary[0].ColumnStart;
    maxcol = Boundary[0].ColumnEnd;
    for(int i=0; i<Boundary.size(); i++){
        if(Boundary[i].ColumnStart < mincol){
            mincol = Boundary[i].ColumnStart;
        }
        if(Boundary[i].ColumnEnd > maxcol){
            maxcol = Boundary[i].ColumnEnd;
        }
    }
    float smth = (area * perimeter)/(2 * (maxcol - mincol + 1) + 2 * (maxline - minline + 1));
    return smth;
}


bool operator < (const BoundaryLine& First, const BoundaryLine& Second){
	 if(First.Line != Second.Line){
  	     return (First.Line < Second.Line);
	 }
	 else{
	 	  return (First.ColumnStart < Second.ColumnStart);
     }
}

bool operator > (const BoundaryLine& First, const BoundaryLine& Second){
	 if(First.Line != Second.Line){
  	     return (First.Line > Second.Line);
	 }
	 else{
	 	  return (First.ColumnStart > Second.ColumnStart);
     }
}

bool operator == (const Neighbor& First, const Neighbor& Second){
	 if(First.id == Second.id){
	     return TRUE;
	 }else{
  	     return FALSE;
	 }
}

bool operator < (const Neighbor& First, const Neighbor& Second){
	 if(First.id < Second.id){
	     return TRUE;
	 }else{
  	     return FALSE;
	 }
}

bool operator > (const Neighbor& First, const Neighbor& Second){
	 if(First.id > Second.id){
	     return TRUE;
	 }else{
  	     return FALSE;
	 }
}

Level::Level(){
    HierarchyID = 0.0;
    Lines = 0;
    Columns = 0;
    ScaleParameter = 0.0;
    Color = 0.0;
    Compactness = 0.0;
    Edge = 0.0;
    cycles = 0;
    ExistenceOfSuperlevel = FALSE;
    ExistenceOfSublevel = FALSE;
    Level* Superlevel = NULL;
    Level* Sublevel = NULL;
}

Level::Level(const Level& L){
    Objects = L.Objects;
    HierarchyID = L.HierarchyID;
    for(int i=0; i<L.raster.size(); i++){
        raster[i] = L.raster[i];
    }
    Lines = L.Lines;
    Columns = L.Columns;
    ScaleParameter = L.ScaleParameter;
    Color = L.Color;
    Compactness = L.Compactness;
    Edge = L.Edge;
    cycles = 0;
    ExistenceOfSuperlevel = L.ExistenceOfSuperlevel;
    ExistenceOfSublevel = L.ExistenceOfSublevel;
    Level* Superlevel = L.Superlevel;
    Level* Sublevel = L.Sublevel;
}    

Level::~Level(){ //or empty????
    Objects.clear();
    raster.clear();
    //cout << "Inside the destructor" << endl;
}

void Level::UpdateRaster(int z){
    int l;
    for(int i=0; i<Objects[z].Boundary.size(); i++){
        l = Objects[z].Boundary[i].Line;
        for(int c = Objects[z].Boundary[i].ColumnStart; c<=Objects[z].Boundary[i].ColumnEnd; c++){
            raster[D1(Lines, Columns, l, c)] = z;
        }
    }        
}

void Level::CalculatePerimeter(int id){
	 int per = 0;
	 int l, c;
	 for(int i=0; i<Objects[id].Boundary.size(); i++){
	     l = Objects[id].Boundary[i].Line;
	     for(c=Objects[id].Boundary[i].ColumnStart; c<=Objects[id].Boundary[i].ColumnEnd; c++){
  		     if(l>1){
			     if(raster[D1(Lines, Columns,l-1, c)]!= id){
				     per++;
 				 }
  			 }
			 else if (l==1){
			 	  per++;
		     }
		     else{
 	  		 }
 	  		 
 	  		 if(l<Lines){
	   		     if(raster[D1(Lines, Columns,l+1, c)]!= id){
			 	     per++;
			     }
			 }
			 else if(l==Lines){
			 	  per++;
 	         }
 	         else{
 	  		 }
 	  		 
 	  		 if(c>1){
			     if(raster[D1(Lines, Columns, l, c-1)]!= id){
				     per++;
 				 }
  			 }
			 else if (c==1){
			 	  per++;
		     }
		     else{
 	  		 }
			 
			 if(l<Columns){
	   		     if(raster[D1(Lines, Columns, l, c+1)]!= id){
			 	     per++;
			     }
			 }
			 else if(c==Columns){
			 	  per++;
 	         }
 	         else{
 	  		 }
		 }
	 }
	Objects[id].perimeter = per;
}

int Level::CalculatePerimeter(int id1, int id2){
	 int per = 0;
	 int l, c;
	 int val;
	 for(int i=0; i<Objects[id1].Boundary.size(); i++){
	     l = Objects[id1].Boundary[i].Line;
	     for(c=Objects[id1].Boundary[i].ColumnStart; c<=Objects[id1].Boundary[i].ColumnEnd; c++){
  		     if(l>1){
			     val = raster[D1(Lines, Columns,l-1, c)];
				 if((val != id1)&&(val != id2)){
				     per++;
 				 }
  			 }
			 else if (l==1){
			 	  per++;
		     }
		     else{
 	  		 }
 	  		 
 	  		 if(l<Lines){
	   		     val = raster[D1(Lines, Columns,l+1, c)];
				 if((val != id1)&&(val != id2)){
			 	     per++;
			     }
			 }
			 else if(l==Lines){
			 	  per++;
 	         }
 	         else{
 	  		 }
 	  		 
 	  		 if(c>1){
			     val = raster[D1(Lines, Columns, l, c-1)];
				 if((val!= id1)&&(val != id2)){
				     per++;
 				 }
  			 }
			 else if (c==1){
			 	  per++;
		     }
		     else{
 	  		 }
			 
			 if(l<Columns){
	   		     val = raster[D1(Lines, Columns, l, c+1)];
				 if((val != id1)&&(val != id2)){
			 	     per++;
			     }
			 }
			 else if(c==Columns){
			 	  per++;
 	         }
 	         else{
 	  		 }
		 }
	 }
	 
	 for(int i=0; i<Objects[id2].Boundary.size(); i++){
	     l = Objects[id2].Boundary[i].Line;
	     for(c=Objects[id2].Boundary[i].ColumnStart; c<=Objects[id2].Boundary[i].ColumnEnd; c++){
  		     if(l>1){
			     val = raster[D1(Lines, Columns,l-1, c)];
				 if((val != id1)&&(val != id2)){
				     per++;
 				 }
  			 }
			 else if (l==1){
			 	  per++;
		     }
		     else{
 	  		 }
 	  		 
 	  		 if(l<Lines){
	   		     val = raster[D1(Lines, Columns,l+1, c)];
				 if((val != id1)&&(val != id2)){
			 	     per++;
			     }
			 }
			 else if(l==Lines){
			 	  per++;
 	         }
 	         else{
 	  		 }
 	  		 
 	  		 if(c>1){
			     val = raster[D1(Lines, Columns, l, c-1)];
				 if((val!= id1)&&(val != id2)){
				     per++;
 				 }
  			 }
			 else if (c==1){
			 	  per++;
		     }
		     else{
 	  		 }
			 
			 if(l<Columns){
	   		     val = raster[D1(Lines, Columns, l, c+1)];
				 if((val != id1)&&(val != id2)){
			 	     per++;
			     }
			 }
			 else if(c==Columns){
			 	  per++;
 	         }
 	         else{
 	  		 }
		 }
	 }
	return per;
}

void Level::SaveRaster(string RasterFile){
    string HeadFile = RasterFile + ".ers";
    
    ofstream out(RasterFile.c_str(), ios::binary);
    ofstream head_out(HeadFile.c_str());
    
    //save binary data
    for(int i=0; i<raster.size(); i++){
        out.write(reinterpret_cast<char *> (&raster[i]), sizeof(int));
    }
    out.close();
    
    //save header
    head_out << "DatasetHeader Begin" << endl;
    head_out << "\tVersion \t= \"" << "6.4" << "\"" << endl;
    //head_out << "\tName\t\t= \"" << HeaderName << "\"" << endl;
    head_out << "\tLastUpdated\t= " << "Tue Mar 09 17:41:16 GMT 2004" << endl;
    head_out << "\tDataSetType\t= ERStorage" << endl;
    head_out << "\tDataType\t= " << "Raster" << endl;
    head_out << "\tByteOrder\t= " << "LSBFirst" << endl;
    head_out << "\tCoordinateSpace Begin" << endl;
    head_out << "\t\tDatum\t\t= \"" << "RAW" << "\"" << endl;
    head_out << "\t\tProjection\t= \"" << "RAW" << "\"" << endl;
    head_out << "\t\tCoordinateType\t= " << "RAW" << endl;
    head_out << "\t\tUnits\t\t= \"" << "natural" << "\"" << endl;
    head_out << "\t\tRotation\t= " << "0:0:0.0" << endl;
    head_out << "\tCoordinateSpace End" << endl;
    head_out << "\tRasterInfo Begin" << endl;
    head_out << "\t\tCellType\t= " << "Unsigned32BitInteger" << endl;
    head_out << "\t\tNrOfLines\t= " << Lines << endl;
    head_out << "\t\tNrOfCellsPerLine\t= " << Columns << endl;
    head_out << "\t\tNrOfBands\t= " << "1" << endl;
    head_out << "\t\tBandId Begin" << endl;
    head_out << "\t\t\tValue\t\t= " << "\"Raster\"" << endl;
    head_out << "\t\tBandId End" << endl;
    
    head_out << "\t\tRegionInfo Begin" << endl;
    head_out << "\t\t\tType\t\t= Polygon" << endl;
    head_out << "\t\t\tRegionName\t= \"All\"" << endl;
    //head_out << "\t\t\tSourceDataset\t= \"" << FileName << "\"" << endl;
    head_out << "\t\t\tRGBcolour Begin" << endl;
    head_out << "\t\t\t\tRed\t\t= 65535" << endl;
    head_out << "\t\t\t\tGreen\t\t= 65535" << endl;
    head_out << "\t\t\t\tBlue\t\t= 65535" << endl;
    head_out << "\t\t\tRGBcolour End" << endl;
    head_out << "\t\t\tSubRegion\t= {" << endl;
    head_out << "\t\t\t\t0\t0" << endl;
    head_out << "\t\t\t\t0\t" << Lines << endl;
    head_out << "\t\t\t\t" << Columns << "\t" << Lines << endl;
    head_out << "\t\t\t\t" << Columns << "\t0" << endl;
    head_out << "\t\t\t}" << endl;
    head_out << "\t\tRegionInfo End" << endl;
    head_out << "\tRasterInfo End" << endl;
    head_out << "DatasetHeader End" << endl;
    
    head_out.close();
    
}    

/*void Level::CreateERS(string ImageFileName, ERS_Image& Img, int BandNumber){
    
}*/    

void Level::SaveAsERS(string ImageFileName, ERS_Image& Img){
    string out_cell_type = "Unsigned8BitInteger";
    string tmp;
    
    ERS_Image Img2(ImageFileName, Img.Bands, Img.Lines, Img.Columns, out_cell_type);
    Img2.Buffer.reserve(Img2.Lines*Img2.Columns*Img2.Bands);
    for(int i=0; i<Img.BandID.size(); i++){
	    tmp = Img.BandID[i];
	    Img2.BandID.push_back(tmp);
    }
    
    int cur_id;
    CELL_TYPE val;
    float m;
    
    for(int l=1; l<=Img2.Lines; l++){
        for(int b=1; b<=Img2.Bands; b++){
            for(int c=1; c<=Img2.Columns; c++){
            	cur_id = raster[D1(Lines, Columns, l, c)];
                //cout << "Reached here" << endl;
                if(!Objects[cur_id].MeanValue.size()){
		  		    Objects[cur_id].MeanValue.reserve(Img.Bands);
		  		    Objects[cur_id].Calculate_MeanValue(Img);
				}
				//cout << "Mean value of x Object size is: " <<  << endl;
				m = Objects[cur_id].MeanValue[b-1];
                //cout << "Mean value is: " << m << endl;
                if((m-floor(m))>=0.5){
                    val = (CELL_TYPE)m + 1;
                }else{
                    val = (CELL_TYPE)m;
                }
                //val = (CELL_TYPE)m+(m-floor(m))>=0.5)?1:0;//Chiossif AnaCoding...
                //cout << "Rounded Value is: " << (int)val << endl;
				Img2.Buffer.push_back(val);
				//cout << "Current buffer size is: " << Img2.Buffer.size() << endl;       
            }
        }
    }
	cout << "Output image buffer size is: " << Img2.Buffer.size() << endl;            
    Img2.SaveAs(ImageFileName);
}    

void Level::CreateBoundaryMap(){
    BoundaryMap.clear();
    BoundaryMap.reserve(Lines*Columns);
    
    //Initialize Boundary map with zeros
    for(int i=0;i<(Lines*Columns);i++)
    	BoundaryMap.push_back(0);
    
    //Set active boundaries for the first line
    for(int c=1;c<=Columns;c++)
		BoundaryMap[D1(Lines, Columns, 1, c)]=1;
    
    //Calculate the rest
    for(int r=2;r<=(Lines-1);r++){
        BoundaryMap[D1(Lines, Columns, r, 1)]=1;//First pixel of every line
        for(int c=2;c<=(Columns-1);c++){
            BoundaryMap[D1(Lines, Columns, r, c)] = ((raster[D1(Lines, Columns, r, c)]!=raster[D1(Lines, Columns, r, c-1)])
            ||(raster[D1(Lines, Columns, r, c)]!=raster[D1(Lines, Columns, r-1, c)]))?1:0;
        }
        BoundaryMap[D1(Lines, Columns, r, Columns)]=1;//Last pixel of every line
    }
    
    //Set active boundaries for the last line
    for(int c=1;c<=Columns;c++)
    	BoundaryMap[D1(Lines, Columns, Lines, c)]=1;

}    

void Level::SaveBoundaryMapERS(string ImageFileName){
    string out_cell_type = "Unsigned8BitInteger";
    string tmp;
    tmp = "\"Band1\"";
    
    ERS_Image Img2(ImageFileName, 1, Lines, Columns, out_cell_type);
    Img2.Buffer.reserve(Img2.Lines*Img2.Columns*Img2.Bands);
    Img2.BandID.push_back(tmp);
    
    CELL_TYPE val;
    
    for(int i=0; i<(Lines*Columns); i++){
        val = (CELL_TYPE)BoundaryMap[i];
        Img2.Buffer.push_back(val);
    }
    
    Img2.SaveAs(ImageFileName, 0);
}

void Level::clear(){
    HierarchyID = 0.0;
    Lines = 0;
    Columns = 0;
    ScaleParameter = 0.0;
    Color = 0.0;
    Compactness = 0.0;
    Edge = 0.0;
    cycles = 0;
    ExistenceOfSuperlevel = FALSE;
    ExistenceOfSublevel = FALSE;
    Level* Superlevel = NULL;
    Level* Sublevel = NULL;
    
    raster.clear();
    Objects.clear();
    BoundaryMap.clear();
}    

PropertyTable::PropertyTable(int ID, int BANDS){
    id = ID;
    PropertyRow tmp;
    tmp.mean_value = 0.0;
    tmp.stdev = 0.0;
    for(int i=0; i<BANDS; i++){
        BandProperties.push_back(tmp);
    }    
}

PropertyTable::~PropertyTable(){

}

PriorityObject::PriorityObject(){
    
}    

PriorityObject::PriorityObject(int pri, int sec, int x){
    primary = pri;
    secondary = sec;
    id = x;
}

PriorityObject::~PriorityObject(){
								  
}

PriorityObject::PriorityObject(const PriorityObject& PO){
    primary = PO.primary;
    secondary = PO.secondary;
    id = PO.id;
}    

bool operator < (const PriorityObject& x, const PriorityObject& y){
    if (x.primary == y.primary){
        return (x.secondary > y.secondary);
    }    
    else {
        return (x.primary > y.primary);
    }    
}

bool operator > (const PriorityObject& x, const PriorityObject& y){
    if (x.primary == y.primary){
        return (x.secondary < y.secondary);
    }    
    else {
        return (x.primary < y.primary);
    }    
}

MSEG_Param::MSEG_Param(){
    Scale = 10.0;
    Color = 0.8;
    Compact = 0.8;
    Edge = 0.0;
    GHH = 0;
    EC = FALSE;
    MA = TRUE;    
}

MSEG_Param::MSEG_Param(float Sc, float Clr, float Cmp, float Edg, bool ec, bool ma, short ghh){
    Scale = Sc;
    Color = Clr;
    Compact = Cmp;
    Edge = Edg;
    GHH = ghh;
    MA = ma;
    EC = ec;
}

MSEG_Param::MSEG_Param(float Sc, float Clr, float Cmp){
    Scale = Sc;
    Color = Clr;
    Compact = Cmp;
    Edge = 0.0;
    GHH = 0;
    EC = FALSE;
    MA = TRUE;
}

MSEG_Param::~MSEG_Param(){    
}    

void MSEG_Init::SaveSPE(string ImageFileName, ERS_Image& Img){
    
    pixel tmp;
    string out_cell_type = "Unsigned8BitInteger";
    string tmp1;
    tmp1 = "Band1";
    
    ERS_Image Img2(ImageFileName, 1, Img.Lines, Img.Columns, out_cell_type);
    Img2.Buffer.reserve(Img2.Lines*Img2.Columns*Img2.Bands);
    Img2.BandID.push_back(tmp1);
    
    for(int i=0; i<(Img2.Lines*Img2.Columns); i++){
        Img2.Buffer.push_back(0);
    }
    
    for(int j=0; j<StartingPoints.size(); j++){
	    tmp = StartingPoints[j];
	    Img2.Buffer[D1(Img2.Lines, Img2.Columns, tmp.line, tmp.column)] = 1;
	}    
    
    Img2.SaveAs(ImageFileName, 0);
    
}    

MSEG_Init::MSEG_Init(){
    last_id = 0;
}

MSEG_Init::~MSEG_Init(){

}

void update_neighbors(Object& O, MSEG_Init& mseg){
	 int tmp, tmp2;
	 for(int i=0; i<O.Neighbors.size(); i++){
	     tmp = O.Neighbors[i].id;
	     if(mseg.Topology[tmp].merged == TRUE){
  		     tmp2 = mseg.Topology[tmp].new_id;
  		     O.Neighbors[i].id = tmp2;
	     }
	 }
}                                                                      

int firstpass(Level& pass1, ERS_Image& Img, MSEG_Init& mseg, float level_ID, MSEG_Param& Param){
    
    //cout << endl << "Inside first pass..." << endl;//Debug mode
    int mergecounter = 0;
	pass1.raster.reserve(Img.Lines*Img.Columns);
    pass1.cycles = 1;
    pass1.Lines = Img.Lines;
    pass1.Columns = Img.Columns;
    
    mseg.last_id = (Img.Lines*Img.Columns); //Set the first object id
    
    //Let all the raster become -1
    for (int i=0; i<(Img.Lines*Img.Columns); i++){
        pass1.raster.push_back(-1); //Mporei na mn xreiazetai afou exoume kanei reserve
        //pass1.raster[i] = -1;//Solution 2, without reallocating memory into the vector
    }
    //cout << "raster size = " << pass1.raster.size() << endl;//Debug mode
        
    //Copy the id number
    pass1.HierarchyID = level_ID;
    
    //Copy segmentation parameters into new level
    pass1.ScaleParameter = Param.Scale;
    pass1.Color = Param.Color;
    pass1.Compactness = Param.Compact;
    pass1.Edge = Param.Edge;
    //cout << "Parameters are: " << pass1.ScaleParameter << ", " << pass1.Color 
    //<< ", " << pass1.Compactness << endl;//Debug mode
    
    //Declare a temporary object
    Object tmp;
    //cout << "Temp Object created" << endl;//Debug mode
    //Generic Declarations
    BoundaryLine tmpBoundaryLine;
	TopologyObject tmpTopologyObject;
 	PriorityObject tmp_priority_object; 
    pixel cur, up, down, right, left, best, best_up, best_down, best_right, best_left, best_best;
    int cur_id, up_id, down_id, right_id, left_id, best_id, best_best_id;
    int best_up_id, best_down_id, best_right_id, best_left_id;
    int i, pri, sec, matches;
    float best_h, up_h, down_h, left_h, right_h, best_up_h, best_down_h, best_right_h, best_left_h, best_best_h;
    float std;
    CELL_TYPE cur_val, near_val;    
    
    //Check the SPE result number
    if (mseg.StartingPoints.size() != Img.Macroblocks.size()){
        cout << "Error in SPE result" << endl << "Quiting at first pass..." << endl;
        exit(1);
    }    
    int msize = Img.Macroblocks.size();
    
    //For every macroblock
    for(i=0;i<(Img.MacroLines * Img.MacroColumns);i++){
        //Assume each pixel is an object and assign them false ids in order to use the priority_queue
        //The pixel false id is the 1D [] representation inside the raster vector
        //cout << "First pass completed: " << ((i+1)/msize)*100 << "%" << endl;//"/r";//Debug mode
        pri = 1;//During first pass only primary is usefull.(inside macroblock)
        sec = 1;
        
        //Load the Starting Point
        cur.line = mseg.StartingPoints[i].line;
        cur.column = mseg.StartingPoints[i].column;
        
        //cout << "Starting point: " << cur.line << "," << cur.column << endl;//Debug mode
        cur_id = D1(Img.Lines, Img.Columns, cur.line, cur.column) + 1;
        //cout << "Corresponding ID for starting point is: " << cur_id << endl;//Debug mode
        
        //Check the boundary for the starting point
        if ((cur.line<Img.Macroblocks[i].LineStart || cur.line>Img.Macroblocks[i].LineEnd)
        ||(cur.column<Img.Macroblocks[i].ColumnStart || cur.column>Img.Macroblocks[i].ColumnEnd)){
            cout << "Starting point" << i << "out of bounds" << endl 
            << "Quiting at first pass..." << endl;
            exit(1);
        }
            
        //Push the first pixel to the priority_queue
        tmp_priority_object.primary = pri;
        tmp_priority_object.secondary = sec;
        tmp_priority_object.id = cur_id;
        //cout << "Priority object created with members:" << tmp_priority_object.primary
        //<< " " << tmp_priority_object.secondary << " " << tmp_priority_object.id << endl;//Debug mode
        mseg.PriorityList.push(tmp_priority_object);
        pass1.raster[cur_id - 1] = 0; //Flag that it is loaded into priority list
        //cout << "Priority object for macroblock created..." << endl;//Debug mode
        //cout << "Priority list size at the begining is " << mseg.PriorityList.size() << endl;//Debug mode
        
        //While the priority list is not empty
        while (!mseg.PriorityList.empty()){
            //cout << "Priority list currently holds " << mseg.PriorityList.size() << "id's" 
			//<< " and " << mergecounter << " merges have occured" << endl;//Debug mode
            //Load the next priority object
            tmp_priority_object = mseg.PriorityList.top();
            pri = tmp_priority_object.primary;
            sec = tmp_priority_object.secondary;
            cur_id = tmp_priority_object.id;
            //cout << "Top-ing object with id: " << cur_id << endl;//Debug mode
            //cout << "Priority object loaded with members:" << tmp_priority_object.primary
            //<< " " << tmp_priority_object.secondary << " " << tmp_priority_object.id << endl;//Debug mode 
            mseg.PriorityList.pop();//Perhaps to the end????
            //If pixel has been merged...goto end of while loop
            if(pass1.raster[cur_id - 1] > 0) continue;
            //cout << "Pop-ing object" << endl;//Debug mode
            //cout << "After pop, the priority list holds: " << mseg.PriorityList.size() << "objects" << endl;//Debug mode
			cur = D2(cur_id, Img.Columns);
			//cout << "made 2d transform" << endl;//Debug mode
            //cout << "Current pixel's x,y is: " << cur.line << " " << cur.column << endl;//Debug mode
            
            best_h = Param.Scale; //The maximum heterogeneity allowed is the scale parameter
            matches = 0;
            
            //Load neigbors
            up.line = cur.line - 1;
            up.column = cur.column;
            up_id = (up.line >= 1) ? (D1(Img.Lines, Img.Columns, up.line, up.column) + 1) : 1;
            down.line = cur.line + 1;
            down.column = cur.column;
            down_id = D1(Img.Lines, Img.Columns, down.line, down.column) + 1;
            right.line = cur.line;
            right.column = cur.column + 1;
            right_id = D1(Img.Lines, Img.Columns, right.line, right.column) + 1;
            left.line = cur.line;
            left.column = cur.column - 1;
            left_id = (left.column >= 1) ? (D1(Img.Lines, Img.Columns, left.line, left.column) + 1) : 1;
            //cout << "loading neighbors..." << endl;//Debug mode
            
            //Check neighbors
            //TODO: Put everything into arrays and loop 4 times...

            if((up.line>=Img.Macroblocks[i].LineStart && up.line<=Img.Macroblocks[i].LineEnd)
            && (up.column>=Img.Macroblocks[i].ColumnStart && up.column<=Img.Macroblocks[i].ColumnEnd)){
                if(pass1.raster[up_id - 1] <= 0){
                      //Calculate heterogeneity
                      up_h = 0.0;
                      for (int b=1; b<=Img.Bands; b++){
                          cur_val = Img.Buffer[Img.D1(cur.line, cur.column, b)];
                          near_val = Img.Buffer[Img.D1(up.line, up.column, b)];
                          std = abs(cur_val - near_val)/2;
                          up_h += Img.BandWeight[b]* 2 * std;
                      }    
                      //cout << "Heterogeneity calculated for up object " << up_h << endl;//Debug mode
                      //Is it a match? --> matches++ --> compare to best_h
                      if (up_h < Param.Scale){
                          matches++;
                          if (up_h < best_h){
                              best_h = up_h;
                              best_id = up_id;//We will copy the best outside this loop in order to spare copies
                          }
                      }        
                      //cout << "Match compared" << endl;//Debug mode
                      //if raster value <0 --> push pixel in queue
                      if (pass1.raster[up_id - 1] == -1){
                          tmp_priority_object.primary = ++pri;//pri + 1
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = up_id;
                          //cout << "Priority object for up neighbor added with members: "
                          //<< tmp_priority_object.primary << " " << tmp_priority_object.secondary
                          //<< " " << tmp_priority_object.id << endl;//Debug mode
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[up_id - 1] = 0;
                          //cout << "Neighbor pushed in priority list" << endl;//Debug mode
                          //cout << "Now priority list holds: " << mseg.PriorityList.size() << " objects" << endl;//Debug mode
                      }
                }
            }

            if((down.line>=Img.Macroblocks[i].LineStart && down.line<=Img.Macroblocks[i].LineEnd)
            && (down.column>=Img.Macroblocks[i].ColumnStart && down.column<=Img.Macroblocks[i].ColumnEnd)){
                if(pass1.raster[down_id - 1] <= 0){
                      //Calculate heterogeneity
                      down_h = 0.0;
                      
                      for (int b=1; b<=Img.Bands; b++){
                          cur_val = Img.Buffer[Img.D1(cur.line, cur.column, b)];
                          near_val = Img.Buffer[Img.D1(down.line, down.column, b)];
                          std = abs(cur_val - near_val)/2;
                          down_h += Img.BandWeight[b]* 2 * std;
                      }    
                      //cout << "Heterogeneity calculated for down object " << down_h << endl;//Debug mode
                      //Is it a match? --> matches++ --> compare to best_h
                      if (down_h < Param.Scale){
                          matches++;
                          if (down_h < best_h){
                              best_h = down_h;
                              best_id = down_id;//We will copy the best outside this loop in order to spare copies
                          }
                      }        
                      //cout << "Match compared" << endl;//Debug mode
                      //if raster value <0 --> push pixel in queue
                      if (pass1.raster[down_id - 1] == -1){
                          tmp_priority_object.primary = ++pri;//+1;
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = down_id;
                          //cout << "Priority object for down neighbor added with members: "
                          //<< tmp_priority_object.primary << " " << tmp_priority_object.secondary
                          //<< " " << tmp_priority_object.id << endl;//Debug mode
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[down_id - 1] = 0;
                          //cout << "Neighbor pushed in priority list" << endl;//Debug mode
                          //cout << "Now priority list holds: " << mseg.PriorityList.size() << " objects" << endl;//Debug mode
                      }
                }
            }
            
            if((right.line>=Img.Macroblocks[i].LineStart && right.line<=Img.Macroblocks[i].LineEnd)
            && (right.column>=Img.Macroblocks[i].ColumnStart && right.column<=Img.Macroblocks[i].ColumnEnd)){
                if(pass1.raster[right_id - 1] <= 0){
                      //Calculate heterogeneity
                      right_h = 0.0;
                      
                      for (int b=1; b<=Img.Bands; b++){
                          cur_val = Img.Buffer[Img.D1(cur.line, cur.column, b)];
                          near_val = Img.Buffer[Img.D1(right.line, right.column, b)];
                          std = abs(cur_val - near_val)/2;
                          right_h += Img.BandWeight[b]* 2 * std;
                      }    
                      //cout << "Heterogeneity calculated for right object " << right_h << endl;//Debug mode
                      //Is it a match? --> matches++ --> compare to best_h
                      if (right_h < Param.Scale){
                          matches++;
                          if (right_h < best_h){
                              best_h = right_h;
                              best_id = right_id;//We will copy the best outside this loop in order to spare copies
                          }
                      }        
                      //cout << "Match compared" << endl;//Debug mode
                      //if raster value <0 --> push pixel in queue
                      if (pass1.raster[right_id - 1] == -1){
                          tmp_priority_object.primary = ++pri;//+1;
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = right_id;
                          //cout << "Priority object for right neighbor added with members: "
                          //<< tmp_priority_object.primary << " " << tmp_priority_object.secondary
                          //<< " " << tmp_priority_object.id << endl;//Debug mode
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[right_id - 1] = 0;
                          //cout << "Neighbor pushed in priority list" << endl;//Debug mode
                          //cout << "Now priority list holds: " << mseg.PriorityList.size() << " objects" << endl;//Debug mode
                      }
                }
            }
            
            if((left.line>=Img.Macroblocks[i].LineStart && left.line<=Img.Macroblocks[i].LineEnd)
            && (left.column>=Img.Macroblocks[i].ColumnStart && left.column<=Img.Macroblocks[i].ColumnEnd)){
                if(pass1.raster[left_id - 1] <= 0){
                      //Calculate heterogeneity
                      left_h = 0.0;
                      
                      for (int b=1; b<=Img.Bands; b++){
                          cur_val = Img.Buffer[Img.D1(cur.line, cur.column, b)];
                          near_val = Img.Buffer[Img.D1(left.line, left.column, b)];
                          std = abs(cur_val - near_val)/2;
                          left_h += Img.BandWeight[b]* 2 * std;
                      }    
                      //cout << "Heterogeneity calculated for left object " << left_h << endl;//Debug mode
                      //Is it a match? --> matches++ --> compare to best_h
                      if (left_h < Param.Scale){
                          matches++;
                          if (left_h < best_h){
                              best_h = left_h;
                              best_id = left_id;//We will copy the best outside this loop in order to spare copies
                          }
                      }        
                      //cout << "Match compared" << endl;//Debug mode
                      //if raster value <0 --> push pixel in queue
                      if (pass1.raster[left_id - 1] == -1){
                          tmp_priority_object.primary = pri+1;
                          tmp_priority_object.secondary = sec;
                          tmp_priority_object.id = left_id;
                          //cout << "Priority object for left neighbor added with members: "
                          //<< tmp_priority_object.primary << " " << tmp_priority_object.secondary
                          //<< " " << tmp_priority_object.id << endl;//Debug mode
                          mseg.PriorityList.push(tmp_priority_object);
                          pass1.raster[left_id - 1] = 0;
                          //cout << "Neighbor pushed in priority list" << endl;//Debug mode
                          //cout << "Now priority list holds: " << mseg.PriorityList.size() << " objects" << endl;//Debug mode
                      }
                }
            }
            
            //Push it's neigbors into the priority list if they are not merged (done above)
            //cout << "Done with neighbors..." << " Best match id is: " << best_id << endl;//Debug mode
            //cout << "matches at this point are: " << matches << endl;//Debug mode
            //Check for match
            if(matches){
                //Check for the best match
                //cout << "Inside scope for best match..." << endl;//Debug mode
                best = D2(best_id, Img.Columns);
                //cout << "Best match x,y: " << best.line << " " << best.column << endl;//Debug mode
                best_best_h = Param.Scale; //The maximum heterogeneity allowed is the scale parameter
                
               	//Load neigbors
               	best_up.line = best.line-1;
               	best_up.column = best.column;
               	best_up_id = (best_up.line >= 1) ? (D1(Img.Lines, Img.Columns, best_up.line, best_up.column) + 1) : 1;
                best_down.line = best.line + 1;
               	best_down.column = best.column;
               	best_down_id = D1(Img.Lines, Img.Columns, best_down.line, best_down.column) + 1;
               	best_right.line = best.line;
               	best_right.column = best.column + 1;
               	best_right_id = D1(Img.Lines, Img.Columns, best_right.line, best_right.column) + 1;
               	best_left.line = best.line;
               	best_left.column = best.column -1;
               	best_left_id = (best_left.column >= 1) ? (D1(Img.Lines, Img.Columns, best_left.line, best_left.column) + 1) : 1;
                //cout << "Neighbors loaded again..." << endl;//Debug mode
                //Check neighbors
               	if((best_up.line>=Img.Macroblocks[i].LineStart && best_up.line<=Img.Macroblocks[i].LineEnd)
               	&& (best_up.column>=Img.Macroblocks[i].ColumnStart && best_up.column<=Img.Macroblocks[i].ColumnEnd)){
                	if(pass1.raster[best_up_id - 1] <= 0){
                      	//Calculate heterogeneity
                      	best_up_h = 0.0;
                      
                     	for (int b=1; b<=Img.Bands; b++){
                          	cur_val = Img.Buffer[Img.D1(best.line, best.column, b)];
                          	near_val = Img.Buffer[Img.D1(best_up.line, best_up.column, b)];
                          	std = abs(cur_val - near_val)/2;
                          	best_up_h += Img.BandWeight[b]* 2 * std;
                     	}    
                      	//cout << "Inside best_up... Heterogeneity: " << best_up_h << endl;//Debug mode
                        //Is it a match? --> compare to best_h
                      	if (best_up_h < best_best_h){
                          	best_best_h = best_up_h;
                          	best_best_id = best_up_id;//We will copy the best outside this loop in order to spare copies
                        }
                        //cout << "calculated best_up heterogeneity" << endl;//Debug mode
                        //if raster value <0 --> push pixel in queue
                        if (pass1.raster[best_up_id - 1] == -1){
                            tmp_priority_object.primary = ++pri;//pri + 1
                            tmp_priority_object.secondary = sec;
                            tmp_priority_object.id = best_up_id;
                            //cout << "Priority object for up neighbor added with members: "
                            //<< tmp_priority_object.primary << " " << tmp_priority_object.secondary
                            //<< " " << tmp_priority_object.id << endl;//Debug mode
                            mseg.PriorityList.push(tmp_priority_object);
                            pass1.raster[best_up_id - 1] = 0;
                            //cout << "Neighbor pushed in priority list" << endl;//Debug mode
                            //cout << "Now priority list holds: " << mseg.PriorityList.size() << " objects" << endl;//Debug mode
                        }
                    }
                 }
                 
                 if((best_down.line>=Img.Macroblocks[i].LineStart && best_down.line<=Img.Macroblocks[i].LineEnd)
               	 && (best_down.column>=Img.Macroblocks[i].ColumnStart && best_down.column<=Img.Macroblocks[i].ColumnEnd)){
                	if(pass1.raster[best_down_id - 1] <= 0){
                      	//Calculate heterogeneity
                      	best_down_h = 0.0;
                      
                     	for (int b=1; b<=Img.Bands; b++){
                          	cur_val = Img.Buffer[Img.D1(best.line, best.column, b)];
                          	near_val = Img.Buffer[Img.D1(best_down.line, best_down.column, b)];
                          	std = abs(cur_val - near_val)/2;
                          	best_down_h += Img.BandWeight[b]* 2 * std;
                     	}    
                      	//cout << "Inside best_down... Heterogeneity: " << best_down_h << endl;//Debug mode
                        //Is it a match? --> compare to best_h
                      	if (best_down_h < best_best_h){
                          	best_best_h = best_down_h;
                          	best_best_id = best_down_id;//We will copy the best outside this loop in order to spare copies
                        }
                        //cout << "calculated best_down heterogeneity" << endl;//Debug mode
                        //if raster value <0 --> push pixel in queue
                        if (pass1.raster[best_down_id - 1] == -1){
                            tmp_priority_object.primary = ++pri;//pri + 1
                            tmp_priority_object.secondary = sec;
                            tmp_priority_object.id = best_down_id;
                            //cout << "Priority object for down neighbor added with members: "
                            //<< tmp_priority_object.primary << " " << tmp_priority_object.secondary
                            //<< " " << tmp_priority_object.id << endl;//Debug mode
                            mseg.PriorityList.push(tmp_priority_object);
                            pass1.raster[best_down_id - 1] = 0;
                            //cout << "Neighbor pushed in priority list" << endl;//Debug mode
                            //cout << "Now priority list holds: " << mseg.PriorityList.size() << " objects" << endl;//Debug mode
                        }
                    }
                 }
                 
                 if((best_right.line>=Img.Macroblocks[i].LineStart && best_right.line<=Img.Macroblocks[i].LineEnd)
               	 && (best_right.column>=Img.Macroblocks[i].ColumnStart && best_right.column<=Img.Macroblocks[i].ColumnEnd)){
                	if(pass1.raster[best_right_id - 1] <= 0){
                      	//Calculate heterogeneity
                      	best_right_h = 0.0;
                      
                     	for (int b=1; b<=Img.Bands; b++){
                          	cur_val = Img.Buffer[Img.D1(best.line, best.column, b)];
                          	near_val = Img.Buffer[Img.D1(best_right.line, best_right.column, b)];
                          	std = abs(cur_val - near_val)/2;
                          	best_right_h += Img.BandWeight[b]* 2 * std;
                     	}    
                      	//cout << "Inside best_right... Heterogeneity: " << best_right_h << endl;//Debug mode
                        //Is it a match? --> compare to best_h
                      	if (best_right_h < best_best_h){
                          	best_best_h = best_right_h;
                          	best_best_id = best_right_id;//We will copy the best outside this loop in order to spare copies
                        }
                        //cout << "calculated best_right heterogeneity" << endl;//Debug mode
                        //if raster value <0 --> push pixel in queue
                        if (pass1.raster[best_right_id - 1] == -1){
                            tmp_priority_object.primary = ++pri;//pri + 1
                            tmp_priority_object.secondary = sec;
                            tmp_priority_object.id = best_right_id;
                            //cout << "Priority object for right neighbor added with members: "
                            //<< tmp_priority_object.primary << " " << tmp_priority_object.secondary
                            //<< " " << tmp_priority_object.id << endl;//Debug mode
                            mseg.PriorityList.push(tmp_priority_object);
                            pass1.raster[best_right_id - 1] = 0;
                            //cout << "Neighbor pushed in priority list" << endl;//Debug mode
                            //cout << "Now priority list holds: " << mseg.PriorityList.size() << " objects" << endl;//Debug mode
                        }
                    }
                 }
                 
                 if((best_left.line>=Img.Macroblocks[i].LineStart && best_left.line<=Img.Macroblocks[i].LineEnd)
               	 && (best_left.column>=Img.Macroblocks[i].ColumnStart && best_left.column<=Img.Macroblocks[i].ColumnEnd)){
                	if(pass1.raster[best_left_id - 1] <= 0){
                      	//Calculate heterogeneity
                      	best_left_h = 0.0;
                      
                     	for (int b=1; b<=Img.Bands; b++){
                          	cur_val = Img.Buffer[Img.D1(best.line, best.column, b)];
                          	near_val = Img.Buffer[Img.D1(best_left.line, best_left.column, b)];
                          	std = abs(cur_val - near_val)/2;
                          	best_left_h += Img.BandWeight[b]* 2 * std;
                     	}    
                      	//cout << "Inside best_left... Heterogeneity: " << best_left_h << endl;//Debug mode
                        //Is it a match? --> compare to best_h
                      	if (best_left_h < best_best_h){
                          	best_best_h = best_left_h;
                          	best_best_id = best_left_id;//We will copy the best outside this loop in order to spare copies
                        }
                        //cout << "calculated best_left heterogeneity" << endl;//Debug mode
                        //if raster value <0 --> push pixel in queue
                        if (pass1.raster[best_left_id - 1] == -1){
                            tmp_priority_object.primary = ++pri;//pri + 1
                            tmp_priority_object.secondary = sec;
                            tmp_priority_object.id = best_left_id;
                            //cout << "Priority object for left neighbor added with members: "
                            //<< tmp_priority_object.primary << " " << tmp_priority_object.secondary
                            //<< " " << tmp_priority_object.id << endl;//Debug mode
                            mseg.PriorityList.push(tmp_priority_object);
                            pass1.raster[best_left_id - 1] = 0;
                            //cout << "Neighbor pushed in priority list" << endl;//Debug mode
                            //cout << "Now priority list holds: " << mseg.PriorityList.size() << " objects" << endl;//Debug mode
                        }
                    }
                 }
                
                //Return best best match --> DONE(best_best_id)
                
                //Is it mutual best match??
                if(cur_id == best_best_id){
               		//If yes create object, return value --> Update raster, clear temp object
               		//cout << "Inside mutual best match..." << endl;//Debug mode
                 	mseg.last_id++;//add 1 to last id to use to the current object
               		tmp.id = mseg.last_id;
               		tmp.area = 2;
               		mergecounter++;
	   				//tmp.cycle = 1;
               		if(cur.line == best.line){//the object will have one BoundaryLine object
               		    tmpBoundaryLine.Line = cur.line;
                     	if(cur.column<best.column){
                     		tmpBoundaryLine.ColumnStart = cur.column;
                     		tmpBoundaryLine.ColumnEnd = best.column;
                   		}
                     	else{
                      		tmpBoundaryLine.ColumnStart = best.column;
                      		tmpBoundaryLine.ColumnEnd = cur.column;
                       	}
                       	tmp.Boundary.push_back(tmpBoundaryLine);
                	}
                 	else if(cur.line < best.line){//the object will have two BoundaryLine objects
                      	tmpBoundaryLine.Line = cur.line;
                      	tmpBoundaryLine.ColumnStart = cur.column;
                      	tmpBoundaryLine.ColumnEnd = cur.column;
                      	tmp.Boundary.push_back(tmpBoundaryLine);
                      	
                      	tmpBoundaryLine.Line = best.line;
                      	tmpBoundaryLine.ColumnStart = best.column;
                      	tmpBoundaryLine.ColumnEnd = best.column;
                      	tmp.Boundary.push_back(tmpBoundaryLine);
                    }
                    else {
                        tmpBoundaryLine.Line = best.line;
                      	tmpBoundaryLine.ColumnStart = best.column;
                      	tmpBoundaryLine.ColumnEnd = best.column;
                      	tmp.Boundary.push_back(tmpBoundaryLine);
                      	
                      	tmpBoundaryLine.Line = cur.line;
                      	tmpBoundaryLine.ColumnStart = cur.column;
                      	tmpBoundaryLine.ColumnEnd = cur.column;
                      	tmp.Boundary.push_back(tmpBoundaryLine);
                    }
                    pass1.Objects.insert(make_pair(tmp.id, tmp));//insert will work better here...
                    //cout << "Object created succesfully..." << endl;//Debug mode
                    tmpTopologyObject.id = tmp.id;
                    tmpTopologyObject.merged = FALSE;
                    tmpTopologyObject.new_id = 0;
                    mseg.Topology.insert(make_pair(tmp.id, tmpTopologyObject));//insert here to
                    //cout << "Topology objects created..." << endl;//Debug mode
                    pass1.raster[cur_id - 1] = tmp.id;
                    pass1.raster[best_id - 1] = tmp.id;
                    tmp.clear();
                    //cout << "Temp Object cleared..." << endl;//Debug mode
          		}  		
               	else{
                    //If no mutual best match, push back in queue in priority +10 and move to next
                    //cout << "No mutual best match found, pushing back for later..." << endl;//Debug mode
                    tmp_priority_object.primary = pri+20;
                    tmp_priority_object.secondary = sec;
                    tmp_priority_object.id = cur_id;
                    mseg.PriorityList.push(tmp_priority_object);
                }    
            }
            else{//if cur has no matches
                //cout << "No matches found, creating one pixel object..." << endl;//Debug mode
                //create single pixel object --> Update raster, clear temp object
                Object tmp(cur_id);
                tmp.area = 1;
                //tmp.cycle = 1;
                tmpBoundaryLine.Line = cur.line;
                tmpBoundaryLine.ColumnStart = cur.column;
                tmpBoundaryLine.ColumnEnd = cur.column;
                tmp.Boundary.push_back(tmpBoundaryLine);
                pass1.Objects.insert(make_pair(tmp.id, tmp));//insert here to
                tmpTopologyObject.id = tmp.id;
                tmpTopologyObject.merged = FALSE;
                tmpTopologyObject.new_id = 0;
                mseg.Topology.insert(make_pair(tmp.id, tmpTopologyObject));//insert here to
                pass1.raster[cur_id - 1] = cur_id;
                tmp.clear();
            }      		
        }    
        //TODO:Check if some pixels have not been treated and create single pixel objects --> Update raster
    }        
    //Create topology --> Load Topology
    cout << endl << "Calculating overall topology now..." << endl;
	Neighbor n;
    n.merging_heterogeneity = 0.0;
    for(int row=1; row<=Img.Lines; row++){
        for(int col=1; col<=Img.Columns; col++){
            cur.line = row;
            cur.column = col;
            up.line = cur.line-1;
            up.column = cur.column;
            down.line = cur.line + 1;
            down.column = down.column;
            right.line = cur.line;
            right.column = cur.column + 1;
            left.line = cur.line;
            left.column = cur.column -1;
            
            cur_id = pass1.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
            if ((up.line>0 && up.line<=Img.Lines) && (up.column>0 && up.column<=Img.Columns)){
                up_id = pass1.raster[D1(Img.Lines, Img.Columns, up.line, up.column)];
                if(cur_id != up_id){
                    n.id = up_id;
                    if(pass1.Objects[cur_id].check_neighbor(up_id) == FALSE){
                        pass1.Objects[cur_id].Neighbors.push_back(n);
                    }    
                }    
            }
            if ((down.line>0 && down.line<=Img.Lines) && (down.column>0 && down.column<=Img.Columns)){
                down_id = pass1.raster[D1(Img.Lines, Img.Columns, down.line, down.column)];
                if(cur_id != down_id){
                    n.id = down_id;
                    if(pass1.Objects[cur_id].check_neighbor(down_id) == FALSE){//Maybe disable and use fix() for speed...???
                        pass1.Objects[cur_id].Neighbors.push_back(n);
                    }
                }    
            }
            if ((right.line>0 && right.line<=Img.Lines) && (right.column>0 && right.column<=Img.Columns)){
                right_id = pass1.raster[D1(Img.Lines, Img.Columns, right.line, right.column)];
                if(cur_id != right_id){
                    n.id = right_id;
                    if(pass1.Objects[cur_id].check_neighbor(right_id) == FALSE){
                        pass1.Objects[cur_id].Neighbors.push_back(n);
                    }
                }    
            }
            if ((left.line>0 && left.line<=Img.Lines) && (left.column>0 && left.column<=Img.Columns)){
                left_id = pass1.raster[D1(Img.Lines, Img.Columns, left.line, left.column)];
                if(cur_id != left_id){
                    n.id = left_id;
                    if(pass1.Objects[cur_id].check_neighbor(left_id) == FALSE){
                        pass1.Objects[cur_id].Neighbors.push_back(n);
                    }
                }    
            }    
        }
    }
	cout << "First pass completed succesfully!" << endl;
	cout << "Objects allocated in first pass Level: " << pass1.Objects.size() << endl;
	return mergecounter;        
}

/*void firstpass(Level& pass1, ERS_Image& Img, MSEG_Init& mseg, float level_ID, MSEG_Param& Param, Level& super_level){
    
} */   


int secondpass(Level& pass2, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass){
    
	int mergecounter = 0;
	pass2.raster.reserve(Img.Lines*Img.Columns);
	//cout << "Raster space reserved" << endl;//Debug mode
	for(int j=0; j<(Img.Lines*Img.Columns); j++){
	    pass2.raster.push_back(-1);
	}    
	pass2.cycles = 2;
    pass2.Lines = Img.Lines;
    pass2.Columns = Img.Columns;
    
    //Copy the id number
    pass2.HierarchyID = previous_pass.HierarchyID;
    
    //Copy segmentation parameters into new level
    pass2.ScaleParameter = previous_pass.ScaleParameter;
    pass2.Color = previous_pass.Color;
    pass2.Compactness = previous_pass.Compactness;
    pass2.Edge = previous_pass.Edge;
    //cout << "Segmentation parameters: " << pass2.ScaleParameter << " " << pass2.Color
    //<< " " << pass2.Compactness << endl;//Debug mode
    
    //Declare a temporary object
    Object tmp;
    pixel cur;
    int pri, sec, matches;
    int cur_id, near_id, best_id, best_best_id;
    float best_f, best_best_f, h_color, h_shape, h_smooth, h_cmpct, f;
    PriorityObject tmp_priority_object;
    TopologyObject tmpTopologyObject;
    //Load the Starting Points into Priority List
    pri = 1;
    sec = 0;
    for(int i=0; i<mseg.StartingPoints.size(); i++){
  		cur = mseg.StartingPoints[i];
  		sec = i+1;
  		cur_id = previous_pass.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
  		tmp_priority_object.id = cur_id;
  		tmp_priority_object.primary = pri;
  		tmp_priority_object.secondary = sec;
  		mseg.PriorityList.push(tmp_priority_object);
  		previous_pass.Objects[cur_id].priority_loaded = TRUE;//Flag that it is loaded into priority list
    }
    //cout << "SP loaded into Priority list..." << endl << "Priority list holds: "
    //<< mseg.PriorityList.size() << " objects" << endl;//Debug mode
    
    //While the priority list is not empty
    while (!mseg.PriorityList.empty()){
	    
		//Load the next priority object
     	tmp_priority_object = mseg.PriorityList.top();
      	pri = tmp_priority_object.primary;
       	sec = tmp_priority_object.secondary;
        cur_id = tmp_priority_object.id;
        mseg.PriorityList.pop();
        //cout << "Priority list holds " << mseg.PriorityList.size() << " objects after pop"
        //<< " and merges occured = " << mergecounter << endl;//Debug mode
       	//Reload Priority object if cur already merged
        if(previous_pass.Objects[cur_id].merged == TRUE) continue;
            
        //If not calculated, calculate object's perimeter
		
		best_f = pass2.ScaleParameter; //The maximum heterogeneity allowed is the scale parameter
		matches = 0;
		
        //Load neigbors --> DONE: dahh, we have objects now...
		
        //For every neighbor
		for(int i=0; i<previous_pass.Objects[cur_id].Neighbors.size(); i++){
            near_id = previous_pass.Objects[cur_id].Neighbors[i].id;
            //cout << "Processing neighbor " << i+1 << endl;//Debug mode
			//if neighbor not pushed to priority --> Push
			if(!previous_pass.Objects[near_id].priority_loaded){
			    tmp_priority_object.id = near_id;
			    tmp_priority_object.primary = pri + 1;
			    tmp_priority_object.secondary = sec;
			    mseg.PriorityList.push(tmp_priority_object);
			    previous_pass.Objects[near_id].priority_loaded = TRUE;
			    //cout << "Neighbor pushed into priority list" << endl;//Debug mode
			}    
			if(!previous_pass.Objects[near_id].merged){//if Neighbor not merged
				 //cout << "Neighbor not merged, inside if..." << endl;//Debug mode
				 if (!previous_pass.Objects[cur_id].virtual_merged){//if neighbor merging heterogeneity not calculated
	 	             //Virtual merge into temp object
                     tmp.clear();
                     tmp.merge(previous_pass.Objects[cur_id], previous_pass.Objects[near_id]);
                     //cout << "Objects virtual merged..." << endl;//Debug mode
                     h_color = tmp.Color_Heterogeneity(Img) - previous_pass.Objects[cur_id].Color_Heterogeneity(Img) - previous_pass.Objects[near_id].Color_Heterogeneity(Img);
                     //h_smooth = tmp.Smoothness() - previous_pass.Objects[cur_id].Smoothness() - previous_pass.Objects[near_id].Smoothness();
                     //h_cmpct = tmp.Compactness() - previous_pass.Objects[cur_id].Compactness() - previous_pass.Objects[near_id].Compactness();
                     //h_shape = (pass2.Compactness * h_cmpct) + ((1-pass2.Compactness)*h_smooth);
                     //f = (pass2.Color*h_color) + ((1-pass2.Color)*h_shape);
                     f = h_color;
                     //cout << "Heterogeneity for neighbor " << i+1 << " calculated: " << f << endl;//Debug mode
                 }
                 else{
                     //get the heterogeneity already calculated
                     //This is not working properly yet, so we don't realy use it here...
                     //cout << "This message should not be printed..." << endl;//Debug mode
                     f = previous_pass.Objects[cur_id].Neighbors[i].merging_heterogeneity;
                 }        
   			     if(f<best_f){
      		     	matches++;
      		     	//cout << "It's a match..." << endl;//Debug mode
  			   		best_f = f;
  			   		best_id = near_id;
   		         }
   		     }    
	 	}
	 	//cout << " Total matches: " << matches << endl;//Debug mode
        //Check for match
	 	if(matches){
	 	    //Load best match
	 	    best_best_f = pass2.ScaleParameter;
	 	    //For every neighbor
	 	    for(int i=0; i<previous_pass.Objects[best_id].Neighbors.size(); i++){
	 	        near_id = previous_pass.Objects[best_id].Neighbors[i].id;
	 	        //cout << "Processing best neighbor " << i+1 << endl;//Debug mode
	 	        //if neighbor not pushed to priority --> Push
	 	        if(!previous_pass.Objects[near_id].priority_loaded){
	 	        	tmp_priority_object.id = near_id;
	 	        	tmp_priority_object.primary = pri + 1;
	 	        	tmp_priority_object.secondary = sec;
	 	        	mseg.PriorityList.push(tmp_priority_object);
	 	        	previous_pass.Objects[near_id].priority_loaded = TRUE;
	 	        	//cout << "Neighbor pushed into priority list" << endl;//Debug mode
 	        	}
	 	        if(!previous_pass.Objects[near_id].merged){//if Neighbor not merged
                    //cout << "Neighbor not merged, inside if..." << endl;//Debug mode
                    if (!previous_pass.Objects[best_id].virtual_merged){//if neighbor merging heterogeneity not calculated
                        //Virtual merge into temp object
                        tmp.clear();
                        tmp.merge(previous_pass.Objects[best_id], previous_pass.Objects[near_id]);
                        h_color = tmp.Color_Heterogeneity(Img) - previous_pass.Objects[best_id].Color_Heterogeneity(Img) - previous_pass.Objects[near_id].Color_Heterogeneity(Img);
                        //h_smooth = tmp.Smoothness() - previous_pass.Objects[best_id].Smoothness() - previous_pass.Objects[near_id].Smoothness();
                        //h_cmpct = tmp.Compactness() - previous_pass.Objects[best_id].Compactness() - previous_pass.Objects[near_id].Compactness();
                        //h_shape = (pass2.Compactness * h_cmpct) + ((1-pass2.Compactness)*h_smooth);
                        //f = (pass2.Color*h_color) + ((1-pass2.Color)*h_shape);
                        f = h_color;
                        //cout << "Heterogeneity for neighbor " << i+1 << " calculated: " << f << endl;//Debug mode
                    }
                    else{
                        //get the heterogeneity already calculated
                        //cout << "This message should not be printed..." << endl;//Debug mode
                        f = previous_pass.Objects[best_id].Neighbors[i].merging_heterogeneity;
                    }
                    if(f<best_best_f){
                    	best_best_f = f;
                    	best_best_id = near_id;
                   	}
                }    
	 	    }
	 	    //Is it mutual best match???
      		if(best_best_id == cur_id){
      		    mark1:
            	//cout << "Inside mutual best match..." << endl;//Debug mode
            	tmp.clear();
            	//merge objects
            	tmp.merge(previous_pass.Objects[best_id], previous_pass.Objects[cur_id]);
            	//calculate perimeter
            	tmp.perimeter = previous_pass.CalculatePerimeter(best_id, cur_id);
            	tmp.merge_neighbors(previous_pass.Objects[best_id], previous_pass.Objects[cur_id]);
            	//cout << "A merge has just occured..." << endl;//Debug mode
				mergecounter++;
				mseg.last_id++;
             	tmp.id = mseg.last_id;
             	tmp.priority_loaded = FALSE;
             	pass2.Objects.insert(make_pair(tmp.id, tmp));
             	
                tmpTopologyObject.id = tmp.id;
                tmpTopologyObject.merged = FALSE;
                tmpTopologyObject.new_id = 0;
                mseg.Topology.insert(make_pair(tmp.id, tmpTopologyObject));
                mseg.Topology[best_id].merged = TRUE;
                mseg.Topology[best_id].new_id = tmp.id;
                mseg.Topology[cur_id].merged = TRUE;
                mseg.Topology[cur_id].new_id = tmp.id;
                
                previous_pass.Objects[cur_id].merged = TRUE;
                previous_pass.Objects[best_id].merged = TRUE;
                
                //Update raster
                pass2.UpdateRaster(tmp.id);
                
                tmp.clear();
                //cout << "Now level holds " << pass2.Objects.size() << " objects" << endl;//Debug mode
      		}
        	else{
             	//push back to priority list
             	//cout << "not mutual best match, pushing back to priority list..." << endl;//Debug mode 
             	if(mseg.Multi_Pushed.find(cur_id)!= mseg.Multi_Pushed.end()){
             	    mseg.Multi_Pushed[cur_id]++;
             	    if(mseg.Multi_Pushed[cur_id] >=5) goto mark1;
             	}else{
                  	mseg.Multi_Pushed.insert(make_pair(cur_id, 1));
                }        
              	tmp_priority_object.id = cur_id;
			    tmp_priority_object.primary = pri + 20;
			    //cout << "Primary priority is now: " << tmp_priority_object.primary << endl;//Debug mode
			    tmp_priority_object.secondary = sec;
			    mseg.PriorityList.push(tmp_priority_object);
        	}        
        }
        else{//no matches
            //push the same object to the next level
            pass2.Objects.insert(make_pair(cur_id, previous_pass.Objects[cur_id]));
            previous_pass.Objects[cur_id].merged = TRUE;
            pass2.Objects[cur_id].priority_loaded = FALSE;
            //cout << "No matches, object copied to the next level" << endl;//Debug mode
            //Update raster
            pass2.UpdateRaster(cur_id);
        }    
        
    }
    return mergecounter; 
}

/*void secondpass(Level& pass2, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass, Level& super_level){
     
}*/


int nthpass(Level& passn, ERS_Image& Img, float level_ID, MSEG_Init& mseg, MSEG_Param& Param, Level& sub_level){
    
	int mergecounter = 0;
	passn.raster.reserve(Img.Lines*Img.Columns);
    passn.cycles = 1;
    
    //Copy the id number
    passn.HierarchyID = level_ID;
    passn.Lines = Img.Lines;
    passn.Columns = Img.Columns;
    
    //Copy segmentation parameters into new level
    passn.ScaleParameter = Param.Scale;
    passn.Color = Param.Color;
    passn.Compactness = Param.Compact;
    passn.Edge = Param.Edge;
    
    //Declare a temporary object
    Object tmp;
    pixel cur;
    int pri, sec, matches;
    int cur_id, near_id, best_id, best_best_id;
    float best_f, best_best_f, h_color, h_shape, h_smooth, h_cmpct, f;
    PriorityObject tmp_priority_object(0,0,0);
    TopologyObject tmpTopologyObject;
    TopologyPair tmpTopologyPair;
    
    //Load the Starting Points into Priority List
    pri = 1;
    sec = 0;
    for(int i=0; i<mseg.StartingPoints.size(); i++){
  		cur = mseg.StartingPoints[i];
  		sec = i+1;
  		cur_id = sub_level.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
  		tmp_priority_object.id = cur_id;
  		tmp_priority_object.primary = pri;
  		tmp_priority_object.secondary = sec;
  		mseg.PriorityList.push(tmp_priority_object);
  		sub_level.Objects[cur_id].priority_loaded = TRUE;//Flag that it is loaded into priority list
    }
    
    //While the priority list is not empty
    while (!mseg.PriorityList.empty()){
	    
		//Load the next priority object
     	tmp_priority_object = mseg.PriorityList.top();
      	pri = tmp_priority_object.primary;
       	sec = tmp_priority_object.secondary;
        cur_id = tmp_priority_object.id;
        mseg.PriorityList.pop();
       	
        //Reload Priority object if cur already merged
        if(sub_level.Objects[cur_id].merged == TRUE) continue;
        
        //If not calculated, calculate object's perimeter
        if(sub_level.Objects[cur_id].perimeter == 0){
  		    sub_level.CalculatePerimeter(cur_id);
		}
		//Update and Fix the neighbors for the previous level in order to call once bellow
		update_neighbors(sub_level.Objects[cur_id], mseg);
		sub_level.Objects[cur_id].fix_neighbors();
		
		best_f = passn.ScaleParameter; //The maximum heterogeneity allowed is the scale parameter
		matches = 0;
		
        //Load neigbors --> DONE: dahh, we have objects now...
		
        //For every neighbor
		for(int i=0; i<sub_level.Objects[cur_id].Neighbors.size(); i++){
            near_id = sub_level.Objects[cur_id].Neighbors[i].id;
            
            //If not calculated, calculate neighbor object's perimeter
            if(sub_level.Objects[near_id].perimeter == 0){
         	    sub_level.CalculatePerimeter(near_id);
		    }
		    update_neighbors(sub_level.Objects[near_id], mseg);
		    sub_level.Objects[near_id].fix_neighbors();
		    
			//if neighbor not pushed to priority --> Push
			if(!sub_level.Objects[near_id].priority_loaded){
			    tmp_priority_object.id = near_id;
			    tmp_priority_object.primary = pri + 1;
			    tmp_priority_object.secondary = sec;
			    mseg.PriorityList.push(tmp_priority_object);
			    sub_level.Objects[near_id].priority_loaded = TRUE;
			}    
			if(!sub_level.Objects[near_id].merged){//if Neighbor not merged
				 if (!sub_level.Objects[cur_id].virtual_merged){//if neighbor merging heterogeneity not calculated
	 	             //Virtual merge into temp object
                     tmp.clear();
                     tmp.merge(sub_level.Objects[cur_id], sub_level.Objects[near_id]);
                     //Calculate merge object's perimeter
                     tmp.perimeter = sub_level.CalculatePerimeter(cur_id, near_id);
                     
                     h_color = tmp.Color_Heterogeneity(Img) - sub_level.Objects[cur_id].Color_Heterogeneity(Img) - sub_level.Objects[near_id].Color_Heterogeneity(Img);
                     h_smooth = tmp.Smoothness() - sub_level.Objects[cur_id].Smoothness() - sub_level.Objects[near_id].Smoothness();
                     h_cmpct = tmp.Compactness() - sub_level.Objects[cur_id].Compactness() - sub_level.Objects[near_id].Compactness();
                     h_shape = (passn.Compactness * h_cmpct) + ((1-passn.Compactness)*h_smooth);
                     f = (passn.Color*h_color) + ((1-passn.Color)*h_shape);
                 }
                 else{
                     //get the heterogeneity already calculated
                     f = sub_level.Objects[cur_id].Neighbors[i].merging_heterogeneity;
                 }        
   			     if(f<best_f){
      		     	matches++;
  			   		best_f = f;
  			   		best_id = near_id;
   		         }
   		     }    
	 	}
	 	
        //Check for match
	 	if(matches){
	 	    //Load best match
	 	    best_best_f = passn.ScaleParameter;
	 	    //For every neighbor
	 	    for(int i=0; i<sub_level.Objects[best_id].Neighbors.size(); i++){
	 	        near_id = sub_level.Objects[best_id].Neighbors[i].id;
	 	        //If not calculated, calculate neighbor object's perimeter
            	if(sub_level.Objects[near_id].perimeter == 0){
			         sub_level.CalculatePerimeter(near_id);
	            }
	            update_neighbors(sub_level.Objects[near_id], mseg);
	            sub_level.Objects[near_id].fix_neighbors();
		    
		        //if neighbor not pushed to priority --> Push
		        if(!sub_level.Objects[near_id].priority_loaded){
		        	tmp_priority_object.id = near_id;
		        	tmp_priority_object.primary = pri + 1;
		        	tmp_priority_object.secondary = sec;
		        	mseg.PriorityList.push(tmp_priority_object);
		        	sub_level.Objects[near_id].priority_loaded = TRUE;
	        	}    
	 	        if(!sub_level.Objects[near_id].merged){//if Neighbor not merged
                    if (!sub_level.Objects[best_id].virtual_merged){//if neighbor merging heterogeneity not calculated
                        //Virtual merge into temp object
                        tmp.clear();
                        tmp.merge(sub_level.Objects[best_id], sub_level.Objects[near_id]);
                        //If not calculated, calculate neighbor object's perimeter
	 	        		tmp.perimeter = sub_level.CalculatePerimeter(best_id, near_id);
	 	        
						h_color = tmp.Color_Heterogeneity(Img) - sub_level.Objects[best_id].Color_Heterogeneity(Img) - sub_level.Objects[near_id].Color_Heterogeneity(Img);
                        h_smooth = tmp.Smoothness() - sub_level.Objects[best_id].Smoothness() - sub_level.Objects[near_id].Smoothness();
                        h_cmpct = tmp.Compactness() - sub_level.Objects[best_id].Compactness() - sub_level.Objects[near_id].Compactness();
                        h_shape = (passn.Compactness * h_cmpct) + ((1-passn.Compactness)*h_smooth);
                        f = (passn.Color*h_color) + ((1-passn.Color)*h_shape);
                    }
                    else{
                        //get the heterogeneity already calculated
                        f = sub_level.Objects[best_id].Neighbors[i].merging_heterogeneity;
                    }
                    if(f<best_best_f){
                    	best_best_f = f;
                    	best_best_id = near_id;
                   	}
                }    
	 	    }
	 	    //Is it mutual best match???
      		if(best_best_id == cur_id){
      		    mark3:
            	tmp.clear();
            	//merge objects
            	tmp.merge(sub_level.Objects[best_id], sub_level.Objects[cur_id]);
            	//Calculate perimeter
            	tmp.perimeter = sub_level.CalculatePerimeter(best_id, cur_id);
            	//Merge neighbors
            	tmp.merge_neighbors(sub_level.Objects[best_id], sub_level.Objects[cur_id]);
            	
				mergecounter++;
				mseg.last_id++;
             	tmp.id = mseg.last_id;
             	passn.Objects.insert(make_pair(tmp.id, tmp));
                tmpTopologyObject.id = tmp.id;
                tmpTopologyObject.merged = FALSE;
                tmpTopologyObject.new_id = 0;
                mseg.Topology.insert(make_pair(tmp.id, tmpTopologyObject));
                
                tmpTopologyPair.from = best_id;
                tmpTopologyPair.to = tmp.id;
                mseg.Temp_Topology.push_back(tmpTopologyPair);
                tmpTopologyPair.from = cur_id;
                tmpTopologyPair.to = tmp.id;
                mseg.Temp_Topology.push_back(tmpTopologyPair);
                
                sub_level.Objects[cur_id].merged = TRUE;
                sub_level.Objects[best_id].merged = TRUE;
                
                //Update raster
                passn.UpdateRaster(tmp.id);
                
                tmp.clear();
      		}
        	else{
             	//push back to priority list
             	if(mseg.Multi_Pushed.find(cur_id)!= mseg.Multi_Pushed.end()){
             	    mseg.Multi_Pushed[cur_id]++;
             	    if(mseg.Multi_Pushed[cur_id] >=5) goto mark3;
             	}else{
                  	mseg.Multi_Pushed.insert(make_pair(cur_id, 1));
                }        
              	tmp_priority_object.id = cur_id;
			    tmp_priority_object.primary = pri + 10;
			    tmp_priority_object.secondary = sec;
			    mseg.PriorityList.push(tmp_priority_object);
        	}        
        }
        else{//no matches
            //push the same object to the next level
            passn.Objects.insert(make_pair(cur_id, sub_level.Objects[cur_id]));
            sub_level.Objects[cur_id].merged = TRUE;
            passn.Objects[cur_id].priority_loaded = FALSE;
            
            //Update raster
            passn.UpdateRaster(cur_id);
        }    
        
    }
    //Update the Topology map
    int t;
    for(int i=0; i<mseg.Temp_Topology.size(); i++){
        t = mseg.Temp_Topology[i].from;
        mseg.Topology[t].merged = TRUE;
        mseg.Topology[t].new_id = mseg.Temp_Topology[i].to;
    }
    mseg.Temp_Topology.clear();    
    
    return mergecounter; 
}

/*void nthpass(Level& passn, ERS_Image& Img, float level_ID, MSEG_Init& mseg, MSEG_Param& Param, Level& sub_level, Level& super_level){
     
}*/

int nthpass(Level& passn, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass){
    
	int mergecounter = 0;
	passn.raster.reserve(Img.Lines*Img.Columns);
    passn.cycles = previous_pass.cycles++;
   	for(int j=0; j<(Img.Lines*Img.Columns); j++){
	    passn.raster.push_back(-1);
	}    
    //Copy the id number
    passn.HierarchyID = previous_pass.HierarchyID;
    passn.Lines = Img.Lines;
    passn.Columns = Img.Columns;
    
    //Copy segmentation parameters into new level
    passn.ScaleParameter = previous_pass.ScaleParameter;
    passn.Color = previous_pass.Color;
    passn.Compactness = previous_pass.Compactness;
    passn.Edge = previous_pass.Edge;
    //cout << "Segmentation parameters: " << passn.ScaleParameter << " " << passn.Color
    //<< " " << passn.Compactness << endl;//Debug mode
    
    //Declare a temporary object
    Object tmp;
    pixel cur;
    int pri, sec, matches;
    int cur_id, near_id, best_id, best_best_id;
    float best_f, best_best_f, h_color, h_shape, h_smooth, h_cmpct, f;
    PriorityObject tmp_priority_object;
    TopologyObject tmpTopologyObject;
    TopologyPair tmpTopologyPair;
    
    //Load the Starting Points into Priority List
    pri = 1;
    sec = 0;
    for(int i=0; i<mseg.StartingPoints.size(); i++){
  		cur = mseg.StartingPoints[i];
  		sec = i+1;
  		cur_id = previous_pass.raster[D1(Img.Lines, Img.Columns, cur.line, cur.column)];
  		tmp_priority_object.id = cur_id;
  		tmp_priority_object.primary = pri;
  		tmp_priority_object.secondary = sec;
  		mseg.PriorityList.push(tmp_priority_object);
  		previous_pass.Objects[cur_id].priority_loaded = TRUE;//Flag that it is loaded into priority list
    }
    //cout << "SP loaded into Priority list..." << endl << "Priority list holds: "
    //<< mseg.PriorityList.size() << " objects" << endl;//Debug mode
    
    //While the priority list is not empty
    while (!mseg.PriorityList.empty()){
	    
		//Load the next priority object
     	tmp_priority_object = mseg.PriorityList.top();
      	pri = tmp_priority_object.primary;
       	sec = tmp_priority_object.secondary;
        cur_id = tmp_priority_object.id;
        mseg.PriorityList.pop();
        //cout << "Priority list holds " << mseg.PriorityList.size() << " objects after pop"
        //<< " and merges occured = " << mergecounter << endl;//Debug mode
       	
		//Reload Priority object if cur already merged
        if(previous_pass.Objects[cur_id].merged == TRUE) continue;
        
		//If not calculated, calculate object's perimeter
        if(previous_pass.Objects[cur_id].perimeter == 0){
  		    previous_pass.CalculatePerimeter(cur_id);
		}
		//Update and Fix the neighbors for the previous level in order to call once bellow
		update_neighbors(previous_pass.Objects[cur_id], mseg);
		previous_pass.Objects[cur_id].fix_neighbors();
		
		best_f = passn.ScaleParameter; //The maximum heterogeneity allowed is the scale parameter
		matches = 0;
		
        //Load neigbors --> DONE: dahh, we have objects now...
		
        //For every neighbor
		for(int i=0; i<previous_pass.Objects[cur_id].Neighbors.size(); i++){
            near_id = previous_pass.Objects[cur_id].Neighbors[i].id;
            //cout << "Processing neighbor " << i+1 << endl;//Debug mode
            //cout << "Neighbor id is " << near_id << endl;
            //If not calculated, calculate neighbor object's perimeter
            if(previous_pass.Objects[near_id].perimeter == 0){
         	    previous_pass.CalculatePerimeter(near_id);
		    }
		    update_neighbors(previous_pass.Objects[near_id], mseg);
		    previous_pass.Objects[near_id].fix_neighbors();
		
			//if neighbor not pushed to priority --> Push
			if(!previous_pass.Objects[near_id].priority_loaded){
			    tmp_priority_object.id = near_id;
			    tmp_priority_object.primary = pri + 1;
			    tmp_priority_object.secondary = sec;
			    mseg.PriorityList.push(tmp_priority_object);
			    previous_pass.Objects[near_id].priority_loaded = TRUE;
			    //cout << "Neighbor pushed into priority list" << endl;//Debug mode
			}    
			if(!previous_pass.Objects[near_id].merged){//if Neighbor not merged
				//cout << "Neighbor not merged, inside if..." << endl;//Debug mode 
     			if (!previous_pass.Objects[cur_id].virtual_merged){//if neighbor merging heterogeneity not calculated
	 	             //Virtual merge into temp object
                     tmp.clear();
                     tmp.merge(previous_pass.Objects[cur_id], previous_pass.Objects[near_id]);
                     //Calculate merge object's perimeter
                     tmp.perimeter = previous_pass.CalculatePerimeter(cur_id, near_id);
                     //cout << "Perimeter calculated " << tmp.perimeter << endl;
                     //Calculate heterogeneity
                     h_color = tmp.Color_Heterogeneity(Img) - previous_pass.Objects[cur_id].Color_Heterogeneity(Img) - previous_pass.Objects[near_id].Color_Heterogeneity(Img);
                     //cout << "Color heterogeneity " << h_color << endl;
                     h_smooth = tmp.Smoothness() - previous_pass.Objects[cur_id].Smoothness() - previous_pass.Objects[near_id].Smoothness();
                     //cout << "Smoothness " << h_smooth << endl;
                     h_cmpct = tmp.Compactness() - previous_pass.Objects[cur_id].Compactness() - previous_pass.Objects[near_id].Compactness();
                     //cout << "Compactness " << h_cmpct << endl;
                     h_shape = (passn.Compactness * h_cmpct) + ((1-passn.Compactness)*h_smooth);
                     f = (passn.Color*h_color) + ((1-passn.Color)*h_shape);
                     //cout << "Heterogeneity for neighbor " << i+1 << " calculated: " << f << endl;//Debug mode
                     //exit(1);//Debug mode
                 }
                 else{
                     //get the heterogeneity already calculated
                    //cout << "This message should not be printed..." << endl;//Debug mode
                     f = previous_pass.Objects[cur_id].Neighbors[i].merging_heterogeneity;
                 }        
   			     if(f<best_f){
      		     	matches++;
      		     	//cout << "It's a match..." << endl;//Debug mode
  			   		best_f = f;
  			   		best_id = near_id;
   		         }
   		     }    
	 	}
	 	//cout << " Total matches: " << matches << endl;//Debug mode
        //Check for match
	 	if(matches){
	 	    //Load best match
	 	    best_best_f = passn.ScaleParameter;
	 	    //For every neighbor
	 	    for(int i=0; i<previous_pass.Objects[best_id].Neighbors.size(); i++){
	 	        near_id = previous_pass.Objects[best_id].Neighbors[i].id;
	 	        //cout << "Processing best neighbor " << i+1 << endl;//Debug mode
	 	        //If not calculated, calculate neighbor object's perimeter
            	if(previous_pass.Objects[near_id].perimeter == 0){
			         previous_pass.CalculatePerimeter(near_id);
 		 		}
 		 		update_neighbors(previous_pass.Objects[near_id], mseg);
		    	previous_pass.Objects[near_id].fix_neighbors();
		    	
 		 		//cout << "Best neighbor perimeter: " << previous_pass.Objects[near_id].perimeter << endl;
 		 		//cout << " Best neighbor area: " << previous_pass.Objects[near_id].area << endl;
 		 		//cout << "near_id " << near_id << endl;
		        //if neighbor not pushed to priority --> Push
		        if(!previous_pass.Objects[near_id].priority_loaded){
		        	tmp_priority_object.id = near_id;
		        	tmp_priority_object.primary = pri + 1;
		        	tmp_priority_object.secondary = sec;
		        	mseg.PriorityList.push(tmp_priority_object);
		        	previous_pass.Objects[near_id].priority_loaded = TRUE;
		        	//cout << "Neighbor pushed into priority list" << endl;//Debug mode
	        	}    
		        if(!previous_pass.Objects[near_id].merged){//if Neighbor not merged
		        	//cout << "Neighbor not merged, inside if..." << endl;//Debug mode
                    if (!previous_pass.Objects[best_id].virtual_merged){//if neighbor merging heterogeneity not calculated
                        //Virtual merge into temp object
                        tmp.clear();
                        tmp.merge(previous_pass.Objects[best_id], previous_pass.Objects[near_id]);
                        //Calculate merge object's perimeter
                        tmp.perimeter = previous_pass.CalculatePerimeter(best_id, near_id);
                        //cout << "Perimeter calculated " << tmp.perimeter << endl;
						h_color = tmp.Color_Heterogeneity(Img) - previous_pass.Objects[best_id].Color_Heterogeneity(Img) - previous_pass.Objects[near_id].Color_Heterogeneity(Img);
                        //cout << "Color heterogeneity " << h_color << endl;
                        h_smooth = tmp.Smoothness() - previous_pass.Objects[best_id].Smoothness() - previous_pass.Objects[near_id].Smoothness();
                        //cout << "Smoothness " << h_smooth << endl;
                        h_cmpct = tmp.Compactness() - previous_pass.Objects[best_id].Compactness() - previous_pass.Objects[near_id].Compactness();
                        //cout << "Compactness " << h_cmpct << endl;
                        h_shape = (passn.Compactness * h_cmpct) + ((1-passn.Compactness)*h_smooth);
                        f = (passn.Color*h_color) + ((1-passn.Color)*h_shape);
                        //cout << "Heterogeneity for neighbor " << i+1 << " calculated: " << f << endl;//Debug mode
                    }
                    else{
                        //get the heterogeneity already calculated
                        //cout << "This message should not be printed..." << endl;//Debug mode
                        f = previous_pass.Objects[best_id].Neighbors[i].merging_heterogeneity;
                    }
                    if(f<best_best_f){
                    	//cout << "It's a match..." << endl;
                     	best_best_f = f;
                    	best_best_id = near_id;
                   	}
                }    
	 	    }
	 	    //Is it mutual best match???
      		if(best_best_id == cur_id){
      		    mark2:
            	//cout << "Inside mutual best match..." << endl;//Debug mode
            	tmp.clear();
            	//merge objects
            	tmp.merge(previous_pass.Objects[best_id], previous_pass.Objects[cur_id]);
            	//Calculate perimeter
            	tmp.perimeter = previous_pass.CalculatePerimeter(best_id, cur_id);
            	//Merge neighbors
            	tmp.merge_neighbors(previous_pass.Objects[best_id], previous_pass.Objects[cur_id]);
            	//cout << "A merge has just occured..." << endl;//Debug mode
             	mergecounter++;
				mseg.last_id++;
             	tmp.id = mseg.last_id;
             	passn.Objects.insert(make_pair(tmp.id, tmp));//insert here
             	
              	tmpTopologyObject.id = tmp.id;
                tmpTopologyObject.merged = FALSE;
                tmpTopologyObject.new_id = 0;
                mseg.Topology.insert(make_pair(tmp.id, tmpTopologyObject));
                tmpTopologyPair.from = best_id;
                tmpTopologyPair.to = tmp.id;
                mseg.Temp_Topology.push_back(tmpTopologyPair);
                tmpTopologyPair.from = cur_id;
                tmpTopologyPair.to = tmp.id;
                mseg.Temp_Topology.push_back(tmpTopologyPair);
                
                previous_pass.Objects[cur_id].merged = TRUE;
                previous_pass.Objects[best_id].merged = TRUE;
                
                //Update raster
                passn.UpdateRaster(tmp.id);
                
                tmp.clear();
                //cout << "Now level holds " << passn.Objects.size() << " objects" << endl;//Debug mode
      		}
        	else{
             	//push back to priority list
             	//cout << "not mutual best match, pushing back to priority list..." << endl;//Debug mode 
             	if(mseg.Multi_Pushed.find(cur_id)!= mseg.Multi_Pushed.end()){
             	    mseg.Multi_Pushed[cur_id]++;
             	    if(mseg.Multi_Pushed[cur_id] >=5) goto mark2;
             	}else{
                  	mseg.Multi_Pushed.insert(make_pair(cur_id, 1));
                }        
              	tmp_priority_object.id = cur_id;
			    tmp_priority_object.primary = pri + 20;
			    tmp_priority_object.secondary = sec;
			    mseg.PriorityList.push(tmp_priority_object);
        	}        
        }
        else{//no matches
            //push the same object to the next level
            passn.Objects.insert(make_pair(cur_id, previous_pass.Objects[cur_id]));//insert here to
            previous_pass.Objects[cur_id].merged = TRUE;
            //cout << "No matches, object copied to the next level" << endl;//Debug mode
            passn.Objects[cur_id].priority_loaded = FALSE;
            
            //Update raster
            passn.UpdateRaster(cur_id);
        }    
        
    }
    //Update the Topology map
    int t;
    for(int i=0; i<mseg.Temp_Topology.size(); i++){
        t = mseg.Temp_Topology[i].from;
        mseg.Topology[t].merged = TRUE;
        mseg.Topology[t].new_id = mseg.Temp_Topology[i].to;
    }
    mseg.Temp_Topology.clear();    
    
    return mergecounter; 
}

/*void nthpass(Level& passn, ERS_Image& Img, MSEG_Init& mseg, Level& previous_pass, Level& super_level){
     
}*/


void SPE (ERS_Image& Img, MSEG_Init& mseg, int SPE_mode){

    if(SPE_mode == 0){
	    cout << "HSI method selected for start point estimation" << endl;
	    vector<float> Intensity_band;
      	Intensity_band.reserve(Img.Lines * Img.Columns);
                
                if(Img.Bands == 1){
                    Img.FullBuffer();
                    for(int i=0; i<(Img.Lines*Img.Columns); i++){
                        Intensity_band[i] = float(Img.Buffer[i]);
                    }
                }    
                else if(Img.Bands == 2){
                    Img.FullBuffer();
                    for(int l=1; l<=Img.Lines; l++){
                        for(int c=1; c<=Img.Columns; c++){
                            Intensity_band[((l-1)*Img.Columns + (c-1))] = float((Img.Buffer[((l-1)*Img.Columns*Img.Bands)+c-1] + Img.Buffer[((l-1)*Img.Columns*Img.Bands)+c+Img.Columns-1])/2);
                        }    
                    }
                }
                else if(Img.Bands > 2){
                    Img.FullBuffer();
                    for(int l=1; l<=Img.Lines; l++){
                        for(int c=1; c<=Img.Columns; c++){
                            Intensity_band[((l-1)*Img.Columns + (c-1))] = float((Img.Buffer[((l-1)*Img.Columns*Img.Bands)+c-1] + Img.Buffer[((l-1)*Img.Columns*Img.Bands)+c+Img.Columns-1] + Img.Buffer[((l-1)*Img.Columns*Img.Bands)+c+Img.Columns+Img.Columns-1])/3);
                        }
                    }
                }    
                
                //Main SPE routine
                int max_line, max_col, min_line, min_col, pixelcount;
                double ss, s, std, mean;
                float max, min, val;
                pixel tmp;
                
                printf("\n Macroblock|Act. Pixels| Min | Max | Aver.|Std.Dev.| SP (x,y)");
                printf("\n-----------+-----------+-----+-----+------+--------|---------");
                
                for(int i=0;i<(Img.MacroLines * Img.MacroColumns);i++){ //For each macroblock
                    
                    ss = 0.0;
                    s = 0.0;
                    std = 0.0;
                    mean = 0.0;
                    pixelcount = 0;
                    
                    //Setting macroblock first value as min and max
                    max = Intensity_band[(Img.Macroblocks[i].LineStart - 1)*Img.Columns + Img.Macroblocks[i].ColumnStart - 1];
                    min = Intensity_band[(Img.Macroblocks[i].LineStart - 1)*Img.Columns + Img.Macroblocks[i].ColumnStart - 1];
                    max_line = Img.Macroblocks[i].LineStart;
                    min_line = Img.Macroblocks[i].LineStart;
                    max_col = Img.Macroblocks[i].ColumnStart;
                    min_col = Img.Macroblocks[i].ColumnStart;
                    //
                    
                    for(int l=Img.Macroblocks[i].LineStart; l<=Img.Macroblocks[i].LineEnd; l++){ //Scaning macroblock
                        for(int c=Img.Macroblocks[i].ColumnStart; c<=Img.Macroblocks[i].ColumnEnd; c++){
                            val = Intensity_band[(l-1)*Img.Columns + c -1];
                            if (val >= 0){//TODO: Import NULL value to Image Library and use it here
                                pixelcount++;
                                s += val;
                                ss += (val*val);
                                if (val > max){
                                    max = val;
                                    max_line = l;
                                    max_col = c;
                                }
                                if (val < min){
                                    min = val;
                                    min_line = l;
                                    min_col = c;
                                }
                            }                                                            
                        }
                    }
                    mean = s/pixelcount;
                    std = sqrt(ss/pixelcount-s*s/pixelcount/pixelcount);
                    
                    if((max-mean)>(mean-min)){
                        tmp.line = max_line;
                        tmp.column = max_col;
                    }
                    else {
                        tmp.line = min_line;
                        tmp.column = min_col;
                    }
                    
                    mseg.StartingPoints.push_back(tmp);
                    //mseg.StartingPoints[i] = tmp;
                    
                    printf("\n %4d | %4d | %3f | %3f | %4.3f | %4.3f | %5d |  %5d",
                                i+1, pixelcount, min, max, mean, std, tmp.line, tmp.column );
                    
                }
                delete &Intensity_band;
        }
        else if(SPE_mode == 1){
        	 cout << "PCA method selected for start point estimation" << endl;
		}
		else if(SPE_mode == 2){
			 cout << "Dithering matrix method selected for start point estimation" << endl;
                Img.FullBuffer();
                FIBITMAP *image = FreeImage_Allocate(Img.Columns, Img.Lines, 8);
                
                BYTE *bits[Img.Lines];
                for(int y = 0; y<Img.Lines; y++){
   					bits[y] = (BYTE*)FreeImage_GetScanLine(image, y);
				}
                
                for(int y = 1; y<=Img.Lines; y++){
				    for(int x = 1; x<=Img.Columns; x++){
   					    bits[y-1][x-1] = Img.Buffer[Img.D1(y, x, 1)];
   					}
				}
                
                RGBQUAD *pal = FreeImage_GetPalette(image);
  				for(int i=0;i<256;i++){
 					pal[i].rgbRed = i;
      				pal[i].rgbGreen = i;
      				pal[i].rgbBlue = i;
				}
				
				FIBITMAP *dithered;
				
				dithered = FreeImage_Dither(image, FID_FS);
  				// algorithms available: "FID_FS" Floyd & Steinberg error diffusion algorithm
  				//                       "FID_BAYER4x4" Bayer ordered dispersed dot dithering (order 2 – 4x4 -dithering matrix)
  				//                       "FID_BAYER8x8" Bayer ordered dispersed dot dithering (order 3 – 8x8 -dithering matrix)
  				//                       "FID_CLUSTER6x6" Ordered clustered dot dithering (order 3 - 6x6 matrix)
  				//                       "FID_CLUSTER8x8" Ordered clustered dot dithering (order 4 - 8x8 matrix)
  				//                       "FID_CLUSTER16x16" Ordered clustered dot dithering (order 8 - 16x16 matrix)
  				// References :
  				//        Ulichney, R., Digital Halftoning. The MIT Press, Cambridge, MA, 1987.
  				//        Hawley S., Ordered Dithering. Graphics Gems, Academic Press, 1990. 
  				
  				FreeImage_Unload(image);
  				delete &bits;
  				delete &pal;
  				
  				
  				BYTE *bits_dith[Img.Lines];
				for(int y = 0; y<Img.Lines; y++){
   					bits_dith[y] = (BYTE*)FreeImage_GetScanLine(dithered, y);
				}
				  
  		        pixel pix;
  				
  				for(int l=1; l<=Img.Lines; l++){
				    for(int c=1; c<=Img.Columns; c++){
						if(bits_dith[l-1][c-1] == TRUE){
			   			    pix.line = l;
			   			    pix.column = c;
			   			    mseg.StartingPoints.push_back(pix);
						}
				    }
				}
				
  				FreeImage_Unload(dithered);
		}
  		else{
		 	cout << "Invalid SPE mode" << endl;
		}       
}


