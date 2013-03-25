/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * msegio.cpp: Image base classes	                                  *
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

#include "msegio.h"

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

/*void ERS_Image::FillBuffer (int LineStart, int LineEnd){

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

}*/

/*void ERS_Image::FillBuffer (Macroblock* macro){

}*/

void ERS_Image::FullBuffer (){

    if(BuffLineStart != BuffLineEnd){
        cout << "Buffer is not empty" << endl;
        cout << "Deleting buffer..." << endl;
        Buffer8.clear();
        Buffer16.clear();
        Buffer16s.clear();
        Buffer32.clear();
    }

    uchar tmp8;
    ushort tmp16;
    sshort tmp16s;
    float tmp32;

    file.seekg(0, ios::beg);

    if(Bytes == 1)
    {
        for(int i=0; i<(Lines*Columns*Bands); i++){
            file.read(reinterpret_cast<char *>(&tmp8), sizeof(uchar));
            Buffer8.push_back(tmp8);
        }
    }
    else if(Bytes == 2)
    {
        for(int i=0; i<(Lines*Columns*Bands); i++){
            file.read(reinterpret_cast<char *>(&tmp16), sizeof(ushort));
            Buffer16.push_back(tmp16);
        }
    }
    else if(Bytes == 3)
    {
        for(int i=0; i<(Lines*Columns*Bands); i++){
            file.read(reinterpret_cast<char *>(&tmp16s), sizeof(sshort));
            Buffer16s.push_back(tmp16s);
        }
    }
    else if(Bytes == 4)
    {
        for(int i=0; i<(Lines*Columns*Bands); i++){
            file.read(reinterpret_cast<char *>(&tmp32), sizeof(float));
            Buffer32.push_back(tmp32);
        }
    }

	BuffLineStart = 1;
    BuffLineEnd = Lines;
    file.seekg(0, ios::beg);
    cout << "Buffer Full!" << endl;

}

float ERS_Image::Buffer(int line, int column, int band){

   if(Bytes == 1)
   {
       return (float)(Buffer8[(((line-1)*Bands*Columns) + ((band-1)*Columns) + column -1)]);
   }
   else if(Bytes == 2)
   {
       return (float)(Buffer16[(((line-1)*Bands*Columns) + ((band-1)*Columns) + column -1)]);
   }
   else if(Bytes == 3)
   {
       return (float)(Buffer16s[(((line-1)*Bands*Columns) + ((band-1)*Columns) + column -1)]);
   }
   else if(Bytes == 4)
   {
       return Buffer32[(((line-1)*Bands*Columns) + ((band-1)*Columns) + column -1)];
   }
   else
   {
       cout << "Unknown image type" << endl;
       return 0.0;
   }
}

float ERS_Image::Buffer(int indx){
    if(Bytes == 1)
   {
       return (float)(Buffer8[indx]);
   }
   else if(Bytes == 2)
   {
       return (float)(Buffer16[indx]);
   }
   else if(Bytes == 3)
   {
       return (float)(Buffer16s[indx]);
   }
   else if(Bytes == 4)
   {
       return Buffer32[indx];
   }
   else
   {
       cout << "Unknown image type" << endl;
       return 0.0;
   }
}

int ERS_Image::D1(int line, int column, int band){
    return (((line-1)*Bands*Columns) + ((band-1)*Columns) + column -1);
}

void ERS_Image::SaveAs(string OutputName){
    ofstream head_out, file_out;
    string tmp;
    tmp = OutputName + ".ers";//TODO: Change that to accept the header of the file

    head_out.open(tmp.c_str());
    file_out.open(OutputName.c_str(), ios::binary);

    //save binary data
	if(Bytes == 1)
	{
        for(int i=0; i<Buffer8.size(); i++){
            file_out.write(reinterpret_cast<char *> (&Buffer8[i]), sizeof(uchar));
        }
	}
	else if(Bytes == 2)
	{
        for(int i=0; i<Buffer16.size(); i++){
            //t = Buffer[i];
            file_out.write(reinterpret_cast<char *> (&Buffer16[i]), sizeof(ushort));
        }
	}
	else if(Bytes == 3)
	{
        for(int i=0; i<Buffer16s.size(); i++){
            //t = Buffer[i];
            file_out.write(reinterpret_cast<char *> (&Buffer16s[i]), sizeof(sshort));
        }
	}
	else if(Bytes == 4)
	{
        for(int i=0; i<Buffer32.size(); i++){
            file_out.write(reinterpret_cast<char *> (&Buffer32[i]), sizeof(float));
        }
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
    if(RegistrationPoint.x != -999 || RegistrationPoint.y != -999)
    {
	head_out << "\t\tRegistrationCoord Begin" << endl;
	head_out << "\t\t\tEastings\t= " << RegistrationPoint.x << endl;
	head_out << "\t\t\tNorthings\t= " << RegistrationPoint.y << endl;
	head_out << "\t\tRegistrationCoord End" << endl;
    }
    
    if(RegistrationCell.x != -999 || RegistrationCell.y != -999)
    {
	head_out << "\t\tRegistrationCellX\t= " << RegistrationCell.x << endl;
	head_out << "\t\tRegistrationCellY\t= " << RegistrationCell.y << endl;
    }
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
    if(Bytes == 1)
    {
        head_out << "\t\t\t\tRed\t\t= 255" << endl;
        head_out << "\t\t\t\tGreen\t\t= 255" << endl;
        head_out << "\t\t\t\tBlue\t\t= 255" << endl;
    }
    else
    {
        head_out << "\t\t\t\tRed\t\t= 65535" << endl;
        head_out << "\t\t\t\tGreen\t\t= 65535" << endl;
        head_out << "\t\t\t\tBlue\t\t= 65535" << endl;
    }
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
    if(Bytes == 1)
	{
        for(int i=0; i<Buffer8.size(); i++){
            file_out.write(reinterpret_cast<char *> (&Buffer8[i]), sizeof(uchar));
        }
	}
	else if(Bytes == 2)
	{
        for(int i=0; i<Buffer16.size(); i++){
            //t = Buffer[i];
            file_out.write(reinterpret_cast<char *> (&Buffer16[i]), sizeof(ushort));
        }
	}
	else if(Bytes == 3)
	{
        for(int i=0; i<Buffer16s.size(); i++){
            //t = Buffer[i];
            file_out.write(reinterpret_cast<char *> (&Buffer16s[i]), sizeof(sshort));
        }
	}
	else if(Bytes == 4)
	{
        for(int i=0; i<Buffer32.size(); i++){
            file_out.write(reinterpret_cast<char *> (&Buffer32[i]), sizeof(float));
        }
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
    if(RegistrationPoint.x != -999 || RegistrationPoint.y != -999)
    {
	head_out << "\t\tRegistrationCoord Begin" << endl;
	head_out << "\t\t\tEastings\t= " << RegistrationPoint.x << endl;
	head_out << "\t\t\tNorthings\t= " << RegistrationPoint.y << endl;
	head_out << "\t\tRegistrationCoord End" << endl;
    }
    
    if(RegistrationCell.x != -999 || RegistrationCell.y != -999)
    {
	head_out << "\t\tRegistrationCellX\t= " << RegistrationCell.x << endl;
	head_out << "\t\tRegistrationCellY\t= " << RegistrationCell.y << endl;
    }
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
    if(Bytes == 1)
    {
        head_out << "\t\t\t\tRed\t\t= 255" << endl;
        head_out << "\t\t\t\tGreen\t\t= 255" << endl;
        head_out << "\t\t\t\tBlue\t\t= 255" << endl;
    }
    else
    {
        head_out << "\t\t\t\tRed\t\t= 65535" << endl;
        head_out << "\t\t\t\tGreen\t\t= 65535" << endl;
        head_out << "\t\t\t\tBlue\t\t= 65535" << endl;
    }
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

   	string word;				//temporary string for input
   	
   	GroundSizeX = 1.0;
	GroundSizeY = 1.0;
	RegistrationPoint.x = -999;
	RegistrationPoint.y = -999;
	RegistrationCell.x = -999;
	RegistrationCell.y = -999;

 	//Header file parser

 	while (header){
 	    header >> word;

	if (strncmp(word.c_str(), "Version", 6) == 0){
 	        header >> word;		/*Skips the "="*/
 	        header >> word;
 	        ERSVersion = word.substr(1,word.size()-2);		/*Removes the "___"*/
 	        //string ERSVersion(word, 2, length-1);				/*Alternative*/
	}

       if (strncmp(word.c_str(), "LastUpdated", 11) == 0){
           header >> word;		/*Skips the "="*/
           for (int i=0; i<5; i++){
               header >> word;
               LastUpdated = LastUpdated + word + " ";/*Reads 5 words*/
           }
           header >> word;
           LastUpdated = LastUpdated + word; /*Does not add the space after the sixth word*/
       }

       if (strncmp(word.c_str(), "DataType", 8) == 0){
           header >> word;		/*Skips the "="*/
           header >> DataType;
       }

       if (strncmp(word.c_str(), "ByteOrder", 9) == 0){
           header >> word;		/*Skips the "="*/
           header >> ByteOrder;
       }

       if (strncmp(word.c_str(), "Datum", 5) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           CoordinateSpace.Datum = word.substr(1,word.size()-2);		/*Removes the "___"*/
       }

       if (strncmp(word.c_str(), "Projection", 10) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           CoordinateSpace.Projection = word.substr(1,word.size()-2);	/*Removes the "___"*/
       }

       if (strncmp(word.c_str(), "CoordinateType", 14) == 0){
           header >> word;		/*Skips the "="*/
           header >> CoordinateSpace.CoordinateType;
       }

       if (strncmp(word.c_str(), "Units", 5) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           CoordinateSpace.Units = word.substr(1,word.size()-2);	/*Removes the "___"*/
       }

       if (strncmp(word.c_str(), "Xdimension", 10) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           GroundSizeX = atof(word.c_str());
       }

       if (strncmp(word.c_str(), "Ydimension", 10) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           GroundSizeY = atof(word.c_str());
       }

       if (strncmp(word.c_str(), "Eastings", 8) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           RegistrationPoint.x = atof(word.c_str());
       }

       if (strncmp(word.c_str(), "Northings", 9) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           RegistrationPoint.y = atof(word.c_str());
       }
       
       if (strncmp(word.c_str(), "RegistrationCellX", 17) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           RegistrationCell.x = atof(word.c_str());
       }
       
       if (strncmp(word.c_str(), "RegistrationCellY", 17) == 0){
           header >> word;		/*Skips the "="*/
           header >> word;
           RegistrationCell.y = atof(word.c_str());
       }

       if (strncmp(word.c_str(), "Rotation", 8) == 0){
           header >> word;		/*Skips the "="*/
           header >> CoordinateSpace.Rotation;
       }

       if (strncmp(word.c_str(), "CellType", 8) == 0){
           header >> word;		/*Skips the "="*/
           header >> CellType;
           if (strncmp(CellType.c_str(), "Unsigned16BitInteger", 20) == 0)
           {
               Bytes = 2;
           }
           else if (strncmp(CellType.c_str(), "Signed16BitInteger", 18) == 0)
           {
               Bytes = 3;
           }
           else if (strncmp(CellType.c_str(), "Unsigned8BitInteger", 19) == 0)
           {
               Bytes = 1;
           }
           else if (strncmp(CellType.c_str(), "IEEE4ByteReal", 13) == 0)
           {
               Bytes = 4;
           }
           else  Bytes = 1;
       }

       if (strncmp(word.c_str(), "NullCellValue", 13) == 0){
           header >> word;		/*Skips the "="*/
           header >> NullValue;
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
           //vector<string> BandID;
           for (int i=0; i<5; i++){
               header >> word;
           }
           BandID.push_back(word);
           cout << word << endl;
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

   //Buffer Initialization
   BuffLineStart = 0;
   BuffLineEnd = 0;   // if BuffLineStart == BuffLineEnd ---> Buffer empty
   //system("Pause");

   /*Default Band Weight values can be modified through the main routine from
   an outside *.txt file */
   BandWeight.push_back(0.0);
   for (int i=1; i<=Bands; i++){
       BandWeight.push_back(1.0);
   }

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
   //Macroblock calculation End
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
    
    GroundSizeX = 1.0;
    GroundSizeY = 1.0;
    RegistrationPoint.x = -999;
    RegistrationPoint.y = -999;
    RegistrationCell.x = -999;
    RegistrationCell.y = -999;

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

ERS_Image::ERS_Image(string Name, int BANDS, int LINES, int COLUMNS, int BYTES){

    Bands = BANDS;
    Lines = LINES;
    Columns = COLUMNS;
    Bytes = BYTES;

    FileName = Name;
    if(Bytes == 1) CellType = "Unsigned8BitInteger";
    else if(Bytes == 2) CellType = "Unsigned16BitInteger";
    else if(Bytes == 4) CellType = "IEEE4ByteReal";
    else if(Bytes == 3) CellType = "Signed16BitInteger";

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

    GroundSizeX = 1.0;
    GroundSizeY = 1.0;
    RegistrationPoint.x = -999;
    RegistrationPoint.y = -999;
    RegistrationCell.x = -999;
    RegistrationCell.y = -999;
    
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

Texture_Image::Texture_Image(){

}

Texture_Image::Texture_Image(ERS_Image& Im, int quantizer, int dist){
    bitdepth=256;
    if(quantizer>bitdepth){
        cout << "Texture quantizer larger than image bitdepth... defaulting to bitdepth" << endl;
        quant = bitdepth;
    }else{
        quant = quantizer;
    }
    distance = dist;
    Img = &Im;

    Lines = Img->Lines;
    Columns = Img->Columns;
    CoCubeBuilt = FALSE;

    CoPair tmpcopair;
    tmpcopair.first = -1;
    tmpcopair.second = -1;
    CoDistribution lolos;
    for(int i=0; i<Im.Lines*Im.Columns; i++){
        lolos.push_back(tmpcopair);
    }
    for(int i=0; i<4; i++){
       CoCube.push_back(lolos);
    }
    cout << "Memory allocation for CoCube completed..." << endl;
}

Texture_Image::~Texture_Image(){

}

void Texture_Image::FullYBuffer(){
    YBuffer.reserve(Lines * Columns);
    for(int i=0; i<Lines*Columns; i++){
        YBuffer.push_back(0);
    }
    cout << "YBuffer size = " << YBuffer.size() << endl;
    if(Img->Bands == 1){
    	cout << "Not enough bands to calculate YUV components... Using HSI by default" << endl;
		Img->FullBuffer();
     	for(int i=0; i<(Img->Lines*Img->Columns); i++){
      		if(Img->Bytes == 1)	YBuffer[i] = (int(Img->Buffer8[i]));
      		else if(Img->Bytes == 2)	YBuffer[i] = (int(Img->Buffer16[i]));
      		else if(Img->Bytes == 3)	YBuffer[i] = (int(Img->Buffer16s[i]));
      		else if(Img->Bytes == 4)	YBuffer[i] = (int(Img->Buffer32[i]));
        }
    }
    else if(Img->Bands == 2){
    	cout << "Not enough bands to calculate YUV components... Using HSI by default" << endl;
		Img->FullBuffer();
     	for(int l=1; l<=Img->Lines; l++){
      		for(int c=1; c<=Img->Columns; c++){
        		if(Img->Bytes == 1) YBuffer[((l-1)*Img->Columns + (c-1))] = (int((Img->Buffer8[((l-1)*Img->Columns*Img->Bands)+c-1] + Img->Buffer8[((l-1)*Img->Columns*Img->Bands)+c+Img->Columns-1])/2));
        		else if(Img->Bytes == 2) YBuffer[((l-1)*Img->Columns + (c-1))] = (int((Img->Buffer16[((l-1)*Img->Columns*Img->Bands)+c-1] + Img->Buffer16[((l-1)*Img->Columns*Img->Bands)+c+Img->Columns-1])/2));
        		else if(Img->Bytes == 3) YBuffer[((l-1)*Img->Columns + (c-1))] = (int((Img->Buffer16s[((l-1)*Img->Columns*Img->Bands)+c-1] + Img->Buffer16s[((l-1)*Img->Columns*Img->Bands)+c+Img->Columns-1])/2));
        		else if(Img->Bytes == 4) YBuffer[((l-1)*Img->Columns + (c-1))] = (int((Img->Buffer32[((l-1)*Img->Columns*Img->Bands)+c-1] + Img->Buffer32[((l-1)*Img->Columns*Img->Bands)+c+Img->Columns-1])/2));
          	}
        }
    }
    else if(Img->Bands > 2){
    	Img->FullBuffer();
     	for(int l=1; l<=Img->Lines; l++){
      		for(int c=1; c<=Img->Columns; c++){
        		if(Img->Bytes == 1) YBuffer[((l-1)*Img->Columns + (c-1))] = (int((0.114*Img->Buffer8[((l-1)*Img->Columns*Img->Bands)+c-1]) + (0.587*Img->Buffer8[((l-1)*Img->Columns*Img->Bands)+c+Img->Columns-1]) + (0.299*Img->Buffer8[((l-1)*Img->Columns*Img->Bands)+c+Img->Columns+Img->Columns-1])));
        		else if(Img->Bytes == 2) YBuffer[((l-1)*Img->Columns + (c-1))] = (int((0.114*Img->Buffer16[((l-1)*Img->Columns*Img->Bands)+c-1]) + (0.587*Img->Buffer16[((l-1)*Img->Columns*Img->Bands)+c+Img->Columns-1]) + (0.299*Img->Buffer16[((l-1)*Img->Columns*Img->Bands)+c+Img->Columns+Img->Columns-1])));
        		else if(Img->Bytes == 3) YBuffer[((l-1)*Img->Columns + (c-1))] = (int((0.114*Img->Buffer16s[((l-1)*Img->Columns*Img->Bands)+c-1]) + (0.587*Img->Buffer16s[((l-1)*Img->Columns*Img->Bands)+c+Img->Columns-1]) + (0.299*Img->Buffer16s[((l-1)*Img->Columns*Img->Bands)+c+Img->Columns+Img->Columns-1])));
        		else if(Img->Bytes == 4) YBuffer[((l-1)*Img->Columns + (c-1))] = (int((0.114*Img->Buffer32[((l-1)*Img->Columns*Img->Bands)+c-1]) + (0.587*Img->Buffer32[((l-1)*Img->Columns*Img->Bands)+c+Img->Columns-1]) + (0.299*Img->Buffer32[((l-1)*Img->Columns*Img->Bands)+c+Img->Columns+Img->Columns-1])));
          	}
        }
    }

    //Find Limits
    BuffMin = 65535;//Larger than 32bit integer...
    BuffMax = 0;

    for (int i=0;i<Lines*Columns;i++){
        if(YBuffer[i]>BuffMax) BuffMax=YBuffer[i];
        if(YBuffer[i]<BuffMin) BuffMin=YBuffer[i];
    }
}

int Texture_Image::Transform(int p){
    //Correct method
    /*float val;
    val = (((YBuffer[p]-BuffMin)*bitdepth)/(BuffMax-BuffMin));
    return (val-floor(val))>=0.5 ? (int)val + 1 : (int)val;*/
    //return (int) ((BuffMax-BuffMin)*YBuffer[p]/bitdepth);
    //Optimized method
    return (int) (((YBuffer[p]-BuffMin)*bitdepth)/(BuffMax-BuffMin));
}

void Texture_Image::BuildHistogram(){
    //Move this to the constructor...
    hist.reserve(bitdepth);
    for(int i=0;i<bitdepth;i++)
        hist.push_back(0);

    for(int i=0;i<Lines*Columns;i++){
        hist[YBuffer[i]]++;//TODO: Maybe use Transform here???????
    }
}

void Texture_Image::CalculateLUT(){ //Chiossif AnaCoding... TODO: Check for more that 8bit images
    //Need histogram to be built
    long sum, target;
    register int k,l;

    sum=0L;
    target=(long)Columns*(long)Lines/(long)quant;

    for (k=0,l=0;k<bitdepth;k++){//Changed 256 to bitdepth... check please!
        sum+=(long) hist[k];

        LUT[k]=(unsigned char)l;

        while ( sum > target ){
            l++;
            sum-=target;
        }

        if (sum<0L) sum=0L;
        if (l>=quant) l=quant-1;
    }
}

void Texture_Image::BuildCoCube(){ //TODO: Check for more that 8bit images
    //need LUT to be built
    int drow, dcol, angle, curp, nextp;
    uchar curval, nextval;

    //cout << "CoCube full with nulls" << endl;
    //cout << "CoCube size is " << CoCube.size() << endl;
    //cout << "CoCube table 2 size is " << CoCube[1].size() << endl;
    for(int i=0;i<Lines;i++){
        for(int j=0;j<Columns;j++){
            for(angle=0;angle<4;angle++){
                switch (angle){
                    case 0: drow =         0; dcol =  distance; break;
                    case 1: drow = -distance; dcol =  distance; break;
                    case 2: drow = -distance; dcol =         0; break;
                    case 3: drow = -distance; dcol = -distance; break;
                }
                if(i+drow>=0 && j+dcol>=0 && j+dcol<Columns){
                    curp = D1(Lines, Columns, i+1, j+1);
                    curval = LUT[YBuffer[curp]];
                    nextp = D1(Lines, Columns, i+drow+1, j+dcol+1);
                    nextval = LUT[YBuffer[nextp]];
                    CoCube[angle][curp].first = (int)curval;
                    CoCube[angle][curp].second = (int)nextval;
                }
            }
        }
    }
}

void Texture_Image::CheckCoCube(){ //TODO: Check for more that 8bit images
    int curp;
    for(int i=0;i<Lines;i++){
        for(int j=0;j<Columns;j++){
            for(int angle=0;angle<4;angle++){
                curp = D1(Lines, Columns, i+1, j+1);
                if(CoCube[angle][curp].first>=quant || CoCube[angle][curp].second>=quant){
                    cout << "problem for CoCube at position " << i+1 << " " << j+1 << endl;
                    cout << "first component is " << CoCube[angle][curp].first << " and second is " << CoCube[angle][curp].second << endl;
                }
            }
        }
    }
}


