/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * *
 * msegcore.cpp: MSEG base modules	                                  *
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

#include "msegcore.h"

string ftoa(float a){
    ostringstream ss;
    ss << a;
    return ss.str();
}

string itoa(int a){
    ostringstream ss;
    ss << a;
    return ss.str();
}

double sqr(double x){
    return x*x;
}

Object::Object(int x){
    id = x;
    merged = FALSE;
    edge = FALSE;
    priority_loaded = FALSE;
    virtual_merged = FALSE;
    properties_exported = FALSE;
    area = 0;
    perimeter = 0;
}

Object::Object(){
    id = 0;
    merged = FALSE;
    edge = FALSE;
    priority_loaded = FALSE;
    virtual_merged = FALSE;
    properties_exported = FALSE;
    area = 0;
    perimeter = 0;
}

Object::Object(const Object& O){
    id = O.id;
    merged = O.merged;
    edge = O.edge;
    priority_loaded = O.priority_loaded;
    virtual_merged = O.virtual_merged;
    //properties_exported = O.properties_exported;//We do not need it yet. Let it be by constructor
    area = O.area;
    perimeter = O.perimeter;

    for(int i=0; i<O.Boundary.size(); i++){
        Boundary.push_back(O.Boundary[i]);
    }
    for(int i=0; i<O.Sum.size(); i++){
        Sum.push_back(O.Sum[i]);
    }
    for(int i=0; i<O.SumSq.size(); i++){
        SumSq.push_back(O.SumSq[i]);
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
    edge = FALSE;
    priority_loaded = FALSE;
    virtual_merged = FALSE;
    properties_exported = FALSE;
    area = 0; //area will be 1 or 2 in 1st pass where this is useful
    perimeter = 0;
    //cycle = 0;
    Boundary.clear();
    Sum.clear();
    SumSq.clear();
    Neighbors.clear();
    SDM.clear();
    ASM.clear();
    R.clear();
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
	 float kokovios = 0.0;
	 merged = FALSE;
	 edge = FALSE;
	 priority_loaded = FALSE;
	 virtual_merged = FALSE;
	 perimeter = 0;
	 area = First.area + Second.area;
	 //calculating the new sum and the new sum of squares
	 for(int i=0; i<First.Sum.size(); i++){
	 	kokovios = First.Sum[i] + Second.Sum[i];
	 	Sum.push_back(kokovios);
	 	kokovios = First.SumSq[i] + Second.SumSq[i];
	 	SumSq.push_back(kokovios);
	 }
     //Boundary.clear();
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
	 for(int i=0; i<Boundtemp.size()-1; i++){//Here was the major pre-0.7.0 bug... the -1 was missing!
	     tmp.Line = Boundtemp[i].Line;
         tmp.ColumnStart = Boundtemp[i].ColumnStart;
         tmp.ColumnEnd = Boundtemp[i].ColumnEnd;
         while((Boundtemp[i].Line == Boundtemp[i+1].Line)&&((Boundtemp[i].ColumnEnd + 1) == Boundtemp[i+1].ColumnStart)){
             tmp.ColumnEnd = Boundtemp[i+1].ColumnEnd;
             i++;
             //in = TRUE;
         }
         /*if(Boundtemp[i].Line == -1){
             break;
         }
         */
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
    float h = 0.0;

//Obsolete method for calculating Color Heterogeneity
//Discontinued from version 0.7
/*    float s[Img.Bands];
    float ss[Img.Bands];
    for(int b=0; b<Img.Bands; b++){
        s[b] = 0.0;
        ss[b] = 0.0;
    }
    int l, cur_pos;
    MeanValue.clear();
    Std.clear();

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
*/
    for(int b=0; b<Img.Bands; b++){
        std[b] = (sqrt(SumSq[b]/area-Sum[b]*Sum[b]/area/area));
        //Maybe here include the area factor inside std number to get rid of two divisions
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

/*void Object::Calculate_Sums(){

}
*/

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

void Object::BuildSDM(Texture_Image& TexImg){
    int l,curpos;

    CoAngle tmp;
    for(int i=0;i<TexImg.quant * TexImg.quant;i++){
        tmp.push_back(0);
    }
    //cout << "SDM1" << endl;
    for(int i=0; i<4;i++){
        SDM.push_back(tmp);
        R.push_back(0L);
    }
    //cout << "SDM2" << endl;
    for(int i=0; i<Boundary.size(); i++){
        l = Boundary[i].Line;
        for(int c = Boundary[i].ColumnStart; c<=Boundary[i].ColumnEnd; c++){
            curpos = D1(TexImg.Lines, TexImg.Columns, l, c);
            for(int angle=0;angle<4;angle++){
                if(TexImg.CoCube[angle][curpos].first == -1) {
                }else{
                    SDM[angle][D1(TexImg.quant, TexImg.quant, TexImg.CoCube[angle][curpos].first + 1, TexImg.CoCube[angle][curpos].second + 1)]++;
                    R[angle]++;
                }
            }
        }
    }
}

void Object::ClearSDM(){
    SDM.clear();
}

void Object::CalculateASM(Texture_Image& TexImg){
    double dval;
    ASM.reserve(4);
    for(int angle=0;angle<4;angle++){
        dval = 0.0;
        if(R[angle] == 0L) {
            ASM.push_back(0.0);
        }else{
            for(int k=0;k<TexImg.quant;k++){
                for(int l=0;l<TexImg.quant;l++){
                    dval+=sqr((double)SDM[angle][D1(TexImg.quant, TexImg.quant, k+1, l+1)] / R[angle]);
                }
            }
            ASM.push_back(dval);
        }
    }
    /*
    if (feat[0])  /* Angular Second Moment *//*Chiossif anacoding...
     if (R==0L)
       ind[0][angle]=0.0;
     else
       {
       for (k=0;k<class;k++)
         for (l=0;l<class;l++)
           dval+=sqr((double)SDM[k][l]/R);
       ind[0][angle]=dval;
       }

    */
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
    Texture = 0.0;
    cycles = 0;
    ExistenceOfSuperlevel = FALSE;
    ExistenceOfSublevel = FALSE;
    Level* Superlevel = NULL;
    Level* Sublevel = NULL;
}

Level::Level(const Level& L){
    Objects = L.Objects;
    Properties = L.Properties;
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
    Texture = L.Texture;
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

void Level::UpdateRaster(int z, int nz){
    int l;
    for(int i=0; i<Objects[z].Boundary.size(); i++){
        l = Objects[z].Boundary[i].Line;
        for(int c = Objects[z].Boundary[i].ColumnStart; c<=Objects[z].Boundary[i].ColumnEnd; c++){
            raster[D1(Lines, Columns, l, c)] = nz;
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

void Level::SaveXML(string XMLFile, ERS_Image& Img){
	int cur_id;
	const char* declaration = "\t<?xml version=\"1.0\"?>\n";
    string tmp;
    vector<TiXmlElement> ObjectItem;
    vector<TiXmlElement> IdItem;
    vector<TiXmlElement> BandIdItem;
    vector<TiXmlElement> MeanItem;
    vector<TiXmlElement> StdItem;
    vector<TiXmlElement> BandItem;
    vector<TiXmlElement> AreaItem;
    vector<TiXmlElement> PerimeterItem;

    vector<string> tmpBandID;

	ObjectItem.reserve(3);
    BandItem.reserve(3);

    for(int b=0;b<Img.BandID.size();b++){
  		tmp = Img.BandID[b].substr(1,Img.BandID[b].size()-2);	//Removes the "___"
  		tmpBandID.push_back(tmp);
	}

	//cout << "Inside xml output routine" << endl;

	TiXmlDocument doc;
    doc.Parse(declaration);
    TiXmlHandle docHandle(&doc);
    assert(docHandle.Node());

	TiXmlElement Levelitem("Level");
	docHandle.Node()->InsertEndChild(Levelitem);

	TiXmlNode* LevelNode = docHandle.FirstChild("Level").Node();
	assert(LevelNode);

	TiXmlElement IDitem("LevelID");
	IDitem.SetAttribute("type", "float");
	tmp = ftoa(HierarchyID);
	TiXmlText text(tmp.c_str());
	IDitem.InsertEndChild(text);
	LevelNode->InsertEndChild(IDitem);

	TiXmlElement FilenameItem("RasterFilename");
	text = XMLFile.c_str();
	FilenameItem.InsertEndChild(text);
	LevelNode->InsertEndChild(FilenameItem);

	TiXmlElement LinesItem("Lines");
	LinesItem.SetAttribute("type", "int");
	tmp = itoa(Lines);
	text = tmp.c_str();
	LinesItem.InsertEndChild(text);
	LevelNode->InsertEndChild(LinesItem);

	TiXmlElement ColumnsItem("Columns");
	ColumnsItem.SetAttribute("type", "int");
	tmp = itoa(Columns);
	text = tmp.c_str();
	ColumnsItem.InsertEndChild(text);
	LevelNode->InsertEndChild(ColumnsItem);

	TiXmlElement ParametersItem("Parameters");
	LevelNode->InsertEndChild(ParametersItem);

	TiXmlNode* ParametersNode = docHandle.FirstChild("Level").FirstChild("Parameters").Node();
	assert(ParametersNode);

	TiXmlElement ScaleParameterItem("ScaleParameter");
	ScaleParameterItem.SetAttribute("type", "float");
	tmp = ftoa(ScaleParameter);
	text = tmp.c_str();
	ScaleParameterItem.InsertEndChild(text);
	ParametersNode->InsertEndChild(ScaleParameterItem);

	TiXmlElement ColorItem("Color");
	ColorItem.SetAttribute("type", "float");
	tmp = ftoa(Color);
	text = tmp.c_str();
	ColorItem.InsertEndChild(text);
	ParametersNode->InsertEndChild(ColorItem);

	TiXmlElement CompactnessItem("Compactness");
	CompactnessItem.SetAttribute("type", "float");
	tmp = ftoa(Compactness);
	text = tmp.c_str();
	CompactnessItem.InsertEndChild(text);
	ParametersNode->InsertEndChild(CompactnessItem);

	TiXmlElement EdgeItem("Edge");
	EdgeItem.SetAttribute("type", "float");
	tmp = ftoa(Edge);
	text = tmp.c_str();
	EdgeItem.InsertEndChild(text);
	ParametersNode->InsertEndChild(EdgeItem);

	TiXmlElement RasterItem("Raster");
	LevelNode->InsertEndChild(RasterItem);

	TiXmlElement PropertiesItem("Properties");
	LevelNode->InsertEndChild(PropertiesItem);

	TiXmlNode* PropertiesNode = docHandle.FirstChild("Level").FirstChild("Properties").Node();
	assert(PropertiesNode);

	TiXmlElement ObjItem("Object");
	ObjectItem.push_back(ObjItem);
	TiXmlNode* ObjectNode;

	TiXmlElement Id_Item("id");
	IdItem.push_back(Id_Item);

	TiXmlElement BndItem("Band");
	BandItem.push_back(BndItem);
	TiXmlNode* BandNode;

	TiXmlElement BandId_Item("BandID");
	BandIdItem.push_back(BandId_Item);

	TiXmlElement Area_Item("area");
	AreaItem.push_back(Area_Item);

	TiXmlElement Perimeter_Item("perimeter");
	PerimeterItem.push_back(Perimeter_Item);

	TiXmlElement Mean_Item("mean");
	Mean_Item.SetAttribute("type", "float");
	MeanItem.push_back(Mean_Item);

	TiXmlElement Std_Item("std");
	Std_Item.SetAttribute("type", "float");
	StdItem.push_back(Std_Item);

	//cout << "Reached step 1" << endl;

	for(int i=0; i<raster.size(); i++){
	    cur_id = raster[i];
	    if(Objects[cur_id].properties_exported == FALSE){
 		    //cout << "Writing object " << cur_id << endl;
 			PropertiesNode->InsertEndChild(ObjectItem[0]);
	    	ObjectNode = docHandle.FirstChild("Level").FirstChild("Properties").Node()->LastChild("Object");
			assert(ObjectNode);

			tmp = itoa(cur_id);
			text = tmp.c_str();
			IdItem[0].InsertEndChild(text);
			ObjectNode->InsertEndChild(IdItem[0]);
			IdItem.clear();
			IdItem.push_back(Id_Item);

			tmp = itoa(Objects[cur_id].area);
			text = tmp.c_str();
			AreaItem[0].InsertEndChild(text);
			ObjectNode->InsertEndChild(AreaItem[0]);
			AreaItem.clear();
			AreaItem.push_back(Area_Item);

			tmp = itoa(Objects[cur_id].perimeter);
			text = tmp.c_str();
			PerimeterItem[0].InsertEndChild(text);
			ObjectNode->InsertEndChild(PerimeterItem[0]);
			PerimeterItem.clear();
			PerimeterItem.push_back(Perimeter_Item);

 		    for(int b=0; b<Img.BandID.size();b++){
			    ObjectNode->InsertEndChild(BandItem[0]);
				BandNode = docHandle.FirstChild("Level").FirstChild("Properties").Node()->LastChild("Object")->LastChild("Band");
				assert(BandNode);

			    text = tmpBandID[b].c_str();
				BandIdItem[0].InsertEndChild(text);
				BandNode->InsertEndChild(BandIdItem[0]);
				BandIdItem.clear();
				BandIdItem.push_back(BandId_Item);

				text = ftoa(Properties[cur_id].MeanValue[b]).c_str();
				MeanItem[0].InsertEndChild(text);
				BandNode->InsertEndChild(MeanItem[0]);
				MeanItem.clear();
				MeanItem.push_back(Mean_Item);

				text = ftoa(Properties[cur_id].StdDev[b]).c_str();
				StdItem[0].InsertEndChild(text);
				BandNode->InsertEndChild(StdItem[0]);
	   			StdItem.clear();
	   			StdItem.push_back(Std_Item);

				BandNode = 0;
	   			BandItem.clear();
	   			BandItem.push_back(BndItem);
			}
 		    ObjectNode = 0;
			ObjectItem.clear();
			ObjectItem.push_back(ObjItem);
 		    Objects[cur_id].properties_exported = TRUE;
	    }
	}
	//cout << "Objects completed succesfully" << endl;
	tmp = XMLFile + ".xml";
	//cout << "Just before writing file" << endl;
	doc.SaveFile(tmp.c_str());
	cout << "xml writing completed" << endl;

	for(int i=0; i<raster.size(); i++){//Restoring the priorities_exported flags
	    cur_id = raster[i];
	    if(Objects[cur_id].properties_exported == TRUE){
 		    Objects[cur_id].properties_exported = FALSE;
	    }
	}
}

/*
//void Level::LoadXML(string XMLFile, int bands){

	 /* TO USE FOR OBJECTS
	//////////////////////////////////////////////////////
	printf ("\n** Iterators. **\n");

	// Walk all the top level nodes of the document.
	count = 0;
	for( node = doc.FirstChild();
		 node;
		 node = node->NextSibling() )
	{
		count++;
	}
	XmlTest( "Top level nodes, using First / Next.", 3, count );

	count = 0;
	for( node = doc.LastChild();
		 node;
		 node = node->PreviousSibling() )
	{
		count++;
	}
	XmlTest( "Top level nodes, using Last / Previous.", 3, count );

	// Walk all the top level nodes of the document,
	// using a different sytax.
	count = 0;
	for( node = doc.IterateChildren( 0 );
		 node;
		 node = doc.IterateChildren( node ) )
	{
		count++;
	}
	XmlTest( "Top level nodes, using IterateChildren.", 3, count );

	// Walk all the elements in a node.
	count = 0;
	for( element = todoElement->FirstChildElement();
		 element;
		 element = element->NextSiblingElement() )
	{
		count++;
	}
	XmlTest( "Children of the 'ToDo' element, using First / Next.",
		3, count );

	// Walk all the elements in a node by value.
	count = 0;
	for( node = todoElement->FirstChild( "Item" );
		 node;
		 node = node->NextSibling( "Item" ) )
	{
		count++;
	}
	XmlTest( "'Item' children of the 'ToDo' element, using First/Next.", 3, count );

	count = 0;
	for( node = todoElement->LastChild( "Item" );
		 node;
		 node = node->PreviousSibling( "Item" ) )
	{
		count++;
	}
	XmlTest( "'Item' children of the 'ToDo' element, using Last/Previous.", 3, count );
	*/
//}

void Level::SaveMiniXML(string XMLFile){

	const char* declaration = "\t<?xml version=\"1.0\"?>\n";
    string tmp;

	TiXmlDocument doc;
    doc.Parse(declaration);
    TiXmlHandle docHandle(&doc);
    assert(docHandle.Node());

	TiXmlElement Levelitem("Level");
	docHandle.Node()->InsertEndChild(Levelitem);

	TiXmlNode* LevelNode = docHandle.FirstChild("Level").Node();
	assert(LevelNode);

	TiXmlElement IDitem("LevelID");
	IDitem.SetAttribute("type", "float");
	tmp = ftoa(HierarchyID);
	TiXmlText text(tmp.c_str());
	IDitem.InsertEndChild(text);
	LevelNode->InsertEndChild(IDitem);

	TiXmlElement FilenameItem("RasterFilename");
	text = XMLFile.c_str();
	FilenameItem.InsertEndChild(text);
	LevelNode->InsertEndChild(FilenameItem);

	TiXmlElement LinesItem("Lines");
	LinesItem.SetAttribute("type", "int");
	tmp = itoa(Lines);
	text = tmp.c_str();
	LinesItem.InsertEndChild(text);
	LevelNode->InsertEndChild(LinesItem);

	TiXmlElement ColumnsItem("Columns");
	ColumnsItem.SetAttribute("type", "int");
	tmp = itoa(Columns);
	text = tmp.c_str();
	ColumnsItem.InsertEndChild(text);
	LevelNode->InsertEndChild(ColumnsItem);

	TiXmlElement ParametersItem("Parameters");
	LevelNode->InsertEndChild(ParametersItem);

	TiXmlNode* ParametersNode = docHandle.FirstChild("Level").FirstChild("Parameters").Node();
	assert(ParametersNode);

	TiXmlElement ScaleParameterItem("ScaleParameter");
	ScaleParameterItem.SetAttribute("type", "float");
	tmp = ftoa(ScaleParameter);
	text = tmp.c_str();
	ScaleParameterItem.InsertEndChild(text);
	ParametersNode->InsertEndChild(ScaleParameterItem);

	TiXmlElement ColorItem("Color");
	ColorItem.SetAttribute("type", "float");
	tmp = ftoa(Color);
	text = tmp.c_str();
	ColorItem.InsertEndChild(text);
	ParametersNode->InsertEndChild(ColorItem);

	TiXmlElement CompactnessItem("Compactness");
	CompactnessItem.SetAttribute("type", "float");
	tmp = ftoa(Compactness);
	text = tmp.c_str();
	CompactnessItem.InsertEndChild(text);
	ParametersNode->InsertEndChild(CompactnessItem);

	TiXmlElement EdgeItem("Edge");
	EdgeItem.SetAttribute("type", "float");
	tmp = ftoa(Edge);
	text = tmp.c_str();
	EdgeItem.InsertEndChild(text);
	ParametersNode->InsertEndChild(EdgeItem);

	TiXmlElement RasterItem("Raster");
	LevelNode->InsertEndChild(RasterItem);

	TiXmlElement PropertiesItem("Properties");
	LevelNode->InsertEndChild(PropertiesItem);

	tmp = XMLFile + ".xml";
	doc.SaveFile(tmp.c_str());

}

void Level::LoadMiniXML(string XMLFile){
	string tmp;
	TiXmlDocument doc(XMLFile.c_str());
    doc.LoadFile();

	TiXmlHandle docHandle(&doc);
    assert(docHandle.Node());

	TiXmlNode* LevelNode = docHandle.FirstChild("Level").Node();
	assert(LevelNode);

	TiXmlNode* tmpNode = docHandle.FirstChild("Level").FirstChild("LevelID").Node()->FirstChild();
	assert(tmpNode);
	tmp = tmpNode->Value();
    HierarchyID = atof(tmp.c_str());

    tmpNode = docHandle.FirstChild("Level").FirstChild("Lines").Node()->FirstChild();
    assert(tmpNode);
    tmp = tmpNode->Value();
    Lines = atoi(tmp.c_str());

    tmpNode = docHandle.FirstChild("Level").FirstChild("Columns").Node()->FirstChild();
    assert(tmpNode);
    tmp = tmpNode->Value();
    Columns = atoi(tmp.c_str());

    tmpNode = docHandle.FirstChild("Level").FirstChild("Parameters").FirstChild("ScaleParameter").Node()->FirstChild();
    assert(tmpNode);
    tmp = tmpNode->Value();
    ScaleParameter = atof(tmp.c_str());

    tmpNode = docHandle.FirstChild("Level").FirstChild("Parameters").FirstChild("Color").Node()->FirstChild();
    assert(tmpNode);
    tmp = tmpNode->Value();
    Color = atof(tmp.c_str());

    tmpNode = docHandle.FirstChild("Level").FirstChild("Parameters").FirstChild("Compactness").Node()->FirstChild();
    assert(tmpNode);
    tmp = tmpNode->Value();
    Compactness = atof(tmp.c_str());

    tmpNode = docHandle.FirstChild("Level").FirstChild("Parameters").FirstChild("Edge").Node()->FirstChild();
    assert(tmpNode);
    tmp = tmpNode->Value();
    Edge = atof(tmp.c_str());

    cout << "XML file loaded..." << endl;
}

void Level::CalculateProperties(ERS_Image& Img){
    PropertyObject tmp;
    int cur_id;
    float std, mean;

	/*A more delicate way to do this would be to use itterators. This way we
	will not need the extra byte for boolean member properties_exported, used in Objects
	The code to use in the upgrade is:*/

	typedef map<int,Object>::const_iterator CI;
	for(CI p=Objects.begin();p!=Objects.end();++p){
 	    //DO STUFF HERE with reference to p->first, p->second at each pair
 	    cur_id = p->first;
 	    for(int b=0; b<Img.Bands; b++){
 	    	std = (sqrt(p->second.SumSq[b] / p->second.area - p->second.Sum[b] * p->second.Sum[b] / p->second.area / p->second.area));
            mean = p->second.Sum[b] / p->second.area;
            tmp.StdDev.push_back(std);
            tmp.MeanValue.push_back(mean);
 	    }
 	    Properties.insert(make_pair(cur_id, tmp));
 	    tmp.clear();
	}


	/* Old code, before the upgrade
	for(int i=0; i<raster.size(); i++){
	    cur_id = raster[i];
	    if(Objects[cur_id].properties_exported == FALSE){
		    if(!Objects[cur_id].MeanValue.size()){
			    Objects[cur_id].MeanValue.reserve(Img.Bands);
				Objects[cur_id].Calculate_MeanValue(Img);
			}
			if(!Objects[cur_id].Std.size()){
			    Objects[cur_id].Std.reserve(Img.Bands);
				Objects[cur_id].Calculate_Std(Img);
			}
 			for(int b=0; b<Img.Bands;b++){
				tmp.MeanValue.push_back(Objects[cur_id].MeanValue[b]);
				tmp.StdDev.push_back(Objects[cur_id].Std[b]);
			}
 		    Properties.insert(make_pair(cur_id, tmp));
 		    tmp.clear();
 		    Objects[cur_id].properties_exported = TRUE;
 		}
	}

    for(int i=0; i<raster.size(); i++){//Restoring the priorities_exported flags
	    cur_id = raster[i];
	    if(Objects[cur_id].properties_exported == TRUE){
 		    Objects[cur_id].properties_exported = FALSE;
	    }
	}*/
}

void Level::DeleteProperties(){
	Properties.clear();
}

/*
void Level::CalculateTextureProperties(Texture_Image& TexImg){

}
*/

void Level::SaveProperties(string CSVFile){
	 string FileName = CSVFile + ".csv";
	 ofstream out(FileName.c_str());
	 int cur_id;
	 typedef map<int, PropertyObject>::const_iterator PO;

	 for(PO p=Properties.begin();p!=Properties.end();++p){
        cur_id = p->first;
        //case all properties
        //out << cur_id << ", " << Objects[cur_id].area << ", " << Objects[cur_id].perimeter;
        //case only spectral
        out << cur_id;
        for(int b=0; b<p->second.MeanValue.size();b++){
		    out << ", " << p->second.MeanValue[b];
		    //out << ", " << p->second.StdDev[b];
		}
 		out << endl;
	 }
	 out.close();
	 //Old code, before 0.7.0
	 /*for(int i=0; i<raster.size(); i++){
	    cur_id = raster[i];
	    if(Objects[cur_id].properties_exported == FALSE){
 		    out << cur_id << ", " << Objects[cur_id].area << ", " << Objects[cur_id].perimeter;
	 		for(int b=0; b<Properties[cur_id].MeanValue.size();b++){
			    out << ", " << Properties[cur_id].MeanValue[b];
			    out << ", " << Properties[cur_id].StdDev[b];
			}
 		    out << endl;
		 	Objects[cur_id].properties_exported = TRUE;
 		}
	}
    for(int i=0; i<raster.size(); i++){//Restoring the priorities_exported flags
	    cur_id = raster[i];
	    if(Objects[cur_id].properties_exported == TRUE){
 		    Objects[cur_id].properties_exported = FALSE;
	    }
	}
	*/
}

void Level::SaveSampleProperties(string CSVFile){
    string FileName = CSVFile + ".csv";
	ofstream out(FileName.c_str());
	int cur_id;
    for(int i=0; i<raster.size(); i++){
	    cur_id = raster[i];
	    if(SampleMap[i]!=0 && SampleMap[i]!=255 && Objects[cur_id].properties_exported == FALSE){
 		    //case all properties
 		    //out << SampleMap[i] << ", " << Objects[cur_id].area << ", " << Objects[cur_id].perimeter;
	 		//case only spectral
	 		out << SampleMap[i];

	 		for(int b=0; b<Properties[cur_id].MeanValue.size();b++){
			    out << ", " << Properties[cur_id].MeanValue[b];
			    //out << ", " << Properties[cur_id].StdDev[b];
			}
 		    out << endl;
		 	Objects[cur_id].properties_exported = TRUE;
 		}
	}
    for(int i=0; i<raster.size(); i++){//Restoring the priorities_exported flags
	    cur_id = raster[i];
	    if(Objects[cur_id].properties_exported == TRUE){
 		    Objects[cur_id].properties_exported = FALSE;
	    }
	}
	out.close();
}

void Level::SaveSVMTesting(string CSVFile){
	 string FileName = CSVFile + ".test";
	 ofstream out(FileName.c_str());
	 int cur_id;
	 typedef map<int, PropertyObject>::const_iterator PO;

	 for(PO p=Properties.begin();p!=Properties.end();++p){
        cur_id = p->first;
        //case all properties
        //out << cur_id << ", " << Objects[cur_id].area << ", " << Objects[cur_id].perimeter;
        //case only spectral
        out << cur_id << " ";
        for(int b=0; b<p->second.MeanValue.size();b++){
		    out << b+1 << ":" << p->second.MeanValue[b] << " ";
		    //out << ", " << p->second.StdDev[b];
		}
 		out << endl;
	 }
	 out.close();
}

void Level::SaveSVMTraining(string CSVFile){
    string FileName = CSVFile + ".train";
	ofstream out(FileName.c_str());
	int cur_id;
    for(int i=0; i<raster.size(); i++){
	    cur_id = raster[i];
	    if(SampleMap[i]!=0 && SampleMap[i]!=255 && Objects[cur_id].properties_exported == FALSE){
 		    //case all properties
 		    //out << SampleMap[i] << ", " << Objects[cur_id].area << ", " << Objects[cur_id].perimeter;
	 		//case only spectral
	 		out << SampleMap[i] << " ";

	 		for(int b=0; b<Properties[cur_id].MeanValue.size();b++){
			    out << b+1 << ":" << Properties[cur_id].MeanValue[b] << " ";
			    //out << ", " << Properties[cur_id].StdDev[b];
			}
 		    out << endl;
		 	Objects[cur_id].properties_exported = TRUE;
 		}
	}
    for(int i=0; i<raster.size(); i++){//Restoring the priorities_exported flags
	    cur_id = raster[i];
	    if(Objects[cur_id].properties_exported == TRUE){
 		    Objects[cur_id].properties_exported = FALSE;
	    }
	}
	out.close();
}

void Level::SaveSVMTrainingVoting(string CSVFile, double percent){
    string FileName = CSVFile + ".train";
	ofstream out(FileName.c_str());
	int cur_id;
	int sample_id;
	int occurences, c, l;
	for(int i=0; i<raster.size(); i++){
	    cur_id = raster[i];
	    if(SampleMap[i]!=0 && SampleMap[i]!=255 && Objects[cur_id].properties_exported == FALSE){
	        sample_id = SampleMap[i];
	        occurences = 0;
	        for(int j=0;j<Objects[cur_id].Boundary.size();j++){
                l = Objects[cur_id].Boundary[j].Line;
                for(c=Objects[cur_id].Boundary[j].ColumnStart;c<=Objects[cur_id].Boundary[j].ColumnEnd;c++){
                    if(SampleMap[D1(Lines, Columns, l, c)] == sample_id) occurences++;
                }
	        }
	        if((occurences / Objects[cur_id].area) >= percent){
                //case all properties
                //out << sample_id << ", " << Objects[cur_id].area << ", " << Objects[cur_id].perimeter;
                //case only spectral
                out << sample_id << " ";
                for(int b=0; b<Properties[cur_id].MeanValue.size();b++){
			    out << b+1 << ":" << Properties[cur_id].MeanValue[b] << " ";
			    //out << ", " << Properties[cur_id].StdDev[b];
                }
                out << endl;
                Objects[cur_id].properties_exported = TRUE;
	        }
	    }
	}
    for(int i=0; i<raster.size(); i++){//Restoring the priorities_exported flags
	    cur_id = raster[i];
	    if(Objects[cur_id].properties_exported == TRUE){
 		    Objects[cur_id].properties_exported = FALSE;
	    }
	}
	out.close();
}

void Level::SaveRaster(string RasterFile, ERS_Image& Img){
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
    head_out << "\t\tDatum\t\t= \"" << Img.CoordinateSpace.Datum << "\"" << endl;
    head_out << "\t\tProjection\t= \"" << Img.CoordinateSpace.Projection << "\"" << endl;
    head_out << "\t\tCoordinateType\t= " << Img.CoordinateSpace.CoordinateType << endl;
    head_out << "\t\tUnits\t\t= \"" << Img.CoordinateSpace.Units << "\"" << endl;
    head_out << "\t\tRotation\t= " << Img.CoordinateSpace.Rotation << endl;
    head_out << "\tCoordinateSpace End" << endl;
    head_out << "\tRasterInfo Begin" << endl;
    head_out << "\t\tCellType\t= " << "Unsigned32BitInteger" << endl;
    head_out << "\t\tNrOfLines\t= " << Lines << endl;
    head_out << "\t\tNrOfCellsPerLine\t= " << Columns << endl;
    if(Img.RegistrationPoint.x != -999 || Img.RegistrationPoint.y != -999)
    {
	head_out << "\t\tRegistrationCoord Begin" << endl;
	head_out << "\t\t\tEastings\t= " << Img.RegistrationPoint.x << endl;
	head_out << "\t\t\tNorthings\t= " << Img.RegistrationPoint.y << endl;
	head_out << "\t\tRegistrationCoord End" << endl;
    }
    
    if(Img.RegistrationCell.x != -999 || Img.RegistrationCell.y != -999)
    {
	head_out << "\t\tRegistrationCellX\t= " << Img.RegistrationCell.x << endl;
	head_out << "\t\tRegistrationCellY\t= " << Img.RegistrationCell.y << endl;
    }
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

void Level::LoadRaster(string RasterFile){
    cout << "Loading Level " << RasterFile << " ..." << endl;
    string HeaderFile = RasterFile + ".ers";
    ifstream header;
    header.open(HeaderFile.c_str());
    if (!header){								//Error handling
	    cout << "Error opening file " << RasterFile << endl;
	    //throw FileError;
     	cout << "Press any key to quit..." << endl;
	    cin.get();
	    exit(1);
	}
	string word;
	//Header file parser
 	while (header){
 	    header >> word;
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
           cout << "Level Number Of Lines : " << Lines << endl;
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
           cout << "Level Number Of Pixels per Line : " << Columns << endl;
       }
    }//Header file parser END
    header.close();
 	cout << "Header closed" << endl;
 	ifstream file;
    //Opening binary data file for input
    file.open(RasterFile.c_str(), ios::binary);
    if (!file){								//Error handling
	    cout << "Error opening file " << RasterFile << endl;
	    //throw FileError;

      	cout << "Press any key to quit..." << endl;
	    cin.get();
	    exit(1);
    }
    cout << "Binary Raster File Opened... " << endl;
    raster.clear();
    int tmp;
    file.seekg(0, ios::beg);
    for(int i=0; i<(Lines*Columns); i++){
	    file.read(reinterpret_cast<char *>(&tmp), sizeof(int));//TODO:Check that it works ok!!!!!
	    raster.push_back(tmp);
	}
	cout << "Binary Raster Read!" << endl;
	file.close();
}

void Level::LoadObjects(){
	Object tmp;
	BoundaryLine tmpBL;
	int cur_id;
	map<int,Object>::const_iterator p = Objects.begin();
	for(int l=1;l<=Lines;l++){
	    tmpBL.Line = l;
		for(int c=1;c<=Columns;c++){
		    cur_id = raster[D1(Lines, Columns, l, c)];
		    p=Objects.find(cur_id);
			if(p==Objects.end()){//If Object id is not found inside the map
			    tmp.id = cur_id;
			    Objects.insert(make_pair(tmp.id, tmp));
			    //tmp.clear();//Do not need
			}
			Objects[cur_id].area++;
			tmpBL.ColumnStart = c;
			while(c<Columns || raster[D1(Lines, Columns, l, c+1)]==cur_id){
		 	    c++;
		 	    Objects[cur_id].area++;
			}
			/*If a bug is present, take a look at the following...
			while(){
			    if(c>=Columns) break;
			    if(raster[D1(Lines, Columns, l, c+1)]!=cur_id) break;
			    c++;
			}
			*/
		    tmpBL.ColumnEnd = c;
		    Objects[cur_id].Boundary.push_back(tmpBL);
	    }
	}
}

void Level::LoadObjectSums(ERS_Image& Img){

    int cur_id, b;
    double cur_val;
    for(int l=1;l<=Img.Lines;l++){
        for(int c=1;c<=Img.Columns;c++){
            cur_id = raster[D1(Lines, Columns, l, c)];
            if(Objects[cur_id].Sum.size() == 0){
                for(b=1;b<=Img.Bands;b++){
                    Objects[cur_id].Sum.push_back(0.0);
                    Objects[cur_id].SumSq.push_back(0.0);
                }
            }
            for(b=1;b<=Img.Bands;b++){
                cur_val = Img.Buffer(l, c, b);
                Objects[cur_id].Sum[b-1] += cur_val;
                Objects[cur_id].SumSq[b-1] += (cur_val*cur_val);
            }
        }
    }
}

void Level::CalculateNeighbors(){

    Object tmp;
    BoundaryLine tmpBoundaryLine;
	int cur_id, up_id, down_id, right_id, left_id;
    pixel cur, up, down, right, left;

    Neighbor n;
    n.merging_heterogeneity = 0.0;
    for(int row=1; row<=Lines; row++){
        for(int col=1; col<=Columns; col++){
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

            cur_id = raster[D1(Lines, Columns, cur.line, cur.column)];
            if ((up.line>0 && up.line<=Lines) && (up.column>0 && up.column<=Columns)){
                up_id = raster[D1(Lines, Columns, up.line, up.column)];
                if(cur_id != up_id){
                    n.id = up_id;
                    if(Objects[cur_id].check_neighbor(up_id) == FALSE){
                        Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((down.line>0 && down.line<=Lines) && (down.column>0 && down.column<=Columns)){
                down_id = raster[D1(Lines, Columns, down.line, down.column)];
                if(cur_id != down_id){
                    n.id = down_id;
                    if(Objects[cur_id].check_neighbor(down_id) == FALSE){
                        Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((right.line>0 && right.line<=Lines) && (right.column>0 && right.column<=Columns)){
                right_id = raster[D1(Lines, Columns, right.line, right.column)];
                if(cur_id != right_id){
                    n.id = right_id;
                    if(Objects[cur_id].check_neighbor(right_id) == FALSE){
                        Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
            if ((left.line>0 && left.line<=Lines) && (left.column>0 && left.column<=Columns)){
                left_id = raster[D1(Lines, Columns, left.line, left.column)];
                if(cur_id != left_id){
                    n.id = left_id;
                    if(Objects[cur_id].check_neighbor(left_id) == FALSE){
                        Objects[cur_id].Neighbors.push_back(n);
                    }
                }
            }
        }
    }
}

void Level::DeleteObjects(){
	Objects.clear();
}

/*void Level::CreateERS(string ImageFileName, ERS_Image& Img, int BandNumber){

}*/

void Level::LoadTTA(string TTAName){
    string tta;
    tta = TTAName + ".asc";
    ifstream ttafile;
    ttafile.open(tta.c_str());
    string word;
    int a;
    ttafile >> word;
    ttafile >> a;
    if (a != Columns){
        cout << "Error on column number" << endl;
    }
    ttafile >> word;
    ttafile >> a;
    if (a != Lines){
        cout << "Error on row number" << endl;
    }
    for (int i=0;i<5;i++){
        ttafile >> word;
        ttafile >> a;
    }
    SampleMap.reserve(Lines*Columns);
    for(int i=0;i<Lines*Columns;i++){
        ttafile >> a;
	//cout << "Reading TTA Value: " << a << endl;
        SampleMap.push_back(a);
    }
    ttafile.close();
}

void Level::LoadAttributes(string TTAName){
    string attr;
    vector<string> kokos;
    ClassAttribute tmp;
    attr = TTAName + ".txt";
    ifstream attrfile;
    attrfile.open(attr.c_str());

    string word;
    getline(attrfile, word);//skips the first line
    while(attrfile){
        getline(attrfile, word);
        int last = word.size();
        // The beginning of the current word:
        size_t current = word.rfind(',');
        // Walk backward through the string:
        while(current != string::npos) {
            // Push each word into the vector.
            // Current is incremented before copying
            // to avoid copying the delimiter:
            ++current;
            kokos.push_back(word.substr(current, last - current));
            // Back over the delimiter we just found,
            // and set last to the end of the next word:
            current -= 2;
            last = current + 1;
            // Find the next delimiter:
            current = word.rfind(',', current);
        }
        // Pick up the first word -- it's not
        // preceded by a delimiter:
        kokos.push_back(word.substr(0, last));
    }
    attrfile.close();
    for(int i=0;i<kokos.size()-5;i+=5){
        tmp.Name = kokos[i];
        //cout << tmp.Name << endl;
        tmp.B = atoi(kokos[i+1].c_str());
        //cout << tmp.B << endl;
        tmp.G = atoi(kokos[i+2].c_str());
        //cout << tmp.G << endl;
        tmp.R = atoi(kokos[i+3].c_str());
        //cout << tmp.R << endl;
        tmp.id = atoi(kokos[i+4].c_str());
        //cout << tmp.id << endl;
        SampleAttributes.insert(make_pair(tmp.id, tmp));
    }
}

void Level::SaveMeanAsERS(string ImageFileName, ERS_Image& Img){
    string tmp;

    ERS_Image Img2(ImageFileName, Img.Bands, Img.Lines, Img.Columns, Img.Bytes);
    if(Img.Bytes == 1) Img2.Buffer8.reserve(Img2.Lines*Img2.Columns*Img2.Bands);
    else if(Img.Bytes == 2) Img2.Buffer16.reserve(Img2.Lines*Img2.Columns*Img2.Bands);
    else if(Img.Bytes == 3) Img2.Buffer16s.reserve(Img2.Lines*Img2.Columns*Img2.Bands);
    else if(Img.Bytes == 4) Img2.Buffer32.reserve(Img2.Lines*Img2.Columns*Img2.Bands);

    for(int i=0; i<Img.BandID.size(); i++){
	    tmp = Img.BandID[i];
	    Img2.BandID.push_back(tmp);
    }

    int cur_id;
    uchar val8;
    ushort val16;
    sshort val16s;

    float val32;

    for(int l=1; l<=Img2.Lines; l++){
        for(int b=1; b<=Img2.Bands; b++){
            for(int c=1; c<=Img2.Columns; c++){
            	cur_id = raster[D1(Lines, Columns, l, c)];
                val32 = Properties[cur_id].MeanValue[b-1];
                if(Img.Bytes == 1)
                {
                    val32 += 0.5;
                    val8 = (uchar)val32;
                    Img2.Buffer8.push_back(val8);
                }
                else if(Img.Bytes == 2)
                {
                    val32 += 0.5;
                    val16 = (ushort)val32;
                    Img2.Buffer16.push_back(val16);
                }
                else if(Img.Bytes == 3)
                {
                    val32 += 0.5;
                    val16s = (sshort)val32;
                    Img2.Buffer16s.push_back(val16s);
                }
                else if(Img.Bytes == 4)
                {
                    Img2.Buffer32.push_back(val32);
                }
           }
        }
    }
	Img2.SaveAs(ImageFileName);
}

void Level::SVM2ClassificationMap(string TestFile, string PredictionFile){

    Properties.clear();//SOS SOS SOS Must clear here!!! This routine is used when level is created back...
    ClassificationMap.clear();
    ClassificationMap.reserve(Lines*Columns);

    ifstream tfile;
    tfile.open(TestFile.c_str());
    ifstream pfile;
    pfile.open(PredictionFile.c_str());

    string t;
    int p, id;
    PropertyObject tmp;
    string word;

    while(pfile){
        pfile >> p;
        getline(tfile, word);
        size_t current = word.find(' ');
        t = word.substr(0, current);
        id = atoi(t.c_str());
        tmp.classid = p;
        if(id!=0){
            Properties.insert(make_pair(id, tmp));
        }
    }

    for(int i=0;i<Lines*Columns;i++){
        ClassificationMap.push_back(Properties[raster[i]].classid);
    }
}

void Level::SaveClassificationMap(string ClassFile){
    string HeadFile = ClassFile + ".ers";

    ofstream out(ClassFile.c_str(), ios::binary);
    ofstream head_out(HeadFile.c_str());

    //save binary data
    for(int i=0; i<ClassificationMap.size(); i++){
        out.write(reinterpret_cast<char *> (&ClassificationMap[i]), sizeof(int));
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

void Level::SaveClassificationAsERS(string ClassFile){
    string out_cell_type = "Unsigned8BitInteger";
    string tmp;

    ERS_Image Img(ClassFile, 3, Lines, Columns, 1);
    Img.Buffer8.reserve(Lines*Columns*3);
    tmp="\"Blue\"";
    Img.BandID.push_back(tmp);
    tmp="\"Green\"";
    Img.BandID.push_back(tmp);
    tmp="\"Red\"";
    Img.BandID.push_back(tmp);


    int cur_id;
    uchar val;
    float m;

    for(int l=1; l<=Lines; l++){
        for(int b=1; b<=3; b++){
            for(int c=1; c<=Columns; c++){
            	cur_id = raster[D1(Lines, Columns, l, c)];
                if(b==1){
                    val = (unsigned char)(SampleAttributes[ClassificationMap[D1(Lines, Columns, l, c)]].B);
                }else if(b==2){
                    val = (unsigned char)(SampleAttributes[ClassificationMap[D1(Lines, Columns, l, c)]].G);
                }else if(b==3){
                    val = (unsigned char)(SampleAttributes[ClassificationMap[D1(Lines, Columns, l, c)]].R);
                }
                Img.Buffer8.push_back(val);
            }
        }
    }
	//cout << "Output image buffer size is: " << Img2.Buffer.size() << endl;
    Img.SaveAs(ClassFile);
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
    Img2.Buffer8.reserve(Img2.Lines*Img2.Columns*Img2.Bands);
    Img2.BandID.push_back(tmp);

    uchar val;

    for(int i=0; i<(Lines*Columns); i++){
        val = (uchar)BoundaryMap[i];
        Img2.Buffer8.push_back(val);
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
    SampleMap.clear();
    ClassificationMap.clear();
    SampleAttributes.clear();
}

PropertyObject::PropertyObject(){

}

PropertyObject::~PropertyObject(){

}

void PropertyObject::clear(){
	 classid = 0;
	 MeanValue.clear();
	 StdDev.clear();
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
    Texture = 0.0;
    GHH = 0;
    EC = FALSE;
    MA = TRUE;
    TH = FALSE;
}

MSEG_Param::MSEG_Param(float Sc, float Clr, float Cmp, float Edg,float Tex, int quant, int dist, bool ec, bool ma, short ghh, bool tx){
    Scale = Sc;
    Color = Clr;
    Compact = Cmp;
    Edge = Edg;
    Texture = Tex;
    quantizer = quant;
    distance = dist;
    GHH = ghh;
    MA = ma;
    EC = ec;
    TH = tx;
}

MSEG_Param::MSEG_Param(float Sc, float Clr, float Cmp){
    Scale = Sc;
    Color = Clr;
    Compact = Cmp;
    Edge = 0.0;
    Texture = 0.0;
    quantizer = 32;
    distance = 1;
    GHH = 0;
    EC = FALSE;
    MA = TRUE;
    TH = FALSE;
}

MSEG_Param::~MSEG_Param(){
}

void MSEG_Init::SaveSPE(string ImageFileName, ERS_Image& Img){

    pixel tmp;
    string out_cell_type = "Unsigned8BitInteger";
    string tmp1;
    tmp1 = "Band1";

    ERS_Image Img2(ImageFileName, 1, Img.Lines, Img.Columns, out_cell_type);
    Img2.Buffer8.reserve(Img2.Lines*Img2.Columns*Img2.Bands);
    Img2.BandID.push_back(tmp1);

    for(int i=0; i<(Img2.Lines*Img2.Columns); i++){
        Img2.Buffer8.push_back(0);
    }

    for(int j=0; j<StartingPoints.size(); j++){
	    tmp = StartingPoints[j];
	    Img2.Buffer8[D1(Img2.Lines, Img2.Columns, tmp.line, tmp.column)] = 1;
	}

    Img2.SaveAs(ImageFileName, 0);

}

MSEG_Init::MSEG_Init(){
    last_id = 0;
}

MSEG_Init::~MSEG_Init(){

}

/*
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
*/

Starting_Points_Estimation::Starting_Points_Estimation(int mode){
    SPE_mode = mode;
}

Starting_Points_Estimation::~Starting_Points_Estimation(){

}

void Starting_Points_Estimation::HSI(ERS_Image& Img, MSEG_Init& mseg){
    cout << "HSI method selected for start point estimation" << endl;
    vector<float> Intensity_band;
    Intensity_band.reserve(Img.Lines * Img.Columns);

    if(Img.Bands == 1){
    	Img.FullBuffer();
     	for(int i=0; i<(Img.Lines*Img.Columns); i++){
      		Intensity_band[i] = Img.Buffer(i);
        }
    }
    else if(Img.Bands == 2){
    	Img.FullBuffer();
     	for(int l=1; l<=Img.Lines; l++){
      		for(int c=1; c<=Img.Columns; c++){
        		Intensity_band[((l-1)*Img.Columns + (c-1))] = float((Img.Buffer(((l-1)*Img.Columns*Img.Bands)+c-1) + Img.Buffer(((l-1)*Img.Columns*Img.Bands)+c+Img.Columns-1))/2);
          	}
        }
    }
    else if(Img.Bands > 2){
    	Img.FullBuffer();
    	for(int l=1; l<=Img.Lines; l++){
      		for(int c=1; c<=Img.Columns; c++){
        		Intensity_band[((l-1)*Img.Columns + (c-1))] = float((Img.Buffer(((l-1)*Img.Columns*Img.Bands)+c-1) + Img.Buffer(((l-1)*Img.Columns*Img.Bands)+c+Img.Columns-1) + Img.Buffer(((l-1)*Img.Columns*Img.Bands)+c+Img.Columns+Img.Columns-1))/3);
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
    //delete &Intensity_band; //Doesn't work on Linux...
}

void Starting_Points_Estimation::YUV(ERS_Image& Img, MSEG_Init& mseg){
    cout << "YUV method selected for start point estimation" << endl;
    vector<float> Intensity_band;
    Intensity_band.reserve(Img.Lines * Img.Columns);

    if(Img.Bands == 1){
    	cout << "Not enough bands to calculate YUV components... Using HSI by default" << endl;
		Img.FullBuffer();
     	for(int i=0; i<(Img.Lines*Img.Columns); i++){
      		Intensity_band[i] = Img.Buffer(i);
        }
    }
    else if(Img.Bands == 2){
    	cout << "Not enough bands to calculate YUV components... Using HSI by default" << endl;
		Img.FullBuffer();
     	for(int l=1; l<=Img.Lines; l++){
      		for(int c=1; c<=Img.Columns; c++){
        		Intensity_band[((l-1)*Img.Columns + (c-1))] = float((Img.Buffer(((l-1)*Img.Columns*Img.Bands)+c-1) + Img.Buffer(((l-1)*Img.Columns*Img.Bands)+c+Img.Columns-1))/2);
          	}
        }
    }
    else if(Img.Bands > 2){
    	Img.FullBuffer();
     	for(int l=1; l<=Img.Lines; l++){
      		for(int c=1; c<=Img.Columns; c++){
        		Intensity_band[((l-1)*Img.Columns + (c-1))] = float((0.114*Img.Buffer(((l-1)*Img.Columns*Img.Bands)+c-1)) + (0.587*Img.Buffer(((l-1)*Img.Columns*Img.Bands)+c+Img.Columns-1)) + (0.299*Img.Buffer(((l-1)*Img.Columns*Img.Bands)+c+Img.Columns+Img.Columns-1)));
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
    //delete &Intensity_band;
}

void Starting_Points_Estimation::Dithering(ERS_Image& Img, MSEG_Init& mseg, string Method){
	float Buffmax = -2000000000.0;
	float Buffmin = 2000000000.0;
	uchar val;
	int val16;
	float val32;

	cout << "Dithering matrix method selected for start point estimation" << endl;
 	Img.FullBuffer();
 	string tmp("test");
 	//return (int) (((YBuffer[p]-BuffMin)*bitdepth)/(BuffMax-BuffMin));
 	ERS_Image Img2(tmp, 1, Img.Lines, Img.Columns, 1);
 	Img2.Buffer8.reserve(Img2.Lines*Img2.Columns);

 	if(Img.Bytes == 1)
 	{
 	    for(int l=1; l<=Img.Lines; l++)
 	    {
 	        for(int c=1; c<=Img.Columns; c++)
 	        {
 	            val = Img.Buffer8[Img2.D1(l, c, 1)];
                Img2.Buffer8.push_back(val);
 	        }
 	    }
 	}
 	else
 	{
 	    for(int l=1; l<=Img.Lines; l++)
 	    {
 	        for(int c=1; c<=Img.Columns; c++)
 	        {
 	            val32 = Img.Buffer(l, c, 1);
 	            if(val32 > Buffmax) Buffmax = val32;
 	            if(val32 < Buffmin) Buffmin = val32;
 	        }
 	    }

 	    for(int l=1; l<=Img.Lines; l++)
 	    {
 	        for(int c=1; c<=Img.Columns; c++)
 	        {
                val32 = Img.Buffer(l, c, 1);
                val16 = (int) (((val32-Buffmin)*255)/(Buffmax-Buffmin));
                val = (uchar)val16;
                Img2.Buffer8.push_back(val);
 	        }
 	    }

 	}

    FIBITMAP *image = FreeImage_Allocate(Img2.Columns, Img2.Lines, 8);

    BYTE *bits[Img2.Lines];
    for(int y = 0; y<Img2.Lines; y++){
		bits[y] = (BYTE*)FreeImage_GetScanLine(image, y);
	}

  	for(int y = 1; y<=Img2.Lines; y++){
		for(int x = 1; x<=Img2.Columns; x++){
		    bits[y-1][x-1] = Img.Buffer8[Img2.D1(y, x, 1)];
		}
	}

  	RGBQUAD *pal = FreeImage_GetPalette(image);
	for(int i=0;i<256;i++){
	    pal[i].rgbRed = i;
		pal[i].rgbGreen = i;
		pal[i].rgbBlue = i;
	}

	FIBITMAP *dithered;
	//Add a case here to choose the method to initialize the FIBITMAP	--> See next comment
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


	BYTE *bits_dith[Img2.Lines];
	for(int y = 0; y<Img2.Lines; y++){
	    bits_dith[y] = (BYTE*)FreeImage_GetScanLine(dithered, y);
	}

  	pixel pix;

	for(int l=1; l<=Img2.Lines; l++){
	    for(int c=1; c<=Img2.Columns; c++){
		    if(bits_dith[l-1][c-1] == TRUE){
			    pix.line = l;
			    pix.column = c;
			    mseg.StartingPoints.push_back(pix);
	        }
	    }
	}

	FreeImage_Unload(dithered);

}

