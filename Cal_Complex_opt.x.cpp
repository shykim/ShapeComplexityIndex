/***************************************************************************
 *   Copyright (C) 2015 by Sun Hyung Kim and Martin Styner		 					   *
 *   NeuroImage Analysis and Research Lab				  												 *
 *   Dept. of Psychiatry, University of North Carolina at Chapel Hill      *
 *   shykim@email.unc.edu	                             			               *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Mesh.h>
#include <Util/SurfaceUtil.h>
//#include "Geodesic/Geodesic.h"
#include "Calculators.h"
#include "ShapeComplexityIndexCLP.h"

#define BOOL bool
#define TRUE true
#define FALSE false
#define TEMPLATE_SA 74291 //SA of surf_reg_model_left.obj
//#define TEMPLATE_SA 48410.9 //SA of AVGmid_LEFT.obj, 12month IBIS data set.

#include <ctime>

using namespace std;

int main(int argc, char *argv[])
{
	
	PARSE_ARGS;
	if(argc<2)
	{
		cout << "Cal_Complex --help" << endl;
		return false;
	}else{
		cout << "Input: " << InputSurface.c_str() << endl;	
		cout << "# of bin: " << pow(static_cast<double>(2),StepSize) + 1 << endl;
	}
	
	if(InputTemplateSurface.size() == 0 )
	{
		//InputTemplateSurface = "@CMAKE_CURRENT_BINARY_DIR@/template.obj";
		InputTemplateSurface = "@CMAKE_BINARY_DIR@/template.obj";
	}
	
	//cout << InputTemplateSurface << endl;
	

	CCalculators* Cal = new CCalculators;
	Mesh *MNI_GetMesh = new Mesh();
	MNI_GetMesh->openFile(InputSurface.c_str());

	
	
	//template surface
	Mesh *Template = new Mesh();
	Template->openFile(InputTemplateSurface.c_str());
	int tn_vertex = Template->nVertex();
	int tn_mesh = Template -> nFace();
	float *T_SA = new float[tn_vertex];
	float SUM_T_SA = 0.0f;
	float MEAN_T_SA = 0.0f;
	SUM_T_SA = Cal->Cal_SA(Template,T_SA);
	MEAN_T_SA = SUM_T_SA/tn_vertex;

	// NN smoothing
	SurfaceUtil::smoothing(MNI_GetMesh, 2);
	
	// principal curvature
	int n_vertex = MNI_GetMesh->nVertex();
	int n_mesh = MNI_GetMesh->nFace();

	float *cmin = new float[n_vertex];
	float *cmax = new float[n_vertex];
	SurfaceUtil::curvature(MNI_GetMesh, cmin, cmax, NULL, NULL);	// if principal directions are unnecessary
	
	// Surface Area
	float* SA = new float[n_vertex];
	float SUM_SA = 0.0f;
	float MEAN_SA = 0.0f;
	SUM_SA = Cal->Cal_SA(MNI_GetMesh,SA);
	MEAN_SA = SUM_SA/n_mesh;  

	// calculate shape index
	float* SI = new float[n_vertex];
	Cal->Cal_ShapeIndex(cmin, cmax, SI, n_vertex);
	
	// make basic shape histogram
	int STEP = StepSize;
	int n_BasicShape =  pow(static_cast<double>(2),STEP) + 1;
  	float step_size = 1/pow(static_cast<double>(2),(STEP-1));
  
  	std::vector<float> Basis_X(n_BasicShape,0);
  	for(int i=0; i<n_BasicShape; i++)
  	{
  		Basis_X[i] = (-1) + i*step_size;
  	}
  
	// set up, geodesic distance
	//Geodesic geo(MNI_GetMesh);
	double dmax;
	//const double *dist = geo.dist();
	std::vector<float> Testing_SI;
	
	float* EMD = new float[n_vertex];
	float* temp_EMD_Ldegree = new float[n_BasicShape];
	float* temp_EMD_Hdegree = new float[5];	


	// Initialize EMD 
	for(int i=0;i<n_vertex; i++)
	{
		EMD[i] = 0;
	}
	

	//set kernel size
	if(scaleLSA == false && scaleGSA == false)
	{
		dmax = KernelSize;
		//std::cout << "dmax:" << dmax << std::endl;
	}
	if(scaleGSA == true)
	{
		dmax = KernelSize*(SUM_SA/SUM_T_SA);
		//std::cout << "dmax:"<< dmax << std::endl;
	}
	
	// lookup table
	std::vector< std::vector<int> > FACE_LOOKUP;
	for(int i=0; i<n_vertex; i++)
	{
		std::vector<int> face_list;
		FACE_LOOKUP.push_back(face_list);
	}
	for(int i=0; i<n_mesh; i++)	
	{
		const int *f = MNI_GetMesh->face(i)->list();
		for (int j=0; j<3; j++)
			FACE_LOOKUP[f[j]].push_back(i);
	}
	
	//int tmp_current_point;
	int _N_nbr=0;
	int START_POINT;
	std::vector<int> temp_nbr;
	std::vector<int> temp_pre_nbr;
	std::vector< std::vector<int> > INPUT_TRI;
	std::vector< std::vector<float> > INPUT_COORD;
	std::vector<float> Geo_Dist;
	
	//int maxiter = 100;


	int _nDivide = nDivide;
	int _nPART = nPart;
	int sVertex = 0;
	int eVertex = 0;

	float temp1 = (float)(_nPART-1)/_nDivide;
        float temp2 = (float)(_nPART)/_nDivide;

	sVertex = (int)(n_vertex*temp1);
	eVertex = (int)(n_vertex*temp2);

	if( _nPART >_nDivide)
	{
		cout << "nPart is lager than nDivide, SET  1<= nPart <= nDivide" << endl;
		exit(1);
	
	}


	//clock_t begin = clock();
	for(int i=sVertex; i<eVertex; i++)
	//for(int i=0; i<maxiter; i++)
	{	
		
	  /////////////////// Geodesic diatance -- by Sun Hyung //////////////////////////
	  // 1.Search neighbour for candiated geodesic distaance
	  temp_nbr = Cal->Find_Candidate_nbr(MNI_GetMesh, i, 7);
	  temp_pre_nbr = Cal->Find_Candidate_nbr(MNI_GetMesh, i, 6); 
	 
          		  
	  // 2. Extract patch
	  _N_nbr = MNI_GetMesh->vertex(i)->nNeighbor();
	  //std::cout << "elapsed_secs:" << double (clock() - begin)/ CLOCKS_PER_SEC<< std::endl;
	  //INPUT_TRI = Cal->Extract_Patch_ReturnInputTri(i, temp_nbr, temp_pre_nbr, MNI_GetMesh,n_vertex,n_mesh, _N_nbr);
	  INPUT_TRI = Cal->Extract_Patch_ReturnInputTri(i, temp_nbr, temp_pre_nbr, MNI_GetMesh,n_vertex,n_mesh, _N_nbr, FACE_LOOKUP);
	  //std::cout << "elapsed_secs:" << double (clock() - begin)/ CLOCKS_PER_SEC<< std::endl;
	  INPUT_COORD = Cal->Extract_Patch_ReturnInputCoord(temp_nbr, MNI_GetMesh);
	  START_POINT = Cal->Extract_Patch_ReturnStartPoint(i, temp_nbr);
	  
	  // 3. Calculate geodesic distance
	  //Geo_Dist.clear();
	  Geo_Dist = Cal->Cal_GeoDist(INPUT_TRI, INPUT_COORD, START_POINT);
	  	
	  //return true;
	  
	  
		////////////////////////////////////////////////////////////////////////////////	 
	 							
		if(scaleLSA == true)
		{
			dmax = KernelSize*(SA[i]/T_SA[i]);
			//std::cout << "dmax:"<< dmax << std::endl;
		}
		if(scaleGSA == true &&  scaleLSA== true)
		{
			dmax = KernelSize * (SUM_SA/SUM_T_SA) * (SA[i]/T_SA[i]);
		}
					
		Testing_SI.clear();
		// max distance: compute geodesic distances within a specific range
		//geo.perform_front_propagation(i, dmax);
		//for (int j = 0; j < n_vertex; j++)
		for (int j = 0; j < Geo_Dist.size(); j++)
		{
			//if (dist[j] <= dmax && j != i) //--> Ilwoo's Distance
			if (Geo_Dist[j] <= dmax && j != i)	//-->Sun's Distance
				{
					//Testing_SI.push_back(SI[j]+(step_size/2));
					Testing_SI.push_back(SI[temp_nbr[j]]+(step_size/2));
				}
		}
		
		
		int Test_N = Testing_SI.size();
		//cout<< "Vertex_ID: " << i << endl;
		std::vector<float> SI_HIST=Cal->Make_Histogram(&Testing_SI[0], step_size, Test_N);
	  	//std::cout << "elapsed_secs:" << double (clock() - begin) / CLOCKS_PER_SEC<< std::endl;
		
		if (STEP<3)
		{
			if(Test_N ==0){
				EMD[i]=0;
			}else{
				Cal->CAL_EMD(n_BasicShape,Test_N, Basis_X, SI_HIST, temp_EMD_Ldegree );
				EMD[i] = 2*Cal->Minimum(temp_EMD_Ldegree, n_BasicShape);
			}
		}else{
			if(Test_N==0){
				EMD[i]=0;
			}else{
				Cal->GoldenSection(n_BasicShape,Test_N, Basis_X, SI_HIST, temp_EMD_Hdegree, STEP);
				EMD[i] = 2*Cal->Minimum(temp_EMD_Hdegree, 5);
			}
		}
	  	//std::cout << "elapsed_secs final :" << double (clock() - begin) / CLOCKS_PER_SEC << std::endl;
		
	}
	//clock_t end = clock();
	//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//std::cout << "elapsed_secs:" << elapsed_secs << std::endl;
	//std::cout << "expected full run time:" << elapsed_secs * n_vertex /3600 / maxiter << std::endl;
	
	if(format == "ASCII")
	{
		Cal->Save_Ascii(EMD,output.c_str(),n_vertex);
	}
	if(format == "KWM")
	{
		Cal->Save_KWMValue(EMD,output.c_str(),n_vertex);
	} 
	
	delete Cal; 
	delete MNI_GetMesh; 
	delete Template; 
	delete [] cmin; 
	delete [] cmax; 
	delete [] SI; 
	delete [] EMD; 
	delete [] temp_EMD_Ldegree; 
	delete [] temp_EMD_Hdegree; 
	delete [] SA; 	
	delete [] T_SA; 
	
	return true;
}

