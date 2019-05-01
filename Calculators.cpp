#include "Calculators.h"
#include "transportSimplex.h"

using namespace t_simplex;
using namespace std;


CCalculators::CCalculators()
{
	
	
}

CCalculators::~CCalculators()
{


}

double DISTANCE(float A,float B) 
{
	//Formula for approximating the distance between a factory and a warehouse;
	return abs(A-B);
};



bool CCalculators::Cal_ShapeIndex(float* cmin, float* cmax, float* SI, int n_vertex)
{
	
	for(int i=0; i<n_vertex; i++)
	{
		if(cmax[i] == cmin[i])
    		{
    			if(cmin[i]>0 && cmax[i]>0)
      			{
      				SI[i] = 1;
      			}
      			if(cmin[i]<0 && cmax[i]<0)
      			{
      				SI[i] = -1;
      			}
			if(cmin[i]==0 && cmax[i]==0)
			{
				SI[i] = 1;
			}
     		}else{
     			SI[i] = (2/PI)*atan( (cmax[i]+cmin[i])/(cmax[i]-cmin[i]) );
		}
	}

	for(int i=0; i<n_vertex; i++)
	{
		if(SI[i]<-1){
			SI[i]=-1;
		}
		if(SI[i]>1){
			SI[i]=1;
		}
		
	}
	return true;
		
}

vector<float> CCalculators::Make_Histogram(float* In, float bucket_size, float TEST_N)
{
	vector<float> out_hist;
	
	
	//initial shift, range 0 to 2
	for(int i=0; i<TEST_N; i++){
		In[i] = In[i] + 1;
	}
	
	int number_of_buckets = (int)ceil(1 / bucket_size)*2; 
	vector<int> histogram(number_of_buckets+1);
	
	for(int i=0; i<TEST_N; i++){
	    int bucket = (int)(In[i] / bucket_size);
	    histogram[bucket] += 1;
	}
	
	for (int i=0; i<number_of_buckets+1; i++)
	{
		out_hist.push_back(histogram[i]);
	}
		
	return out_hist;	
}


float CCalculators::Minimum(float* A, int size)
{ 
	float temp = 0.5f;
	for(int i=0;i<size; i++){
		if(A[i]<temp){
			temp = A[i];
		}
	}	
		
	return temp;
}
	
bool CCalculators::Save_Ascii(float* mWrite,const char* filename, int n_vertex)
{
	FILE* writefile;
  
	writefile = fopen( filename, "w" );
	
	for(int i=0; i<n_vertex; i++){
		fprintf( writefile, "%f\n", mWrite[i]);
	}
			
	fclose(writefile);
	return true;
}		

bool CCalculators::Save_KWMValue(float* mWrite,const char* filename, int n_vertex)
{
	FILE* writefile;
  
	writefile = fopen( filename, "w" );
	
	fprintf(writefile, "NUMBER_OF_POINTS=");
	fprintf(writefile, "%i \n", n_vertex);
	fprintf(writefile, "DIMENSION=1 \n"); 
	fprintf(writefile, "TYPE=Scalar \n");
	
	for(int i=0; i<n_vertex; i++){
		fprintf( writefile, "%f\n", mWrite[i]);
	}
			
	fclose(writefile);
	return true;
}		


bool CCalculators::CAL_EMD(int n_BasicShape,int Test_N, vector<float> Basis_X, vector<float> SI_HIST, float* temp_EMD )
{
		
	for(int j=0;j<n_BasicShape;j++)
		{	
  		vector<nBin> IDEAL_SI_HIST(n_BasicShape, nBin(n_BasicShape,0));	                                 
			for(int i=0;i<n_BasicShape;i++)
			{	
				IDEAL_SI_HIST[i][i]=Test_N;	
			}
						
			TsSignature<string> srcSig(n_BasicShape, Basis_X, SI_HIST);
			TsSignature<string> snkSig(n_BasicShape, Basis_X, IDEAL_SI_HIST[j]);
		
			TsFlow flow[n_BasicShape];
			int flowVars = 0;
			double result = transportSimplex(&srcSig, &snkSig, DISTANCE, flow, &flowVars);
			
			//cout << "Total cost: " << result << endl;
			//cout << "Flows:" << endl;
		
			float temp1=0.0f;
			float temp2=0.0f;
			for (int k = 0; k < flowVars; k++)
			{ 
				//cout << basis_x[flow[i].from] << " to " << basis_x[flow[i].to] << " : " << flow[i].amount << endl;
				temp1 = temp1 + abs(Basis_X[flow[k].from]-Basis_X[flow[k].to])*flow[k].amount;
				temp2 = temp2 +flow[k].amount;	
			}
			
			temp_EMD[j] = temp1/temp2;
		}
	return true;
}


bool CCalculators::GoldenSection(int n_BasicShape,int Test_N, vector<float> Basis_X, vector<float> SI_HIST, float* temp_EMD, int STEP)
{ 
	typedef vector<float> IN;
	vector<IN> INIT_TEST_BIN(5,IN(n_BasicShape,0));

	vector<float> SLOPE(4,0);
	vector<nBin> IDEAL_SI_HIST(n_BasicShape, nBin(n_BasicShape,0));	                                 
	for(int i=0;i<n_BasicShape;i++)
	{	
		IDEAL_SI_HIST[i][i]=Test_N;	
	}
	
	// Initial Step
	for(int i=0; i<5; i++)
	{
		INIT_TEST_BIN[i] = IDEAL_SI_HIST[i*((n_BasicShape-1)/4)];
		
		TsSignature<string> srcSig(n_BasicShape, Basis_X, SI_HIST);
		TsSignature<string> snkSig(n_BasicShape, Basis_X, INIT_TEST_BIN[i]);
		
		TsFlow flow[n_BasicShape];
		int flowVars = 0;
		double result = transportSimplex(&srcSig, &snkSig, DISTANCE, flow, &flowVars);
		
		float temp1=0.0f;
		float temp2=0.0f;
		for (int k = 0; k < flowVars; k++)
		{ 
			temp1 = temp1 + abs(Basis_X[flow[k].from]-Basis_X[flow[k].to])*flow[k].amount;
			temp2 = temp2 +flow[k].amount;	
		}
		
		temp_EMD[i] = temp1/temp2;
	}
	
	for(int i=0;i<4; i++)
	{
		SLOPE[i] = temp_EMD[i+1] - temp_EMD[i];
	}
	int N_BETWEEN =100;
	for(int i=0;i<3; i++)
	{
		if(SLOPE[i] * SLOPE[i+1] <= 0)
			{
				N_BETWEEN = i;
			}
	}
	if(N_BETWEEN ==100)
	{
    		if(SLOPE[0]<0 && SLOPE[1]<0 && SLOPE[2]<0 && SLOPE[3]<0)
    		{
     			N_BETWEEN = 2;
    		}
    		if(SLOPE[0]>0 && SLOPE[1]>0 && SLOPE[2]>0 && SLOPE[3]>0)
    		{
      			N_BETWEEN = 0;
    		}
  	}
  
  	if(N_BETWEEN ==100)
	{
  		cout << "N_BETWEEN = 100 !!!!!" << endl;	
  	}
  
	// Loop Step
  	int tmp_START_POINT = 0;
  	for(int j=0;j<(STEP-2);j++)
  	{
  		int NUM_BIN = ((n_BasicShape-1)/pow(static_cast<double>(2),j+1));
      		tmp_START_POINT = tmp_START_POINT + (NUM_BIN/2)*(N_BETWEEN);
     
      		for(int i=0; i<5; i++)
      		{
          		INIT_TEST_BIN[i] = IDEAL_SI_HIST[tmp_START_POINT+(i)*(NUM_BIN/4)];
      		}
      		N_BETWEEN = OPT_EMD(temp_EMD,N_BETWEEN, Basis_X,INIT_TEST_BIN[1],INIT_TEST_BIN[3],SI_HIST,n_BasicShape);
  	}

	return true;
}

int CCalculators::OPT_EMD(float* temp_EMD, int START, vector<float> BASE, vector<float> HIST_A2, vector<float> HIST_A4, vector<float> HIST_B, int n_BasicShape)
{
	vector<float> SLOPE(4,0);
	
	float* mEMD = new float[5];	
	mEMD[0] = temp_EMD[START];
	mEMD[2] = temp_EMD[START+1];
	mEMD[4] = temp_EMD[START+2];
	
  TsSignature<string> srcSig(n_BasicShape, BASE, HIST_B);
	TsSignature<string> snkSigA2(n_BasicShape, BASE, HIST_A2);
		
	TsFlow flow[n_BasicShape];
	int flowVars = 0;
	double result = transportSimplex(&srcSig, &snkSigA2, DISTANCE, flow, &flowVars);
		
	float temp1=0.0f;
	float temp2=0.0f;
	for (int k = 0; k < flowVars; k++)
	{ 
		temp1 = temp1 + abs(BASE[flow[k].from]-BASE[flow[k].to])*flow[k].amount;
		temp2 = temp2 +flow[k].amount;	
	}
		
	mEMD[1] = temp1/temp2;
	

	TsSignature<string> snkSigA4(n_BasicShape, BASE, HIST_A4);
	flowVars = 0;
	result = transportSimplex(&srcSig, &snkSigA4, DISTANCE, flow, &flowVars);
		
	temp1=0.0f;
	temp2=0.0f;
	for (int k = 0; k < flowVars; k++)
	{ 
		temp1 = temp1 + abs(BASE[flow[k].from]-BASE[flow[k].to])*flow[k].amount;
		temp2 = temp2 +flow[k].amount;	
	}
		
	mEMD[3] = temp1/temp2;	
	
	for(int i=0;i<4; i++)
	{
		SLOPE[i] = mEMD[i+1] - mEMD[i];
	}
	int N_BETWEEN =100;
	for(int i=0;i<3; i++)
	{
		if(SLOPE[i] * SLOPE[i+1] <= 0)
			{
				N_BETWEEN = i;
			}
	}
	
	if(N_BETWEEN ==100)
	{
    if(SLOPE[0]<0 && SLOPE[1]<0 && SLOPE[2]<0 && SLOPE[3]<0)
    {
     	N_BETWEEN = 2;
    }
    if(SLOPE[0]>0 && SLOPE[1]>0 && SLOPE[2]>0 && SLOPE[3]>0)
    {
      N_BETWEEN = 0;
    }
  }
  
  for(int i=0; i<5; i++)
  {
  	temp_EMD[i] = mEMD[i];
  }
  
  delete [] mEMD;
	mEMD=NULL;

	return N_BETWEEN;	
}	

vector<int> CCalculators::Find_Candidate_nbr(Mesh* mesh, int current_point_index, int n_ring)
{
	vector<int> tmp_NBR;
	tmp_NBR.clear();	
	const int n_vertex = mesh->nVertex();
	
	int temp_sum_n_nbr=0;
	int cp = current_point_index;
	
	const int *neighbor = mesh->vertex(cp)->list();
	int n_nbr = mesh->vertex(cp)->nNeighbor();
	//temp_sum_n_nbr += n_nbr;
	temp_sum_n_nbr = n_nbr;
		
	for(int i=0; i<n_nbr;i++)
	{
		tmp_NBR.push_back(neighbor[i]);
	}
	
	//cout << "1st n_nbr:" << n_nbr << endl;
	
	int loop_condition=1;
	int sub_n_nbr;
	const int *tmp_neighbor;
	while(loop_condition<n_ring)
	{
		for(int i=0;i<n_nbr; i++)
		{	
			tmp_neighbor = mesh->vertex(tmp_NBR[i])->list();
			sub_n_nbr = mesh->vertex(tmp_NBR[i])->nNeighbor();
			for(int j=0; j<sub_n_nbr;j++)
			{
				tmp_NBR.push_back(tmp_neighbor[j]);
			}
			temp_sum_n_nbr += sub_n_nbr;
		}
		
		n_nbr = temp_sum_n_nbr;
		loop_condition = loop_condition + 1;
	}
	
	vector<int> NBR(n_nbr,0);
	
	for(int i=0; i<n_nbr; i++)
	{
		NBR[i] = tmp_NBR[i];
	}
	
	return cUNIQUE(NBR);
	
}

vector< vector<int> > CCalculators::Extract_Patch_ReturnInputTri(int Vertex_Orig_ID, vector<int> Uniq_NBR, vector<int> Uniq_NBR_preStep, Mesh* mesh, int n_vertex, int nFace, int _N_nbr, std::vector<std::vector<int> > Face_LookupTable)
{
	
	int N_search_Face = int(Uniq_NBR_preStep.size());
	int N_NBR = int(Uniq_NBR.size());

	int FACE_size = N_search_Face*_N_nbr;
	
	vector<int> FACE;
	
	//Need Optimize///////////////////////////////////////////
	/*int ver_id;
	for(int i=0; i<N_search_Face; i++)
	{
		ver_id = Uniq_NBR_preStep[i];
		for(int j=0;j<nFace;j++)
		{
			const int *f = mesh->face(j)->list();
			if(f[0] == ver_id || f[1] == ver_id || f[2] == ver_id)
			{
				FACE.push_back(j);
			} 
		}
	}*/

	int ver_id;
	for(int i=0; i<N_search_Face; i++)
	{
		ver_id = Uniq_NBR_preStep[i];
		FACE.insert(FACE.end(), Face_LookupTable[ver_id].begin(), Face_LookupTable[ver_id].end());
	}
	///////////////////////////////////////////////////////////
	
	for(int i=0;i<FACE_size;i++)
	{
		if(FACE[i]==0)
		{
			FACE[i]=FACE[0];	
		}
		//cout << FACE[i] << endl;
	} 
	
  sort(FACE.begin(), FACE.end());
  FACE.erase(unique(FACE.begin(), FACE.end()),FACE.end()); 	 
	
	int N_face = FACE.size();

	//vector<int> Patch_Face_DB01(N_face);
	//vector<int> Patch_Face_DB02(N_face);
	//vector<int> Patch_Face_DB03(N_face);
	vector<int> Patch_Face_DB04(N_face);	
	vector<int> Patch_Face_DB05(N_face);
	vector<int> Patch_Face_DB06(N_face);	
	
	for(int i=0; i<N_face; i++)
	{
		int p1, p2, p3;
		p1 = mesh->face(FACE[i])->list()[0];
		p2 = mesh->face(FACE[i])->list()[1];
		p3 = mesh->face(FACE[i])->list()[2];
		//cout << p1 << endl;
		for(int j=0;j<N_NBR;j++)
		{
			if(Uniq_NBR[j]==p1)
				{
					//Patch_Face_DB01[i] = p1;
					Patch_Face_DB04[i] = j;
				}
			if(Uniq_NBR[j]==p2)
			{
				//Patch_Face_DB02[i] = p2;
				Patch_Face_DB05[i] = j;
			}
			if(Uniq_NBR[j]==p3)
			{
				//Patch_Face_DB03[i] = p3;
				Patch_Face_DB06[i] = j;
			}
		}
		//cout << Patch_Face_DB01[i] << ',' << Patch_Face_DB02[i]<< ',' << Patch_Face_DB03[i]<< ',' << Patch_Face_DB04[i]<< ',' << Patch_Face_DB05[i]<< ',' << Patch_Face_DB06[i] << endl;
		
	}

	vector< vector<int> > Input_Tri(N_face);
	//vector< vector<int>* > Input_Tri(N_face);
	
	int tmp_a, tmp_b, tmp_c;
	for(int i=0;i<N_face;i++)
	{
		tmp_a = Patch_Face_DB04[i];
		tmp_b = Patch_Face_DB05[i];
		tmp_c = Patch_Face_DB06[i];
		Input_Tri[i].push_back(tmp_a);
		Input_Tri[i].push_back(tmp_b);
		Input_Tri[i].push_back(tmp_c);
		
		//cout << Input_Tri[i][0] << ',' << Input_Tri[i][1] << ',' << Input_Tri[i][2] << endl;
	}	 

	return Input_Tri;
		
}

vector<int> CCalculators::cUNIQUE(vector<int> Matrix)
{
	int count=0;	
  
  sort(Matrix.begin(), Matrix.end());
  vector<int>::iterator last;
  
  last = unique(Matrix.begin(), Matrix.end()); // 1 2 3 2 1 3 2 2 1
  Matrix.erase(last, Matrix.end());
  
	return Matrix;
}


vector< vector<float> > CCalculators::Extract_Patch_ReturnInputCoord( vector<int> Uniq_NBR, Mesh* mesh )
{
	int N_NBR = int(Uniq_NBR.size());
	vector< vector<float> > Input_Coord(N_NBR);
	
	for(int i=0;i<N_NBR;i++)
	{
		Input_Coord[i].push_back(mesh->vertex(Uniq_NBR[i])->fv()[0]);
		Input_Coord[i].push_back(mesh->vertex(Uniq_NBR[i])->fv()[1]);
		Input_Coord[i].push_back(mesh->vertex(Uniq_NBR[i])->fv()[2]);
		//cout << Input_Coord[i][0] << ',' << Input_Coord[i][1] << ',' << Input_Coord[i][2] << endl;
	}	 
	
	return Input_Coord;
 	
}

int CCalculators::Extract_Patch_ReturnStartPoint(int Vertex_Orig_ID, vector<int> Uniq_NBR)
{
	int Start_Point;
	int nSize = Uniq_NBR.size();
	for(int i=0; i< nSize; i++)
	{
		if(Uniq_NBR[i]==Vertex_Orig_ID)
			{
				Start_Point = i;
			}
	}	
	
	return Start_Point;
}	

vector<float> CCalculators::Cal_GeoDist(vector< vector<int> > faces, vector< vector<float> > vertex, int INIT_POINT_ID)
{
	vector<float> GDist;
	
	//n=size(POINTS,2);
	int n = vertex.size();
	//W = ones(n,1);
	vector<int> W(n,1);
	//I = [INIT_POINT_ID];
	int I = INIT_POINT_ID;
	//U = zeros(n,1);
	vector<float> U(n,0);
	//i = [faces(1,:) faces(2,:) faces(3,:) ];
	//j = [faces(2,:) faces(3,:) faces(1,:) ];
	//k = [faces(3,:) faces(1,:) faces(2,:) ];
	int nfaces = faces.size();
	vector<int> xi(nfaces*3);
	vector<int> xj(nfaces*3);
	vector<int> xk(nfaces*3);
	for(int i=0;i<nfaces;i++)
	{
		xi[i]= (faces[i][0]);
		xi[i+(nfaces)] = (faces[i][1]);
		xi[i+(2*nfaces)] = (faces[i][2]);
		
		xj[i]= (faces[i][1]);
		xj[i+(nfaces)] = (faces[i][2]);
		xj[i+(2*nfaces)] = (faces[i][0]);
		
		xk[i]= (int(faces[i][2]));
		xk[i+(nfaces)] = (faces[i][0]);
		xk[i+(2*nfaces)] = (faces[i][1]);
	}	
	
	//x  = vertex(:,i);
	//x1 = vertex(:,j) - x;
	//x2 = vertex(:,k) - x;
	
	int NSIZE = nfaces*3;
	//cout << NSIZE << endl;
	
	TEMP* x = new TEMP[nfaces*3];	
	TEMP* x1 = new TEMP[nfaces*3];	
	TEMP* x2 = new TEMP[nfaces*3];	
	
	
	for(int i=0; i<nfaces*3;i++)
	{
		float xia = x[i].a = vertex[xi[i]][0];
		float xib = x[i].b = vertex[xi[i]][1];
		float xic = x[i].c = vertex[xi[i]][2];
	
		x1[i].a = vertex[xj[i]][0] - xia;
		x1[i].b = vertex[xj[i]][1] - xib;
		x1[i].c = vertex[xj[i]][2] - xic;
		
		x2[i].a = vertex[xk[i]][0] - xia;
		x2[i].b = vertex[xk[i]][1] - xib;
		x2[i].c = vertex[xk[i]][2] - xic;
		
		//cout << x2[i].a << ',' << x2[i].b << ',' << x2[i].c << endl;
	}	
		
	
	float terminate = 0.01;
	float condition = 1;
	
	int Nsize = nfaces*3;	
	
	  
	//int iter=0;
		
    vector<float> C11(Nsize);
		vector<float> C12(Nsize);
		vector<float> C21(Nsize);
		vector<float> C22(Nsize);//_Nsize = 162
		
    vector<float> d(Nsize);     
    vector<int> J; 
    vector<float> U1(vertex.size());
    vector<float> d1(Nsize); 
    vector<float> d2(Nsize); 
 		
	while(condition >=	terminate)
	{
		//uj = U(j);
    //uk = U(k);
	  //u = [R(uj); R(uk)];
    //w = R( W(i) );
    //C = [R(dotp(x1,x1)) R(dotp(x1,x2)); R(dotp(x2,x1)) R(dotp(x2,x2))];
		
    //S = Inv(C);
    //Inv1 = @(M,d)[M(2,2,:)./d -M(1,2,:)./d; -M(2,1,:)./d M(1,1,:)./d];
		//Inv  = @(M)Inv1(M, M(1,1,:).*M(2,2,:) - M(1,2,:).*M(2,1,:));
		//a = sum(sum(S));
 		//b = dotp( sum(S,2), u );
 		//Mult = @(M,u)[M(1,1,:).*u(1,1,:) + M(1,2,:).*u(2,1,:);  M(2,1,:).*u(1,1,:) + M(2,2,:).*u(2,1,:)];
 		//c = dotp( Mult(S,u), u ) - w.^2;
    //delta = max( b.^2 - a.*c, 0);
    //d = (b + sqrt(delta) )./a;
    //alpha = Mult( S, u - repmat(d, 2, 1) );
    //Mult = @(M,u)[M(1,1,:).*u(1,1,:) + M(1,2,:).*u(2,1,:);  M(2,1,:).*u(1,1,:) + M(2,2,:).*u(2,1,:)];
    //J = find( alpha(1,1,:)>0 | alpha(2,1,:)>0 );
    //d1 = sqrt(dotp(x1,x1)); 
    //d1 = d1(:).*w(:) + uj(:);
    //d2 = sqrt(dotp(x2,x2)); 
    //d2 = d2(:).*w(:) + uk(:);
    //d = d(:);
    //d(J) = min(d1(J), d2(J));
    
		C11=dotp(x1,x1,Nsize); 
		C12=dotp(x1,x2,Nsize); 
		C21=dotp(x2,x1,Nsize); 
		C22=dotp(x2,x2,Nsize);
		
		for(int i=0; i<Nsize;i++)
		{
			float u11 = U[xj[i]];
			float u12 = U[xk[i]];
			
			float w = W[xi[i]];
			
			float det = C11[i]*C22[i] - C12[i]*C21[i];
			
			float S11i =		   C22[i]/det; 
			float S12i =  (-1)*C12[i]/det;
			float S21i =  (-1)*C21[i]/det; 
			float S22i = 		   C11[i]/det;
			float a = S11i + S12i + S21i + S22i;
 			float b = (S11i+S12i)*u11 + (S21i+S22i)*u12;
 			float mult11 = S11i*u11 + S12i*u12;
 			float mult12 = S21i*u11 + S22i*u12; 
 			float c = (mult11*u11 + mult12*u12) - w*w;
			//cout << det << endl;
			float delta = 0;
			
			float bsquare = b*b;
    	if( (bsquare - a*c) > 0)
    	{
    		delta = bsquare - a*c;
    	}	 
    	
    	d[i] = (b + sqrt(delta))/a;
    	//cout << d[i] << endl; 
    	
    	float repmat11 = u11-d[i];
    	float repmat12 = u12-d[i];
    	float alpha11 = S11i*repmat11 + S12i*repmat12;
    	float alpha12 = S21i*repmat11 + S22i*repmat12;
    	//cout<< alpha11[i] << endl;
    	
    	if(alpha11>0 || alpha12>0)
    		{
    			J.push_back(i);
    		}
    		
    	d1[i] = sqrt(dotp(x1, x1,Nsize)[i]);
    	d1[i] = d1[i]*w + u11;
    	d2[i] = sqrt(dotp(x2, x2,Nsize)[i]);
    	d2[i] = d2[i]*w + u12;
    	
    }
    
    for(int i=0; i<J.size(); i++)
    {
    	d[J[i]] = min(d1[J[i]], d2[J[i]]);
    	//cout<< d[i] << endl;
    }
    
    //U1 = accumarray(i', d, [n 1], @min);  U1(U1==0) = Inf;
    //U1(I) = 0;
    U1 = accumarray(xi,d);
    U1[I]=0;
   
    //condition = mean(abs(U1 - U),1);
    //U = U1;
    float mean_tmp =0.0f;
    for(int i=0;i<U1.size();i++)
    {
    	mean_tmp = mean_tmp + abs(U1[i]-U[i]);
    	U[i] = U1[i];
    }
    condition = mean_tmp/U1.size();
    
		//cout << "iter:" << iter << endl;
	}
	
	
	for(int i=0;i<U.size();i++)
	{
		GDist.push_back(U[i]);
	}
	
	delete [] x;
	x=NULL;	
	delete [] x1;
	x1=NULL;
	delete [] x2;
	x2=NULL;
	
	return GDist;
}	

vector<float> CCalculators::accumarray(vector<int> Index, vector<float> Val)
{
	vector<float> Out;Out.clear();
	
	for(int j=0; j<Index.size(); j++)
	{
		float temp = 100000.0f;
		for(int i=0; i<Index.size(); i++)
		{
			if( j == Index[i])
			{
				temp = min(Val[i], temp);
			}
		}
		if(temp < 100000.0f)
		{
			Out.push_back(temp);
		}
	}	
	
	return Out;
}	

vector<float> CCalculators::dotp(TEMP* x1, TEMP* x2, int size)
{
	vector<float> dot_out(size);
	//dotp = @(u,v)sum(u.*v,1);
	TEMP* tmp = new TEMP[size];
	for(int i=0;i<size;i++)
	{
		tmp[i].a = x1[i].a * x2[i].a;
		tmp[i].b = x1[i].b * x2[i].b;
		tmp[i].c = x1[i].c * x2[i].c;
		dot_out[i] = tmp[i].a + tmp[i].b + tmp[i].c;
		//cout << dot_out[i] << endl;
	}

	delete [] tmp;
	tmp = NULL;

	return dot_out;
	
}
	
double CCalculators::Cal_SA(Mesh* mesh, float* SA)
{
	double SUM_SA;
	const int n_vertex = mesh->nVertex();
	vector<float> n1;
	vector<float> n2;
	vector<float> pt;
		
	float edge1,edge2,edge3;
			
	for (int i = 0; i < n_vertex; i++)
	{
		pt.clear();
		const int *neighbor = mesh->vertex(i)->list();
		const int n_nbr = mesh->vertex(i)->nNeighbor();
		float temp =0.0f;
		//cout << n_nbr << endl;	
	
		for(int j=0; j<3; j++)
		{	
			pt.push_back(mesh->vertex(i)->fv()[j]);
			//cout << pt[j] << endl;
		}	
		
		for (int j = 0; j <  n_nbr-1; j++)
		{
			n1.clear();
			n2.clear();
			for(int k=0;k<3;k++)
			{
				n1.push_back(mesh->vertex(neighbor[j])->fv()[k]);
				n2.push_back(mesh->vertex(neighbor[j+1])->fv()[k]);
				//cout << n1[k] << endl;
			}
			
			edge1 = Distance_TwoPoint(pt,n1);
			edge2 = Distance_TwoPoint(pt,n2);
			edge3 = Distance_TwoPoint(n1,n2);
			temp = temp + Cal_TRI(edge1, edge2, edge3);
		}
		SA[i] = temp/n_nbr;
		SUM_SA += SA[i];
	}
	
	return SUM_SA;	
}

float CCalculators::Distance_TwoPoint(vector<float> pt1, vector<float> pt2)
{
	float Dist=0.0;
	Dist = sqrt( (pt1[0]-pt2[0])*(pt1[0]-pt2[0]) + (pt1[1]-pt2[1])*(pt1[1]-pt2[1]) + (pt1[2]-pt2[2])*(pt1[2]-pt2[2]) );
	
	return Dist;
}

float CCalculators::Cal_TRI(float edge1,float edge2,float edge3)
{
	float S = 0.0f;
	float m = 0.0f;
	m= (edge1+edge2+edge3)/2;
  S = sqrt(m*(m-edge1)+m*(m-edge2) + m*(m-edge3));
	return S;
}	         

        
