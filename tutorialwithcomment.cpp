
/*


VIEW FACTOR CALCULATION
SOURCE:
1) Point Source
2) Multi-dimensional object


METHODS:
-------
1) Area Integration method
2) MOnte Carlo Method

GEOMETRIES:
----------
Select geometry for source and object
1 for circular disk
2 for rectangular plate
3 for tiangular plate
4 for sphere
5 for cylinder

There are 12 tutorials in Optix
All the methods are defined in trace function
*/

#include "random.h"
#include <GLUTDisplay.h>
#include <OptiXMesh.h>
#include <sutil.h>
#include <optixu/optixpp_namespace.h>
#include <optixu/optixu_math_namespace.h>
#include <iostream>
#include <GLUTDisplay.h>
#include <ImageLoader.h>
#include "commonStructs.h"
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <math.h>
#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>
#include <OptiXMesh.h>
#include <Mouse.h>
#include <random>
#include <string>



float total = 0.0;
int icount = -1;
double time_final = GetTickCount();

using namespace optix;
using namespace std;

std::string  m_ptx_path;

std::ofstream myfile;		//To store data of viewfactor, normals, points
std::ofstream myfile3;		//we have "myfile"  output buffers
std::ofstream myfile4;
std::ofstream myfile2;
bool flag = true;           // source sampling- true

int total_points;
int total_triangles;

GeometryGroup geometrygroup;	//a feature of optix for geometry instance creation

//float eyex[1250] = { 0 };
//float eyey[1250] = { 0 };
//float eyez[1250] = { 0 };

float* eyex = new float[100000000];		// The eye is center of triangles. To Find it,we take average of 3 corner points.
float* eyey = new float[100000000];		//They are seperately accumulated in x,y,z float nine dimension arrays .
float* eyez = new float[100000000];		//

float* normx = new float[100000000];	//We have normals of triangles which are perpendicular to surface. They are at the eye
float* normy = new float[100000000];	//and aim other obstacles in 3D space domain.
float* normz = new float[100000000];

int count_s = 0;


//variables for calculating view factor
float enr[7056] = { 0 };		// inital energy
float total_p[7056] = { 0 };	//total primitives
float hitting_p[7056] = { 0 };	//hitting primitives
float vol[7056] = { 0 };		//volume?   //not used
int energy_size = 324;			//  //not used
float initial_cond[9830] = { 0 };  //not used


int Faces[7056][3] = { 0 }; //9830
float3 centers[7056]; //324 //720 for sphere  //7056
float3 normal[7056]; //69984

//random generators
std::random_device rd;	//Random number generator
std::minstd_rand gen(rd());	//
std::uniform_int_distribution<> dis(0, 720);	//uniform distribution between 0 and 720


//function definition

float area(double a[], double b[], double c[], int length);	//

void material(Material box_matl, optix::Context m_context);		// Create geometry, material type can be chosen in Optix
void material_parabola(Material box_matl, optix::Context m_context);	//All depend on closest hit or any hit or miss hit
void material_object(Material box_matl, optix::Context m_context);

void CreateBox(optix::Context m_context, Geometry box, std::string box_ptx);

void CreateTriangle(optix::Context m_context, Geometry triangle, std::string ptx_path);
//void CreateTriangle2(optix::Context m_context, Geometry triangle, std::string ptx_path);

void Face_center(void);	//Face center function gives the required triangle primitive`s centers







static float rand_range(float min, float max)
{
	return min + (max - min) * (float)rand() / (float)RAND_MAX;		//return any average random number
}






class Tutorial : public SampleScene		//Define scene basic, area sizes, triangulars, centers, normals, hitting&total primitives
{
public:
	Tutorial(int tutnum, const std::string& texture_path)
		: SampleScene(), m_tutnum(tutnum), m_width(1000u), m_height(1000u), texture_path(texture_path)
	{}

	// From SampleScene
	void   initScene(InitialCameraData& camera_data);
	void   trace(const RayGenCameraData& camera_data);

	void   setDimensions(const unsigned int w, const unsigned int h) { m_width = w; m_height = h; }
	Buffer getOutputBuffer();
	Buffer getOutputBufferr();


	void genRndSeeds(unsigned int width, unsigned int height);

private:
	std::string texpath(const std::string& base);
	void createGeometry();

	void InsertModel(const std::string& name, GeometryGroup ggroup, Material mat, Program mesh_intersect, const optix::Matrix4x4 matx, bool suppressErrors);

	unsigned int m_tutnum;
	unsigned int m_width;
	unsigned int m_height;
	std::string   texture_path;

	float3 radiuss;



};


void Tutorial::initScene(InitialCameraData& camera_data)
{
	myfile2.open("VF.csv");		//View Factor file with coordinates
	// set up path to ptx file associated with tutorial number
	std::stringstream ss;
	ss << "tutorial" << m_tutnum << ".cu";
	m_ptx_path = ptxpath("tutorial", ss.str());	//Making string to load file name to a variable

	// context 
	m_context->setRayTypeCount(2);		//
	m_context->setEntryPointCount(2);	//
	m_context->setStackSize(4640); //4640 this is arbitrary

	m_context["max_depth"]->setInt(100); //100
	m_context["radiance_ray_type"]->setUint(0);//0
	m_context["shadow_ray_type"]->setUint(1); // 1
	m_context["visual_ray_type"]->setUint(2); // 1
	m_context["frame_number"]->setUint(0u);
	m_context["scene_epsilon"]->setFloat(1.e-3f);
	m_context["importance_cutoff"]->setFloat(0.01f);
	m_context["ambient_light_color"]->setFloat(0.1f, 0.1f, 0.1f); //red,green,blue 
	m_context->setPrintEnabled(1);//
	m_context->setExceptionEnabled(RT_EXCEPTION_ALL, true);






	m_context["output_buffer"]->set(createOutputBuffer(RT_FORMAT_UNSIGNED_BYTE4, m_width, m_height));

	m_context["output_bufferl"]->set(createOutputBuffer(RT_FORMAT_UNSIGNED_BYTE4, m_width, m_height));// not used




	// Ray gen program
	std::string camera_name, camera_name1;
	
    {
		camera_name = "pinhole_camera";  //for source,calculation
		camera_name1 = "random_camera"; // for visualization
	}

	Program ray_gen_program_source = m_context->createProgramFromPTXFile(m_ptx_path, camera_name1);//Optix library provides us ray generation toolbox.
	m_context->setRayGenerationProgram(0, ray_gen_program_source);	//camera represents ray source center.,type of rays calculation or visualization

	Program ray_gen_program = m_context->createProgramFromPTXFile(m_ptx_path, camera_name);

	m_context->setRayGenerationProgram(1, ray_gen_program);


	// Exception / miss programs
	Program exception_program = m_context->createProgramFromPTXFile(m_ptx_path, "exception");
	m_context->setExceptionProgram(0, exception_program);	//if there is no obstacle in ray direction, it is miss otherwise there is error



	m_context["bad_color"]->setFloat(0.0f, 50.0f, 0.0f);

	std::string miss_name;
	if (m_tutnum >= 5)				//tutorials btw 0 and 11 shapes the view.For instance, we have a miss until 4.tutorial, called it "miss"
		miss_name = "envmap_miss";	//Otherwise it is environmental blue?
	else
		miss_name = "miss";

	m_context->setMissProgram(0, m_context->createProgramFromPTXFile(m_ptx_path, miss_name));
	const float3 default_color = make_float3(0.0f, 0.55f, 0.85f);
	m_context["envmap"]->setTextureSampler(loadTexture(m_context, texpath("CedarCity.hdr"), default_color));



	m_context["bg_color"]->setFloat(make_float3(255.0f, 255.0f, 255.0f)); //white background


	// prims
	BasicLight lights[] = {
		{ make_float3(-10.0f, 15.0f, -20.0f), make_float3(1.0f, 1.0f, 1.0f), 1 },                //which color?
		// { make_float3(10.0f, 9.0f, 12.0f), make_float3(1.0f, 0.0f, 0.0f), 1 },
		// { make_float3(5.0f, 9.0f, 12.0f), make_float3(0.0f, 1.0f, 0.0f), 1 }
	};





	Buffer light_buffer = m_context->createBuffer(RT_BUFFER_INPUT);
	light_buffer->setFormat(RT_FORMAT_USER);
	light_buffer->setElementSize(sizeof(BasicLight));
	light_buffer->setSize(sizeof(lights) / sizeof(lights[0]));	//ligts over total lights
	memcpy(light_buffer->map(), lights, sizeof(lights));  //lights(as much as length of sizeof(lights)) are copied to light_buffer->map() 
	light_buffer->unmap();  // clear buffer and set new variable to "lights"

	m_context["lights"]->set(light_buffer);


	// Buffer points = m_context["points"]->getBuffer();


	Buffer Energy_init = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);	// After this point, some variable setups, definitions are similar.Therefore,
	Energy_init->setFormat(RT_FORMAT_FLOAT);								//Set Format, size, push and pull buffers
	// float valsf[1120] = { 0.0 };
	Energy_init->setSize(sizeof(enr));											//enr is Initial energy					
	memcpy(Energy_init->map(), enr, sizeof(enr));
	Energy_init->unmap();
	m_context["Energy_init"]->set(Energy_init);

	Buffer Extinction = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);		// Extinction
	Extinction->setFormat(RT_FORMAT_FLOAT);
	// float valsf[1120] = { 0.0 };
	Extinction->setSize(sizeof(enr));
	memcpy(Extinction->map(), enr, sizeof(enr));
	Extinction->unmap();
	m_context["Extinction"]->set(Extinction);

	Buffer hitting_prim = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);		//hitting primitive // buffer for counting number of rays hitting primitives
	hitting_prim->setFormat(RT_FORMAT_INT);
	// float valsf[1120] = { 0.0 };
	hitting_prim->setSize(sizeof(hitting_p));
	memcpy(hitting_prim->map(), hitting_p, sizeof(hitting_p));
	hitting_prim->unmap();
	m_context["hitting_prim"]->set(hitting_prim);

	Buffer total_prim = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);		// total_primitive  // buffer for counting number of rays emmiting from primitives
	total_prim->setFormat(RT_FORMAT_INT);
	// float valsf[1120] = { 0.0 };
	total_prim->setSize(sizeof(total_p));
	memcpy(total_prim->map(), total_p, sizeof(total_p));
	total_prim->unmap();
	m_context["total_prim"]->set(total_prim);





	Buffer normals = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);			//Normals
	normals->setFormat(RT_FORMAT_FLOAT);

	normals->setSize(sizeof(normal));
	memcpy(normals->map(), normal, sizeof(normal));
	normals->unmap();
	m_context["normals"]->set(normals);


	Buffer centersb = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);   // centers buffer  //location of emission from primitves
	centersb->setFormat(RT_FORMAT_FLOAT);
	centersb->setSize(sizeof(centers));
	memcpy(centersb->map(), centers, sizeof(centers));
	centersb->unmap();
	m_context["centers"]->set(centersb);


	Buffer b = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);			// total number of rays generated
	b->setFormat(RT_FORMAT_INT);
	int vals[] = { 0 };
	b->setSize(sizeof(vals));
	memcpy(b->map(), vals, sizeof(vals));
	b->unmap();
	int* buffers = static_cast<int*>(b->map());
	//std::cout << buffers[0] << std::endl; // double check to see if we inserted the correct value, it comes out correct
	b->unmap();
	m_context["out"]->set(b);

	Buffer b2 = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);
	b2->setFormat(RT_FORMAT_INT);
	int vals2[] = { 0 };
	b2->setSize(sizeof(vals2));
	// memcpy(b2->map(), vals2, sizeof(vals2));
	// b2->unmap();
	int* buffers2 = static_cast<int*>(b2->map());
	//std::cout << buffers2[0] << std::endl; // double check to see if we inserted the correct value, it comes out correct
	b2->unmap();
	m_context["out2"]->set(b2);


	Buffer b3 = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);  // total number of rays hitted
	b3->setFormat(RT_FORMAT_INT);

	b3->setSize(sizeof(vals2));
	// memcpy(b3->map(), vals2, sizeof(vals2));
	//b3->unmap();


	buffers2 = static_cast<int*>(b3->map());
	//std::cout << buffers2[0] << std::endl; // double check to see if we inserted the correct value, it comes out correct
	b3->unmap();
	m_context["out3"]->set(b3);



	Buffer count = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);		//count  // number of simulation (trials)
	count->setFormat(RT_FORMAT_INT);
	count->setSize(sizeof(vals2));
	memcpy(count->map(), vals2, sizeof(vals2));
	count->unmap();

	//std::cout << buffers2[0] << std::endl; // double check to see if we inserted the correct value, it comes out correct

	m_context["countbuffer"]->set(count);


	// Set up camera

	float3 eyee = make_float3(-18.0f, 14.0f, 18.0f);	//to create a vector, we determine camera  as a start point
	float3 lookkat = make_float3(0.0f, 6.0f, 0.0f);		//target is the lookat point. This is for initialization.
								//Because we need to specify a starting point.Target can be anywhere.
	camera_data = InitialCameraData(eyee,   // eye 
		lookkat,  // lookat
		make_float3(0.0f, 1.0f, 0.0f),//  up
		60.0f);

	radiuss = (lookkat - eyee);

	//m_context["new_eye"]->setFloat(eyee);
	m_context["count"]->setInt(0);					//these are CUDA functions and store initial data as zeros.
	m_context["eye"]->setFloat(make_float3(0.0f, 0.0f, 0.0f));	//count is counting number of simulations.
	m_context["U"]->setFloat(make_float3(0.0f, 0.0f, 0.0f));	//U,V,W 
	m_context["V"]->setFloat(make_float3(0.0f, 0.0f, 0.0f));
	m_context["W"]->setFloat(make_float3(0.0f, 0.0f, 0.0f));

	// 3D solid noise buffer, 1 float channel, all entries in the range [0.0, 1.0].
	srand(0); // Make sure the pseudo random numbers are the same every run.

	// Populate scene hierarch
	createGeometry();

	// Prepare to run
	m_context->validate();
	m_context->compile();
}


Buffer Tutorial::getOutputBuffer()
{
	return m_context["output_buffer"]->getBuffer();
}




void Tutorial::trace(const RayGenCameraData& camera_data)	//most important function for visualization
{



	if (flag)

	{
		for (int i = 0; i <= total_triangles - 1; i++)
		{
			enr[i] = 0.0; //rand();
			vol[i] = 0; // inital energy setting   
			hitting_p[i] = 0.0;
			total_p[i] = 0.0;


			//initial_cond[i] = 300;
		}



		myfile.open("Eye_direction.csv");	//Excel files to  store data
		//myfile3.open("normal.csv");
		myfile4.open("time.csv");


		m_context["eye"]->setFloat(camera_data.eye);	//Collect data from initialScene and set as a new camera for visualization
		m_context["U"]->setFloat(camera_data.U);	//Initial camera data from Sample Scene assigned to trace
		m_context["V"]->setFloat(camera_data.V);
		m_context["W"]->setFloat(camera_data.W);

		Buffer buffer = m_context["output_buffer"]->getBuffer();
		RTsize buffer_width, buffer_height;
		buffer->getSize(buffer_width, buffer_height);






		int vals2[] = { 0 };
		//Buffer hittingb = m_context["hitting_cord"]->getBuffer();

		//memcpy(hittingb->map(), vals2, sizeof(vals2));
		//hittingb->unmap();



		Buffer Vertex = m_context["vertex_buffer"]->getBuffer();	//It keeps the huge data of all points from GPU in a buffer
		//	Buffer Vertex_n = m_context["vertex_normal"]->getBuffer();

		//Buffer points = m_context["points"]->getBuffer();

		//Buffer hittingn = m_context["hitting_normal"]->getBuffer();

		Buffer b2 = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);
		b2->setFormat(RT_FORMAT_INT);

		b2->setSize(sizeof(vals2));
		memcpy(b2->map(), vals2, sizeof(vals2));
		b2->unmap();
		m_context["out2"]->set(b2);



		m_context["totaltraingles"]->setInt(total_triangles);

		//	m_context->launch(0, static_cast<unsigned int>(buffer_width),
		//	static_cast<unsigned int>(buffer_height));

		//		float3* hitting_cord = static_cast<float3*>(hittingb->map());

		//	hittingb->unmap();


		float3* vert = static_cast<float3*>(Vertex->map());

		Vertex->unmap();

		

		int* buffers1 = static_cast<int*>(b2->map());  // grab the resulting buffer. // number of rays hitting
		b2->unmap();

		int hitting_rays = buffers1[0];




		for (int i = 0; i <= total_points - 1; i++)
		{
			eyex[i] = vert[i].x;//hitting_cord[i].x;
			eyey[i] = vert[i].y;//hitting_cord[i].y;
			eyez[i] = vert[i].z;//hitting_cord[i].z;



		}


		//printf("\n\n yahan!");

		Face_center();  // calculate emission points and normals


		Buffer centersb = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);
		centersb->setFormat(RT_FORMAT_FLOAT);
		centersb->setSize(sizeof(centers));
		memcpy(centersb->map(), centers, sizeof(centers));
		centersb->unmap();
		m_context["centers"]->set(centersb);


		Buffer normals = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);
		normals->setFormat(RT_FORMAT_FLOAT);

		normals->setSize(sizeof(normal));
		memcpy(normals->map(), normal, sizeof(normal));
		normals->unmap();
		m_context["normals"]->set(normals);

		//m_context["vertex_buffer"]->setBuffer(centers);
		//m_context["vertex_buffer"]->setBuffer(normal);
		myfile3.open("normal.csv");
		for (int i = 0; i <= total_triangles - 1; i++)
		{
			myfile3 << normal[i].x << ";" << normal[i].y << ";" << normal[i].z << std::endl;
		}


		printf("\n\n source done");

		m_context["top_object"]->set(geometrygroup);		//determine priority here.We assign most important objects
		m_context["top_shadower"]->set(geometrygroup);
		flag = false;



	}

	else



	{

		//myfile.open("sphere_rectangle.csv");
		//RTresult RTAPI rtGeometryDestroy(GeometryGroup gemoetry_group 	geometry);

		int flags = 0;						// ALmost every buffer, variable should be created for else part


		m_context["flag"]->setInt(flags);




		double time = GetTickCount();	//store the starting time
		double b = 0;
		float final_end;
		int flag = 1;


		//	m_context["flag"]->setInt(flag);

		m_context["eye"]->setFloat(camera_data.eye);	//Previous camera data also should be assign here for else
		m_context["U"]->setFloat(camera_data.U);	//
		m_context["V"]->setFloat(camera_data.V);
		m_context["W"]->setFloat(camera_data.W);








		Buffer buffer = m_context["output_buffer"]->getBuffer();
		RTsize buffer_width, buffer_height;
		buffer->getSize(buffer_width, buffer_height);



		Buffer bb = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);
		bb->setFormat(RT_FORMAT_INT);
		int vals[] = { 0 };
		bb->setSize(sizeof(vals));
		memcpy(bb->map(), vals, sizeof(vals));
		bb->unmap();
		m_context["out"]->set(bb);



	

		Buffer Extinction = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);
		Extinction->setFormat(RT_FORMAT_FLOAT);
		Extinction->setSize(sizeof(enr)); // inital energy
		memcpy(Extinction->map(), enr, sizeof(enr));
		Extinction->unmap();
		m_context["Extinction"]->set(Extinction);

		Buffer hitting_prim = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);
		hitting_prim->setFormat(RT_FORMAT_INT);
		hitting_prim->setSize(sizeof(hitting_p)); // inital energy
		memcpy(hitting_prim->map(), hitting_p, sizeof(hitting_p));
		hitting_prim->unmap();
		m_context["hitting_prim"]->set(hitting_prim);

		Buffer total_prim = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);
		total_prim->setFormat(RT_FORMAT_INT);
		total_prim->setSize(sizeof(total_p)); // inital energy
		memcpy(total_prim->map(), total_p, sizeof(total_p));
		total_prim->unmap();
		m_context["total_prim"]->set(total_prim);






		Buffer b3 = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);
		b3->setFormat(RT_FORMAT_INT);
		int vals2[] = { 0 };
		b3->setSize(sizeof(vals2));
		memcpy(b3->map(), vals2, sizeof(vals2));
		b3->unmap();
		m_context["out3"]->set(b3);

		Buffer count = m_context->createBuffer(RT_BUFFER_INPUT_OUTPUT);
		count->setFormat(RT_FORMAT_INT);

		count->setSize(sizeof(vals2));
		//	memcpy(count->map(), vals2, sizeof(vals2));
		//	count->unmap();

		m_context["countbuffer"]->set(count);

		int max = total_triangles;


		m_context->launch(1, static_cast<unsigned int>(buffer_width),  
			static_cast<unsigned int>(buffer_height));



		int i = 0;

		int* countt = static_cast<int*>(count->map());

		count->unmap();





		int* buffers = static_cast<int*>(bb->map());   // grab the resulting buffer. ..  generating rays
		bb->unmap();

		int* buffers1 = static_cast<int*>(b3->map());  // grab the resulting buffer. // number of rays hitting
		b3->unmap();

		int hitting_rays = buffers1[0];
		int total_rays = buffers[0];



		float* total_pr = static_cast<float*>(total_prim->map());

		total_prim->unmap();

		float* hitting_pr = static_cast<float*>(hitting_prim->map());

		hitting_prim->unmap();














		//if (count_s < 7056)
		{

			float rad = length(radiuss);
			//float view_factor1 = (effected_area) / (M_PI*(rad*rad)); //(Total_pixel[250][250]*hitting_rays) // divide by 4 pi for sphere

			icount = countt[0];

			printf("\n\n count: %d", icount);

			float viewfactor;

			double VF = (float(hitting_rays)) / total_rays;   // for sphere


			//	double Flux = VF / (0.111111*0.111111);

			b = GetTickCount() - time;




			//if (!(VF < 1))
			//{

			//	VF = 0;

			//}


			printf(" \n\n View Factor is %f \n hitting rays is %d \n generated rays are : %d", VF, hitting_rays, buffers[0]);

			printf("\n\n\n time taken per sample in miliseconds is: %f  ", b);


			if (icount >= 0)
				//{

				myfile2 << VF << "," << centers[count_s].x << "," << centers[count_s].y << "," << centers[count_s].z << "," << b << "," << std::endl;
			//}
			float volume[5832] = { 0 }; //5832
			float ext[5832] = { 0 };



		}
		//if (count_s > 100)
		{
			//memcpy(Extinction->map(), enr, sizeof(enr));
			//Extinction->unmap();
			//m_context["Extinction"]->set(Extinction);
			m_context->launch(0, static_cast<unsigned int>(buffer_width),	//Everything(Calculation,view factor,geometries) is ready to launch
				static_cast<unsigned int>(buffer_height));		//All reaction happens,starts with this "launch" function

		}



		for (i = 0; i <= total_triangles - 1; i++) //total_triangles
		{

			total_p[i] = total_pr[i];
			hitting_p[i] = hitting_pr[i];

			float v = hitting_pr[i] / total_pr[i];



			//{
			//	myfile2 << hitting_p[i] << ";";
			//}

			//	myfile4 << v << ";";

		}

		//myfile4 << std::endl;
		myfile4 << b << std::endl;

	}

	count_s++;






}




std::string Tutorial::texpath(const std::string& base)
{
	return texture_path + "/" + base;
}

float4 make_plane(float3 n, float3 p)
{
	n = normalize(n);
	float d = -dot(n, p);
	return make_float4(n, d);
}




void Tutorial::InsertModel(const std::string& name, GeometryGroup ggroup, Material mat, Program mint, const optix::Matrix4x4 matx, bool suppressErrors)
{
	try {


		OptiXMesh mesh(m_context, ggroup, mat, m_accel_desc);


		mesh.setDefaultIntersectionProgram(mint);
		mesh.setLoadingTransform(matx);
		mesh.loadBegin_Geometry(name);

		total_points = mesh.getNumVertices();
		total_triangles = mesh.getNumTriangles();

		printf("\n no of triangles is : %d", total_triangles);

		printf("\n\n number of vertices are :  %d", total_points);

		mesh.loadFinish_Materials();

		//extracting face centers from obj file



		std::ifstream myfile;
		myfile.open("C:/ProgramData/NVIDIA Corporation/OptiX SDK 3.9.1/SDK/tutorial/cylinder.obj", std::ios::in);
		if (!(myfile.is_open()))
		{
			std::cout << "Error Opening File";
			getchar();
			std::exit(0);
		}

		std::string firstline;
		char my_character;
		int test;

		int i = 0;

		while (!myfile.eof())
		{

			//std::getline(myfile, firstline);
			getline(myfile, firstline);

			if ((firstline[0] == 'f'))
			{

				const char * c = firstline.c_str();	//Object files have #of surfaces, normal directions,locations,view factors,vectors
									//Every line begins with name of information in that line
				int i1, i2, i3, i4, i5, i6, i7, i8, i9;	//For example, faces line has "f" in first column of line or "vf" for view factor
				if (6 == sscanf(c, "%*[^0123456789]%d%*[^0123456789]%d%*[^0123456789]%d%*[^0123456789]%d%*[^0123456789]%d%*[^0123456789]%d",
					&i1,
					&i2,
					&i3, &i4, &i5, &i6))
				{

					Faces[i][0] = i1;
					Faces[i][1] = i3;
					Faces[i][2] = i5;

				}

				//if (9 == sscanf(c, "%*[^0123456789]%d%*[^0123456789]%d%*[^0123456789]%d%*[^0123456789]%d%*[^0123456789]%d%*[^0123456789]%d%*[^0123456789]%d%*[^0123456789]%d%*[^0123456789]%d",
				//&i1,
				//&i2,
				//&i3, &i4, &i5, &i6, &i7, &i8, &i9))
				//{
				//Faces[i][0] = i1;
				//Faces[i][1] = i4;
				//Faces[i][2] = i7;
				//}

				i++;

			}







		}

		myfile.close();



	}
	catch (...) {
		if (suppressErrors)
			return;
		else
			throw;
	}

}

void Tutorial::createGeometry()
{


	const float matrix_0[4 * 4] = { 1.0, 0, 0, 0,	//create matrixes as much as objects 
		0, 1.0, 0, 4.5,				//to place them a place in our space
		0, 0, 1.0, 0.0,
		0, 0, 0, 1 };

	const optix::Matrix4x4 m0(matrix_0);


	const float matrix_1[4 * 4] = { 1, 0, 0, 0.0,
		0, 1, 0, 0.0,
		0, 0, 1, 0.0,
		0, 0, 0, 1 };

	const optix::Matrix4x4 m1(matrix_1);




	std::string prog_path(ptxpath("glass", "triangle_mesh_iterative.cu"));
	Program mesh_intersect = m_context->createProgramFromPTXFile(prog_path, "mesh_intersect");






	// TRIANGLE GEOMETRY
	//Create Triangle

	std::string ptx_path = ptxpath("primitiveIndexOffsets", "triangle_mesh.cu");
	Geometry triangle = m_context->createGeometry();
	CreateTriangle(m_context, triangle, ptx_path);



	//Geometry triangle2 = m_context->createGeometry();
	//CreateTriangle2(m_context, triangle2, ptx_path);




	// Create box

	std::string box_ptx(ptxpath("tutorial", "box.cu"));
	Geometry box = m_context->createGeometry();
	CreateBox(m_context, box, box_ptx);

	Program box_bounds = m_context->createProgramFromPTXFile(box_ptx, "box_bounds");
	Program box_intersect = m_context->createProgramFromPTXFile(box_ptx, "box_intersect");
	Geometry box1 = m_context->createGeometry();
	box1->setPrimitiveCount(1u);
	box1->setBoundingBoxProgram(box_bounds);

	box1->setIntersectionProgram(box_intersect);	// to be able to know whether there is an object or not, 
							//Optix intersectionprogram will do the job

	box1["boxmin"]->setFloat(-1.0f, -1.0f, 2.25f);    //0.0f, 40000.0f, -600000.0f //0.3888889
	box1["boxmax"]->setFloat(1.0f, 1.0f, 2.25f);


	Geometry box2 = m_context->createGeometry();

	box2->setPrimitiveCount(1u);
	box2->setBoundingBoxProgram(box_bounds);
	box2->setIntersectionProgram(box_intersect);
	box2["boxmin"]->setFloat(7.0f, -15.0f, -15.0f);    //0.0f, 40000.0f, -600000.0f //0.3888889
	box2["boxmax"]->setFloat(7.0f, 15.0f, 15.0f);


	// Materials

	Material box_matl = m_context->createMaterial();	//assign material type of object 
	material(box_matl, m_context);


	//Material 
	//parabola_matl = m_context->createMaterial();
	// material_parabola(parabola_matl, m_context);

	Material object_matl = m_context->createMaterial();		
	material_object(object_matl, m_context);

	// Create GIs for each piece of geometry
	std::vector<GeometryInstance> gis;



	// which geometry you want ? 
	{
		//gis.push_back(m_context->createGeometryInstance(triangle, &box_matl, &box_matl + 1));
		//gis.push_back(m_context->createGeometryInstance(triangle2, &object_matl, &object_matl + 1));
		//gis.push_back(m_context->createGeometryInstance(box, &object_matl, &object_matl + 1));
		//	gis.push_back(m_context->createGeometryInstance(box1, &object_matl, &object_matl + 1));
		//	gis.push_back(m_context->createGeometryInstance(box2, &object_matl, &object_matl + 1));
		//	gis.push_back(m_context->createGeometryInstance(cylinder, &object_matl, &object_matl + 1));
		// gis.push_back(m_context->createGeometryInstance(cylinder2, &object_matl, &object_matl + 1));
		//  gis.push_back(m_context->createGeometryInstance(boxa, &box_matl, &box_matl + 1));
		//gis.push_back(m_context->createGeometryInstance(disk, &box_matl, &box_matl + 1));
		//	gis.push_back(m_context->createGeometryInstance(disk2, &object_matl, &object_matl + 1));
		//   gis.push_back(m_context->createGeometryInstance(sphere, &box_matl, &box_matl + 1));
		//gis.push_back( m_context->createGeometryInstance( parallelogram, &floor_matl, &floor_matl+1 ) );
	}


	//printf("\n\n yahan tk!!");

	// Place all in group
	geometrygroup = m_context->createGeometryGroup();
	geometrygroup->setChildCount(static_cast<unsigned int>(gis.size()));

	GeometryGroup geometrygroupobj = m_context->createGeometryGroup();
	//geometrygroup->setChild(0, gis[0]);
	//geometrygroup->setChild(1, gis[1]);
	//geometrygroup->setChild(2, gis[2]);

	//InsertModel(prog_path, geometrygroup, box_matl, mesh_intersect, m0, true);

	InsertModel("C:/ProgramData/NVIDIA Corporation/OptiX SDK 3.9.1/SDK/tutorial/sphere.obj", geometrygroup, object_matl, mesh_intersect, m0, false); //C:/ProgramData/NVIDIA Corporation/OptiX SDK 3.8.0/SDK/tutorial
	InsertModel("C:/ProgramData/NVIDIA Corporation/OptiX SDK 3.9.1/SDK/tutorial/cylinder.obj", geometrygroup, box_matl, mesh_intersect, m1, false);//geometry type, material, intersection program, location should be defined








	// geometrygroup->setChild(1, gis[1]);
	// geometrygroup->setChild(1, gis[1]);  //mainbox



	std::string path(std::string(sutilSamplesDir()) + "/tutorial"); //for texture path
	//	geometrygroup->setAcceleration(m_context->createAcceleration("NoAccel", "NoAccel"));  //NoAccel
	//acceleration accel;
	//rtAccelerationCreate(m_context, &accel);
	//rtAccelerationSetBuilder(accel, "Trbvh");
	//rtAccelerationSetTraverser(accel, "Bvh");
	geometrygroup->setAcceleration(m_context->createAcceleration("Bvh", "Bvh"));  	//Bounding Volume Hierarchy
											//there are bvh,sah,NoAcc but fastest is bvh.


	m_context["top_object"]->set(geometrygroup);
	m_context["top_shadower"]->set(geometrygroup);




}




//-----------------------------------------------------------------------------
//
// Main driver
//gis
//-----------------------------------------------------------------------------

void printUsageAndExit(const std::string& argv0, bool doExit = true)
{
	std::cerr
		<< "Usage  : " << argv0 << " [options]\n"
		<< "App options:\n"
		<< "  -h  | --help                               Print this usage message\n"
		<< "  -T  | --tutorial-number <num>              Specify tutorial number\n"
		<< "  -t  | --texture-path <path>                Specify path to texture directory\n"
		<< "        --dim=<width>x<height>               Set image dimensions\n"
		<< std::endl;
	GLUTDisplay::printUsage();

	if (doExit) exit(1);
}


int main(int argc, char** argv)
{

	GLUTDisplay::init(argc, argv);

	unsigned int width = 500u, height = 500u; //1366x768  // Number of dots



	std::string texture_path;
	int tutnum = 1;
	for (int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		if (arg == "--help" || arg == "-h") {
			printUsageAndExit(argv[0]);
		}
		else if (arg.substr(0, 6) == "--dim=") {
			std::string dims_arg = arg.substr(6);
			if (sutilParseImageDimensions(dims_arg.c_str(), &width, &height) != RT_SUCCESS) {
				std::cerr << "Invalid window dimensions: '" << dims_arg << "'" << std::endl;
				printUsageAndExit(argv[0]);
			}
		}
		else if (arg == "-t" || arg == "--texture-path") {
			if (i == argc - 1) {
				printUsageAndExit(argv[0]);
			}
			texture_path = argv[++i];
		}
		else if (arg == "-T" || arg == "--tutorial-number") {
			if (i == argc - 1) {
				printUsageAndExit(argv[0]);
			}
			tutnum = atoi(argv[++i]);
		}
		else {
			std::cerr << "Unknown option: '" << arg << "'\n";
			printUsageAndExit(argv[0]);
		}
	}

	if (!GLUTDisplay::isBenchmark()) printUsageAndExit(argv[0], false);

	if (tutnum < 0 || tutnum > 11) {
		std::cerr << "Tutorial number (" << tutnum << ") is out of range [0..11]\n";
		exit(1);
	}

	if (texture_path.empty()) {
		texture_path = std::string(sutilSamplesDir()) + "/tutorial/data";
	}

	std::stringstream title;
	title << "Tutorial " << tutnum;
	try {

		Tutorial scene(tutnum, texture_path);


		scene.setDimensions(width, height);
		GLUTDisplay::run(title.str(), &scene);

	}
	catch (Exception& e){
		sutilReportError(e.getErrorString().c_str());
		exit(1);
	}

	return 0;
	getchar();
	getchar();
}


float area(double a[], double b[], double c[], int length)
{


	int i = 0;
	double * x;
	x = new double[length];
	double* y;
	y = new double[length];
	double * z;
	z = new double[length];
	int l = length;
	for (i = 0; i <= length - 1; i++)
	{
		x[i] = a[i];
		y[i] = b[i];
		z[i] = c[i];
	}
	double edge0[] = { x[1] - x[0], y[1] - y[0], z[1] - z[0] };
	double edge1[] = { x[2] - x[0], y[2] - y[0], z[2] - z[0] };
	double  nor3[] = { edge0[1] * edge1[2] - edge0[2] * edge1[1], edge0[2] * edge1[0] - edge0[0] * edge1[2], edge0[0] * edge1[1] - edge0[1] * edge1[0] };
	float inveln = 1 / (sqrt(nor3[0] * nor3[0] + nor3[1] * nor3[1] + nor3[2] * nor3[2]));
	for (i = 0; i <= 2; i++)
	{
		nor3[i] = inveln*nor3[i];
	}
	float csumx[4] = { 0 }, csumy[4] = { 0 }, csumz[4] = { 0 };
	i = 0;
	while (i <= (l - 2))
	{
		csumx[i] = y[i] * z[i + 1] - z[i] * y[i + 1];
		csumy[i] = z[i] * x[i + 1] - x[i] * z[i + 1];
		csumz[i] = x[i] * y[i + 1] - y[i] * x[i + 1];
		i++;
	}

	float aa = 0.0, bb = 0.0, cc = 0.0;
	for (i = 0; i <= (l - 2); i++)
	{
		aa = aa + csumx[i];
		bb = bb + csumy[i];
		cc = cc + csumz[i];
	}

	double sumx = aa + y[l - 1] * z[0] - z[l - 1] * y[0];
	double sumy = bb + z[l - 1] * x[0] - x[l - 1] * z[0];
	double sumz = cc + x[l - 1] * y[0] - y[l - 1] * x[0];
	double area = (abs(nor3[0] * sumx + nor3[1] * sumy + nor3[2] * sumz)) / 2;
	if (!(area < 1))
	{
		area = 0;
	}
	return area;

}


void material(Material box_matl, Context m_context)
{


	std::string box_chname;
	std::string box_ahname;
	box_chname = "closest_hit_radiance1";    //closest_hit_radiance1

	//	box_ahname = "any_hit2";
	Program box_ch = m_context->createProgramFromPTXFile(m_ptx_path, box_chname);
	//Program box_ah = m_context->createProgramFromPTXFile(m_ptx_path, box_ahname);


	box_matl->setClosestHitProgram(0, box_ch);	//determine material type while creating instance with intersection
	//	box_matl->setAnyHitProgram(0, box_ah);

	Material matl = m_context->createMaterial();

	box_matl["Ka"]->setFloat(0.1f, 0.1f, 0.1f);//0.3
	box_matl["Kd"]->setFloat(0.3f, 0.3f, 0.3f);// 0.6 0.7 0.8
	box_matl["Ks"]->setFloat(0.8f, 0.9f, 0.8f);// 0.8 0.9 0.8
	box_matl["phong_exp"]->setFloat(88);//88
	box_matl["reflectivity_n"]->setFloat(0.0f, 0.0f, 0.0f); //0.2


}


void material_parabola(Material parabola_matl, Context m_context)
{


	std::string par_chname;
	std::string par_ahname;
	par_chname = "closest_hit_radiance_parabola";
	par_ahname = "any_hit_shadow";
	Program box_ch = m_context->createProgramFromPTXFile(m_ptx_path, par_chname);
	Program box_ah = m_context->createProgramFromPTXFile(m_ptx_path, par_ahname);


	parabola_matl->setClosestHitProgram(0, box_ch);
	parabola_matl->setAnyHitProgram(1, box_ah);

	Material matl = m_context->createMaterial();

	parabola_matl["Ka"]->setFloat(0.3f, 0.3f, 0.3f);//0.3
	parabola_matl["Kd"]->setFloat(0.3f, 0.3f, 0.3f);// 0.6 0.7 0.8
	parabola_matl["Ks"]->setFloat(0.8f, 0.9f, 0.8f);// 0.8 0.9 0.8
	parabola_matl["phong_exp"]->setFloat(88);//88
	parabola_matl["reflectivity_n"]->setFloat(0.0f, 0.0f, 0.0f); //0.2


}



void CreateBox(optix::Context m_context, Geometry box, std::string box_ptx){

	Program box_bounds = m_context->createProgramFromPTXFile(box_ptx, "box_bounds");
	Program box_intersect = m_context->createProgramFromPTXFile(box_ptx, "box_intersect");

	box->setPrimitiveCount(1u);
	box->setBoundingBoxProgram(box_bounds);
	box->setIntersectionProgram(box_intersect);
	box["boxmin"]->setFloat(-1.0f, -1.0f, 1.5f);    //0.0f, 40000.0f, -600000.0f //0.3888889
	box["boxmax"]->setFloat(1.0f, 1.0f, 1.5f);    //100000.0f, 80000.0f, -400000.0     //0.5

}




void material_object(Material obj_matl, Context m_context)
{


	std::string par_chname;
	std::string par_ahname;
	par_chname = "closest_hit_radiance_object"; //closest_hit_radiance_object
	//par_ahname = "any_hit";
	Program box_ch = m_context->createProgramFromPTXFile(m_ptx_path, par_chname);
	//Program box_ah = m_context->createProgramFromPTXFile(m_ptx_path, par_ahname);


	obj_matl->setClosestHitProgram(0, box_ch);
	//obj_matl->setAnyHitProgram(0, box_ah);

	Material matl = m_context->createMaterial();

	obj_matl["Ka"]->setFloat(0.1f, 0.1f, 0.1f);//0.3
	obj_matl["Kd"]->setFloat(0.3f, 0.3f, 0.3f);// 0.6 0.7 0.8
	obj_matl["Ks"]->setFloat(0.8f, 0.9f, 0.8f);// 0.8 0.9 0.8
	obj_matl["phong_exp"]->setFloat(88);//88
	obj_matl["reflectivity_n"]->setFloat(0.0f, 0.0f, 0.0f); //0.2


}





void CreateTriangle(optix::Context m_context, Geometry triangle, std::string ptx_path){



	Program mesh_intersect = m_context->createProgramFromPTXFile(ptx_path, "mesh_intersect");
	Program mesh_bounds = m_context->createProgramFromPTXFile(ptx_path, "mesh_bounds");
	triangle->setPrimitiveCount(1u);
	triangle->setBoundingBoxProgram(mesh_bounds);
	triangle->setIntersectionProgram(mesh_intersect);
	Buffer verts_buffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT3, 4 + 4);
	Buffer index_buffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_INT3, 4 + 2);
	Buffer matl_buffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_UNSIGNED_INT, 4 + 2);

	Buffer tindex_buffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_INT3, 4 + 2);
	Buffer nindex_buffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_INT3, 4 + 2);
	Buffer nbuffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT3, 0);
	Buffer tbuffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT2, 0);


	float3*   verts = reinterpret_cast<float3*>(verts_buffer->map());
	int3*     indices = reinterpret_cast<int3*>(index_buffer->map());
	unsigned* mats = static_cast<unsigned*>(matl_buffer->map());


	triangle["vertex_buffer"]->setBuffer(verts_buffer);
	triangle["vindex_buffer"]->setBuffer(index_buffer);
	m_context["material_buffer"]->setBuffer(matl_buffer);


	triangle["normal_buffer"]->setBuffer(nbuffer);
	triangle["nindex_buffer"]->setBuffer(nindex_buffer);
	triangle["texcoord_buffer"]->setBuffer(tbuffer);
	triangle["tindex_buffer"]->setBuffer(tindex_buffer);

	mats[0] = 0u;
	mats[1] = 0u;
	mats[2] = 0u;
	mats[3] = 0u;

	verts[0] = make_float3(-1.0f, -1.0f, 3.25f);//   /(sqrtf(2.0f))
	verts[1] = make_float3(0.9f, -1.0f, 3.25f);
	verts[2] = make_float3(0.9f, 0.9f, 3.25f);
	verts[3] = make_float3(0.9f, 0.9f, 3.25f);

	indices[0] = make_int3(0, 1, 2); // make_int3( i*4+0, i*4+1, i*4+2 );
	indices[1] = make_int3(1, 3, 2);
	indices[2] = make_int3(3, 0, 2);
	indices[3] = make_int3(1, 0, 3);



	verts_buffer->unmap();
	index_buffer->unmap();
	matl_buffer->unmap();
}


/*void CreateTriangle2(optix::Context m_context, Geometry triangle, std::string ptx_path){



Program mesh_intersect = m_context->createProgramFromPTXFile(ptx_path, "mesh_intersect2");
Program mesh_bounds = m_context->createProgramFromPTXFile(ptx_path, "mesh_bounds");
triangle->setPrimitiveCount(1u);
triangle->setBoundingBoxProgram(mesh_bounds);
triangle->setIntersectionProgram(mesh_intersect);
Buffer verts_buffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT3, 4 + 4);
Buffer index_buffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_INT3, 4 + 2);
Buffer matl_buffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_UNSIGNED_INT, 4 + 2);

Buffer tindex_buffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_INT3, 4 + 2);
Buffer nindex_buffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_INT3, 4 + 2);
Buffer nbuffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT3, 0);
Buffer tbuffer = m_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT2, 0);


float3*   verts = reinterpret_cast<float3*>(verts_buffer->map());
int3*     indices = reinterpret_cast<int3*>(index_buffer->map());
unsigned* mats = static_cast<unsigned*>(matl_buffer->map());


triangle["vertex_buffer"]->setBuffer(verts_buffer);
triangle["vindex_buffer"]->setBuffer(index_buffer);
m_context["material_buffer"]->setBuffer(matl_buffer);


triangle["normal_buffer"]->setBuffer(nbuffer);
triangle["nindex_buffer"]->setBuffer(nindex_buffer);
triangle["texcoord_buffer"]->setBuffer(tbuffer);
triangle["tindex_buffer"]->setBuffer(tindex_buffer);

mats[0] = 0u;
mats[1] = 0u;
mats[2] = 0u;
mats[3] = 0u;

verts[0] = make_float3(1.0f, 0.0f, 2.0f);//   /(sqrtf(2.0f))
verts[1] = make_float3(1.0f, 0.0f, 9.0f);
verts[2] = make_float3(1.0f, 7.0f, 2.0f);
verts[3] = make_float3(1.0f, 7.0f, 2.0f);

indices[0] = make_int3(0, 1, 2); // make_int3( i*4+0, i*4+1, i*4+2 );
indices[1] = make_int3(1, 3, 2);
indices[2] = make_int3(3, 0, 2);
indices[3] = make_int3(1, 0, 3);


verts_buffer->unmap();
index_buffer->unmap();
matl_buffer->unmap();
}
*/
void Face_center()
{




	for (int i = 0; i < total_triangles; i++)

	{

		int a = Faces[i][0];
		int b = Faces[i][1];
		int c = Faces[i][2];

		//	float r1 = ((float)rand() / (float)RAND_MAX);
		//	float r2 = ((float)rand() / (float)RAND_MAX);
		//	float r3 = ((float)rand() / (float)RAND_MAX);


		double bx, by, bz;

		double cx = ((eyex[a - 1] + eyex[b - 1] + eyex[c - 1]) / 3);

		double cy = ((eyey[a - 1] + eyey[b - 1] + eyey[c - 1]) / 3);
		double cz = ((eyez[a - 1] + eyez[b - 1] + eyez[c - 1]) / 3);

		/*	double cx = (((eyex[a - 1] + eyex[c - 1]) / 2) + eyex[b - 1]) / 2 ;

		double cy = (((eyey[a - 1] + eyey[c - 1]) / 2) + eyey[b - 1]) / 2 ;
		double cz = (((eyez[a - 1] + eyez[c - 1]) / 2) + eyez[b - 1]) / 2 ;*/


		//	printf("\n\n centers is : %f, %f ,%f", cx, cy, cz);
		//getchar();

		centers[i].x = cx;
		centers[i].y = cy;
		centers[i].z = cz;


		float3 p1 = make_float3(eyex[a - 1], eyey[a - 1], eyez[a - 1]);
		float3 p2 = make_float3(eyex[b - 1], eyey[b - 1], eyez[b - 1]);
		float3 p3 = make_float3(eyex[c - 1], eyey[c - 1], eyez[c - 1]);

		normal[i] = normalize(cross(p2 - p1, p3 - p1));  //p2 - p1, p3 - p1

		if (normal[i].x == 0)
		{
			normal[i].x = 0.0000000001;
		}

		if (normal[i].y == 0)
		{
			normal[i].y = 0.0000000001;

		}
		if (normal[i].z == 0)
		{
			normal[i].z = 0.0000000001;
		}








	}





}

/*for (i = 0; i <= buffers1[0] - 2; i++)
{

int indxi = buffer_launch[i].x;
int indxj = buffer_launch[i].y;



//printf("\n\ indxi is : %d", indxi);
//	double xx[] = { node_coordinates[indxi][indxj].x, node_coordinates[indxi][indxj + 1].x, node_coordinates[indxi + 1][indxj + 1].x, node_coordinates[indxi + 1][indxj].x, node_coordinates[indxi][indxj].x };
//	double yy[] = { node_coordinates[indxi][indxj].y, node_coordinates[indxi][indxj + 1].y, node_coordinates[indxi + 1][indxj + 1].y, node_coordinates[indxi + 1][indxj].y, node_coordinates[indxi][indxj].y };
//	double zz[] = { node_coordinates[indxi][indxj].z, node_coordinates[indxi][indxj + 1].z, node_coordinates[indxi + 1][indxj + 1].z, node_coordinates[indxi + 1][indxj].z, node_coordinates[indxi][indxj].z };

//Areapixel[i] = area(xx, yy, zz, 5);

//if (!(Areapixel[i] > 0))

{
// Areapixel[i] = 0;
} */a
