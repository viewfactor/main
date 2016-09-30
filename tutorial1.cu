


#include "tutorial.h"
//#include <optix_device.h>
#include "random.h"
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>
#include <stdio.h>
//#include<random>


rtDeclareVariable(float3, shading_normal, attribute shading_normal, );	//all basic functions to create geometry occurs here,
rtDeclareVariable(float3, geometric_normal, attribute geometric_normal, );//Defined functions can be used in everywhere  
rtDeclareVariable(int, primid, attribute primid, );

rtDeclareVariable(PerRayData_radiance, prd_radiance, rtPayload, );

rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(float, t_hit, rtIntersectionDistance, );
rtDeclareVariable(uint2, launch_index, rtLaunchIndex, );
rtDeclareVariable(unsigned int, radiance_ray_type, , );
rtDeclareVariable(float, scene_epsilon, , );
rtDeclareVariable(rtObject, top_object, , );
rtDeclareVariable(uint, launch_index1, rtLaunchIndex, );


static __device__ __inline__ float fold(const float value)
{
	return fminf(value, 1.0f - value) * 2.0f;
}



//
// Pinhole camera implementation
//
rtDeclareVariable(float3, eye, , );		//we use them in main code for calculation camera and visualization camera
rtDeclareVariable(float3, U, , );
rtDeclareVariable(float3, V, , );
rtDeclareVariable(float3, W, , );



rtDeclareVariable(float3, bad_color, , );
rtBuffer<uchar4, 2>              output_buffer;


rtBuffer<float3, 1>            normals;		//buffers provide temporary memory for us to be able to pass them to other variables
rtBuffer<float3, 1>            centers;
rtBuffer<float3, 1>            vertex_buffer;
//rtBuffer<float3, 1>            points;
//rtBuffer<float3, 1>            vertex_buffer2;
//rtBuffer<float, 1>             Energy;
rtBuffer<float, 1>             Energy_init;
rtBuffer<float, 1>             Extinction;
rtBuffer<float, 1>             hitting_prim;
rtBuffer<float, 1>             total_prim;
//rtBuffer<float, 1>             Energy_volume;
rtBuffer<float3, 1>            vertex_normal;
//rtBuffer<float3, 1>            vertex_normal2;
//rtBuffer<float3, 1>            hitting_normal;
//rtBuffer<float3, 1>             hitting_cord;
//rtBuffer<float, 2>              output_bufferx;
//rtBuffer<float, 2>              output_buffery;
//rtBuffer<float, 2>              output_bufferz;
rtBuffer<float, 2>              output_bufferl;
//rtBuffer<uint2, 1>              output_launch_index;
rtBuffer<int, 1>                  out;
rtBuffer<int, 1>                  out2;
rtBuffer<int, 1>                  out3;



rtBuffer<float2, 1>                  a;
rtBuffer<uint, 1>                  countbuffer;
rtDeclareVariable(float3, new_eye, , );
rtDeclareVariable(float3, new_normal, , );
rtDeclareVariable(int, ID, , );
rtDeclareVariable(int, count, , );
rtDeclareVariable(int, maxcount, , );
rtDeclareVariable(float, viewfactor, , );
rtDeclareVariable(int, totaltraingles, , );
rtDeclareVariable(int, flag, , );

//random yahan



RT_PROGRAM void random_camera()
{

	size_t2 screen = output_buffer.size();

	float2 d = make_float2(launch_index) / make_float2(screen) * 2.f - 1.f;
	float3 ray_origin = eye;							//ray source point
	float3 ray_direction = normalize(d.x*U + d.y*V + W);				// make unit position to real position

	optix::Ray ray(ray_origin, ray_direction, radiance_ray_type, scene_epsilon);	//its generation function

	PerRayData_radiance prd;	
	prd.importance = 1.f;
	prd.depth = 0;

	rtTrace(top_object, ray, prd);							//important function to trace and find intersection point

	output_buffer[launch_index] = make_color(prd.result);				//demonstrate result as different colors

}



RT_PROGRAM void pinhole_camera()
{

	PerRayData_radiance prd;


	// double coef = (pow(a, e));

	size_t2 screen = output_buffer.size();

	uint2 seed;
	seed.x = tea<16>(count + screen.x*launch_index.y + launch_index.x, (count + 100));
	seed.y = tea<16>(count + screen.y*launch_index.x + launch_index.y, (count + 100));
	uint seedz = tea<16>(2 + +screen.y*launch_index.x + launch_index.y, 112);

	uint seez = seed.x + seed.y;
	double r1 = rnd(seed.x);
	double r2 = rnd(seed.y);
	double r3 = rnd((seez));



	//prd.random = r1;

	int ran = 0;

	//curandState state;

	//curand_init(seed.x, 0, 0, &state);   // generate random numbers




	if (screen.x*launch_index.y + launch_index.x < totaltraingles)
	{
		ran = screen.x*launch_index.y + launch_index.x;//curand(&state) % totaltraingles;    //screen.x*launch_index.y + launch_index.x; //curand(&state) % totaltraingles; 



	}
	else
	{
		ran = (screen.x*launch_index.y + launch_index.x) % (totaltraingles - 1);//curand(&state) % ((totaltraingles));//(screen.x*launch_index.y + launch_index.x) % (totaltraingles - 1); //curand(&state) % ((totaltraingles - 1));// //(screen.x*launch_index.y + launch_index.x) % (totaltraingles);//  //// / 

	}
	//rtPrintf("rand is %d", ran);


	float ph = 0 * M_PI / 180;    //180
	float th = (90)*M_PI / 180; //60



	float3 W_normal = normals[ran]; //make_float3(0.000001, 1, 0.0000001);//normals[ran]; //make_float3(sin(th)*sin(ph), cos(th), -sin(th)*cos(ph));//normals[screen.x*launch_index.y + launch_index.x]; // normals[count];// normals[ran];// vertex_normal[count]; //normals[count];////normals[ran];//new_normal;//normals[count];//new_normal;//vertex_normal[count]; // new_normal //new_normal



	float3 up_new = make_float3(0, 1, 0);
	up_new = normalize(up_new);

	float3 U_new = cross(W_normal, up_new);
	float3 V_new = cross(U_new, W_normal);




	//random




	double rad = length(W_normal);   //W_normal

	float3 rx = (((rad*normalize((U_new))))); //U_new  
	float3 ry = (rad*normalize(V_new));   // V_new
	float3 rz = (rad*normalize(W_normal)); //W_new

	// mapping into hemisphere

	float ee = 0;
	float cos_phi = cos((2 * M_PI)*r1);                 // for whole sphere cos( (M_PI)*d.x); where d.x is random from -1 to 1
	float sin_phi = sin((2 * M_PI)*r1);                 // // for hemisohere sphere cos( (2*M_PI)*r1); where 1 is random from 0 to 1
	//float cos_theta = (((1 - 2 * r2)));                       // remove 2 fro original method
	double sin_theta = sqrt(r2);
	float cos_theta = sqrt(1 - (sin_theta*sin_theta));  // from pinar hoca's book
	// for new try
	//float sin_theta = sqrt(1 - (cos_theta*cos_theta));  // cosine distribtion
	float px = sin_theta*cos_phi;//sin_theta*sin_phi;
	float py = sin_theta*sin_phi; //cos_theta;
	float pz = cos_theta;//sin_theta*cos_phi;


	float3 d_mapping = normalize(px*(rx)+py*(ry)+pz*(rz));



	//d_f

	//selecting random position within cube



	float3 ray_origin = centers[ran];// point;//make_float3(0.01f, 0.02f, 0.2f);//point;//centers[screen.x*launch_index.y + launch_index.x];// centers[count];// centers[count];// point;// //point;// centers[ran]; //point;// // centers[ran];  //vertex_buffer[count];// centers[ran];// new_eye;//vertex_buffer[count];  centers[count];// //new_eye for sampling



	int a = atomicAdd(&out[0], 1);  // counting generated rays

	atomicAdd(&total_prim[ran], 1);


	//hitting_prim[ran] = 0;

	prd.prim = ran;



	optix::Ray ray(ray_origin, d_mapping, radiance_ray_type, scene_epsilon); //ray_direction //d_mapping


	rtTrace(top_object, ray, prd);





	//	output_bufferl[launch_index] = ray_dot;

	countbuffer[0] = count;

	countbuffer[1] = maxcount;




}


//
// Returns solid color for miss rays
//


rtDeclareVariable(float3, bg_color, , );
RT_PROGRAM void miss()
{
	prd_radiance.result = bg_color;
	// atomicAdd(&Energy_init[output_buffer.size().x*launch_index.y + launch_index.x], (1.81 * 1000 / 24)*abs(dot(prd_radiance.normal, prd_radiance.direction)));
	//atomicAdd(&Energy_init[output_buffer.size().x*launch_index.y + launch_index.x], 1);

}

//
// (UPDATED)
// Implements basic lambertian surface shading model 
//
rtDeclareVariable(float3, Ka, , );
rtDeclareVariable(float3, Kd, , );

rtDeclareVariable(float3, ambient_light_color, , );
rtBuffer<BasicLight> lights;


RT_PROGRAM void closest_hit_radiance1()
{


	float3 hit_point = ray.origin + t_hit * ray.direction;
	float vf = 0.5;// float(hitting_prim[primid]) / float(total_prim[primid]);

	//rtPrintf("\n\n vf of prim id %f is %f", hitting_prim[6000], total_prim[6000]);

	float r = min(max(0.0f, 1.5 - abs(1 - 4 * (vf - 0.4))), 1.0f);
	float g = min(max(0.0f, 1.5 - abs(1 - 4 * (vf - 0.15))), 1.0f);
	float b = min(max(0.0f, 1.5 - abs(1 - 4 * vf)), 1.0f);


	prd_radiance.result = make_float3(r, g, b);  //Extinction[prim_index]
	//prd_radiance.result = make_float3(1.0, 0.0, 0.0);

}

RT_PROGRAM void closest_hit_radiance_object()
{

	//hitting_normal[count] = ffnormal;
	int countt = atomicAdd(&out3[0], 1);

	int id = prd_radiance.prim;


	atomicAdd(&hitting_prim[id], 1);

	prd_radiance.result = make_float3(1, 0, 0);

}


//
// Set pixel to solid color upon failure
//
RT_PROGRAM void exception()
{
	output_buffer[launch_index] = make_color(bad_color);
}

RT_PROGRAM void any_hit() // object
{




}


RT_PROGRAM void any_hit2()  // box_matl (source)
{


}

RT_PROGRAM void closest_hit_radiance_test()
{
	float3 world_geo_normal = normalize(rtTransformNormal(RT_OBJECT_TO_WORLD, geometric_normal));
	float3 world_shade_normal = normalize(rtTransformNormal(RT_OBJECT_TO_WORLD, shading_normal));
	float3 ffnormal = faceforward(world_shade_normal, -ray.direction, world_geo_normal);
	float3 color = Ka * ambient_light_color;

	float3 hit_point = ray.origin + t_hit * ray.direction;

	for (int i = 0; i < lights.size(); ++i) {
		BasicLight light = lights[i];
		float3 L = normalize(light.pos - hit_point);
		float nDl = dot(ffnormal, L);

		if (nDl > 0.0f){
			// cast shadow ray
			PerRayData_shadow shadow_prd;
			shadow_prd.attenuation = make_float3(1.0f);
			float Ldist = length(light.pos - hit_point);
			//			optix::Ray shadow_ray(hit_point, L, shadow_ray_type, scene_epsilon, Ldist);
			//			rtTrace(top_shadower, shadow_ray, shadow_prd);
			float3 light_attenuation = shadow_prd.attenuation;

			if (fmaxf(light_attenuation) > 0.0f) {
				float3 Lc = light.color * light_attenuation;
				color += Kd * nDl * Lc;

				float3 H = normalize(L - ray.direction);
				float nDh = dot(ffnormal, H);
				//if (nDh > 0)
				//	color += 0.3 * Lc * pow(nDh, 0.3);  //phong // ks 0.3
			}
		}
	}
	prd_radiance.result = color;
}



RT_PROGRAM void closest_hit_radiance_test2()
{
	float3 world_geo_normal = normalize(rtTransformNormal(RT_OBJECT_TO_WORLD, geometric_normal));
	float3 world_shade_normal = normalize(rtTransformNormal(RT_OBJECT_TO_WORLD, shading_normal));
	float3 ffnormal = faceforward(world_shade_normal, -ray.direction, world_geo_normal);
	float3 color = Ka * make_float3(0.5, 0.5, 4.0);
	float3 hit_point = ray.origin + t_hit * ray.direction;

	for (int i = 0; i < lights.size(); ++i) {
		BasicLight light = lights[i];
		float3 L = normalize(light.pos - hit_point);
		float nDl = dot(ffnormal, L);

		if (nDl > 0.0f){
			PerRayData_shadow shadow_prd;
			shadow_prd.attenuation = make_float3(1.0f);
			float Ldist = length(light.pos - hit_point);
			float3 light_attenuation = shadow_prd.attenuation;
			if (fmaxf(light_attenuation) > 0.0f)
			{
				float3 Lc = light.color * light_attenuation;
				color += Kd * nDl * Lc;
				float3 H = normalize(L - ray.direction);
				float nDh = dot(ffnormal, H);

			}
		}
	}

	prd_radiance.result = color;

}


// camera spherical
// spherical camera
/*'
double rt = sqrt(d_centers.x* d_centers.x + d_centers.y* d_centers.y);
//if ((rt*rt) <= 1)
{


double SI = rt*(M_PIf / 2);  //si max should be pi/2//
double sa = d_centers.y / rt;
double ca = d_centers.x / rt;
float3 d_spherical = normalize(sin(SI)*ca* (rx)+sin(SI)*sa* (ry)+cos(SI)*(rz));

float rt_nodes = sqrt(d_nodes.x* d_nodes.x + d_nodes.y* d_nodes.y);

// nodes
float SI_nodes = rt_nodes*(M_PI / 2);
float sa_cord = d_nodes.y / rt_nodes;
float ca_cord = d_nodes.x / rt_nodes;

float3 d_spherical_cord = eye + (sin(SI_nodes)*ca_cord* (rx)+sin(SI_nodes)*sa_cord* (ry)+cos(SI_nodes)*(rz));


//camera spherical

double lambda = d_centers.x* (M_PIf / 2);
double si = d_centers.y* (M_PIf / 2);


double phi = M_PIf - lambda;
double theta = (M_PIf / 2) - si;

float3 d_2 = normalize(sin(theta)*sin(phi)*(rx)+cos(theta)*(ry)-sin(theta)*cos(phi)*(rz)); */


