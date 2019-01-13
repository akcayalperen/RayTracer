#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include "parser.h"
#include "tinyxml2.h"
#include "tinyply.h"
#include "sstream"
#include "rtVector.h"
#include "sphere.h"
#include <fcntl.h>
#include <stdio.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "tinyexr.h"

#define RETURN_IF_ERROR(condition, error_string) \
if(condition){ \
	std::cout << error_string << std::endl; \
	return -1; \
}


using namespace std;
using namespace tinyxml2;
using namespace tinyply;

std::chrono::high_resolution_clock c;

inline std::chrono::time_point<std::chrono::high_resolution_clock> now()
{
	return c.now();
}

inline double difference_millis(std::chrono::time_point<std::chrono::high_resolution_clock> start, std::chrono::time_point<std::chrono::high_resolution_clock> end)
{
	return (double)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

void Parser::read_faces_from_ply_file(const char* file_name, Mesh* mesh, int vertex_offset, Scene* scene)
{
	try
	{
		std::ifstream ss(file_name, std::ios::binary);

		if (ss.fail())
		{
			throw std::runtime_error("failed to open ply file");
		}

		PlyFile file;

		file.parse_header(ss);

		std::cout << "================================================================\n";

		for (auto c : file.get_comments()) std::cout << "Comment: " << c << std::endl;

		for (auto e : file.get_elements())
		{
			std::cout << "element - " << e.name << " (" << e.size << ")" << std::endl;
			for (auto p : e.properties)
			{
				std::cout << "\tproperty - " << p.name << " (" << tinyply::PropertyTable[p.propertyType].str << ")" << std::endl;
			}
		}

		std::cout << "================================================================\n";

		std::shared_ptr<PlyData> vertices, uvs, normals, colors, faces, faces1, faces2, texcoords;

		try { vertices = file.request_properties_from_element("vertex", { "x", "y", "z" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		try { normals = file.request_properties_from_element("vertex", { "nx", "ny", "nz" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		try { uvs = file.request_properties_from_element("vertex", { "u", "v" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		try { colors = file.request_properties_from_element("vertex", { "red", "green", "blue", "alpha" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		try { faces1 = file.request_properties_from_element("face", { "vertex_index" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		try { faces2 = file.request_properties_from_element("face", { "vertex_indices" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		try { texcoords = file.request_properties_from_element("face", { "texcoord" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		std::chrono::time_point<std::chrono::high_resolution_clock> before = now();
		file.read(ss);
		std::chrono::time_point<std::chrono::high_resolution_clock> after = now();

		std::cout << "Parsing took " << difference_millis(before, after) << " ms: " << std::endl;
		if (vertices) std::cout << "\tRead " << vertices->count << " total vertices " << std::endl;
		if (normals) std::cout << "\tRead " << normals->count << " total vertex normals " << std::endl;
		if (uvs) std::cout << "\tRead " << uvs->count << " total vertex uvs " << std::endl;
		if (colors) std::cout << "\tRead " << colors->count << " total vertex colors " << std::endl;
		if (faces1) std::cout << "\tRead " << faces1->count << " total faces (triangles) " << std::endl;
		if (faces2) std::cout << "\tRead " << faces2->count << " total faces (triangles) " << std::endl;
		if (texcoords) std::cout << "\tRead " << texcoords->count << " total texcoords " << std::endl;

		if (faces1) faces = faces1;
		if (faces2) faces = faces2;

		if (uvs)
		{
			for (size_t i = 0; i < scene->vertex_data.size(); i++)
			{
				scene->texture_coordinate_data.push_back(Vector(0.0f));
			}

			struct float3 { float x, y, z; };
			struct float2 { float u, v; };
			const size_t numVerticesBytes = vertices->buffer.size_bytes();
			std::vector<float3> verts(vertices->count);
			std::memcpy(verts.data(), vertices->buffer.get(), numVerticesBytes);

			const size_t numNormalsBytes = normals->buffer.size_bytes();
			std::vector<float3> norms(normals->count);
			std::memcpy(norms.data(), normals->buffer.get(), numNormalsBytes);

			const size_t numUVBytes = uvs->buffer.size_bytes();
			std::vector<float2> uv(uvs->count);
			std::memcpy(uv.data(), uvs->buffer.get(), numUVBytes);

			for (int i = 0; i < verts.size(); i++)
			{
				Vector vertex = Vector(verts[i].x, verts[i].y, verts[i].z);
				scene->vertex_data.push_back(vertex);

				if (normals)
				{
					Vector normal = Vector(norms[i].x, norms[i].y, norms[i].z);
					scene->normal_data.push_back(normal);
				}
				else
				{
					scene->normal_data.push_back(Vector(0.0f));
				}

				Vector texture_coordinate = Vector(uv[i].u, uv[i].v, 0.0f);

				scene->texture_coordinate_data.push_back(texture_coordinate);
			}
		}
		else
		{
			const size_t numVerticesBytes = vertices->buffer.size_bytes();
			struct float3 { float x, y, z; };
			std::vector<float3> verts(vertices->count);
			std::memcpy(verts.data(), vertices->buffer.get(), numVerticesBytes);

			for (auto it = verts.begin(); it != verts.end(); it++)
			{
				Vector vertex = Vector((*it).x, (*it).y, (*it).z);
				scene->vertex_data.push_back(vertex);
				scene->normal_data.push_back(Vector(0.0f));
			}
		}

		if (faces1)
		{
			const size_t numFacesBytes = faces->buffer.size_bytes();
			struct int4 { int x, y, z, w; };
			std::vector<int4> faces4(faces->count);
			std::memcpy(faces4.data(), faces->buffer.get(), numFacesBytes);

			for (auto it = faces4.begin(); it != faces4.end(); it++)
			{
				Triangle triangle1;
				triangle1.v0_id = (*it).x;
				triangle1.v1_id = (*it).y;
				triangle1.v2_id = (*it).z;

				Triangle* face1 = new Triangle(triangle1.v0_id, triangle1.v1_id, triangle1.v2_id, vertex_offset);
				face1->shading_mode = mesh->shading_mode;
				face1->vertex_offset = vertex_offset;

				Vector a = scene->vertex_data[face1->v0_id + vertex_offset];
				Vector b = scene->vertex_data[face1->v1_id + vertex_offset];
				Vector c = scene->vertex_data[face1->v2_id + vertex_offset];

				Vector normal = ((b - a).cross(c - a)).normalized();

				scene->normal_data[face1->v0_id + vertex_offset] = scene->normal_data[face1->v0_id + vertex_offset] + normal;
				scene->normal_data[face1->v1_id + vertex_offset] = scene->normal_data[face1->v1_id + vertex_offset] + normal;
				scene->normal_data[face1->v1_id + vertex_offset] = scene->normal_data[face1->v1_id + vertex_offset] + normal;

				mesh->faces.push_back(face1);

				Triangle triangle2;
				triangle2.v0_id = (*it).x;
				triangle2.v1_id = (*it).z;
				triangle2.v2_id = (*it).w;

				Triangle* face2 = new Triangle(triangle2.v0_id, triangle2.v1_id, triangle2.v2_id, vertex_offset);
				face2->shading_mode = mesh->shading_mode;
				face2->vertex_offset = vertex_offset;

				a = scene->vertex_data[face2->v0_id + vertex_offset];
				b = scene->vertex_data[face2->v1_id + vertex_offset];
				c = scene->vertex_data[face2->v2_id + vertex_offset];

				normal = ((b - a).cross(c - a)).normalized();

				scene->normal_data[face2->v0_id + vertex_offset] = scene->normal_data[face2->v0_id + vertex_offset] + normal;
				scene->normal_data[face2->v1_id + vertex_offset] = scene->normal_data[face2->v1_id + vertex_offset] + normal;
				scene->normal_data[face2->v1_id + vertex_offset] = scene->normal_data[face2->v1_id + vertex_offset] + normal;

				mesh->faces.push_back(face2);
			}
		}
		
		if (faces2)
		{
			const size_t numFacesBytes = faces->buffer.size_bytes();
			struct int3 { int x, y, z; };
			std::vector<int3> faces3(faces->count);
			std::memcpy(faces3.data(), faces->buffer.get(), numFacesBytes);

			for (auto it = faces3.begin(); it != faces3.end(); it++)
			{
				Triangle triangle1;
				triangle1.v0_id = (*it).x;
				triangle1.v1_id = (*it).y;
				triangle1.v2_id = (*it).z;

				Triangle* face = new Triangle(triangle1.v0_id, triangle1.v1_id, triangle1.v2_id, vertex_offset);
				face->shading_mode = mesh->shading_mode;
				face->vertex_offset = vertex_offset;
				face->texture_offset = vertex_offset;

				if (normals)
				{

				}
				else
				{
					Vector a = scene->vertex_data[face->v0_id + vertex_offset];
					Vector b = scene->vertex_data[face->v1_id + vertex_offset];
					Vector c = scene->vertex_data[face->v2_id + vertex_offset];

					Vector normal = ((b - a).cross(c - a)).normalized();

					scene->normal_data[face->v0_id + vertex_offset] = scene->normal_data[face->v0_id + vertex_offset] + normal;
					scene->normal_data[face->v1_id + vertex_offset] = scene->normal_data[face->v1_id + vertex_offset] + normal;
					scene->normal_data[face->v1_id + vertex_offset] = scene->normal_data[face->v1_id + vertex_offset] + normal;
				}

				mesh->faces.push_back(face);
			}
		}

		mesh->construct_bvh();

	}
	catch (const std::exception & e)
	{
		std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
	}
}

int Parser::loadXmlintoScene(const char* file_name, Scene* scene)
{
	XMLDocument file;

	XMLError error = file.LoadFile(file_name);
	RETURN_IF_ERROR(error, "Cannot load XML file.");

	auto root = file.FirstChild();
	RETURN_IF_ERROR(root == nullptr, "Root node not found.");

	//Background Color
	stringstream s_background;
	auto element = root->FirstChildElement("BackgroundColor");
	if (element)
	{
		s_background << element->GetText() << endl;
		s_background >> scene->background_color.x >> scene->background_color.y >> scene->background_color.z;
	}

	//Background Texture
	stringstream s_backgroundtexture;
	element = root->FirstChildElement("BackgroundTexture");
	if (element)
	{
		scene->background_texture = new Texture();
		scene->background_texture->interpolation = "bilinear";

		s_backgroundtexture << element->GetText() << endl;
		s_backgroundtexture >> scene->background_texture->image_name;

		scene->background_texture->image_data = stbi_load(scene->background_texture->image_name.c_str(), 
			&(scene->background_texture->width), &(scene->background_texture->height), &(scene->background_texture->number_of_components), 0);
	}

	//Max Recursion Depth
	stringstream s_recursiondepth;
	element = root->FirstChildElement("MaxRecursionDepth");
	if (element)
	{
		s_recursiondepth << element->GetText() << endl;
		s_recursiondepth >> scene->max_recursion_depth;
	}

	//Shadow Ray Epsilon
	stringstream s_rayepsilon;
	element = root->FirstChildElement("ShadowRayEpsilon");
	if (element)
	{
		s_rayepsilon << element->GetText() << endl;
		s_rayepsilon >> scene->shadow_ray_epsilon;
	}

	//Intersection Test Epsilon
	stringstream s_intersectionepsilon;
	element = root->FirstChildElement("IntersectionTestEpsilon");
	if (element)
	{
		s_intersectionepsilon << element->GetText() << endl;
		s_intersectionepsilon >> scene->intersection_test_epsilon;
	}

	//Zero Based Indexing
	stringstream s_zeroindex;
	element = root->FirstChildElement("ZeroBasedIndexing");
	if(element)
	{
		std::string zero_index;

		s_zeroindex << element->GetText() << endl;
		s_zeroindex >> zero_index;

		if (zero_index == "true")
		{
			rt_scene->zero_index = true;
		}
	}

	//Integrator
	stringstream s_integrator;
	element = root->FirstChildElement("Integrator");
	if (element)
	{
		std::string integrator;

		s_integrator << element->GetText() << endl;
		s_integrator >> integrator;

		if (integrator == "PathTracing")
		{
			scene->integrator = Scene::Integrator::PathTracing;
		}
	}

	//Integrator Parameters
	stringstream s_integratorparams;
	element = root->FirstChildElement("IntegratorParams");
	if (element)
	{
		std::string param;

		s_integratorparams << element->GetText() << endl;
		s_integratorparams >> param;

		if (param == "ImportanceSampling")
		{
			scene->importance_sampling = true;
		}
	}

	//Cameras
	element = root->FirstChildElement("Cameras");
	element = element->FirstChildElement("Camera");
	while (element)
	{
		Camera camera;

		//Position
		stringstream s_position;
		auto child = element->FirstChildElement("Position");
		s_position << child->GetText() << endl;
		s_position >> camera.position.x >> camera.position.y >> camera.position.z;

		//Up
		stringstream s_up;
		child = element->FirstChildElement("Up");
		s_up << child->GetText() << endl;
		s_up >> camera.up.x >> camera.up.y >> camera.up.z;

		//Image Resolution
		stringstream s_resolution;
		child = element->FirstChildElement("ImageResolution");
		s_resolution << child->GetText() << endl;
		s_resolution >> camera.image_width >> camera.image_height;

		//Near Distance
		stringstream s_neardistance;
		child = element->FirstChildElement("NearDistance");
		s_neardistance << child->GetText() << endl;
		s_neardistance >> camera.near_distance;

		//Camera Type
		if (element->Attribute("type", "simple"))
		{
			//Gaze Point
			Vector gaze_point = Vector(0.0f);

			stringstream s_gazepoint;
			child = element->FirstChildElement("GazePoint");
			if (child)
			{
				s_gazepoint << child->GetText() << endl;
				s_gazepoint >> gaze_point.x >> gaze_point.y >> gaze_point.z;

				camera.gaze = gaze_point - camera.position;
			}

			child = element->FirstChildElement("Gaze");
			if (child)
			{
				s_gazepoint << child->GetText() << endl;
				s_gazepoint >> gaze_point.x >> gaze_point.y >> gaze_point.z;

				camera.gaze = gaze_point;
			}

			//Fov Y
			float fovy;

			stringstream s_fov;
			child = element->FirstChildElement("FovY");
			s_fov << child->GetText() << endl;
			s_fov >> fovy;

			float aspect_ratio = (float)camera.image_width / camera.image_height;

			camera.near_plane.w = tan(fovy / 2 * (PI / 180)) * camera.near_distance;
			camera.near_plane.z = -camera.near_plane.w;

			camera.near_plane.y = camera.near_plane.w * aspect_ratio;
			camera.near_plane.x = -camera.near_plane.y;
		}
		else
		{
			//Gaze
			stringstream s_gaze;
			child = element->FirstChildElement("Gaze");
			s_gaze << child->GetText() << endl;
			s_gaze >> camera.gaze.x >> camera.gaze.y >> camera.gaze.z;

			//Near Plane
			stringstream s_nearplane;
			child = element->FirstChildElement("NearPlane");
			s_nearplane << child->GetText() << endl;
			s_nearplane >> camera.near_plane.x >> camera.near_plane.y >> camera.near_plane.z >> camera.near_plane.w;
		}

		if (element->Attribute("handedness", "left"))
		{
			camera.handedness = Camera::Handedness::Left;
		}

		//Focus Distance
		stringstream s_focusdistance;
		child = element->FirstChildElement("FocusDistance");
		if (child)
		{
			s_focusdistance << child->GetText() << endl;
			s_focusdistance >> camera.focus_distance;
		}

		//Aperture Size
		stringstream s_aperturesize;
		child = element->FirstChildElement("ApertureSize");
		if (child)
		{
			s_aperturesize << child->GetText() << endl;
			s_aperturesize >> camera.aperture_size;
		}

		//Number of Samples
		stringstream s_numsamples;
		child = element->FirstChildElement("NumSamples");
		if (child)
		{
			s_numsamples << child->GetText() << endl;
			s_numsamples >> camera.num_samples;
		}

		//Image Name
		stringstream s_name;
		child = element->FirstChildElement("ImageName");
		s_name << child->GetText() << endl;
		s_name >> camera.image_name;
		camera.image_name = camera.image_name.substr(0, camera.image_name.find_last_of("."));

		//Tonemap
		auto tonemap = element->FirstChildElement("Tonemap");
		if (tonemap)
		{
			camera.is_hdr = true;

			//TMO
			stringstream s_tmo;
			child = tonemap->FirstChildElement("TMO");
			s_tmo << child->GetText() << endl;
			s_tmo >> camera.tmo;

			//TMO Options
			stringstream s_options;
			child = tonemap->FirstChildElement("TMOOptions");
			s_options << child->GetText() << endl;
			s_options >> camera.image_key >> camera.burn_percentage;

			//Saturation
			stringstream s_saturation;
			child = tonemap->FirstChildElement("Saturation");
			s_saturation << child->GetText() << endl;
			s_saturation >> camera.saturation;

			//Gamma
			stringstream s_gamma;
			child = tonemap->FirstChildElement("Gamma");
			if (child)
			{
				s_gamma << child->GetText() << endl;
				s_gamma >> camera.gamma;
			}
		}

		//Gamma Correction
		stringstream s_gammacorrection;
		child = element->FirstChildElement("GammaCorrection");
		if(child)
		{
			s_gammacorrection << child->GetText() << endl;
			s_gammacorrection >> camera.gamma_correction_method;
		}

		scene->cameras.push_back(camera);
		element = element->NextSiblingElement("Camera");
	}

	//Lights
	element = root->FirstChildElement("Lights");
	if (element)
	{
		//Ambient Light
		stringstream s_ambientlight;
		auto child = element->FirstChildElement("AmbientLight");
		if (child)
		{
			s_ambientlight << child->GetText() << std::endl;
			s_ambientlight >> scene->ambient_light.x >> scene->ambient_light.y >> scene->ambient_light.z;
		}

		//Point Lights
		auto light = element->FirstChildElement("PointLight");
		while (light)
		{
			PointLight *point_light = new PointLight;

			//Position
			stringstream s_position;
			child = light->FirstChildElement("Position");
			s_position << child->GetText() << endl;
			s_position >> point_light->position.x >> point_light->position.y >> point_light->position.z;

			//Intensity
			stringstream s_intensity;
			child = light->FirstChildElement("Intensity");
			s_intensity << child->GetText() << endl;
			s_intensity >> point_light->intensity.x >> point_light->intensity.y >> point_light->intensity.z;

			scene->lights.push_back(point_light);
			light = light->NextSiblingElement("PointLight");
		}

		//Area Lights
		light = element->FirstChildElement("AreaLight");
		while (light)
		{
			AreaLight *area_light = new AreaLight;

			//Position
			stringstream s_position;
			child = light->FirstChildElement("Position");
			s_position << child->GetText() << endl;
			s_position >> area_light->position.x >> area_light->position.y >> area_light->position.z;

			//Intensity
			stringstream s_intensity;
			child = light->FirstChildElement("Intensity");
			s_intensity << child->GetText() << endl;
			s_intensity >> area_light->intensity.x >> area_light->intensity.y >> area_light->intensity.z;
			
			//Edge Vector 1
			stringstream s_edge1;
			child = light->FirstChildElement("EdgeVector1");
			s_edge1 << child->GetText() << endl;
			s_edge1 >> area_light->edge_vector_1.x >> area_light->edge_vector_1.y >> area_light->edge_vector_1.z;

			//Edge Vector 2
			stringstream s_edge2;
			child = light->FirstChildElement("EdgeVector2");
			s_edge2 << child->GetText() << endl;
			s_edge2 >> area_light->edge_vector_2.x >> area_light->edge_vector_2.y >> area_light->edge_vector_2.z;

			//Normal
			area_light->normal = (area_light->edge_vector_1).cross(area_light->edge_vector_2).normalized();

			scene->lights.push_back(area_light);
			light = light->NextSiblingElement("AreaLight");
		}

		//Directional Lights
		light = element->FirstChildElement("DirectionalLight");
		while (light)
		{
			DirectionalLight *directional_light = new DirectionalLight;

			//Direction
			stringstream s_direction;
			child = light->FirstChildElement("Direction");
			s_direction << child->GetText() << endl;
			s_direction >> directional_light->direction.x >> directional_light->direction.y >> directional_light->direction.z;

			//Radiance
			stringstream s_radiance;
			child = light->FirstChildElement("Radiance");
			s_radiance << child->GetText() << endl;
			s_radiance >> directional_light->radiance.x >> directional_light->radiance.y >> directional_light->radiance.z;

			directional_light->position = directional_light->direction * (-100000);

			scene->lights.push_back(directional_light);
			light = light->NextSiblingElement("DirectionalLight");
		}

		//Spot Lights
		light = element->FirstChildElement("SpotLight");
		while (light)
		{
			SpotLight *spot_light = new SpotLight;

			//Position
			stringstream s_position;
			child = light->FirstChildElement("Position");
			s_position << child->GetText() << endl;
			s_position >> spot_light->position.x >> spot_light->position.y >> spot_light->position.z;

			//Direction
			stringstream s_direction;
			child = light->FirstChildElement("Direction");
			s_direction << child->GetText() << endl;
			s_direction >> spot_light->direction.x >> spot_light->direction.y >> spot_light->direction.z;
			spot_light->direction = spot_light->direction.normalized();

			//Intensity
			stringstream s_intensity;
			child = light->FirstChildElement("Intensity");
			s_intensity << child->GetText() << endl;
			s_intensity >> spot_light->intensity.x >> spot_light->intensity.y >> spot_light->intensity.z;

			//Coverage Angle
			stringstream s_coverage;
			child = light->FirstChildElement("CoverageAngle");
			s_coverage << child->GetText() << endl;
			s_coverage >> spot_light->coverage_angle;

			//Falloff Angle
			stringstream s_falloff;
			child = light->FirstChildElement("FalloffAngle");
			s_falloff << child->GetText() << endl;
			s_falloff >> spot_light->falloff_angle;

			scene->lights.push_back(spot_light);
			light = light->NextSiblingElement("SpotLight");
		}

		//Spherical Directional Lights
		light = element->FirstChildElement("SphericalDirectionalLight");
		if(light)
		{
			Texture *environment_map = new Texture;

			//Environment Map Name
			stringstream s_envmap;
			child = light->FirstChildElement("EnvMapName");
			s_envmap << child->GetText() << endl;
			s_envmap >> environment_map->image_name;

			const char* err;

			int ret = LoadEXR(&(environment_map->radiance_data), &environment_map->width, &environment_map->height, environment_map->image_name.c_str(), &err);
			RETURN_IF_ERROR(ret != 0, "Cannot load EXR file.");

			environment_map->number_of_components = 4;
			environment_map->texture_type = Texture::TextureType::Image;
			environment_map->interpolation = "bilinear";

			rt_scene->environment_map = environment_map;
		}
	}

	//BRDFs
	element = root->FirstChildElement("BRDFs");
	if (element)
	{
		//Original Phong
		auto child = element->FirstChildElement("OriginalPhong");
		while (child)
		{
			OriginalPhongBRDF *brdf = new OriginalPhongBRDF;

			//Exponent
			auto exponent = child->FirstChildElement("Exponent");

			stringstream s_exponent;
			s_exponent << exponent->GetText() << endl;
			s_exponent >> brdf->exponent;

			scene->brdfs.push_back(brdf);
			child = child->NextSiblingElement("OriginalPhong");
		}

		//Modified Phong
		child = element->FirstChildElement("ModifiedPhong");
		while (child)
		{
			ModifiedPhongBRDF *brdf = new ModifiedPhongBRDF;

			if (child->Attribute("normalized", "true"))
			{
				brdf->normalizer = 2;
			}

			//Exponent
			auto exponent = child->FirstChildElement("Exponent");

			stringstream s_exponent;
			s_exponent << exponent->GetText() << endl;
			s_exponent >> brdf->exponent;

			scene->brdfs.push_back(brdf);
			child = child->NextSiblingElement("ModifiedPhong");
		}

		//Original Blinn-Phong
		child = element->FirstChildElement("OriginalBlinnPhong");
		while (child)
		{
			OriginalBlinnPhongBRDF *brdf = new OriginalBlinnPhongBRDF;

			//Exponent
			auto exponent = child->FirstChildElement("Exponent");

			stringstream s_exponent;
			s_exponent << exponent->GetText() << endl;
			s_exponent >> brdf->exponent;

			scene->brdfs.push_back(brdf);
			child = child->NextSiblingElement("OriginalBlinnPhong");
		}

		//Modified Blinn-Phong
		child = element->FirstChildElement("ModifiedBlinnPhong");
		while (child)
		{
			ModifiedBlinnPhongBRDF *brdf = new ModifiedBlinnPhongBRDF;

			if (child->Attribute("normalized", "true"))
			{
				brdf->normalizer = 8;
			}

			//Exponent
			auto exponent = child->FirstChildElement("Exponent");

			stringstream s_exponent;
			s_exponent << exponent->GetText() << endl;
			s_exponent >> brdf->exponent;

			scene->brdfs.push_back(brdf);
			child = child->NextSiblingElement("ModifiedBlinnPhong");
		}

		//Modified Blinn-Phong
		child = element->FirstChildElement("TorranceSparrow");
		while (child)
		{
			TorranceSparrowBRDF *brdf = new TorranceSparrowBRDF;

			//Exponent
			auto exponent = child->FirstChildElement("Exponent");

			stringstream s_exponent;
			s_exponent << exponent->GetText() << endl;
			s_exponent >> brdf->exponent;

			//Refractive Index
			auto refractiveindex = child->FirstChildElement("RefractiveIndex");

			stringstream s_index;
			s_index << refractiveindex->GetText() << endl;
			s_index >> brdf->refractive_index;

			scene->brdfs.push_back(brdf);
			child = child->NextSiblingElement("TorranceSparrow");
		}
	}

	//Transformations
	element = root->FirstChildElement("Transformations");
	if (element)
	{
		//Translation
		auto child = element->FirstChildElement("Translation");
		while (child)
		{
			Vector translation_vector = Vector();

			stringstream s_translation;
			s_translation << child->GetText() << endl;
			s_translation >> translation_vector.x >> translation_vector.y >> translation_vector.z;

			scene->translations.push_back(Matrix::TranslationMatrix(translation_vector));
			child = child->NextSiblingElement("Translation");
		}

		//Scaling
		child = element->FirstChildElement("Scaling");
		while (child)
		{
			Vector scaling_vector = Vector();

			stringstream s_scaling;
			s_scaling << child->GetText() << endl;
			s_scaling >> scaling_vector.x >> scaling_vector.y >> scaling_vector.z;

			scene->scalings.push_back(Matrix::ScaleMatrix(scaling_vector));
			child = child->NextSiblingElement("Scaling");
		}

		//Rotation
		child = element->FirstChildElement("Rotation");
		while (child)
		{
			float angle;
			Vector rotation_vector = Vector();

			stringstream s_rotation;
			s_rotation << child->GetText() << endl;
			s_rotation >> angle >> rotation_vector.x >> rotation_vector.y >> rotation_vector.z;

			scene->rotations.push_back(Matrix::RotationMatrix(rotation_vector, angle));
			child = child->NextSiblingElement("Rotation");
		}
	}

	//Materials
	element = root->FirstChildElement("Materials");
	element = element->FirstChildElement("Material");
	while (element)
	{
		Material material;

		material.brdf_id = element->IntAttribute("BRDF") - 1;

		//Ambient Reflectance
		stringstream s_ambient;
		auto child = element->FirstChildElement("AmbientReflectance");
		if (child)
		{
			s_ambient << child->GetText() << endl;
			s_ambient >> material.ambient.x >> material.ambient.y >> material.ambient.z;
		}

		//Diffuse Reflectance
		stringstream s_diffuse;
		child = element->FirstChildElement("DiffuseReflectance");
		if (child)
		{
			s_diffuse << child->GetText() << endl;
			s_diffuse >> material.diffuse.x >> material.diffuse.y >> material.diffuse.z;
		}

		//Specular Reflectance
		stringstream s_specular;
		child = element->FirstChildElement("SpecularReflectance");
		if (child)
		{
			s_specular << child->GetText() << endl;
			s_specular >> material.specular.x >> material.specular.y >> material.specular.z;
		}

		//Mirror Reflectance
		stringstream s_mirror;
		child = element->FirstChildElement("MirrorReflectance");
		if (child)
		{
			s_mirror << child->GetText() << endl;
			s_mirror >> material.mirror.x >> material.mirror.y >> material.mirror.z;
		}

		//Roughness
		stringstream s_roughness;
		child = element->FirstChildElement("Roughness");
		if (child)
		{
			s_roughness << child->GetText() << endl;
			s_roughness >> material.roughness;
		}

		//Transparency
		stringstream s_transparency;
		child = element->FirstChildElement("Transparency");
		if (child)
		{
			s_transparency << child->GetText() << endl;
			s_transparency >> material.transparency.x >> material.transparency.y >> material.transparency.z;
		}

		//Degamma
		if (element->Attribute("degamma", "true"))
		{
			material.ambient.x = std::pow(material.ambient.x, 2.2);
			material.ambient.y = std::pow(material.ambient.y, 2.2);
			material.ambient.z = std::pow(material.ambient.z, 2.2);

			material.diffuse.x = std::pow(material.diffuse.x, 2.2);
			material.diffuse.y = std::pow(material.diffuse.y, 2.2);
			material.diffuse.z = std::pow(material.diffuse.z, 2.2);

			material.specular.x = std::pow(material.specular.x, 2.2);
			material.specular.y = std::pow(material.specular.y, 2.2);
			material.specular.z = std::pow(material.specular.z, 2.2);
		}

		stringstream s_refractionindex;
		child = element->FirstChildElement("RefractionIndex");
		if (child)
		{
			s_refractionindex << child->GetText() << endl;
			s_refractionindex >> material.refraction_index;
		}

		//Phong Exponent
		stringstream s_phong;
		child = element->FirstChildElement("PhongExponent");
		if (child)
		{
			s_phong << child->GetText() << endl;
			s_phong >> material.phong_exponent;
		}

		scene->materials.push_back(material);
		element = element->NextSiblingElement("Material");
	}

	//Textures
	element = root->FirstChildElement("Textures");
	if (element)
	{
		element = element->FirstChildElement("Texture");
		while (element)
		{
			Texture texture;

			//Image Name
			stringstream s_imagename;
			auto child = element->FirstChildElement("ImageName");
			if (child)
			{
				s_imagename << child->GetText() << endl;
				s_imagename >> texture.image_name;
			}

			//Bump Map
			if (element->Attribute("bumpmap", "true"))
			{
				texture.bump_map = true;
			}

			//Bump Multiplier
			texture.bump_map_multiplier = element->FloatAttribute("bumpmapMultiplier");

			if (texture.bump_map_multiplier == 0.0f)
			{
				texture.bump_map_multiplier = 1.0f;
			}

			//Interpolation
			stringstream s_interpolation;
			child = element->FirstChildElement("Interpolation");
			if (child)
			{
				s_interpolation << child->GetText() << endl;
				s_interpolation >> texture.interpolation;
			}

			//Decal Mode
			stringstream s_decalmode;
			child = element->FirstChildElement("DecalMode");
			if (child)
			{
				s_decalmode << child->GetText() << endl;
				s_decalmode >> texture.decal_mode_string;
				texture.set_decal_mode(texture.decal_mode_string);
			}

			//Normalizer
			stringstream s_normalizer;
			child = element->FirstChildElement("Normalizer");
			if (child)
			{
				s_normalizer << child->GetText() << endl;
				s_normalizer >> texture.normalizer;
			}

			//Appearance
			stringstream s_appearance;
			child = element->FirstChildElement("Appearance");
			if (child)
			{
				s_appearance << child->GetText() << endl;
				s_appearance >> texture.appearance_string;
				texture.set_appearance(texture.appearance_string);
			}

			//Perlin
			if (texture.image_name == "perlin")
			{
				//Scaling Factor
				stringstream s_scaling;
				child = element->FirstChildElement("ScalingFactor");
				if (child)
				{
					s_scaling << child->GetText() << endl;
					s_scaling >> texture.scaling_factor;
				}
				texture.texture_type = Texture::TextureType::Perlin;
			}
			else
			{
				texture.image_data = stbi_load(texture.image_name.c_str(), &texture.width, &texture.height, &texture.number_of_components, 0);
				texture.texture_type = Texture::TextureType::Image;
			}

			if (element->Attribute("degamma", "true"))
			{
				texture.degamma = true;
			}

			scene->textures.push_back(texture);
			element = element->NextSiblingElement("Texture");
		}
	}


	//VertexData
	element = root->FirstChildElement("VertexData");
	if (element)
	{
		if (element->Attribute("binaryFile"))
		{
			const char* file_name = element->Attribute("binaryFile");

			int N = 0;
			ifstream file(file_name, std::ios::binary);

			file.read(reinterpret_cast<char*>(&N), sizeof(int));

			Vector vertex;

			for (int i = 0; i < N; i++)
			{
				file.read(reinterpret_cast<char*>(&(vertex.x)), sizeof(float));
				file.read(reinterpret_cast<char*>(&(vertex.y)), sizeof(float));
				file.read(reinterpret_cast<char*>(&(vertex.z)), sizeof(float));

				scene->vertex_data.push_back(vertex);
				scene->normal_data.push_back(Vector(0.0f));
			}

		}
		else
		{
			stringstream s_vertex;
			s_vertex << element->GetText() << endl;
			Vector vertex;
			while (!(s_vertex >> vertex.x).eof())
			{
				s_vertex >> vertex.y >> vertex.z;
				scene->vertex_data.push_back(vertex);
				scene->normal_data.push_back(Vector(0.0f));
			}
			s_vertex.clear();
		}
	}


	//Texture Coordinates
	element = root->FirstChildElement("TexCoordData");
	if (element)
	{
		if (element->Attribute("binaryFile"))
		{
			const char* file_name = element->Attribute("binaryFile");

			int N = 0;
			ifstream file(file_name, std::ios::binary);

			file.read(reinterpret_cast<char*>(&N), sizeof(int));

			Vector texture_coordinate;

			for (int i = 0; i < N; i++)
			{
				file.read(reinterpret_cast<char*>(&(texture_coordinate.x)), sizeof(float));
				file.read(reinterpret_cast<char*>(&(texture_coordinate.y)), sizeof(float));

				scene->texture_coordinate_data.push_back(texture_coordinate);
			}

		}
		else
		{
			stringstream s_texcoord;
			s_texcoord << element->GetText() << endl;
			Vector texture_coordinate;
			while (!(s_texcoord >> texture_coordinate.x).eof())
			{
				s_texcoord >> texture_coordinate.y;
				scene->texture_coordinate_data.push_back(texture_coordinate);
			}
			s_texcoord.clear();
		}
	}
	

	//Meshes
	element = root->FirstChildElement("Objects");
	element = element->FirstChildElement("Mesh");
	while (element)
	{
		Mesh* mesh = new Mesh;

		//Shading Mode
		if (element->Attribute("shadingMode", "smooth"))
		{
			mesh->shading_mode = Object::ShadingMode::Smooth;
		}

		//Material
		stringstream s_material;
		auto child = element->FirstChildElement("Material");
		s_material << child->GetText() << endl;
		s_material >> mesh->material_id;
		mesh->material_id--;

		//Texture
		stringstream s_texture;
		child = element->FirstChildElement("Texture");
		if (child)
		{
			s_texture << child->GetText() << endl;
			s_texture >> mesh->texture_id;
			mesh->texture_id--;
		}

		//Transformations
		stringstream s_transformation;
		child = element->FirstChildElement("Transformations");
		if (child)
		{
			s_transformation << child->GetText() << std::endl;
		}

		std::string token;
		while (!(s_transformation >> token).eof())
		{
			const char * token_arr = token.c_str();
			char type = token_arr[0];
			token_arr++;
			int idx = atoi(token_arr) - 1;
			switch (type)
			{
			case 'r':
				mesh->base_transformation = scene->rotations[idx] * mesh->base_transformation;
				break;
			case 's':
				mesh->base_transformation = scene->scalings[idx] * mesh->base_transformation;
				break;
			case 't':
				mesh->base_transformation = scene->translations[idx] * mesh->base_transformation;
				break;
			default:
				break;
			}
		}

		//Motion Blur
		stringstream s_blur;
		child = element->FirstChildElement("MotionBlur");
		if (child)
		{
			s_blur << child->GetText() << endl;
			s_blur >> mesh->motion_vector.x >> mesh->motion_vector.y >> mesh->motion_vector.z;
		}

		//Faces
		stringstream s_face;
		child = element->FirstChildElement("Faces");
		int vertex_offset = child->IntAttribute("vertexOffset");
		int texture_offset = child->IntAttribute("textureOffset");

		//Ply File
		auto ply_file = child->Attribute("plyFile");
		auto binary_file = child->Attribute("binaryFile");
		if (ply_file != nullptr)
		{
			vertex_offset = scene->vertex_data.size();
			read_faces_from_ply_file(ply_file, mesh, vertex_offset, scene);
		}
		else if(binary_file != nullptr)
		{
			int N = 0;
			ifstream file(binary_file, std::ios::binary);

			file.read(reinterpret_cast<char*>(&N), sizeof(int));

			Triangle triangle;

			for (int i = 0; i < N; i++)
			{
				file.read(reinterpret_cast<char*>(&(triangle.v0_id)), sizeof(int));
				file.read(reinterpret_cast<char*>(&(triangle.v1_id)), sizeof(int));
				file.read(reinterpret_cast<char*>(&(triangle.v2_id)), sizeof(int));
				/*if (!scene->zero_index)
				{
					triangle.v0_id--;
					triangle.v1_id--;
					triangle.v2_id--;
				}*/

				Triangle* face = new Triangle(triangle.v0_id, triangle.v1_id, triangle.v2_id, vertex_offset);
				face->vertex_offset = vertex_offset;
				face->texture_offset = texture_offset;

				face->shading_mode = mesh->shading_mode;

				Vector a = scene->vertex_data[face->v0_id + vertex_offset];
				Vector b = scene->vertex_data[face->v1_id + vertex_offset];
				Vector c = scene->vertex_data[face->v2_id + vertex_offset];

				Vector normal = ((b - a).cross(c - a)).normalized();

				scene->normal_data[face->v0_id + vertex_offset] = scene->normal_data[face->v0_id + vertex_offset] + normal;
				scene->normal_data[face->v1_id + vertex_offset] = scene->normal_data[face->v1_id + vertex_offset] + normal;
				scene->normal_data[face->v2_id + vertex_offset] = scene->normal_data[face->v2_id + vertex_offset] + normal;

				mesh->faces.push_back(face);
			}

			mesh->construct_bvh();
		}
		else
		{
			s_face << child->GetText() << endl;
			Triangle triangle;
			while (!(s_face >> triangle.v0_id).eof())
			{
				s_face >> triangle.v1_id >> triangle.v2_id;
				triangle.v0_id--;
				triangle.v1_id--;
				triangle.v2_id--;

				Triangle* face = new Triangle(triangle.v0_id, triangle.v1_id, triangle.v2_id, vertex_offset);
				face->vertex_offset = vertex_offset;
				face->texture_offset = texture_offset;

				face->shading_mode = mesh->shading_mode;

				Vector a = scene->vertex_data[face->v0_id + vertex_offset];
				Vector b = scene->vertex_data[face->v1_id + vertex_offset];
				Vector c = scene->vertex_data[face->v2_id + vertex_offset];

				Vector normal = ((b - a).cross(c - a)).normalized();

				scene->normal_data[face->v0_id + vertex_offset] = scene->normal_data[face->v0_id + vertex_offset] + normal;
				scene->normal_data[face->v1_id + vertex_offset] = scene->normal_data[face->v1_id + vertex_offset] + normal;
				scene->normal_data[face->v2_id + vertex_offset] = scene->normal_data[face->v2_id + vertex_offset] + normal;

				mesh->faces.push_back(face);
			}
			s_face.clear();

			mesh->construct_bvh();
		}
		
		if (mesh->texture_id != -1)
		{
			if (scene->textures[mesh->texture_id].bump_map && scene->textures[mesh->texture_id].texture_type == Texture::TextureType::Image)
			{
				for (auto it = mesh->faces.begin(); it != mesh->faces.end(); it++)
				{
					((Triangle*)(*it))->calculate_dp();
				}
			}
		}

		scene->meshes.push_back(mesh);

		MeshInstance* mesh_instance = new MeshInstance;
		mesh_instance->base_mesh_id = scene->meshes.size() - 1;
		mesh_instance->material_id = mesh->material_id;
		mesh_instance->texture_id = mesh->texture_id;
		mesh_instance->shading_mode = mesh->shading_mode;
		mesh_instance->base_transformation = mesh->base_transformation;
		mesh_instance->motion_vector = mesh->motion_vector;
		mesh_instance->bbox = mesh->bbox;
		mesh_instance->bbox.transformBBox(mesh_instance->base_transformation, mesh_instance->motion_vector);

		mesh_instance->inverse_transformation = mesh_instance->base_transformation.inverse();
		mesh_instance->transpose_transformation = (mesh_instance->inverse_transformation).transpose();

		scene->objects.push_back(mesh_instance);


		element = element->NextSiblingElement("Mesh");
	}

	for (auto it = scene->normal_data.begin(); it != scene->normal_data.end(); it++)
	{
		(*it) = (*it).normalized();
	}

	//Mesh Instances
	element = root->FirstChildElement("Objects");
	element = element->FirstChildElement("MeshInstance");
	while (element)
	{
		MeshInstance* mesh_instance = new MeshInstance;

		//Base Mesh
		mesh_instance->base_mesh_id = element->IntAttribute("baseMeshId") - 1;

		//Material
		stringstream s_material;
		auto child = element->FirstChildElement("Material");
		s_material << child->GetText() << endl;
		s_material >> mesh_instance->material_id;
		mesh_instance->material_id--;

		//Texture
		stringstream s_texture;
		child = element->FirstChildElement("Texture");
		if (child)
		{
			s_texture << child->GetText() << endl;
			s_texture >> mesh_instance->texture_id;
			mesh_instance->texture_id--;
		}

		//Motion Blur
		stringstream s_blur;
		child = element->FirstChildElement("MotionBlur");
		if (child)
		{
			s_blur << child->GetText() << endl;
			s_blur >> mesh_instance->motion_vector.x >> mesh_instance->motion_vector.y >> mesh_instance->motion_vector.z;
		}

		//Reset Transform
		if (!element->Attribute("resetTransform", "true"))
		{
			mesh_instance->base_transformation = scene->meshes[mesh_instance->base_mesh_id]->base_transformation;
		}
		else
		{
			mesh_instance->base_transformation = Matrix::Identity();
		}

		//Transformations
		stringstream s_transformation;
		child = element->FirstChildElement("Transformations");
		if (child)
		{
			s_transformation << child->GetText() << std::endl;
		}

		std::string token;
		while (!(s_transformation >> token).eof())
		{
			const char * token_arr = token.c_str();
			char type = token_arr[0];
			token_arr++;
			int idx = atoi(token_arr) - 1;
			switch (type)
			{
			case 'r':
				mesh_instance->base_transformation = scene->rotations[idx] * mesh_instance->base_transformation;
				break;
			case 's':
				mesh_instance->base_transformation = scene->scalings[idx] * mesh_instance->base_transformation;
				break;
			case 't':
				mesh_instance->base_transformation = scene->translations[idx] * mesh_instance->base_transformation;
				break;
			default:
				break;
			}
		}

		mesh_instance->shading_mode = scene->meshes[mesh_instance->base_mesh_id]->shading_mode;
		mesh_instance->motion_vector = mesh_instance->motion_vector;
		mesh_instance->bbox = scene->meshes[mesh_instance->base_mesh_id]->bbox;
		mesh_instance->bbox.transformBBox(mesh_instance->base_transformation, mesh_instance->motion_vector);

		mesh_instance->inverse_transformation = mesh_instance->base_transformation.inverse();
		mesh_instance->transpose_transformation = (mesh_instance->inverse_transformation).transpose();

		scene->objects.push_back(mesh_instance);

		element = element->NextSiblingElement("MeshInstance");
	}

	//Light Meshes
	element = root->FirstChildElement("Objects");
	element = element->FirstChildElement("LightMesh");
	while (element)
	{
		LightMesh* light_mesh = new LightMesh;

		//Material
		stringstream s_material;
		auto child = element->FirstChildElement("Material");
		s_material << child->GetText() << endl;
		s_material >> light_mesh->material_id;
		light_mesh->material_id--;

		//Radiance
		stringstream s_radiance;
		child = element->FirstChildElement("Radiance");
		s_radiance << child->GetText() << endl;
		s_radiance >> light_mesh->radiance.x >> light_mesh->radiance.y >> light_mesh->radiance.z;

		//Transformations
		stringstream s_transformation;
		child = element->FirstChildElement("Transformations");
		if (child)
		{
			s_transformation << child->GetText() << std::endl;
		}

		std::string token;
		while (!(s_transformation >> token).eof())
		{
			const char * token_arr = token.c_str();
			char type = token_arr[0];
			token_arr++;
			int idx = atoi(token_arr) - 1;
			switch (type)
			{
			case 'r':
				light_mesh->base_transformation = scene->rotations[idx] * light_mesh->base_transformation;
				break;
			case 's':
				light_mesh->base_transformation = scene->scalings[idx] * light_mesh->base_transformation;
				break;
			case 't':
				light_mesh->base_transformation = scene->translations[idx] * light_mesh->base_transformation;
				break;
			default:
				break;
			}
		}

		//Faces
		stringstream s_face;
		child = element->FirstChildElement("Faces");
		int vertex_offset = child->IntAttribute("vertexOffset");
		int texture_offset = child->IntAttribute("textureOffset");

		s_face << child->GetText() << endl;
		Triangle triangle;
		while (!(s_face >> triangle.v0_id).eof())
		{
			s_face >> triangle.v1_id >> triangle.v2_id;
			triangle.v0_id--;
			triangle.v1_id--;
			triangle.v2_id--;

			Triangle* face = new Triangle(triangle.v0_id, triangle.v1_id, triangle.v2_id, vertex_offset);
			face->vertex_offset = vertex_offset;
			face->texture_offset = texture_offset;

			face->shading_mode = light_mesh->shading_mode;

			Vector a = scene->vertex_data[face->v0_id + vertex_offset];
			Vector b = scene->vertex_data[face->v1_id + vertex_offset];
			Vector c = scene->vertex_data[face->v2_id + vertex_offset];

			Vector normal = ((b - a).cross(c - a)).normalized();

			scene->normal_data[face->v0_id + vertex_offset] = scene->normal_data[face->v0_id + vertex_offset] + normal;
			scene->normal_data[face->v1_id + vertex_offset] = scene->normal_data[face->v1_id + vertex_offset] + normal;
			scene->normal_data[face->v1_id + vertex_offset] = scene->normal_data[face->v1_id + vertex_offset] + normal;

			light_mesh->faces.push_back(face);
			light_mesh->areas.push_back(0);
		}
		s_face.clear();

		light_mesh->construct_bvh();

		scene->meshes.push_back(light_mesh);

		MeshInstance* mesh_instance = new MeshInstance;
		mesh_instance->base_mesh_id = scene->meshes.size() - 1;
		mesh_instance->material_id = light_mesh->material_id;
		mesh_instance->texture_id = light_mesh->texture_id;
		mesh_instance->shading_mode = light_mesh->shading_mode;
		mesh_instance->base_transformation = light_mesh->base_transformation;
		mesh_instance->motion_vector = light_mesh->motion_vector;
		mesh_instance->bbox = light_mesh->bbox;
		mesh_instance->bbox.transformBBox(mesh_instance->base_transformation, mesh_instance->motion_vector);

		mesh_instance->inverse_transformation = mesh_instance->base_transformation.inverse();
		mesh_instance->transpose_transformation = (mesh_instance->inverse_transformation).transpose();

		scene->objects.push_back(mesh_instance);
		scene->lights.push_back(light_mesh);

		light_mesh->calculate_area();

		element = element->NextSiblingElement("LightMesh");
	}

	//Light Spheres
	element = root->FirstChildElement("Objects");
	element = element->FirstChildElement("LightSphere");
	while (element)
	{
		LightSphere* light_sphere = new LightSphere;

		//Material
		stringstream s_material;
		auto child = element->FirstChildElement("Material");
		s_material << child->GetText() << endl;
		s_material >> light_sphere->material_id;
		light_sphere->material_id--;

		//Transformations
		stringstream s_transformation;
		child = element->FirstChildElement("Transformations");
		if (child)
		{
			s_transformation << child->GetText() << std::endl;
		}

		std::string token;
		while (!(s_transformation >> token).eof())
		{
			const char * token_arr = token.c_str();
			char type = token_arr[0];
			token_arr++;
			int idx = atoi(token_arr) - 1;
			switch (type)
			{
			case 'r':
				light_sphere->base_transformation = scene->rotations[idx] * light_sphere->base_transformation;
				break;
			case 's':
				light_sphere->base_transformation = scene->scalings[idx] * light_sphere->base_transformation;
				break;
			case 't':
				light_sphere->base_transformation = scene->translations[idx] * light_sphere->base_transformation;
				break;
			default:
				break;
			}
		}

		light_sphere->inverse_transformation = light_sphere->base_transformation.inverse();
		light_sphere->transpose_transformation = (light_sphere->inverse_transformation).transpose();

		//Center
		stringstream s_center;
		child = element->FirstChildElement("Center");
		s_center << child->GetText() << endl;
		s_center >> light_sphere->center_vertex_id;
		light_sphere->center_vertex_id--;

		//Radius
		stringstream s_radius;
		child = element->FirstChildElement("Radius");
		s_radius << child->GetText() << endl;
		s_radius >> light_sphere->radius;


		//Radiance
		stringstream s_radiance;
		child = element->FirstChildElement("Radiance");
		s_radiance << child->GetText() << endl;
		s_radiance >> light_sphere->radiance.x >> light_sphere->radiance.y >> light_sphere->radiance.z;


		scene->spheres.push_back(light_sphere);
		scene->lights.push_back(light_sphere);
		element = element->NextSiblingElement("LightSphere");
	}

	//Spheres
	element = root->FirstChildElement("Objects");
	element = element->FirstChildElement("Sphere");
	while (element)
	{
		Sphere* sphere = new Sphere;

		//Material
		stringstream s_material;
		auto child = element->FirstChildElement("Material");
		s_material << child->GetText() << endl;
		s_material >> sphere->material_id;
		sphere->material_id--;

		//Transformations
		stringstream s_transformation;
		child = element->FirstChildElement("Transformations");
		if (child)
		{
			s_transformation << child->GetText() << std::endl;
		}

		std::string token;
		while (!(s_transformation >> token).eof())
		{
			const char * token_arr = token.c_str();
			char type = token_arr[0];
			token_arr++;
			int idx = atoi(token_arr) - 1;
			switch (type)
			{
			case 'r':
				sphere->base_transformation = scene->rotations[idx] * sphere->base_transformation;
				break;
			case 's':
				sphere->base_transformation = scene->scalings[idx] * sphere->base_transformation;
				break;
			case 't':
				sphere->base_transformation = scene->translations[idx] * sphere->base_transformation;
				break;
			default:
				break;
			}
		}

		sphere->inverse_transformation = sphere->base_transformation.inverse();
		sphere->transpose_transformation = (sphere->inverse_transformation).transpose();

		//Texture
		stringstream s_texture;
		child = element->FirstChildElement("Texture");
		if (child)
		{
			s_texture << child->GetText() << endl;
			s_texture >> sphere->texture_id;
			sphere->texture_id--;
		}

		//Center
		stringstream s_center;
		child = element->FirstChildElement("Center");
		s_center << child->GetText() << endl;
		s_center >> sphere->center_vertex_id;
		sphere->center_vertex_id--;

		//Radius
		stringstream s_radius;
		child = element->FirstChildElement("Radius");
		s_radius << child->GetText() << endl;
		s_radius >> sphere->radius;

		//Motion Blur
		stringstream s_blur;
		child = element->FirstChildElement("MotionBlur");
		if (child)
		{
			s_blur << child->GetText() << endl;
			s_blur >> sphere->motion_vector.x >> sphere->motion_vector.y >> sphere->motion_vector.z;
		}

		scene->spheres.push_back(sphere);
		element = element->NextSiblingElement("Sphere");
	}

	
	return 0;
}