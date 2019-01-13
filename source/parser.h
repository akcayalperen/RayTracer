#pragma once
#ifndef __PARSER__H_
#define __PARSER__H_

#include <string>
#include "scene.h"
#include "mesh.h"

using namespace std;

namespace Parser
{
	int loadXmlintoScene(const char *file_name, Scene* scene);
	void read_faces_from_ply_file(const char* file_name, Mesh* mesh, int vertex_offset, Scene* scene);

};


#endif