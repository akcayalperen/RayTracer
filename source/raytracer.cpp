#include <iostream>
#include "parser.h"
#include "scene.h"
#include "mesh.h"
#include "lodepng.h"
#include <chrono>


Scene* rt_scene;

int main(int argc, const char ** argv)
{
	rt_scene = new Scene();

	if (argc > 1)
	{
		int result = Parser::loadXmlintoScene(argv[1], rt_scene);
		if (result == -1)
		{
			std::cout << "Exiting program" << std::endl;
		}
	}
	else
	{
		int result = Parser::loadXmlintoScene("sponza_direct.xml", rt_scene);
		if (result == -1)
		{
			std::cout << "Exiting program" << std::endl;
		}
	}

	std::cout << "Parsing completed." << std::endl;

	std::chrono::time_point<std::chrono::system_clock> start_time;
	start_time = std::chrono::system_clock::now();

	rt_scene->render();

	using second = std::chrono::duration<double, std::ratio <1>>;
	double elapsed_time =  std::chrono::duration_cast<second>(std::chrono::system_clock::now() - start_time).count();

	std::cout << "Elapsed Time: " << elapsed_time << " seconds" << std::endl;

	system("PAUSE");

	return 0;
}