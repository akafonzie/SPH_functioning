#ifndef Container_h
#define Container_h
#include "Vec3f.h"
#include "include\GL\glew.h"
#include "include\glm\glm.hpp"
#include <GLFW\glfw3.h>
#include <vector>

#define HEIGHT .75f
#define WIDTH .75f

struct BB{
	Vec3f min, max;
};

class Container{
private:
	BB containerBB;
	glm::mat4 Model;
public:
	void genVAO();
	void Render();
	Container();
	~Container();
	glm::mat4 GetModel(){return Model;}
	BB getDims(){return containerBB;}
	Vec3f cont[24];
	Vec3f contNormals[24];
	Vec3f fNormals[6];
	GLuint VBO, VAO, CBO, VNO;
};
#endif