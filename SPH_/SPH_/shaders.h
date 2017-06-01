#pragma
#ifndef SHADERS_H
#define SHADERS_H
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "include\GL\glew.h"
#include "include\GLFW\glfw3.h"
#include "include\glm\glm.hpp"
#include "Particle.h"
#include "Camera.h"
#include "ppf.h"

using namespace std;

struct uniformData{
	float camPos[3], lPos[3], lInt[3], l_Ac, l_At, pCol[3];
	float cam[16];
};


class shaders{
private:
	ppf* _pr;
public:
	GLuint loadShaders(const char * vertex_file_path,const char * fragment_file_path);
	//void UseShaderInstanced(GLuint, Camera*, GLuint, vector<glm::vec4>*, vector<glm::mat4>*, Vec3f, Vec3f, float, float, std::vector<Particle>*);
	void UseShaderInstanced(GLuint, Camera*, Vec3f, Vec3f, float, float);
	void UseSBShader(GLuint, Camera*);
	void UseShader(GLuint, glm::mat4, Camera*, Vec3f, Vec3f, float, float);
	void UseTextShader(GLuint, glm::mat4, glm::vec3);
	shaders();
};
#endif