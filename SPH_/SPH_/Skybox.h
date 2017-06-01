#ifndef Skybox_h
#define Skybox_h

#include "Vec3f.h"
#include "include\GL\glew.h"
#include "include\GLFW\glfw3.h"
#include <string>
#include <assert.h>

enum _sbfaces{
	front, 
	back,
	top,
	bottom,
	left, 
	right
};

class Skybox{
private:
	std::string _fn[6];
	GLuint  _sbVBO, _sbVAO;
	GLuint _tex;
public:
	Skybox();
	bool loadSide(GLuint*, GLenum, const char*);
	GLuint GetVBO(){return _sbVBO;}
	GLuint GetVAO(){return _sbVAO;}
	GLuint GetTex(){return _tex;}
	void Render();
};
#endif