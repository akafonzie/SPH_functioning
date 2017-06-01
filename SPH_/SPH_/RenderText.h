#ifndef RenderText_h
#define RenderText_h
//GL includes
#include "include\GL\glew.h"
#include "include\GLFW\glfw3.h"
#include "include\glm\glm.hpp"
#include "include\glm\gtc\matrix_transform.hpp"
//std includes
#include <stdio.h>
#include <string>
#include <map>
// Freetype includes
#include "include\FreeType\include\ft2build.h"
#include FT_FREETYPE_H

struct Character{
	GLuint		TextureID;
	glm::ivec2  Size;
	glm::ivec2  Bearing;
	GLuint		Advance;
};

class RenderText{
private:
	FT_Library	ft;
	FT_Face		fc;
	size_t		size;
	std::map<GLchar, Character> Characters;
	GLuint t_VAO, t_VBO;
public:
	RenderText();
	~RenderText();
	void PrintThis(GLuint, std::string, glm::mat4, GLfloat, GLfloat, GLfloat, glm::vec3);
	void GenFont();
	bool LoadFont(const char*);
	void SetFontSize(size_t);
};
#endif