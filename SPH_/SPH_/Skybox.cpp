#include "Skybox.h"
#define STB_IMAGE_IMPLEMENTATION
#include "include\stb\stb_image.h"



Skybox::Skybox(){
	float points[] = {
  -10.0f,  10.0f, -10.0f,
  -10.0f, -10.0f, -10.0f,
   10.0f, -10.0f, -10.0f,
   10.0f, -10.0f, -10.0f,
   10.0f,  10.0f, -10.0f,
  -10.0f,  10.0f, -10.0f,
  
  -10.0f, -10.0f,  10.0f,
  -10.0f, -10.0f, -10.0f,
  -10.0f,  10.0f, -10.0f,
  -10.0f,  10.0f, -10.0f,
  -10.0f,  10.0f,  10.0f,
  -10.0f, -10.0f,  10.0f,
  
   10.0f, -10.0f, -10.0f,
   10.0f, -10.0f,  10.0f,
   10.0f,  10.0f,  10.0f,
   10.0f,  10.0f,  10.0f,
   10.0f,  10.0f, -10.0f,
   10.0f, -10.0f, -10.0f,
   
  -10.0f, -10.0f,  10.0f,
  -10.0f,  10.0f,  10.0f,
   10.0f,  10.0f,  10.0f,
   10.0f,  10.0f,  10.0f,
   10.0f, -10.0f,  10.0f,
  -10.0f, -10.0f,  10.0f,
  
  -10.0f,  10.0f, -10.0f,
   10.0f,  10.0f, -10.0f,
   10.0f,  10.0f,  10.0f,
   10.0f,  10.0f,  10.0f,
  -10.0f,  10.0f,  10.0f,
  -10.0f,  10.0f, -10.0f,
  
  -10.0f, -10.0f, -10.0f,
  -10.0f, -10.0f,  10.0f,
   10.0f, -10.0f, -10.0f,
   10.0f, -10.0f, -10.0f,
  -10.0f, -10.0f,  10.0f,
   10.0f, -10.0f,  10.0f
};

	//_sbVBO = _sbVAO = 0;
	glGenBuffers(1, &_sbVBO);
	glBindBuffer(GL_ARRAY_BUFFER, _sbVBO);
	glBufferData(GL_ARRAY_BUFFER, 3 * 36 * sizeof(float), &points, GL_STATIC_DRAW);
	
	glGenVertexArrays(1, &_sbVAO);
	glBindVertexArray(_sbVAO);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, _sbVBO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

	glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &_tex);
	//glBindTexture(GL_TEXTURE_CUBE_MAP, _tex);
	
	loadSide(&_tex, GL_TEXTURE_CUBE_MAP_NEGATIVE_Z,	"SB/front.png");
	loadSide(&_tex, GL_TEXTURE_CUBE_MAP_POSITIVE_Z,	"SB/back.png");
	loadSide(&_tex, GL_TEXTURE_CUBE_MAP_POSITIVE_Y,	"SB/top.png");
	loadSide(&_tex, GL_TEXTURE_CUBE_MAP_NEGATIVE_Y,	"SB/bottom.png");
	loadSide(&_tex, GL_TEXTURE_CUBE_MAP_NEGATIVE_X,	"SB/left.png");
	loadSide(&_tex, GL_TEXTURE_CUBE_MAP_POSITIVE_X,	"SB/right.png");


	glTexParameteri ( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri ( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri ( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	glTexParameteri ( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri ( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

}

bool Skybox::loadSide(GLuint* _t, GLenum _side, const char* _fp){
	glBindTexture(GL_TEXTURE_CUBE_MAP, *_t);
	int x, y, n, _fc;
	_fc = 4;
	unsigned char* imgData;
	
	imgData = stbi_load(_fp, &x, &y, &n, _fc);
	if(!imgData){
		printf("Sorry Jabroni, could not load %s\n", _fp);
		return false;
	}
	
	//check  img is !^2
	//if( (x & (x - 1)) != 0 || (y & (y - 1)) != 0 ){
	//	printf("Your image is not ^2 dims sucka!\n");
	//}
	glTexImage2D(_side, 0, GL_RGBA, x, y, 0, GL_RGBA, GL_UNSIGNED_BYTE, imgData);
	free(imgData);
	return true;
}


void Skybox::Render(){


	glActiveTexture (GL_TEXTURE0);
	glBindTexture (GL_TEXTURE_CUBE_MAP, _tex);
	glBindVertexArray (_sbVAO);
	glDrawArrays (GL_TRIANGLES, 0, 36);

}