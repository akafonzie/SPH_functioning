#include "RenderText.h"

RenderText::RenderText(){
	if(FT_Init_FreeType(&ft))
		printf("ERROR: could not init freetype lib\n");

	if(FT_New_Face(ft, "fonts/arial.ttf", 0, &fc))
		printf("ERROR: could not find the font Pendejo!\n");
	FT_Set_Pixel_Sizes(fc, 0, 48);
	GenFont();
}


void RenderText::GenFont(){
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	for(GLubyte i = 0; i < 128; i++){
		if(FT_Load_Char(fc, i, FT_LOAD_RENDER)){
			printf("ERROR: Failed to load glyph!\n");
			continue;
		}
		GLuint tex;
		glGenTextures(1, &tex);
		glBindTexture(GL_TEXTURE_2D, tex);
		glTexImage2D(
			GL_TEXTURE_2D, 
			0, 
			GL_RED, 
			fc->glyph->bitmap.width, 
			fc->glyph->bitmap.rows, 
			0, 
			GL_RED, 
			GL_UNSIGNED_BYTE,
			fc->glyph->bitmap.buffer);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		Character ch = {
			tex, 
			glm::ivec2(fc->glyph->bitmap.width, fc->glyph->bitmap.rows),
			glm::ivec2(fc->glyph->bitmap_left, fc->glyph->bitmap_top), 
			fc->glyph->advance.x
		};
		Characters.insert(std::pair<GLchar, Character>(i, ch));
	}
	glBindTexture(GL_TEXTURE_2D, 0);
	FT_Done_Face(fc);
	FT_Done_FreeType(ft);

	glGenVertexArrays(1, &t_VAO);
	glGenBuffers(1, &t_VBO);
	glBindVertexArray(t_VAO);
	glBindBuffer(GL_ARRAY_BUFFER, t_VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

RenderText::~RenderText(){

}

void RenderText::PrintThis(GLuint shader, std::string text, glm::mat4 p, GLfloat x, GLfloat y, GLfloat sz, glm::vec3 col){
	
	glUseProgram(shader);
	glUniform3f(glGetUniformLocation(shader, "textColor"), col.x, col.y, col.z);
	//glm::mat4 proj = glm::ortho(0.0f, 1280.0f, 0.0f, 720.0f);
	glUniformMatrix4fv(glGetUniformLocation(shader, "projection"), 1, GL_FALSE, &p[0][0]);
	
	glActiveTexture(GL_TEXTURE0);
	glBindVertexArray(t_VAO);

	std::string::const_iterator c;
	for(c = text.begin(); c != text.end(); c++){
		Character ch = Characters[*c];
		GLfloat xpos = x + ch.Bearing.x * sz;
		GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * sz;

		GLfloat w = ch.Size.x * sz;
		GLfloat h = ch.Size.y * sz;

		GLfloat verts[6][4] = {
			{xpos, ypos + h,	0.0f, 0.0f}, 
			{xpos, ypos,		0.0f, 1.0f}, 
			{xpos + w, ypos,	1.0f, 1.0f},

			{xpos, ypos + h,	0.0f, 0.0f}, 
			{xpos + w, ypos,	1.0f, 1.0f},
			{xpos + w, ypos + h, 1.0f, 0.0f}
		};
		glBindTexture(GL_TEXTURE_2D, ch.TextureID);
		glBindBuffer(GL_ARRAY_BUFFER, t_VBO);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(verts), verts);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glDrawArrays(GL_TRIANGLES, 0, 6);
		x+= (ch.Advance >> 6) * sz;
	}
	glBindVertexArray(0);
	glBindTexture(GL_TEXTURE_2D, 0);
}

bool RenderText::LoadFont(const char* fp){
	if(FT_New_Face(ft, fp, 0, &fc)){//returns non-zero if an error occurred 
		return false;
		printf("ERROR: could not load the font Pendejo!\n");
		GenFont();
	}else
		return true;
}

void RenderText::SetFontSize(size_t sz){
	size = sz;
	FT_Set_Pixel_Sizes(fc, 0, size);
}