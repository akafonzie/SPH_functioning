#include<stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include <GL\glew.h>
#define GLM_FORCE_RADIANS
#include <glm\glm.hpp>
#include <glm\gtc\matrix_transform.hpp>
#include <GLFW\glfw3.h>
#include <time.h>
#include <assert.h>
//#include <IL\il.h>
#include <thread>
#include "include\ImGUI\imgui.h"



//User generated includes
#include "shaders.h"
#include "modelLoader.h"
#include "Container.h"
#include "ParticleSystem.h"
#include "Camera.h"
#include "Light.h"
#include "RenderText.h"
#include "Skybox.h"

using namespace std;

/*--------------------------------------------------------------------------------------

	- Main program below this line, ya know

---------------------------------------------------------------------------------------*/



//Global Pointers to user generated classes
shaders sPointer;
//modelLoader mPointer;
GLuint  c_Shader, i_Shader, p_Shader, t_Shader, sb_Shader,  p_Matrix;
GLFWwindow* window;
//Container container;
ParticleSystem* p_System;
Camera* cam;
Light* l_Amb;
RenderText* rt;
Skybox* _sb;

//Global variables
int g_gl_w=1280;
int g_gl_h=720;
int butt = 0;
double xPos=NULL, yPos=NULL;
float cTime, delta, f_Delta;

const float m_Sens = 0.1f;
bool simPause = true, renderContainer = true, renderParticles = true, 
	mFocus = true, mFocFirst = true, _debug = false, renderSB = true,
	mRightButt = false, DemoMode = false;

//GLuint ubo;
//struct ubo_Data_s{//data to be passed to shaders in the uniform buffer
//	float l_Pos[3];
//	float l_Int[3];
//	float cam_Pos[3];
//	float ambCoef;
//	float atten;
//}ubo_Data;

//Global method prototypes
void Demo();
void fps(GLFWwindow*);
float getDT();
void init(void);
void Input(float);
void Keys(GLFWwindow*, int, int, int, int);
//void LeUBO(GLuint);
int main();
void mouseButt(GLFWwindow*, int, int, int);
void PrintDebugInfo(size_t);
void render();
void resize(GLFWwindow*, int, int);
void RunSim();
void setColVec(std::vector<Particle>*);
void setupGLHints(void);
void setupShaders(void);
void setupWindow(void);
void update();

float g_FPS;

//Demo Mode, scroll through types and reset the simulation after a certain time
void Demo(){
	static double _counter = glfwGetTime();
	double _dCounter = glfwGetTime();
	double _elap  = _dCounter - _counter;
	if((int)_elap % 25 == 0){
		p_System->ResetSimulation();
	}
	
}

//Calculate FPS and show in window Title
void fps(GLFWwindow*)
{
	static double previous_seconds = glfwGetTime();
	static int frame_count;
	double current_seconds = glfwGetTime();
	double elapsed_seconds = current_seconds - previous_seconds;
	glm::vec3 t = cam->GetPosition();
	if(elapsed_seconds > 0.25)
	{
		previous_seconds = current_seconds;
		double fps = (double)frame_count / elapsed_seconds;
		char tmp[128];
		sprintf(tmp, "SPH | Delta: %.10f ", f_Delta);
		glfwSetWindowTitle(window, tmp);
		frame_count = 0;
		g_FPS = (float)fps;
	}
	frame_count ++;
}

//Get Delta TIme
float getDT(){
	float nTime = (float)glfwGetTime();
	float fTime = nTime - cTime;
	cTime = nTime;
	return fTime;
}

//Init method, to set everything up!
void init()
{
	setupWindow();
	setupGLHints();
	glewInit();
	glViewport(0, 0, g_gl_w, g_gl_h);
	// Set up scene camera
	cam = new Camera(glm::vec3(0.0f, 0.8f, 1.8f), 45.0f, (GLfloat)g_gl_w / (GLfloat)g_gl_h, 0.1f, 1000.0f);
	cam->SetAngles(0.0f, 28.0f);


	//instantiate any pointer objects
	p_System = new ParticleSystem;

	//	Set up the light for the scene
	l_Amb = new Light(	1,						 //	Type
						Vec3f(0.0f, 2.0f, 0.0f), //	Position
						Vec3f(1.0f, 1.0f, 1.0f), // Intensity
						0.2f,					 // Ambient Coefficient (Default 0.2f)
						0.001f);				 // Attenuation (Default 0.005f)

	rt = new RenderText();

	_sb = new Skybox();

	

}

//Handle mouse input, and keys for movement since the callback method has input lag
void Input(float dt){
	float camMoveSpeed = 0.5f;

	//	Camera interaction keys
	if(glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS){
		cam->OffsetPosition(camMoveSpeed  * dt *  cam->Forward());
	}
	if(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS){
		cam->OffsetPosition(camMoveSpeed  * dt * -cam->Forward());
	}
	if(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS){
		cam->OffsetPosition(camMoveSpeed  * dt * -cam->Right());
	}
	if(glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS){
		cam->OffsetPosition(camMoveSpeed  * dt *cam->Right());
	}
	if(glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS){
		cam->OffsetPosition(camMoveSpeed  * dt *cam->Up());
	}
	if(glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS){
		cam->OffsetPosition(camMoveSpeed  * dt * -cam->Up());
	}



	if(glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS){
		p_System->ResetSimulation();
	}
	
	///	MOUSE INPUT
	glfwGetCursorPos(window, &xPos, &yPos);
	glfwSetCursorPos(window, 0, 0);	
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	//double _xp,  _yp;
	//glfwGetCursorPos(window, &_xp, &_yp);
	if(mRightButt)
		cam->OffsetOrientation(m_Sens * (float)xPos, m_Sens * (float)yPos);
	glfwSetCursorPos(window, 0, 0);	
	//add mouse input to camera
}

void Keys(GLFWwindow* w, int a, int b, int c, int d){ //a == the key | b == scancode | c == action | d == mods
	/* FOR EASE OF COPYING

	if(a == GLFW_KEY_W && c == GLFW_PRESS){

	}

	*/
	
	
	float camMoveSpeed = 0.01f;
	
	if(a == GLFW_KEY_G && c == GLFW_PRESS){
		p_System->CLGetNeighbours(200);
	}

		if(a == GLFW_KEY_H && c == GLFW_PRESS){
		p_System->InitCL();
	}



	///	PROGRAM EXIT
	if(a == GLFW_KEY_ESCAPE && c == GLFW_PRESS){
		glfwSetWindowShouldClose(window, 1);
	}

	///	SIMULATION INTERACTION KEYS
	//Pause and unPause the simulation
	if(a == GLFW_KEY_P && c == GLFW_PRESS){
		if(p_System->GetPauseState())
			p_System->SetPause(false);
		else
			p_System->SetPause(true);
	}

	//Set the container to render, or not
	if(a == GLFW_KEY_C && c == GLFW_PRESS){
		if(renderContainer)
			renderContainer = false;
		else
			renderContainer = true;
	}

	//Set the particles to render, or not
	if(a == GLFW_KEY_X && c == GLFW_PRESS){
		if(!renderParticles)
			renderParticles = true;
		else
			renderParticles = false;
	}

	if(a == GLFW_KEY_Z && c == GLFW_PRESS){
		if(!renderSB)
			renderSB = true;
		else
			renderSB = false;
	}

	// Set distribution to In BB
	if(a == GLFW_KEY_1 && c == GLFW_PRESS){
		p_System->SetDistribution(0);
	}

	// Set distribution to DamBreak
	if(a == GLFW_KEY_2 && c == GLFW_PRESS){
		p_System->SetDistribution(1);
	}

	//turn pressure colouring on and off
	if(a == GLFW_KEY_N && c == GLFW_PRESS){
		bool _b = p_System->GetColPress();
		if (_b == 1)
			p_System->SetPressCol(false);
		else
			p_System->SetPressCol(true);
	}

	//toggle debug - i.e. printing out of info to the screen
	if(a == GLFW_KEY_M && c == GLFW_PRESS){
		if(_debug)
			_debug = false;
		else
			_debug = true;
	}
	
	if(a == GLFW_KEY_Y && c == GLFW_PRESS){
		if(DemoMode)
			DemoMode = false;
		else
			DemoMode = true;
	}
}

//void LeUBO(GLuint program){
//	{///	update the data held in the UBO
//		size_t b_Index = glGetUniformLocation(program, "ubo_Data");
//		GLint blockSize;
//		glGetActiveUniformBlockiv(program, b_Index, GL_UNIFORM_BLOCK_DATA_SIZE, &blockSize);
//		GLubyte* blockBuffer =	(GLubyte*) malloc(blockSize);
//		
//		glBindBuffer(GL_UNIFORM_BUFFER, ubo);
//		GLvoid* p = glMapBuffer(GL_UNIFORM_BUFFER, GL_WRITE_ONLY);
//		memcpy(p, &ubo_Data, sizeof(ubo_Data));
//		glUnmapBuffer(GL_UNIFORM_BUFFER);
//					
//		
//		GLuint bp_Index = 2;
//		glBindBufferBase(GL_UNIFORM_BUFFER, bp_Index, ubo);
//		glUniformBlockBinding(program, b_Index, bp_Index);
//		
//		ubo = 0;
//		glGenBuffers(1, &ubo);
//		glBindBuffer(GL_UNIFORM_BUFFER, ubo);
//		glBufferData(GL_UNIFORM_BUFFER, sizeof(ubo_Data), &ubo_Data, GL_DYNAMIC_DRAW);
//		glBindBuffer(GL_UNIFORM_BUFFER, 0);
//	}
//} //ubo

//Main method containing the main loop
int main(){
	init();
	setupShaders();
	cTime = (float)glfwGetTime();
	while (!glfwWindowShouldClose(window)){
		fps(window);
		update();
		render();
		glfwPollEvents();
		glfwSwapBuffers(window);
	}
	glfwTerminate();
	//system("cls"); //clears the console window
	//printf("\tHAVE NICE DAY xoxo\n\n\n");
	return 0;
}

//This method handles mouse presses for the program!
void mouseButt(GLFWwindow*, int a, int b, int c){
	
	//left mouse button
	if(a == 0 && b == 1){
	}else if(a == 0 && b == 0){
		
	}
	
	//middle mouse button
	if(a==2 && b == 1){printf("middle button pressed @ (%.0f, %.0f)\n", xPos, yPos);}

	//right mouse button
	if(a==1 && b == 1){
		mRightButt = true;
	}else if (a==1 && b == 0){
		mRightButt = false;
	}

}

//This method prints out debug info for the selected particle
void PrintDebugInfo(size_t _p){
		glFrontFace(GL_CCW); //need to disable this prior to rendering else it won't work!
		glm::mat4 p = glm::ortho(0.0f, (float)g_gl_w, 0.0f, (float)g_gl_h);
		//sPointer.UseTextShader(t_Shader, p, glm::vec3(1.0f, 1.0f, 1.0f));
		// rt->PrintThis(t_Shader, "", p, 0.0f, 0.0f, 1.0f, glm::vec3(0.0f, 1.0f, 0.0f));
		string boob, bewb, m, d, pr, v;
		// PRINT FPS
		bewb = std::to_string((int)g_FPS);
		boob = "FPS: " + bewb + " ";
		rt->PrintThis(t_Shader, boob, p, 0.0f, 708.0f, 0.3f, glm::vec3(0.0f, 0.85f, 0.0f));
		// PRINT SIM STATE
		if(p_System->GetPauseState())
			boob = "SIMULATION = PAUSED ";
		if(!p_System->GetPauseState())
			boob = "SIMULATION = RUNNING ";
		rt->PrintThis(t_Shader, boob, p, 0.0f, 690.0f, 0.3f, glm::vec3(0.0f, 0.85f, 0.0f));
		// PRINT RANDOM PARTICLE INFO FOR DEBUG YA KNOW
		if(_debug){
			Particle& tmp = p_System->getParticle(_p);
			m = std::to_string((float)tmp.GetMass());
			d = std::to_string((float)tmp.GetDensity());
			pr = std::to_string((float)tmp.GetPressure());
			v = std::to_string((float)tmp.GetViscosity());
			bewb = std::to_string((int)tmp.GetID());
			boob = "Particle ("+std::to_string(_p)+") - ID = "+bewb+" : MASS = "+ m +" | DENSITY = "+ d +" | PRESSURE = "+pr+" | VISCOSITY = "+v;
			rt->PrintThis(t_Shader, boob, p, 0.0f, 672.0f, 0.3f, glm::vec3(1.0f, 0.0f, 0.0f));
			m = std::to_string((float)tmp.GetPressureForce().x);
			d = std::to_string((float)tmp.GetPressureForce().y);
			pr =std::to_string((float)tmp.GetPressureForce().z);
			boob = "PressureForce = ("+m+", "+d+", "+pr+") ";
			rt->PrintThis(t_Shader, boob, p, 0.0f, 654.0f, 0.3f, glm::vec3(1.0f, 0.0f, 0.0f));
			m = std::to_string((float)tmp.GetViscosityForce().x);
			d = std::to_string((float)tmp.GetViscosityForce().y);
			pr =std::to_string((float)tmp.GetViscosityForce().z);
			boob = "ViscosityForce = ("+m+", "+d+", "+pr+") ";
			rt->PrintThis(t_Shader, boob, p, 0.0f, 636.0f, 0.3f, glm::vec3(1.0f, 0.0f, 0.0f));
			m = std::to_string((float)tmp.GetPosition().x);
			d = std::to_string((float)tmp.GetPosition().y);
			pr =std::to_string((float)tmp.GetPosition().z);
			boob = "Position = ("+m+", "+d+", "+pr+") ";
			rt->PrintThis(t_Shader, boob, p, 0.0f, 618.0f, 0.3f, glm::vec3(1.0f, 0.0f, 0.0f));
			m = std::to_string((float)tmp.GetColour().x);
			d = std::to_string((float)tmp.GetColour().y);
			pr = std::to_string((float)tmp.GetColour().z);
			boob = "Color = ("+m+", "+d+", "+pr+") ";
			rt->PrintThis(t_Shader, boob, p, 0.0f, 600.0f, 0.3f, glm::vec3(1.0f, 0.0f, 0.0f));
		}
		glFrontFace(GL_CW); //then re-enable it again!
}

//Render deals with all instances that need to be drawn
void render(void){
	{///	Skybox
		if(renderSB){
			//glm::mat4 id = glm::mat4(1.0f);
			glDisable(GL_CULL_FACE);		
		//	glm::mat4 sbp = glm::ortho(0.0f, (float)g_gl_w, 0.0f, (float)g_gl_h);
			//use shader
			glDepthMask (GL_FALSE);
			sPointer.UseSBShader(sb_Shader, cam);
			_sb->Render();
			glDepthMask (GL_TRUE);
			glEnable(GL_CULL_FACE);
		}
	}

	{///	Render Container
		if(renderContainer){
			sPointer.UseShader(	c_Shader, p_System->container.GetModel(), cam, l_Amb->GetPos(), 
								l_Amb->GetInt(), l_Amb->GetAC(), l_Amb->GetAtt());
			p_System->container.Render();
		}
	}


	{///	Render Particles
		if(renderParticles){
			sPointer.UseShaderInstanced(i_Shader, cam, l_Amb->GetPos(), l_Amb->GetInt(), l_Amb->GetAC(), l_Amb->GetAtt());
			p_System->m_Loader.RenderIModel(p_System->GetModel(), p_System->getParticleCount(), p_System->GetIPosVec(), 
											p_System->GetIModVec(), p_System->GetIColVec());
		}
	}

	

	{///	Render Screen Text
		PrintDebugInfo(4999);
	}

}

//Resize deals with the resizing of the program window
void resize(GLFWwindow*, int w, int h)
{
	g_gl_w=w;
	g_gl_h=h;
	glViewport(0, 0, g_gl_w, g_gl_h);
}

//Setup the GL hints
void setupGLHints(void){
	// *** SHADING
	glShadeModel	(GL_SMOOTH);
	glEnable		(GL_NORMALIZE);
	glHint			(GL_LINE_SMOOTH_HINT, GL_NICEST);
	
	// *** CULLING
	glPolygonMode	(GL_FRONT_AND_BACK, GL_FILL);
	glEnable		(GL_CULL_FACE);
	glFrontFace		(GL_CW);
	glCullFace		(GL_BACK);
	glEnable		(GL_DEPTH_TEST);
	glDepthFunc		(GL_LESS);
	glDepthMask		(GL_TRUE);
	
	// *** BLENDING
	glEnable		(GL_BLEND);
	glBlendFunc		(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

//Setup the Shaders
void setupShaders(void){
	c_Shader	= sPointer.loadShaders("shaders/vs.glsl","shaders/fs.glsl");
	i_Shader	= sPointer.loadShaders("shaders/i_VS.glsl", "shaders/i_FS.glsl");
	t_Shader	= sPointer.loadShaders("shaders/t_VS.glsl", "shaders/t_FS.glsl");
	sb_Shader	= sPointer.loadShaders("shaders/sb_VS.glsl", "shaders/sb_FS.glsl");
}

//Setup the window, and set the hints to get the best possible look.
void setupWindow(void){
	glfwInit();
	if(!glfwInit()){
		fprintf(stderr, "ERROR | Could not start GLFW3\n");
	}

	window = glfwCreateWindow(1280, 720, "SPH", NULL, NULL);
	if(!window){
		fprintf(stderr, "ERROR | Could not create window with GLFW\n"); 
		glfwTerminate();
		std::exit(0);
	}
	glfwMakeContextCurrent(window);
	glfwSetWindowSizeCallback(window, resize);
	glfwSetMouseButtonCallback(window, mouseButt);
	glfwSetKeyCallback(window, Keys);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SAMPLES, 4); //Multi sample Anti Aliasing
	GLFWmonitor* mon=glfwGetPrimaryMonitor();
	glfwSwapInterval(1);
}

//This method handles per-frame updates, and any methods/vars that need updating
void update(void){
	f_Delta = getDT();
	Input(f_Delta);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	if(!p_System->GetPauseState()){	
		std::thread simThread(RunSim);
		simThread.join();
	}
	if(DemoMode)
		Demo();
	
}

void RunSim(){
	p_System->Run(f_Delta);
}

