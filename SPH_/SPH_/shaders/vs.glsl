#version 400

layout (location=0) in vec3 v_Pos;
layout (location=1) in vec3 v_Color;
layout (location=2) in vec3 v_Normals;

uniform mat4 Model;
uniform mat4 Camera;
uniform vec3 camPos;


//layout (std140) uniform ubo_Data{
//	vec3 l_Pos;
//	vec3 l_Int;
//	vec3 cam_Pos;
//	float ambCoef;
//	float atten;
//}ubo;

out vec4 fragmentColor; 
out vec3 fragmentNormal;
out vec3 fragmentVertex;
out mat4 m_Out;
out vec3 CamPos;


//out vec3 lpos;
//out vec3 lint;
//out float ambco;
//out float att; 


void main()
{
	
	fragmentColor	= vec4(v_Color.rgb, 1.0f);;
	fragmentVertex	= v_Pos;
	fragmentNormal  = v_Normals;
	m_Out			= Model;
	CamPos			= camPos;
	gl_Position		=  Camera * Model * vec4(v_Pos ,1);

	//lpos			= ubo.l_Pos;
	//lint			= ubo.l_Int;
	//ambco			= ubo.ambCoef;
	//att				= ubo.atten;
}