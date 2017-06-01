#ifndef Particle_h
#define Particle_h
#include <vector>
#include "Vec3f.h"
#include "modelLoader.h"

const static float p_Scale = 0.015f;

class Particle{
	private:
	Vec3f pos, dir, colour, acceleration, p_accel, velocity, v_Force, p_Force, grav, netForce;
	size_t p_ID;
	float mass, pressure, viscosity, radius, density;
	bool inBounds, isSurfaceParticle, isFalling;
	model*	p_Model;
	bool LoadModel();
	glm::mat4 ModelView;
	std::vector<unsigned int> n_Indexes; //indexes of neighbour particles

public:
	//Methods
	Particle();
	Particle(Vec3f);
	~Particle();
	void GenParticle(Vec3f);
	void Render();
	void UpdateMV(Vec3f);
	void RemoveNeighbourIndex(unsigned int);
	void AccumulateForces(Vec3f _ff){netForce += _ff;}
	void ResetForces();

	//Getters
	Vec3f GetAcceleration(){return acceleration;}
	bool GetBoundCheck(){return inBounds;}
	Vec3f GetColour(){return colour;};
	glm::vec3 GetColourGLM(){return glm::vec3(colour.x, colour.y, colour.z);}
	float GetDensity(){return density;}
	bool GetFalling(){return isFalling;}
	Vec3f GetForce(){return netForce;}
	Vec3f GetGravityForce(){return grav;}
	size_t GetID(){return p_ID;}
	float GetMass(){return mass;}
	model* GetModel(){return p_Model;}
	glm::mat4 GetMVM(){return ModelView;}
	std::vector<unsigned int> GetNeighbourIndexes(){return n_Indexes;}
	Vec3f GetPAcceleration(){return p_accel;}
	Vec3f GetPosition(){return pos;}
	glm::vec3 GetPositionGLM(){return glm::vec3(pos.x, pos.y, pos.z);}
	float GetPressure(){return pressure;}
	Vec3f GetPressureForce(){return p_Force;}
	float GetRadius(){return radius;}
	float GetScale(){return p_Scale;}
	Vec3f GetVelocity(){ return velocity;}
	float GetViscosity(){return viscosity;}
	Vec3f GetViscosityForce(){return v_Force;}
	
	//Setters
	void SetAccel(Vec3f _a){p_accel = acceleration; acceleration = _a;}
	void SetBounds(bool b){inBounds = b;}
	void SetColour(float r, float g, float b){colour = Vec3f(r, g, b);}
	void SetColour(Vec3f c){colour = c;}
	void SetDensity(float _d){density = _d;}
	void SetFalling(bool b){isFalling = b;}
	void SetGravity(Vec3f _g){grav = _g;}
	void SetID(size_t id){p_ID = id;}
	void SetMass(float _m){mass = _m;}
	void SetModel(model* p){p_Model = p;}
	void SetPosition(Vec3f p){pos = p; UpdateMV(p);}
	void SetPressure(float p){pressure = p;}
	void SetPressureForce(Vec3f _pf){p_Force = _pf;}
	void SetVelocity(Vec3f v){velocity = v;}
	void SetViscosity(float _v){viscosity = _v;}
	void SetViscosityForce(Vec3f _vf){v_Force = _vf;}



	void AddNeighbourIndex(unsigned int i){n_Indexes.push_back(i);}
	
	//Variables
	modelLoader m_Loader;
};
#endif


