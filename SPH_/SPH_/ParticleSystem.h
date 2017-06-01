#ifndef ParticleSystem_h
#define ParticleSystem_h
#include "Particle.h"
#include "Container.h"
#include "modelLoader.h"
#include "nanoflann.hpp"
#include "Integrate.h"
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#include "include\OpenCL\CL\cl.hpp"
#include "include\OpenCL\CL\CLUtil.hpp"



#include <vector>
#include<random>

using namespace std;

//				MAGIC HAPPENS				//

#define MAX_PARTICLES	5000 
#define GRAVITY			1.6f
#define SEARCH_RADIUS	0.3f
#define REST_COEF		0.5f		/// SS: .5f
#define SMOOTHINGLENGTH 0.3f		/// SS: 0.3 / 0.4 
#define GAS_CONSTANT	5.0f		/// SS: 10.0f
#define PI				3.141592653589f
#define GRAVITYII		-9.28f
#define MAT_DENSITY		998.2f
#define VOLUME			100.0f		/// SS: 250.0f
#define VISC_CONST		1000.0f		/// SS:	5000.0f 
#define TIME_STEP		0.00099f	/// SS: 0.0005f 

//				NO MOAR MAGIC =[			//

enum Dist{
	inBB,
	DamBreak,
	Random
};

class ParticleSystem{
private:
	//	Methods
	void CalculatePressure();
	void CalculateForces();
	string ConvertString(const char*);
	void UpdateVecs();
	void InBounds(Particle&);
	std::vector<glm::vec3> GetPositions();
	


	template <typename num_t>
	void GetNeighbours();

	//	Variables / objects
	vector<Particle>	Particles;
	std::vector<std::vector<Particle*>> __nnv;
	BB					dim;
	model*				p_Model;
	Integrate			Integrator;
	bool				simPause, colPressure;
	std::vector<glm::vec4> i_Positions;
	std::vector<glm::mat4> i_Models;
	std::vector<glm::vec3> i_Color;
	std::vector<glm::vec3> p_Positions;
	size_t					DistType;
	//std::vector<std::pair<Particle, std::vector<Particle>>> PNPairs;

	//	Kernel specific methods

	float m_weightPoly, m_weightPolyGrad, m_weightPolyLap, m_weightPressGrad, m_weightViscLap, m_SpikyGrad;
	float KernelPoly(const float);
	Vec3f KernelPolyGrad(const Vec3f);
	float KernelPolyLap(const float);
	Vec3f KernelPresGrad(const Vec3f);
	float kernelViscLap(const float);
	Vec3f KernelSpikyGrad(const Vec3f);


public:
	//	Variables / objects 
	modelLoader			m_Loader;
	Container			container;


	//	Methods/Functions
	ParticleSystem();
	~ParticleSystem();
	void GenParticles();
	void InitCL();
	void ResetSimulation();
	void Run(float);

	//	Getters
	size_t				GetDistType(){return DistType;}
	vector<Particle>*	getParticles(){return &Particles;}
	const int			getParticleCount(){return Particles.size();}
	vector<Particle>*	GetContainer(){return &Particles;}
	bool				GetColPress(){return colPressure;}
	Particle			getParticle(int i){return Particles[i];}
	bool				GetPauseState(){return simPause;}
	model*				GetModel(){return p_Model;}
	vector<glm::vec4>*  GetIPosVec(){return &i_Positions;}
	vector<glm::mat4>*	GetIModVec(){return &i_Models;}
	vector<glm::vec3>*	GetIColVec(){return &i_Color;}
	void				CLGetNeighbours(size_t);
	//	Setters
	void				SetPause(bool b){simPause = b;}
	void				SetParticlePCol();
	void				SetDistribution(int);
	void				SetPressCol(bool _b){colPressure = _b;}
	void				SetCol(size_t, bool);
};


// NANOFLANN USAGE BELOW


using namespace nanoflann;

template <typename T>
struct PointCloud
{
	struct Point
	{
		T  x,y,z;
	};

	std::vector<Point>  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline T kdtree_distance(const T *p1, const size_t idx_p2,size_t /*size*/) const
	{
		const T d0=p1[0]-pts[idx_p2].x;
		const T d1=p1[1]-pts[idx_p2].y;
		const T d2=p1[2]-pts[idx_p2].z;
		return d0*d0+d1*d1+d2*d2;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim==0) return pts[idx].x;
		else if (dim==1) return pts[idx].y;
		else return pts[idx].z;
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

};

template <typename T>
void generateRandomPointCloud(PointCloud<T> &point, const size_t N, const T max_range = 10)
{
	std::cout << "Generating "<< N << " point cloud...";
	point.pts.resize(N);
	for (size_t i=0;i<N;i++)
	{
		point.pts[i].x = max_range * (rand() % 1000) / T(1000);
		point.pts[i].y = max_range * (rand() % 1000) / T(1000);
		point.pts[i].z = max_range * (rand() % 1000) / T(1000);
	}

	std::cout << "done\n";
}

template <typename T>
void addPtoPC(PointCloud<T> &point, std::vector<Particle>& p){
	point.pts.resize(p.size());
	for(size_t i = 0; i < p.size(); i++){
		point.pts[i].x = p.at(i).GetPosition().x;
		point.pts[i].y = p.at(i).GetPosition().y;
		point.pts[i].z = p.at(i).GetPosition().z;
	}
}


#endif

/// UNUSED OLD TING#



	//void CalculatePressure(std::vector<std::pair<Particle, std::vector<Particle>>>*);
	//void CalculateForces(std::vector<std::pair<Particle, std::vector<Particle>>>*);
	//void BoundaryCheck();
	//template <typename num_t>
	//void BuildKD(std::vector<Particle>&);
		//void MoveParticles(float);
	//void ApplyGravity(float);

	//void CollisionResolvePVP(Particle&, std::vector<Particle>&);
	//bool PPCollision(Particle&, Particle&);


	//void MovePFromBounds(Particle*);
	//void MoveOutOfCollision(Particle*);
	//void LoadModel();
	/*KDTree *pTree;
	kdtree *aTree;
	kdIter *aIter;
	void ppCollResponse(OcTree*);
	bool pwCollTest(Particle*, Wall);
	void pwCollResponse(OcTree*);
	void MoveParticles(OcTree*, float);
	void NeighbourSearch(Particle*);
	std::vector<Particle> NSVec(Particle*);
	void PColCheck(Particle*, std::vector<Particle>*);
	void UpdateNeighbours();
	void BuildKD();
	bool ppCollTest(Particle*, Particle*);

	void BuildOT();
	void QueryOTNN(Particle*);
	OcTree* ot;*/