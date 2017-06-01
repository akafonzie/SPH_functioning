#ifndef OcTree_h
#define OcTree_h
#include "Vec3f.h"
#include "Particle.h"
#include <vector>
#include <set>

#define MAX_OCTREE_DEPTH 6
#define MIN_P_PER_OCTREE 3
#define MAX_P_PER_OCTREE 6

enum Wall{
	W_LEFT ,
	W_RIGHT,
	W_BACK,
	W_FRONT,
	W_TOP,
	W_BOTTOM
};

struct pPair{
	Particle* p1; 
	Particle* p2;
};

struct pwPair{
	Particle* p;
	Wall wall;
};

class OcTree{
public:
	Vec3f corner1, corner2, center;
	OcTree *Children[2][2][2];
	bool hasChildren;
	int depth;
	int numParticles;


	set<Particle*> oParticles;
	void FileParticle(Particle*, Vec3f, bool);
	void HaveChildren();
	void CollectParticles(set<Particle*> &);
	void DestroyChildren();
	void Remove(Particle*, Vec3f);
	OcTree(Vec3f, Vec3f, int);
	~OcTree();
	void Add(Particle*);
	void Remove(Particle*);
	void ParticleMoved(Particle*, Vec3f);
	void PotentialPPColl(std::vector<pPair> &);
	void PotentialPPColl(std::vector<pPair> &, std::vector<Particle> &, OcTree*);
	void PotentialPWColl(std::vector<pwPair> &);
	void PotentialPWColl(std::vector<pwPair> &, Wall, char, int);
	//void PotentialPWColl(std::vector<pwPair> &, std::vector<Particle> &, OcTree*);
	Vec3f GetWallDir(Wall);
};


//struct OcPoint{
//	Vec3f Pos;
//	size_t idx;
//	Vec3f GetPos(){return Pos;}
//	void SetPos(Vec3f& p){Pos = p;}
//	void SetIDX(size_t& i){idx = i;}
//};
//
//class OcTree{
//private: 
//	Vec3f Origin, hDim;
//	OcTree* Children[8];
//	OcPoint* Data;
//public:
//	OcTree(const Vec3f&, const Vec3f&);
//	OcTree(const OcTree&);
//	~OcTree();	
//	int GetOctantContainingPoint(const Vec3f&);
//	bool IsLeafNode();
//	void Insert(OcPoint*);
//	void GetPointsInCube(const Vec3f&, const Vec3f&, std::vector<OcPoint*>&);
//	void SetPoints(size_t, std::vector<Particle>*);
//};
#endif