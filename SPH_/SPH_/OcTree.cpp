#include "OctTree.h"

OcTree::OcTree(Vec3f c1, Vec3f c2, int d){
	corner1			= c1;
	corner2			= c2;
	center			= (c1 + c2) / 2;
	depth			= d;
	numParticles	= 0;
	hasChildren		= false;
}

OcTree::~OcTree(){
	if(hasChildren){
		DestroyChildren();
	}
}

void OcTree::Add(Particle* p){
	numParticles++;
	if(!hasChildren && depth < MAX_OCTREE_DEPTH &&	numParticles > MAX_P_PER_OCTREE){
		HaveChildren();
	}

	if(hasChildren){
		FileParticle(p, p->GetPosition(), true);
	}else{
		oParticles.insert(p);
	}

}

void OcTree::Remove(Particle* p){
	Remove(p, p->GetPosition());
}

void OcTree::Remove(Particle* p, Vec3f pos){
	numParticles--;
	if(hasChildren && numParticles < MIN_P_PER_OCTREE){
		DestroyChildren();
	}

	if(hasChildren)
		FileParticle(p, pos, false);
	else
		oParticles.erase(p);
}

void OcTree::FileParticle(Particle* p, Vec3f pos, bool a){
	for(int i = 0; i < 2; i++){
		if(i == 0){
			if(pos.x - p->GetRadius() > center.x)
				continue;
		}
		else if(pos.x + p->GetRadius() < center.x)
			continue;
		
		for(int j = 0; j < 2; j++){
			if(j == 0){
				if(pos.y - p->GetRadius() > center.y)
					continue;
			}
			else if (pos.y + p->GetRadius() < center.y)
				continue;
		

			for(int k = 0; k < 2; k++){
				if(k == 0){
					if(pos.z - p->GetRadius() > center.z)
						continue;
				}
				else if (pos.z + p->GetRadius() < center.z)
					continue;

				if(a)
					Children[i][j][k]->Add(p);
				else
					Children[i][j][k]->Remove(p, pos);
			}
		}
	}
}

void OcTree::HaveChildren(){
	for(int i = 0; i < 2; i++){
		float minx, maxx;
		if(i == 0){
			minx = corner1.x;
			maxx = center.x;
		}else {
			minx = center.x;
			maxx = corner2.x;
		}
		for(int j = 0; j < 2; j++){
			float miny, maxy;
			if(j == 0){
				miny = corner1.y;
				maxy = center.y;
			}else{
				miny = center.y;
				maxy = corner2.y;
			}
			for(int k = 0; k < 2;k++){
				float minz, maxz;
				if(k == 0){
					minz = corner1.z;
					maxz = center.z;
				}else{
					minz = center.z;
					maxz = corner2.z;
				}

				Children[i][j][k] = new OcTree(Vec3f(minx, miny, minz), Vec3f(maxx, maxy, maxz), depth + 1);
			} // end k for
		}// end j for
	} // end i for

	for(set<Particle*>::iterator it = oParticles.begin(); it != oParticles.end(); it++){
		Particle* p = *it;
		FileParticle(p, p->GetPosition(), true);
	}
	oParticles.clear();
	hasChildren = true;
}

void OcTree::CollectParticles(set<Particle*> &ps){
	if(hasChildren){
		for(int i = 0; i < 2; i++){
			for(int j = 0; j < 2; j++){
				for(int k = 0; k < 2; k++){
					Children[i][j][k]->CollectParticles(ps);
				}
			}
		}
	}else{
		for(set<Particle*>::iterator it = oParticles.begin(); it != oParticles.end(); it++){
			Particle* p = *it;
			ps.insert(p);
		}
	}
}

void OcTree::DestroyChildren(){
	CollectParticles(oParticles);

	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			for(int k = 0; k < 2; k++){
				delete Children[i][j][k];
			}
		}
	}
	hasChildren = false;
}

void OcTree::ParticleMoved(Particle* p, Vec3f pos){
	Remove(p, pos);
	Add(p);
}

void OcTree::PotentialPPColl(std::vector<pPair> &coll){
	if(hasChildren){
		for(int i = 0; i < 2; i++){
			for(int j = 0; j < 2; j++){
				for(int k = 0; k < 2; k++){
					Children[i][j][k]->PotentialPPColl(coll);
				}
			}
		}
	}else{
		for(set<Particle*>::iterator it = oParticles.begin(); it != oParticles.end(); it++){
			Particle* p = *it;
			for(set<Particle*>::iterator it2 = oParticles.begin(); it2!=oParticles.end(); it2++){
				Particle* p2 = *it2;
				if(p < p2){
					pPair bp;
					bp.p1 = p;
					bp.p2 = p2;
					coll.push_back(bp);
				}
			}
		}
	}
}

void OcTree::PotentialPPColl(std::vector<pPair> &pColl, std::vector<Particle> &p, OcTree* ot){
	ot->PotentialPPColl(pColl);
}

void OcTree::PotentialPWColl(std::vector<pwPair> &coll, Wall w, char coord, int dir){
	if(hasChildren){
		for(int d2 = 0; d2 < 2; d2++){
			for(int d3 = 0; d3 < 2; d3++){
				OcTree* child;
				switch(coord){
				case 'x':
					child = Children[dir][d2][d3];
					break;
				case 'y':
					child = Children[d2][dir][d3];
					break;
				case 'z':
					child = Children[d2][d3][dir];
					break;
				}
				child->PotentialPWColl(coll, w, coord, dir);
			}
		}
	}else{
		for(set<Particle*>::iterator it = oParticles.begin(); it != oParticles.end(); it++){
			Particle *p = *it;
			pwPair bwp;
			bwp.p = p;
			bwp.wall = w;
			coll.push_back(bwp);
		}
	}
}

void OcTree::PotentialPWColl(std::vector<pwPair> &coll){
	PotentialPWColl(coll, W_LEFT,  'x',  0);
	PotentialPWColl(coll, W_RIGHT, 'x',  1);
	PotentialPWColl(coll, W_BOTTOM,'y',  0);
	PotentialPWColl(coll, W_TOP,   'y',  1);
	PotentialPWColl(coll, W_BACK,  'z',  0);
	PotentialPWColl(coll, W_FRONT, 'z',  1);

}

//void OcTree::PotentialPWColl(std::vector<pwPair> &pColl, std::vector<Particle> &p, OcTree* ot){
//	ot->PotentialPWColl(pColl);
//}

Vec3f OcTree::GetWallDir(Wall wall){
	switch(wall){
	case W_LEFT:
		return Vec3f(-1.0f, 0.0f, 0.0f);
		break;
	case W_RIGHT:
		return Vec3f(1.0f, 0.0f, 0.0f);
		break;
	case W_BACK:
		return Vec3f(0.0f, 0.0f, -1.0f);
		break;
	case W_FRONT:
		return Vec3f(0.0f, 0.0f, 1.0f);
		break;
	case W_TOP:
		return Vec3f(0.0f, 1.0f, 0.0f);
		break;
	case W_BOTTOM:
		return Vec3f(0.0f, -1.0f, 0.0f);
		break;
	default:
		return Vec3f(0.0f, 0.0f, 0.0f);
		break;
	}
}




































////Initialise the object and set all pertinent data.
//OcTree::OcTree(const Vec3f& o, const Vec3f& hd){
//	Origin	= o;
//	hDim	= hd;
//	Data	= NULL;
//	for(size_t i = 0; i < 8; i++){
//		Children[i] = NULL;
//	}
//}
//
//OcTree::OcTree(const OcTree& c){
//	Origin	= c.Origin;
//	hDim	= c.hDim;
//	Data	= c.Data;
//}
//
//OcTree::~OcTree(){
//	for(size_t i = 0; i < 8; i ++){
//		delete Children[i];
//	}
//}
//
//int OcTree::GetOctantContainingPoint(const Vec3f& p){
//	int oct = 0;
//	if(p.x > Origin.x) oct |= 4;
//	if(p.y > Origin.y) oct |= 2;
//	if(p.z > Origin.z) oct |= 1;
//	return oct;
//}
//
//bool OcTree::IsLeafNode(){
//	return Children[0] == NULL;
//}
//
//void OcTree::Insert(OcPoint* ocp){
//
//	if(IsLeafNode()){
//		if(Data == NULL){
//		Data = ocp;
//		return;
//	}else{
//		OcPoint* oldPoint	= Data;
//		Data				= NULL;
//		for(size_t i = 0; i < 8; i++){
//			Vec3f nOrigin	 = Origin;
//			nOrigin.x		+= hDim.x * (i&4 ? .5f : -.5f);
//			nOrigin.y		+= hDim.y * (i&2 ? .5f : -.5f);
//			nOrigin.z		+= hDim.z * (i&1 ? .5f : -.5f);
//			Children[i]		 = new OcTree(nOrigin, hDim * .5f);
//		}
//
//		Children[GetOctantContainingPoint(oldPoint->GetPos())]->Insert(oldPoint);
//		Children[GetOctantContainingPoint(ocp->GetPos())]->Insert(ocp);
//		}
//	}else{
//		int octant = GetOctantContainingPoint(ocp->GetPos());
//		Children[octant]->Insert(ocp);
//	}
//
//}
//
//void OcTree::GetPointsInCube(const Vec3f& bmin, const Vec3f& bmax, std::vector<OcPoint*>& results){
//
//	if(IsLeafNode()){
//		if(Data != NULL){
//			const Vec3f& p = Data->GetPos();
//			if(p.x > bmax.x || p.y > bmax.y || p.z > bmax.z) return;
//			if(p.x < bmin.x || p.y < bmin.y || p.z < bmin.z) return;
//			results.push_back(Data);
//		}
//	} else {
//		for(size_t i = 0; i < 8; i++){
//
//			Vec3f cmax = Children[i]->Origin + Children[i]->hDim;
//			Vec3f cmin = Children[i]->Origin + Children[i]->hDim;
//
//			if(cmax.x < bmin.x || cmax.y < bmin.y || cmax.z < bmin.z) continue;
//			if(cmin.x > bmax.x || cmin.y > bmax.y || cmin.z > bmax.z) continue;
//			
//			Children[i]->GetPointsInCube(bmin, bmax, results);
//		}
//	}
//
//}
//
//
//void OcTree::SetPoints(size_t it, std::vector<Particle>* p){
//	for(size_t i = 0; i < it; i++){
//		OcPoint* temp = new OcPoint();
//		temp->SetPos(p->at(i).GetPosition());
//		temp->SetIDX(i);
//		Insert(temp);
//		delete temp;
//	}
//}