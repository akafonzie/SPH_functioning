#ifndef Integrate_h
#define Integrate_h

#include "Particle.h"
#include "Vec3f.h"

enum iType{
	SemiImplicitEuler,
	LeapFrog
};


class Integrate{
public:
	Integrate();
	Integrate(const iType, const float);
	void IntegrateNext(Particle&);
	iType GetiType(){return m_iType;}
	void SetiType(const iType _it){m_iType = _it;}
	float GetTimeStep(){return m_TimeStep;}
	void SetTimeStep(float ts){m_TimeStep = ts;}
	
private:
	iType m_iType;
	float m_TimeStep;
	void EvaluateSIE(Particle&);
	void EvaluateLF(Particle&);
};
#endif