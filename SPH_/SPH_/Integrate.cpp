#include "Integrate.h"

Integrate::Integrate(){}

Integrate::Integrate(const iType _it, const float ts){
	m_iType		= _it;
	m_TimeStep	= ts;
}

void Integrate::IntegrateNext(Particle& _p){
	switch(m_iType){
	case SemiImplicitEuler: { EvaluateSIE(_p); break;}
	case LeapFrog:			{ EvaluateLF(_p); break;}						
	default:
		break;	
	}
}

void Integrate::EvaluateSIE(Particle& _p){
	Vec3f velo, pos;
	velo	= _p.GetVelocity() + (_p.GetAcceleration() * m_TimeStep);
	pos		= _p.GetPosition() + (_p.GetVelocity() * m_TimeStep);
	_p.SetVelocity(velo);
	_p.SetPosition(pos);
}

void Integrate::EvaluateLF(Particle& _p){
	Vec3f velo, pos;
	velo	= _p.GetVelocity() + (((_p.GetPAcceleration() + _p.GetAcceleration()) / 2.0f ) * m_TimeStep);
	pos		= _p.GetPosition() + (_p.GetVelocity() * m_TimeStep) + ((_p.GetPAcceleration() / 2.0f ) * m_TimeStep *m_TimeStep);
	_p.SetVelocity(velo);
	_p.SetPosition(pos);
}