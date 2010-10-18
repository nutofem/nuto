// $Id$ 
#ifndef VIRTUALFUNCTIONSVSDYNAMICCASTCLASSES_H_
#define VIRTUALFUNCTIONSVSDYNAMICCASTCLASSES_H_

typedef double NuToVec2[2];
typedef double NuToVec3[3];

class BaseA;
class BaseB;
class AB;

class Base
{
public:
	virtual void FunctionA(NuToVec2 rValue)
	{
		throw 1;
	}
	virtual void FunctionB(NuToVec3 rValue)
	{
		throw 1;
	}
	virtual BaseA* asA()
	{
		throw 1;
	}
	virtual BaseB* asB()
	{
	    throw 1;
	}
	virtual AB* asAB()
	{
	    throw 1;
	}
};

class BaseA : public virtual Base
{
public:
	BaseA()
	{
		mA[0]=0.;
		mA[1]=0.;
	}
	virtual void FunctionA(NuToVec2 rValue);
	virtual void FunctionB(NuToVec3 rValue);
	virtual BaseA* asA();

	NuToVec2 mA;
};

class BaseB : public virtual Base
{
public:
	BaseB()
	{
		mB[0]=1.;
		mB[1]=1.;
		mB[2]=1.;
	}
	virtual void FunctionB(NuToVec3 rValue);
	virtual BaseB* asB();

	NuToVec3 mB;
};

class AB : public BaseA, public BaseB
{
public:
	AB* asAB();
	void FunctionB(NuToVec3 rValue);
};

#endif /* VIRTUALFUNCTIONSVSDYNAMICCASTCLASSES_H_ */
