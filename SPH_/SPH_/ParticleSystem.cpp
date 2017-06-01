#include "ParticleSystem.h"

const char *getErrorString(cl_int error)
{
switch(error){
    // run-time and JIT compiler errors
    case 0: return "CL_SUCCESS";
    case -1: return "CL_DEVICE_NOT_FOUND";
    case -2: return "CL_DEVICE_NOT_AVAILABLE";
    case -3: return "CL_COMPILER_NOT_AVAILABLE";
    case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case -5: return "CL_OUT_OF_RESOURCES";
    case -6: return "CL_OUT_OF_HOST_MEMORY";
    case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case -8: return "CL_MEM_COPY_OVERLAP";
    case -9: return "CL_IMAGE_FORMAT_MISMATCH";
    case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case -11: return "CL_BUILD_PROGRAM_FAILURE";
    case -12: return "CL_MAP_FAILURE";
    case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
    case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
    case -15: return "CL_COMPILE_PROGRAM_FAILURE";
    case -16: return "CL_LINKER_NOT_AVAILABLE";
    case -17: return "CL_LINK_PROGRAM_FAILURE";
    case -18: return "CL_DEVICE_PARTITION_FAILED";
    case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

    // compile-time errors
    case -30: return "CL_INVALID_VALUE";
    case -31: return "CL_INVALID_DEVICE_TYPE";
    case -32: return "CL_INVALID_PLATFORM";
    case -33: return "CL_INVALID_DEVICE";
    case -34: return "CL_INVALID_CONTEXT";
    case -35: return "CL_INVALID_QUEUE_PROPERTIES";
    case -36: return "CL_INVALID_COMMAND_QUEUE";
    case -37: return "CL_INVALID_HOST_PTR";
    case -38: return "CL_INVALID_MEM_OBJECT";
    case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case -40: return "CL_INVALID_IMAGE_SIZE";
    case -41: return "CL_INVALID_SAMPLER";
    case -42: return "CL_INVALID_BINARY";
    case -43: return "CL_INVALID_BUILD_OPTIONS";
    case -44: return "CL_INVALID_PROGRAM";
    case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
    case -46: return "CL_INVALID_KERNEL_NAME";
    case -47: return "CL_INVALID_KERNEL_DEFINITION";
    case -48: return "CL_INVALID_KERNEL";
    case -49: return "CL_INVALID_ARG_INDEX";
    case -50: return "CL_INVALID_ARG_VALUE";
    case -51: return "CL_INVALID_ARG_SIZE";
    case -52: return "CL_INVALID_KERNEL_ARGS";
    case -53: return "CL_INVALID_WORK_DIMENSION";
    case -54: return "CL_INVALID_WORK_GROUP_SIZE";
    case -55: return "CL_INVALID_WORK_ITEM_SIZE";
    case -56: return "CL_INVALID_GLOBAL_OFFSET";
    case -57: return "CL_INVALID_EVENT_WAIT_LIST";
    case -58: return "CL_INVALID_EVENT";
    case -59: return "CL_INVALID_OPERATION";
    case -60: return "CL_INVALID_GL_OBJECT";
    case -61: return "CL_INVALID_BUFFER_SIZE";
    case -62: return "CL_INVALID_MIP_LEVEL";
    case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
    case -64: return "CL_INVALID_PROPERTY";
    case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
    case -66: return "CL_INVALID_COMPILER_OPTIONS";
    case -67: return "CL_INVALID_LINKER_OPTIONS";
    case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";

    // extension errors
    case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
    case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
    case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
    case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
    case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
    case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
    default: return "Unknown OpenCL error";
    }
}


///	///	///	///	///	///	///			KERNEL METHODS BEGINNETH		///	///	///	///	///	///	///

float ParticleSystem::KernelPoly(const float _sl){
	float weight = 0.0f;
	weight = m_weightPoly * (pow(pow(SMOOTHINGLENGTH, 2) - pow(_sl, 2), 3));
	return weight;
}

Vec3f ParticleSystem::KernelPolyGrad(const Vec3f _sv){
	Vec3f weight;
	float _sl = _sv.Length();
	weight = m_weightPolyGrad * _sv * pow(pow(SMOOTHINGLENGTH, 2) - pow(_sl, 2), 2);
	return weight;
}

float ParticleSystem::KernelPolyLap(const float _sl){
	float weight = 0.0f;
	weight = m_weightPolyLap * (pow(SMOOTHINGLENGTH, 2) - pow(_sl, 2)) * ((3.0f * pow(SMOOTHINGLENGTH, 2)) - (7.0f * pow(_sl, 2 )));
	return weight;
}

Vec3f ParticleSystem::KernelPresGrad(const Vec3f _sv){
	Vec3f weight;
	float _sl = _sv.Length();
	weight = m_weightPressGrad * pow(SMOOTHINGLENGTH - _sl, 2) * _sv.normalize();
	return weight;
}

float ParticleSystem::kernelViscLap(const float _sl){
	float weight = 0.0f;
	weight = m_weightViscLap * (SMOOTHINGLENGTH - _sl);
	return weight;
}

Vec3f ParticleSystem::KernelSpikyGrad(const Vec3f _sv){
	float _l = _sv.Length();
	float _t = SMOOTHINGLENGTH - _l;
	if (_l >= 0.0001f)
		_sv.normalize();
	Vec3f weight = (_t * _t) * _sv;
	return m_SpikyGrad * weight;

}

///	///	///	///	///	///				KERNEL METHODS ENDETH			///	///	///	///	///	///	///


//	This constructor instantiates some variables and ting
ParticleSystem::ParticleSystem(){
	dim			= container.getDims();
	container.genVAO();
	simPause	= true;
	colPressure	= true;
	
	p_Model = m_Loader.loadModel("models/sphere.obj", true);
	GenParticles();

	m_weightPoly		= (315.0f / (64.0f * PI * pow(SMOOTHINGLENGTH, 9)));
	m_weightPolyGrad	= (-945.0f / (32.0f * PI * pow(SMOOTHINGLENGTH, 9)));
	m_weightPolyLap		= (-945.0f / (32.0f * PI * pow(SMOOTHINGLENGTH, 9)));
	m_weightPressGrad	= (-45.0f / (PI * pow(SMOOTHINGLENGTH, 6)));
	m_SpikyGrad			= (45.0f / (PI * pow(SMOOTHINGLENGTH, 6)));
	m_weightViscLap		= (45.0f / (PI * pow(SMOOTHINGLENGTH, 6)));

	Integrator.SetiType(SemiImplicitEuler);
}

//	This destructor may eventually destruct things
ParticleSystem::~ParticleSystem(){}

//	This method calculates the forces exerted on a particle and updates the data 
void ParticleSystem::CalculateForces(){
	Vec3f _rv = Vec3f(0.0f, 0.0f, 0.0f);
	float _rf = 0.0f;
	float _mpd, _seplen;
	Vec3f _sepvec, _press, _visc;
	
	//Calculate the pressure and viscosity from the neighbors
	for(size_t i = 0; i < Particles.size(); i++){
		//Particle* _p = Particles[i];
		/// RESET ALL THE THINGS
		Particles[i].ResetForces();
		_mpd = _rf; _seplen = _rf;
		_sepvec = _rv; _press = _rv; _visc = _rv;


		for(size_t j = 0; j < __nnv.at(i).size(); j++){
			Particle _n = Particles[__nnv.at(i).at(j)->GetID()];
				_mpd		= _n.GetMass() / _n.GetDensity();
				_sepvec		= Particles[i].GetPosition() - _n.GetPosition();
				_seplen		= _sepvec.Length();
				//if(_n.GetID() != Particles[i].GetID() ){
				_press = _press + 
					(((Particles[i].GetPressure() / pow(Particles[i].GetDensity(), 2)) + 
					(_n.GetPressure() / pow(_n.GetDensity(), 2))) * 
					_n.GetMass() *
					KernelPresGrad(_sepvec));

				_visc = _visc + ( _mpd * 
					(_n.GetVelocity() - Particles[i].GetVelocity() ) * 
					kernelViscLap(_seplen));
				//}

		}

		Particles[i].SetPressureForce(-1.0f * Particles[i].GetDensity() * _press);
		Particles[i].SetViscosityForce(_visc * VISC_CONST);//  Particles[i].GetViscosity());

		//calculating gravity
		Particles[i].SetGravity(Vec3f(0.0f, Particles[i].GetDensity() * GRAVITYII, 0.0f));

		//calculate the acceleration from accumulated forces : force / mass
		Particles[i].AccumulateForces(Particles[i].GetPressureForce());
		Particles[i].AccumulateForces(Particles[i].GetViscosityForce());
		Particles[i].AccumulateForces(Particles[i].GetGravityForce());
		Particles[i].SetAccel(Particles[i].GetForce() / Particles[i].GetDensity());
	}
}

//	This preliminary force calculation method computes the inital force and density for use in the method below
void ParticleSystem::CalculatePressure(){

	//First calculate the density, then the pressure.
	for(size_t i = 0; i < Particles.size(); i ++){
		Particle* _p = &Particles[i];
		float theD = 0.0f; 
		for(size_t j = 0; j < __nnv.at(i).size(); j++){ 
			Particle* _n = &Particles[__nnv.at(i).at(j)->GetID()];
			theD += _n->GetMass() * KernelPoly((_p->GetPosition() - _n->GetPosition()).Length());

			_p->SetDensity(theD);
			//set the pressure as the ideal gas state equation
			float _press = GAS_CONSTANT * (_p->GetDensity() - MAT_DENSITY);
			_p->SetPressure(_press);
		}
	}
	
}

void ParticleSystem::CLGetNeighbours(size_t root){
	///get all platforms and choose available one
	cl_uint			numPlatforms;
	cl_platform_id	platform = NULL;
	cl_int			status = clGetPlatformIDs(0, NULL, &numPlatforms);
	//error check
	if(status != CL_SUCCESS)
		printf("CL_ERROR, could not get platforms\n");
	//choose first available one
	if(numPlatforms > 0){
		cl_platform_id* platforms = (cl_platform_id*)malloc(numPlatforms* sizeof(cl_platform_id));
		status = clGetPlatformIDs(numPlatforms, platforms, NULL);
		platform = platforms[0];
		free(platforms);
	}
	/// query platform and choose first GPU if present, else use CPU
	cl_uint			numDevices;
	cl_device_id	*devices;
	status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
	if(numDevices == 0){ //ie no GPU
		printf("OpenCL - No GPU available\n");
		printf("OpenCL - Choosing CPU as default device\n");
		status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 0, NULL, &numDevices);
		devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
		status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, numDevices, devices, NULL);
	}else{
		devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
		status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numDevices, devices, NULL);
	}
	/// create context
	cl_context context = clCreateContext(NULL, 1, devices, NULL, NULL, NULL);
	/// create commanding queue associated with context
	cl_command_queue commandQueue = clCreateCommandQueue(context, devices[0], 0, NULL);
	///create program object
	string sourceStr = ConvertString("nn.cl");
	//convert kernel file to string
	const char *source = sourceStr.c_str();
	size_t sourceSize[] = {strlen(source)};
	cl_program program = clCreateProgramWithSource(context, 1, &source, sourceSize, NULL);

	/// build the program
	status = clBuildProgram(program, 1, devices, NULL, NULL, NULL);
	if(status != CL_SUCCESS){
		cout<< getErrorString(status)<<endl;
	}
	/// initial input, output for the host and create mem objects for the kernel
	if(p_Positions.size() == 0)
		p_Positions = GetPositions();

	cl_mem p_InputBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(cl_float4) * p_Positions.size(), &p_Positions[0], NULL);
	size_t p_NeighbourRetSize = 500;
	size_t* p_Neighbours = (size_t*)malloc(p_NeighbourRetSize);
	cl_mem p_OutputBuffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, p_NeighbourRetSize * sizeof(size_t), NULL, NULL);
	glm::vec3 p_Root = glm::vec3(Particles[root].GetPosition().x, Particles[root].GetPosition().y, Particles[root].GetPosition().z);
	cl_mem p_RootBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(cl_float4) * sizeof(p_Root), &p_Root, NULL);
	///create kernel object
	cl_kernel kernel = clCreateKernel(program, "GetNearestNeighbour", NULL);

	/// set kernel arguments
	status = clSetKernelArg(kernel, 0, sizeof(size_t), (void*)&p_RootBuffer);
	status = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&p_InputBuffer);
	status = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&p_OutputBuffer);

	/// run the kernel
	size_t global_work_size[1] = {p_Positions.size()};
	status = clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL, global_work_size, NULL, 0, NULL, NULL);

	///read the output back to host mem
	status = clEnqueueReadBuffer(commandQueue, p_OutputBuffer, CL_TRUE, 0, p_NeighbourRetSize * sizeof(size_t), p_Neighbours, 0, NULL, NULL);

	///print out array elements if they are 1
	printf("Root - (%f, %f, %f)\n", p_Root.x, p_Root.y, p_Root.z);
	for(size_t i = 0; i < p_NeighbourRetSize; i++){
		if(p_Neighbours[i] == 1){
			printf("Neighbour %i - (%f, %f, %f)\n",i, Particles[i].GetPosition().x, Particles[i].GetPosition().y, Particles[i].GetPosition().z);
			if(i == root)
				SetCol(i, true);
			else
				SetCol(i, false);
		}
	}

	///clean up
	status = clReleaseKernel(kernel);
	status = clReleaseProgram(program);
	status = clReleaseMemObject(p_InputBuffer);
	status = clReleaseMemObject(p_OutputBuffer);
	status = clReleaseCommandQueue(commandQueue);
	status = clReleaseContext(context);
	if(p_Neighbours != NULL){
		free(p_Neighbours);
		p_Neighbours = NULL;
	}
	if(devices != NULL){
		free(devices);
		devices = NULL;
	}
	printf("finished\n");
}

//OpenCL convert string method
string ParticleSystem::ConvertString(const char* _fn){
	//const char *filename = "nn.cl";
	string sourceStr;
	size_t size; 
	char* str;
	std::fstream f(_fn, (std::fstream::in | std::fstream::binary));
	if(f.is_open()){
		size_t fileSize;
		f.seekg(0, std::fstream::end);
		size = fileSize = (size_t)f.tellg();
		f.seekg(0, std::fstream::beg);
		str = new char[size+1];
		if(!str){
			f.close();
			printf("OpenCL - ERROR - Could not convert .cl to string\n");
		}
		f.read(str, fileSize);
		f.close();
		str[size] = '\0';
		sourceStr = str;
		delete[] str;
	}else{printf("OpenCL - Error - could not open file");}
	return sourceStr;
}

//Generate the particles up to the particle count
void ParticleSystem::GenParticles(){

	//calculate mass
	float _m = MAT_DENSITY * (VOLUME / (float)MAX_PARTICLES);

	for(size_t i = 0; i < MAX_PARTICLES; i++){
		Particle temp; 
		temp.GenParticle(Vec3f(0.0f, 0.0f, 0.0f));
		temp.SetColour(0.0f, 0.0f, 1.0f);
		temp.SetModel(p_Model);
		temp.SetID(i);
		temp.SetVelocity(Vec3f(0.0f, 0.0f, 0.0f));
		temp.SetMass(_m);
		Particles.push_back(temp);
		i_Positions.push_back(glm::vec4(temp.GetPositionGLM(), 1.0f));
		i_Models.push_back(temp.GetMVM() );
		Vec3f t = temp.GetColour();
		i_Color.push_back(glm::vec3(t.x, t.y, t.z));
	}
	SetDistribution(inBB);	
}

//	This method utilises nanoflann to build a kd-tree and return a vector of neighbours
template <typename num_t>
void ParticleSystem::GetNeighbours(){
		PointCloud<num_t> cloud;
	addPtoPC(cloud, Particles);
	typedef KDTreeSingleIndexAdaptor<
		L2_Simple_Adaptor<num_t, PointCloud<num_t> >, 
		PointCloud<num_t>, 
		3
		>my_kd_tree_t;
	my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10) );
	index.buildIndex();
#if 0
	cloud.pts.resize(clouds.pts.size()*0.5);
	index.buildIndex();
#endif 
	/// Search for NN using the Particles vector
	{
		std::vector<Particle*> _tmp;
		__nnv.clear();
		//Vector to hold results
		for(size_t i = 0; i < Particles.size(); i++){
			_tmp.clear();
			Particle* _p = &Particles[i];
			const num_t query_pt[3] = {_p->GetPosition().x, _p->GetPosition().y, _p->GetPosition().z};
		
			const num_t sr = static_cast<num_t>(0.002); // should be equal to SMOOTHINGLENGTH ideally
			std::vector<std::pair<size_t, num_t> > retmatches;
			nanoflann::SearchParams params;
	
			const size_t nMatches = index.radiusSearch(&query_pt[0], sr, retmatches, params);
			
			for(size_t j = 0; j < nMatches; j++){
			
				if( Particles[retmatches[j].first].GetID() != _p->GetID() ) //ensure that the root particle is not added to neighbours
				_tmp.push_back(&Particles[retmatches[j].first]);
				//printf("nnMatches size = %i\n", nnMatches->size());
			}
			if(retmatches.size() > 0)
				__nnv.push_back(_tmp);
		}
	}
}

std::vector<glm::vec3> ParticleSystem::GetPositions(){
	std::vector<glm::vec3> temp;
	for(size_t i = 0; i < Particles.size(); i++){
		glm::vec3 tv = glm::vec3(Particles[i].GetPosition().x, Particles[i].GetPosition().y, Particles[i].GetPosition().z);
		temp.push_back(tv);
	}
	return temp;
}


//	This method checks to see whether a particular particle is within the bounding volume
void ParticleSystem::InBounds(Particle& p){
	Vec3f temp = p.GetPosition();
	float r = 0.03f; //p.GetRadius();
	
	//set the bounds to the BB dims
	float xmin = dim.min.x;
	float xmax = dim.max.x;
	float ymin = dim.min.y;
	float ymax = dim.max.y;
	float zmin = dim.min.z;
	float zmax = dim.max.z;

	// + || -  off the radi of the particles
	xmin += r;
	xmax -= r;
	ymin += r;
	ymax -= r;
	zmin += r;
	zmax -= r;

	if(temp.x  < xmin){ //if the particle x position is < the min boundary x position
		p.SetPosition(Vec3f(xmin, temp.y, temp.z));
		p.SetVelocity(Vec3f(REST_COEF * -p.GetVelocity().x,  p.GetVelocity().y, p.GetVelocity().z));
	}

	if(temp.x > xmax){ //if the particle x position is > the max boundary x position
		p.SetPosition(Vec3f(xmax, temp.y, temp.z));
		p.SetVelocity(Vec3f(REST_COEF * -p.GetVelocity().x,  p.GetVelocity().y, p.GetVelocity().z));
	}


	if(temp.y < ymin ){ //if the particle y position is < the min boundary y position
		p.SetPosition(Vec3f(temp.x, ymin, temp.z));
		p.SetVelocity(Vec3f(p.GetVelocity().x, REST_COEF * -p.GetVelocity().y, p.GetVelocity().z));
	}

	if(temp.y > ymax){ //if the particle y position is > the min boundary y position
		p.SetPosition(Vec3f(temp.x, ymax, temp.z));
		p.SetVelocity(Vec3f(p.GetVelocity().x, REST_COEF * -p.GetVelocity().y, p.GetVelocity().z));
	}


	if(temp.z < zmin){ //if the particle z position is < the min boundary z position
		p.SetPosition(Vec3f(temp.x, temp.y, zmin));
		p.SetVelocity(Vec3f(p.GetVelocity().x, p.GetVelocity().y, REST_COEF * -p.GetVelocity().z));
	}

	if(temp.z > zmax){ //if the particle z position is > the min boundary z position
		p.SetPosition(Vec3f(temp.x, temp.y, zmax));
		p.SetVelocity(Vec3f(p.GetVelocity().x, p.GetVelocity().y, REST_COEF * -p.GetVelocity().z));
	}

}

void ParticleSystem::InitCL(){
	///get all platforms and choose available one
	cl_uint			numPlatforms;
	cl_platform_id	platform = NULL;
	cl_int			status = clGetPlatformIDs(0, NULL, &numPlatforms);
	//error check
	if(status != CL_SUCCESS)
		printf("CL_ERROR, could not get platforms\n");
	//choose first available one
	if(numPlatforms > 0){
		cl_platform_id* platforms = (cl_platform_id*)malloc(numPlatforms* sizeof(cl_platform_id));
		status = clGetPlatformIDs(numPlatforms, platforms, NULL);
		platform = platforms[0];
		free(platforms);
	}
	/// query platform and choose first GPU if present, else use CPU
	cl_uint			numDevices;
	cl_device_id	*devices;
	status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
	if(numDevices == 0){ //ie no GPU
		printf("OpenCL - No GPU available\n");
		printf("OpenCL - Choosing CPU as default device\n");
		status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 0, NULL, &numDevices);
		devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
		status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, numDevices, devices, NULL);
	}else{
		devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
		status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numDevices, devices, NULL);
	}
	/// create context
	cl_context context = clCreateContext(NULL, 1, devices, NULL, NULL, NULL);
	/// create commanding queue associated with context
	cl_command_queue commandQueue = clCreateCommandQueue(context, devices[0], 0, NULL);
	///create program object
	string sourceStr = ConvertString("hello.cl");
	//convert kernel file to string
	const char *source = sourceStr.c_str();
	size_t sourceSize[] = {strlen(source)};
	cl_program program = clCreateProgramWithSource(context, 1, &source, sourceSize, NULL);

	/// build the program
	status = clBuildProgram(program, 1, devices, NULL, NULL, NULL);
	if(status != CL_SUCCESS){
		cout<< getErrorString(status)<<endl;
	}
	/// initial input, output for the host and create mem objects for the kernel
	const char* input = "GdkknVnqkc";
	size_t strLength = strlen(input);
	printf("OpenCL - Input string: %s\n", input);
	char* output = (char*) malloc(strLength + 1 );

	cl_mem inputBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, (strLength + 1)*sizeof(char), (void*) input, NULL);
	cl_mem outputBuffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, (strLength + 1) * sizeof(char), NULL, NULL);

	///creeate kernel object
	cl_kernel kernel = clCreateKernel(program, "helloworld", NULL);

	/// set kernel arguments
	status = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&inputBuffer);
	status = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&outputBuffer);

	/// run the kernel
	size_t global_work_size[1] = {strLength};
	status = clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL, global_work_size, NULL, 0, NULL, NULL);

	///read the output back to host mem
	status = clEnqueueReadBuffer(commandQueue, outputBuffer, CL_TRUE, 0, strLength * sizeof(char), output, 0, NULL, NULL);
	output[strLength] = '\0';
	printf("OpenCL - Output string: %s \n", output);

	///clean up
	status = clReleaseKernel(kernel);
	status = clReleaseProgram(program);
	status = clReleaseMemObject(inputBuffer);
	status = clReleaseMemObject(outputBuffer);
	status = clReleaseCommandQueue(commandQueue);
	status = clReleaseContext(context);
	if(output != NULL){
		free(output);
		output = NULL;
	}
	if(devices != NULL){
		free(devices);
		devices = NULL;
	}

}

//	This method resets the particle data, and the simulation
void ParticleSystem::ResetSimulation(){
	float _m = MAT_DENSITY * (VOLUME / (float)MAX_PARTICLES);
	
	SetDistribution(DistType);

	for(size_t i = 0; i < Particles.size(); i ++){
		Particles[i].SetColour(Vec3f(0.0f, 0.0f, 1.0f));
		Particles[i].SetVelocity(Vec3f(0.0f, 0.0f, 0.0f));
		Particles[i].SetMass(_m);
		Particles[i].ResetForces();
		Particles[i].SetDensity(0.0f);
		Particles[i].SetPressure(0.0f);
	}
	UpdateVecs();
}

//	This method is called when the simulation should run and is essentially the main fluid sim loop
void ParticleSystem::Run(float delta){
	/// Set the time step (may be able to be done outside of sim loop)
	Integrator.SetTimeStep(TIME_STEP);

	///		Build the kd-tree and get dems neighbours
	GetNeighbours<float>();

	///		First caclculate the pressure 
	CalculatePressure();

	///		Then calculate the forces
	CalculateForces();

	
	///		Integrate the movement
	for(size_t i = 0; i < Particles.size(); i++){
		Integrator.IntegrateNext(Particles[i]);
	}

	///		Check the particles are within the bounding structure
	for(size_t i = 0; i < Particles.size();i++)
		InBounds(Particles[i]);


	///		This will change the particle colours based upon pressure (differing shades of blue)
	SetParticlePCol();

	//finally update the vectors for the buffers
	UpdateVecs();
}

void ParticleSystem::SetCol(size_t _r, bool _is){
	if(_is)
		Particles[_r].SetColour(Vec3f(1.0f, 1.0f, 1.0f));
	else
		Particles[_r].SetColour(Vec3f(0.0f, 1.0f, 0.0f));

	UpdateVecs();
}

//This method sets the distribution of the particles, e.g. dam break, random etc
void ParticleSystem::SetDistribution(int type){
	float xmin = dim.min.x, xmax = dim.max.x, ymin = dim.min.y, ymax = dim.max.y, zmin = dim.min.z, zmax = dim.max.z;
	float PADDING		= 0.03f;	// particle scale * 2
	float PADDING_INIT	= 0.02f;	// particle scale


	if(type == 0){ //inBB
		DistType = 0;
		size_t titties = 0;
		float p2 = 0.04f;
		for(float i = ymax - p2; i >= ymin; i -= 0.03f){
			for(float j = zmin + p2; j <= (zmax - p2); j += 0.03f){
				for(float k = xmin + p2; k <= (xmax - p2); k+= 0.03f){
					if(titties <= Particles.size()){
						Particles[titties].SetPosition(Vec3f(k, i, j));
						if(k >= (xmax - 0.06f)){
							Particles[titties].SetPosition(Vec3f(k + .022f, i, j));
							if(j >= (zmax - .06f))
								Particles[titties].SetPosition(Vec3f(k + .022f, i, j +.022f));
						}
					}
					titties++;
				}				
			}
		}
	}
	if(type == 1){ //DamBreak
		DistType = 1;
		size_t titties = 0;
		float p2 = 0.04f;
		for(float k = xmin + p2; k <= (xmax - p2); k+= 0.03f){
			for(float j = zmin + p2; j <= (zmax - p2); j += 0.03f){
				for(float i = ymax - p2; i >= ymin + .5f; i -= 0.03f){
					if(titties <= Particles.size()){
						Particles[titties].SetPosition(Vec3f(k, i, j));
						if(k >= (xmax - 0.06f)){
							Particles[titties].SetPosition(Vec3f(k + .022f, i, j));
							if(j >= (zmax - .06f))
								Particles[titties].SetPosition(Vec3f(k + .022f, i, j +.022f));
						}
					}
					titties++;
				}				
			}
		}
	}

	if(type == 2){ //Random
		for(size_t i = 0; i < Particles.size(); i++){
			//float x, y, z;
			
		}
	}
	UpdateVecs();
}

//This method is used to set particle colours when pressure force exceeds a threshold
void ParticleSystem::SetParticlePCol(){
	if(colPressure){

		for(size_t i = 0; i < Particles.size(); i ++){
			float _tp = Particles[i].GetPressure();

			if(_tp > 10000 && _tp < 20000)
				Particles[i].SetColour(Vec3f(0.0f, 0.0f, 0.5f));
			else if(_tp > 20000)
				Particles[i].SetColour(Vec3f(0.8f, 0.0f, 0.0f));
			else
				Particles[i].SetColour(Vec3f(0.0f, 0.0f, 1.0f));
		}
	}else{
		for(size_t i = 0; i < Particles.size(); i ++){
			Particles[i].SetColour(Vec3f(0.0f, 0.0f, 1.0f));
		}
	}

}

//	This method updates the relevant vectors ready to be sent to the GPU buffers.
void ParticleSystem::UpdateVecs(){
	for(size_t i = 0; i < MAX_PARTICLES; i++){
		i_Positions.at(i) = glm::vec4(Particles[i].GetPositionGLM(), 1.0f);
		i_Models.at(i) = Particles[i].GetMVM();
		i_Color.at(i) = Particles[i].GetColourGLM();
	}
}
































///			OLD METHODS AND TING BELOW
////This method
//void ParticleSystem::MovePFromBounds(Particle* p){
//	Vec3f temp	= p->GetPosition();
//	float r		= p->GetRadius();
//	float m		= r/2;
//
//	if(temp.x + r > dim.max.x)
//		p->SetPosition(Vec3f(temp.x - m, temp.y, temp.z));
//	else if (temp.y + r*2 > dim.max.y)
//		p->SetPosition(Vec3f(temp.x, temp.y - m, temp.z));
//	else if (temp.z + r*2 > dim.max.z)
//		p->SetPosition(Vec3f(temp.x, temp.y, temp.z - m));
//	else if (temp.x - r*2 < dim.min.x)
//		p->SetPosition(Vec3f(temp.x + m, temp.y, temp.z));
//	else if (temp.y - r*2 < dim.min.y)
//		p->SetPosition(Vec3f(temp.x, temp.y + m, temp.z));
//	else if (temp.z - r*2 < dim.min.z)
//		p->SetPosition(Vec3f(temp.x, temp.y, temp.z + m));
//}


////Check the postion of the particle + or - the radius is within
////the bounding box of the container. If the particle is outside
////then move it and check again until it is inside the bounds
//void ParticleSystem::BoundaryCheck(){
//	for(size_t i = 0; i < Particles.size();i++){
//		Vec3f temp	= Particles[i].GetPosition();
//		float r		= Particles[i].GetRadius();
//		float p		= 0.03f;
//		if( !((temp.x < dim.max.x && temp.x > dim.min.x) &&
//			  (temp.y < dim.max.y - p && temp.y > dim.min.y ) &&
//			  (temp.z < dim.max.z && temp.z > dim.min.z)) ){
//				Vec3f dir;
//			if(temp.x < dim.min.x + p)
//				dir = Vec3f(-1.0f, 0.0f, 0.0f);
//			if(temp.x > dim.max.x - p)
//				dir = Vec3f(1.0f, 0.0f, 0.0f);
//			if(temp.y < dim.min.y + p)
//				dir = Vec3f(0.0f, -1.0f, 0.0f);
//			if(temp.y > dim.max.y - p)
//				dir = Vec3f(0.0f, 1.0f, 0.0f);
//			if(temp.z < dim.min.z + p)
//				dir = Vec3f(0.0f, 0.0f, -1.0f);
//			if(temp.z > dim.max.z - p)
//				dir = Vec3f(0.0f, 0.0f, 1.0f);
//
//			dir.normalize();
//			Vec3f temp2 = Particles[i].GetVelocity();
//
//
//			temp2 -=   dir *2.0f * Particles[i].GetVelocity().dot(dir);
//			Particles[i].SetVelocity(temp2);
//			Particles[i].SetFalling(false);
//
//			Particles[i].SetBounds(false);
//		}else{
//			Particles[i].SetBounds(true);
//		}
//	}
//}



////This method updates the position of the particle after gravity has been applied etc
//void ParticleSystem::MoveParticles(float dt){
//	for(size_t i = 0; i < Particles.size(); i++){
//		Particle* p = &Particles[i];
//		Vec3f op = p->GetPosition();
//		Vec3f temp = p->GetPosition() -= p->GetVelocity() * dt;
//		p->SetPosition(temp);
//	}
//}



////This method applies the gravity to the particles in the simulation
//void ParticleSystem::ApplyGravity(float dt){
//	for(size_t i = 0; i < Particles.size(); i++){
//		Particle* p = &Particles[i];
//		if(p->GetFalling())
//			p->SetVelocity(Vec3f(0.0f, GRAVITY * dt ,0.0f));
//	}
//}


//void ParticleSystem::CollisionResolvePVP(Particle& p, std::vector<Particle>& tp){
//	//loop over each instance in the vector and test it against the test particle for collisions
//	for(size_t i = 0; i <tp.size(); i++){
//		if(PPCollision(p, tp[i])){
//			Vec3f disp = (p.GetPosition() - tp[i].GetPosition()).normalize();
//			Vec3f t = p.GetVelocity();
//			t -= 2 * disp * p.GetVelocity().dot(disp);
//			p.SetVelocity(  t);
//
//			t = tp[i].GetVelocity();
//			t -= 2 * disp * tp[i].GetVelocity().dot(disp);
//			tp[i].SetVelocity(t);
//		}
//	}
//
//
//	//float r = p.GetRadius() + tp[i].GetRadius();
//	//	if( (p.GetPosition() - tp[i].GetPosition()).magnitudeSquared() < r * r ){
//	//		Vec3f netVel = p.GetVelocity() - tp[i].GetVelocity();
//	//		Vec3f disp   = p.GetPosition() - tp[i].GetPosition();
//	//		//if true there's been a collision
//	//		if(netVel.dot(disp) < 0){
//	//			Vec3f disp2 = (p.GetPosition() - tp[i].GetPosition()).normalize();
//	//			Vec3f p1 = p.GetVelocity();
//	//			p1 -=  disp2 * 2.0f * p.GetVelocity().dot(disp2);
//	//			p.SetVelocity(p1);
//	//			p1 = tp[i].GetVelocity();
//	//			p1 -=   disp2 * 2.0f * tp[i].GetVelocity().dot(disp2);
//	//			tp[i].SetVelocity(p1);
//	//		}
//	//	}
//}
//
//bool ParticleSystem::PPCollision(Particle& p1, Particle& p2){
//	Vec3f p1p = p1.GetPosition();
//	Vec3f p2p = p2.GetPosition();
//	float r = p1.GetRadius() + p2.GetRadius();
//
//	if( (p1p - p2p).magnitudeSquared() < r * r){
//		Vec3f netVel = p1.GetVelocity() - p2.GetVelocity();
//		Vec3f disp   = p1.GetPosition() - p2.GetPosition();
//		return netVel.dot(disp) < 0;
//	}else
//		return false;
//}

//template <typename num_t>
//void ParticleSystem::BuildKD(std::vector<Particle>& p){
//	
//	PointCloud<num_t> cloud;
//	addPtoPC(cloud, p);
//	typedef KDTreeSingleIndexAdaptor<
//		L2_Simple_Adaptor<num_t, PointCloud<num_t> >, 
//		PointCloud<num_t>, 
//		3
//		>my_kd_tree_t;
//	my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10) );
//	index.buildIndex();
//#if 0
//	cloud.pts.resize(clouds.pts.size()*0.5);
//	index.buildIndex();
//#endif 
//
//
//
//	/// Search for NN using the Particles vector
//	{
//		//Vector to hold results
//		std::vector<Particle> nnMatches;
//		PNPairs.clear();
//		for(size_t i = 0; i < p.size(); i++){
//			nnMatches.clear();
//			const num_t query_pt[3] = {p.at(i).GetPosition().x, p.at(i).GetPosition().y, p.at(i).GetPosition().z};
//			const num_t sr = static_cast<num_t>(0.002);
//			std::vector<std::pair<size_t, num_t> > retmatches;
//			nanoflann::SearchParams params;
//			const size_t nMatches = index.radiusSearch(&query_pt[0], sr, retmatches, params);
//			for(size_t j = 0; j < nMatches; j++){
//				if( retmatches[j].first != i ) //ensure that the root particle is not added to neighbours
//					nnMatches.push_back(p.at(retmatches[j].first));
//			}
//			//make a pair of the current particle and all it's neighbours!
//			if(nnMatches.size()>0){
//				std::pair <Particle, std::vector<Particle>> temp;
//				temp = std::make_pair (p.at(i), nnMatches);
//				PNPairs.push_back(temp);
//			}
//		}
//
//	}
//	
//		//{   //Tings used in search
//	//	size_t nR = 20;
//	//	size_t rng = rand() %MAX_PARTICLES;
//	//	const num_t query_pt[3] = {p.at(rng).GetPosition().x, p.at(rng).GetPosition().y, p.at(rng).GetPosition().z};
//	//	/// find n neighbours of a particle
//	//	//{
//	//	//	
//	//	//	std::vector<size_t> ret_index(nR);
//	//	//	std::vector<num_t> out_dist_sq(nR);
//	//	//	index.knnSearch(&query_pt[0], nR, &ret_index[0], &out_dist_sq[0]);
//	//	//	cout<< "knnSearch(): num results="<<nR<<"\n";	
//	//	//	for(size_t i = 0 ; i < nR;i++){
//	//	//		cout<<"TestPoint = "<<rng<<" | idx["<<i<<"]="<<ret_index[i]<<" dist["<<i<<"]"<<out_dist_sq[i]<<endl;
//	//	//		p.at(ret_index[i]).SetColour(Vec3f(1.0f, 1.0f, 1.0f));
//	//	//	}
//	//	//	cout<<"\n";
//	//	//}
//	//
//
//		/// test within a radius
//		//{
//		//	size_t nR = 20;
//		//	size_t rng = rand() %MAX_PARTICLES;
//		//	const num_t query_pt[3] = {p.at(rng).GetPosition().x, p.at(rng).GetPosition().y, p.at(rng).GetPosition().z};
//		//	const num_t sr = static_cast<num_t>(0.002);
//		//	std::vector<std::pair<size_t, num_t> > retmatches;
//		//	nanoflann::SearchParams params;
//		//	const size_t nMatches = index.radiusSearch(&query_pt[0], sr, retmatches, params);
//		//	cout<<"radiusSearch(): radius ="<<sr<< " -> "<<nMatches<<" matches\n";
//		//	for(size_t i =0; i < nMatches;i++){
//		//		cout<<"TestPoint = "<<rng<<" | idx["<<i<<"]="<<retmatches[i].first<<" dist["<<i<<"]="<< retmatches[i].second<<endl;
//		//		p.at(retmatches[i].first).SetColour(Vec3f(1.0f, 1.0f, 0.0f));
//		//		//p.at(retmatches[i].second).SetColour(Vec3f(1.0f, 1.0f, 0.0f));
//		//	}
//		//	cout<<"\n";
//		//
//
//		////Colour the seed particles GREEN
//		//p.at(rng).SetColour(Vec3f(0.0f, 1.0f, 0.0f));
//		//}
//	
//}





//
//void ParticleSystem::CalculateForces(std::vector<std::pair<Particle, std::vector<Particle>>>* PandN){
//	Vec3f _rv = Vec3f(0.0f, 0.0f, 0.0f);
//	float _rf = 0.0f;
//	float _mpd, _seplen, _surflap;
//	Vec3f _sepvec, _press, _visc;
//	
//	//Calculate the pressure and viscosity from the neighbors
//	for(size_t i = 0; i < PandN->size(); i++){
//		Particle _t = PandN->at(i).first;
//		_t.ResetForces();
//		//reset the tings
//		_mpd		= _rf;
//		_seplen		= _rf;
//		_surflap	= _rf;
//		_sepvec		= _rv;
//		_press		= _rv;
//		_visc		= _rv;
//
//		for(size_t j = 0; j < PandN->at(i).second.size(); j++){
//			Particle _n = PandN->at(i).second.at(j);
//
//			_mpd		= _n.GetMass() / _n.GetDensity();
//			_sepvec		= _t.GetPosition() - _n.GetPosition();
//			_seplen		= _sepvec.Length();
//			if(_t.GetID() != _n.GetID() ){
//				_press = _press + (((_t.GetPressure() / pow(_t.GetDensity(), 2)) + 
//					(_n.GetPressure() / pow(_n.GetDensity(), 2))) * _n.GetMass() *KernelPresGrad(_sepvec));
//
//				_visc = _visc + ( _mpd * (_n.GetVelocity() - _t.GetVelocity() ) * kernelViscLap(_seplen));
//			}
//
//		}
//		
//		Particles[i].SetPressureForce(-1.0f * Particles[i].GetDensity() * _press);
//		Particles[i].SetViscosityForce(_visc * VISC_CONST);
//
//		//calculating gravity - density * grav
//		Particles[i].SetGravity(Vec3f(0.0f, Particles[i].GetDensity() * GRAVITYII, 0.0f));
//
//		//calculate the acceleration from accumulated forces : force / mass
//		Particles[i].AccumulateForces(Particles[i].GetPressureForce());
//		Particles[i].AccumulateForces(Particles[i].GetViscosityForce());
//		Particles[i].AccumulateForces(Particles[i].GetGravityForce());
//		Particles[i].SetAccel(Particles[i].GetForce() / Particles[i].GetDensity());
//	}
//}


//void ParticleSystem::CalculatePressure(std::vector<std::pair<Particle, std::vector<Particle>>>* PandN){
//
//	//First calculate the density, then the pressure.
//	for(size_t i = 0; i < PandN->size(); i ++){
//		float theD = 0.0f; 
//		//Particle* t = &PandN->at(i).first;
//		
//		for(size_t j = 0; j < PandN->at(i).second.size(); j++){
//			theD += PandN->at(i).second.at(j).GetMass() * 
//				KernelPoly((PandN->at(i).first.GetPosition() - PandN->at(i).second.at(j).GetPosition()).Length());
//		}
//		//PandN->at(i).first.SetDensity(theD);
//		Particles[i].SetDensity(theD);
//		//set the pressure as the ideal gas state equation
//		float _press = GAS_CONSTANT * (Particles[i].GetDensity() - MAT_DENSITY);///   REST_COEF);
//		//PandN->at(i).first.SetPressure(_press);
//		Particles[i].SetPressure(_press);
//	}
//}




/*
bool ParticleSystem::ppCollTest(Particle* p1, Particle* p2){
	float r = p1->GetRadius() + p2->GetRadius();
	if( (p1->GetPosition() - p2->GetPosition()).magnitudeSquared() < r * r ){
		Vec3f netVel = p1->GetVelocity() - p2->GetVelocity();
		Vec3f disp = p1->GetPosition() - p2->GetPosition();
		return netVel.dot(disp) < 0;
	}else
		return false;
}

void ParticleSystem::ppCollResponse(OcTree* ot){
	vector<pPair> bps;
	ot->PotentialPPColl(bps, Particles, ot);
	for(size_t i = 0; i < bps.size();i++){
		pPair bp = bps[i];
		Particle* p1 = bp.p1;
		Particle* p2 = bp.p2;
		if(ppCollTest(p1, p2)){
			Vec3f disp = (p1->GetPosition() - p2->GetPosition()).normalize();
			Vec3f temp = p1->GetVelocity();
			temp -= disp * 2.0f * p1->GetVelocity().dot(disp);
			p1->SetVelocity(temp);
			temp = p2->GetVelocity();
			temp -= disp * 2.0f * p2->GetVelocity().dot(disp);
			p2->SetVelocity(temp);
		}
	}
}



void ParticleSystem::PColCheck(Particle* p, std::vector<Particle>* pv){

	for(size_t i = 0; i < pv->size(); i ++){
		Vec3f coll = pv->at(i).GetPosition() - p->GetPosition();
		float dist = coll.Length();
		float radSum = pv->at(i).GetRadius() + p->GetRadius();
		if(dist < radSum){ //there's been a collision
			printf("collision\n");
		}
	}
}






void ParticleSystem::MoveParticles(OcTree* oct, float dt){
	for(size_t i = 0; i < Particles.size(); i++){
		Particle* p = &Particles[i];
		Vec3f op = p->GetPosition();
		Vec3f temp = p->GetPosition() -= p->GetVelocity() * dt;
		p->SetPosition(temp);
		//ot->ParticleMoved(p, op);
	}
}

//Same as neighbour search above but return the vector so that the neighbours can be dealt with
std::vector<Particle> ParticleSystem::NSVec(Particle* p){

	std::vector<Particle> Neighbours;
	if(Neighbours.size() > 1) //to make sure it's clear (so we don't use up loads of memory)
		Neighbours.clear();

	size_t c = 0;
	Vec3f ss1, ss2, pp;
	size_t j;
	float r = SEARCH_RADIUS; //the search radius
	pp = p->GetPosition();
	
	ss1 = Vec3f(pp.x - r, pp.y -r, pp.z - r); //space search MIN dims
	ss2 = Vec3f(pp.x + r, pp.y + r, pp.z + r);//space search MAX dims
	pTree->KDTreeSearchSpace(aTree, &aIter, ss1.x, ss2.x, ss1.y, ss2.y, ss1.z, ss2.z);


	j = pTree->kdIterGetNext(aIter);

	while(j != KD_TREE_END){
		//do something here
		//Particles.at(j).SetColour(Vec3f(1.0f, 1.0f, 1.0f));
		Neighbours.push_back(Particles.at(j));
		c++;
		j = pTree->kdIterGetNext(aIter);
	}
	return Neighbours;
}



void ParticleSystem::BuildKD(){
	aTree = new kdtree();
	pTree->BuildKDTree(&Particles, Particles.size(), &aTree);
}



void ParticleSystem::NeighbourSearch(Particle* p){
	//aTree = new kdtree();
	//aTree = pTree->BuildKDTree(&Particles, Particles.size(), &aTree);

	//std::vector<Particle> Neighbours;
	size_t c = 0;
	Vec3f ss1, ss2, pp;
	size_t j;
	float r = 0.03f; //the search radius
	pp = p->GetPosition();
	
	ss1 = Vec3f(pp.x - r, pp.y -r, pp.z - r); //space search MIN dims
	ss2 = Vec3f(pp.x + r, pp.y + r, pp.z + r);//space search MAX dims
	pTree->KDTreeSearchSpace(aTree, &aIter, ss1.x, ss2.x, ss1.y, ss2.y, ss1.z, ss2.z);
	j = pTree->kdIterGetNext(aIter);

	while(j != KD_TREE_END){
		//do something here
		//Particles.at(j).SetColour(Vec3f(1.0f, 1.0f, 1.0f));
		//Neighbours.push_back(Particles.at(j));
		c++;
		j = pTree->kdIterGetNext(aIter);
	}
	
	printf(" Neighbours? %i\n", c);

}





void ParticleSystem::UpdateNeighbours(){}




bool ParticleSystem::pwCollTest(Particle* p, Wall w){
	Vec3f dir = ot->GetWallDir(w);
	return p->GetPosition().dot(dir) + p->GetRadius() > (dim.max.y - dim.min.y)  && p->GetVelocity().dot(dir) > 0;
}

void ParticleSystem::pwCollResponse(OcTree* oc){
	vector<pwPair> bwps;
	ot->PotentialPWColl(bwps);
	for(size_t i = 0; i < bwps.size(); i++){
		pwPair bwp = bwps[i];
		Particle* p = bwp.p;
		Wall w		= bwp.wall;
		if(pwCollTest(p, w)){
			Vec3f dir = (ot->GetWallDir(w)).normalize();
			Vec3f temp = p->GetVelocity();
			temp -= dir * 2.0f * p->GetVelocity().dot(dir);
			p->SetVelocity(temp);
		}
	}
}






void ParticleSystem::BuildOT(){
	//if(ot == NULL){
	//	ot = new OcTree(Vec3f(0.0f, 0.0f, 0.0f), Vec3f(1.0f, 1.0f, 1.0f));
	//	ot->SetPoints(Particles.size(), &Particles);
	//}else{
	//	ot = NULL;
	//	ot = new OcTree(Vec3f(0.0f, 0.0f, 0.0f), Vec3f(1.0f, 1.0f, 1.0f));
	//	ot->SetPoints(Particles.size(), &Particles);
	//}

	ot = new OcTree(dim.min, dim.max, 1);

}

void ParticleSystem::QueryOTNN(Particle* p){
	Vec3f pp  = p->GetPosition();
	Vec3f min = Vec3f(pp.x - (-.2f), pp.y -(-.2f), pp.z - (-.2f));
	Vec3f max = Vec3f(pp.x + .2f, pp.y + .2f, pp.z +.2f);
	std::vector<OcPoint*> results;
	ot->GetPointsInCube(min, max, results);
	printf("found %i neighbours\n", results.size());

}

*/
