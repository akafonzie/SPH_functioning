#pragma OPENCL EXTENSION cl_amd_printf :enable

__kernel void 
GetNearestNeighbour(__global float4* root, __global float4* p_Array, __global unsigned int* n_Array)
{
	int gID = get_global_id(0);
	float4 testP = p_Array[gID];

	float dist = sqrt( ( 
						 pow((testP.s0 - root[0].s0), 2) +
						 pow((testP.s1 - root[0].s1), 2) + 
						 pow((testP.s2 - root[0].s2), 2) 
					   ) 
					 );
	if(dist < 0.3){
		n_Array[gID] = 1;
	}
}

