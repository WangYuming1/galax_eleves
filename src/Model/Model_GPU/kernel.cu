#ifdef GALAX_MODEL_GPU

#include "cuda.h"
#include "kernel.cuh"
#define DIFF_T (0.1f)
#define EPS (1.0f)

__global__ void compute_acc(float3 * positionsGPU, float3 * velocitiesGPU, float3 * accelerationsGPU, float* massesGPU, int n_particles)
{
	 unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	 for(int j=0; j<n_particles;j++){
		if (i==j)
			continue;
		else{
			float3 distVector={0.0f, 0.0f,0.0f};
			distVector.x=positionsGPU[j].x-positionsGPU[i].x;
			distVector.y=positionsGPU[j].y-positionsGPU[i].y;
			distVector.z=positionsGPU[j].z-positionsGPU[i].z;
			float dijSqr=distVector.x*distVector.x+distVector.y*distVector.y+distVector.z*distVector.z;
			float dij=min(10.0f/(dijSqr * std::sqrt(dijSqr)),10.0f); 
			accelerationsGPU[i].x+=distVector.x*dij*massesGPU[j];
			accelerationsGPU[i].y+=distVector.y*dij*massesGPU[j];
			accelerationsGPU[i].z+=distVector.z*dij*massesGPU[j];
		}
	 }
}

__global__ void maj_pos(float3 * positionsGPU, float3 * velocitiesGPU, float3 * accelerationsGPU)
{
	 unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	 velocitiesGPU[i].x	+=	accelerationsGPU[i].x*2.0f;
	 velocitiesGPU[i].y	+=	accelerationsGPU[i].y*2.0f;
	 velocitiesGPU[i].z	+=	accelerationsGPU[i].z*2.0f;
	 positionsGPU[i].x	+=	velocitiesGPU[i].x*0.1f;
	 positionsGPU[i].y	+=	velocitiesGPU[i].y*0.1f;
	 positionsGPU[i].z	+=	velocitiesGPU[i].z*0.1f;
	
}

void update_position_cu(float3* positionsGPU, float3* velocitiesGPU, float3* accelerationsGPU, float* massesGPU, int n_particles)
{
	int nthreads = 128;
	int nblocks =  (n_particles + (nthreads -1)) / nthreads;

	compute_acc<<<nblocks, nthreads>>>(positionsGPU, velocitiesGPU, accelerationsGPU, massesGPU, n_particles);
	maj_pos    <<<nblocks, nthreads>>>(positionsGPU, velocitiesGPU, accelerationsGPU);
}


#endif // GALAX_MODEL_GPU