#ifdef GALAX_MODEL_CPU_FAST

#include <cmath>

#include "Model_CPU_fast.hpp"

#include <mipp.h>
#include <omp.h>

Model_CPU_fast
::Model_CPU_fast(const Initstate& initstate, Particles& particles)
: Model_CPU(initstate, particles)
{
}

void Model_CPU_fast
::step()
{
    std::fill(accelerationsx.begin(), accelerationsx.end(), 0);
    std::fill(accelerationsy.begin(), accelerationsy.end(), 0);
    std::fill(accelerationsz.begin(), accelerationsz.end(), 0);

/*
 #pragma omp parallel
 {
    #pragma omp for
        for (int i = 0; i < n_particles; i += mipp::N<float>())
        {
            for (int j = 0; j < n_particles; j += mipp::N<float>())
                    {
                        if (i==j) 
                            continue;
                        
                        
                        const float diffx = particles.x[j] - particles.x[i];
                        const float diffy = particles.y[j] - particles.y[i];
                        const float diffz = particles.z[j] - particles.z[i];

                        float dij = diffx * diffx + diffy * diffy + diffz * diffz;

                        if (dij < 1.0)
                        {
                            dij = 10.0;
                        }
                        else
                        {
                            
                            dij = 10.0 / (dij * std::sqrt(dij));
                        }

                        accelerationsx[i] += diffx * dij * initstate.masses[j];
                        accelerationsy[i] += diffy * dij * initstate.masses[j];
                        accelerationsz[i] += diffz * dij * initstate.masses[j];
                    }
            }	
		
	
	
		#pragma omp for
			for (int i = 0; i < n_particles; i++)
			{
				velocitiesx[i] += accelerationsx[i] * 2.0f;
				velocitiesy[i] += accelerationsy[i] * 2.0f;
				velocitiesz[i] += accelerationsz[i] * 2.0f;
				particles.x[i] += velocitiesx   [i] * 0.1f;
				particles.y[i] += velocitiesy   [i] * 0.1f;
				particles.z[i] += velocitiesz   [i] * 0.1f;
			}
}	


	
*/
  
    const mipp::Reg<float> rG=10.0f;
    auto vecLoopSize = (n_particles / mipp::N<float>()) * mipp::N<float>();
#pragma omp parallel
 {
    #pragma omp for
    for (int i = 0; i < vecLoopSize; i += mipp::N<float>())
    {
            // load registers body i
            const mipp::Reg<float> rposx_i = &particles.x[i];
            const mipp::Reg<float> rposy_i = &particles.y[i];
            const mipp::Reg<float> rposz_i = &particles.z[i];
                mipp::Reg<float> raccx_i = &accelerationsx[i];
                mipp::Reg<float> raccy_i = &accelerationsy[i];
                mipp::Reg<float> raccz_i = &accelerationsz[i];
                mipp::Reg<float> rvcx_i ;
                mipp::Reg<float> rvcy_i ;
                mipp::Reg<float> rvcz_i ;
  
        //for (int j = 0; j < vecLoopSize; j += mipp::N<float>()){
        for (int j = 0; j < vecLoopSize; j += mipp::N<float>()){
            if (i==j) 
                continue;
            const mipp::Reg<float> rposx_j = &particles.x[j];
            const mipp::Reg<float> rposy_j = &particles.y[j];
            const mipp::Reg<float> rposz_j = &particles.z[j];
            const mipp::Reg<float> diffx=rposx_j-rposx_i;
            const mipp::Reg<float> diffy=rposy_j-rposy_i;
            const mipp::Reg<float> diffz=rposz_j-rposz_i;   
            const mipp::Reg<float> rmass_j = &initstate.masses[j];
            mipp::Reg<float> dij;
            dij=diffx * diffx + diffy * diffy + diffz * diffz;
        
            dij=mipp::min(rG/(dij * mipp::sqrt(dij)),rG); 
           // dij =  rG/ (dij * mipp::sqrt(dij));
                
                // r1.store(&myVector[(i+1)*mipp::N<float>()]);

                raccx_i+=diffx*dij*rmass_j;
                raccy_i+=diffy*dij*rmass_j;
                raccz_i+=diffz*dij*rmass_j;
                raccx_i.store(&accelerationsx[i]);
                raccy_i.store(&accelerationsy[i]);
                raccz_i.store(&accelerationsz[i]);
        }

        
    }
    #pragma omp for
    for (int i = 0; i < n_particles; i++)
			{
				velocitiesx[i] += accelerationsx[i] * 2.0f;
				velocitiesy[i] += accelerationsy[i] * 2.0f;
				velocitiesz[i] += accelerationsz[i] * 2.0f;
				particles.x[i] += velocitiesx   [i] * 0.1f;
				particles.y[i] += velocitiesy   [i] * 0.1f;
				particles.z[i] += velocitiesz   [i] * 0.1f;
			}

 }




// OMP + MIPP version
 /*   #pragma omp parallel {
    #pragma omp for
        for (int i = 0; i < n_particles; i += mipp::N<float>())
        {
                // load registers body i
                const mipp::Reg<float> rposx_i = &particles.x[i];
                const mipp::Reg<float> rposy_i = &particles.y[i];
                const mipp::Reg<float> rposz_i = &particles.z[i];
                    mipp::Reg<float> raccx_i = &accelerationsx[i];
                    mipp::Reg<float> raccy_i = &accelerationsy[i];
                    mipp::Reg<float> raccz_i = &accelerationsz[i];
                    mipp::Reg<float> rvcx_i 
                    mipp::Reg<float> rvcy_i 
                    mipp::Reg<float> rvcz_i 

     
            for (int j = 0; j < n_particles; j += mipp::N<float>()){
                if (i==j) 
                    continue;
                const mipp::Reg<float> diffx=rposx_j-ropsx_i;
                const mipp::Reg<float> diffy=rposx_j-ropsx_i;
                const mipp::Reg<float> diffz=rposx_j-ropsx_i;   
                const mipp::Reg<float> rmass_j = &initstate.masses[j];
                mipp::Reg<float> dij=diffx * diffx + diffy * diffy + diffz * diffz;

                if (dij < 1.0)
                    {
                        dij = 10.0;
                    }
                    else
                    {
                        
                        dij = 10.0 / (dij * mipp::sqrt(dij));
                    }
                   // r1.store(&myVector[(i+1)*mipp::N<float>()]);

                   raccx_i+=diffx*dij*rmass_j;
                   raccy_i+=diffy*dij*rmass_j;
                   raccz_i+=diffz*dij*rmass_j;

            }

            	rvcx_i += raccx_i * 2.0f;
				rvcy_i += raccy_i * 2.0f;
				rvcz_i += raccz_i * 2.0f;
				rposx_i += rvcx_i * 0.1f;
				rposy_i += rvcy_i * 0.1f;
				rposz_i += rvcz_i  * 0.1f;
            


        }



}*/

}
#endif // GALAX_MODEL_CPU_FAST
