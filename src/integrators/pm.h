#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PM_H
#define PBRT_INTEGRATORS_PM_H

// integrators/sppm.h*
#include "pbrt.h"
#include "integrator.h"
#include "camera.h"
#include "film.h"

namespace pbrt {
	// SPPM Declarations
class PMIntegrator : public Integrator {
  public:
    // SPPMIntegrator Public Methods
    PMIntegrator(std::shared_ptr<const Camera> &camera, int nPhotons, int K, int spp, int maxDepth)
        : camera(camera),
          maxDepth(maxDepth),
          nPhotons(nPhotons),
          K(K),
          spp(spp)
          {};
    void Render(const Scene &scene);

  private:
    // SPPMIntegrator Private Data
    std::shared_ptr<const Camera> camera;
    const int maxDepth;
    const int nPhotons;
    const int K; //KNN search
    const int spp; //spp for aa, motion blur, dof
};

Integrator *CreatePMIntegrator(const ParamSet &params,
                                 std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif