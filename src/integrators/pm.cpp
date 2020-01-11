#include "integrators/pm.h"
#include "parallel.h"
#include "scene.h"
#include "imageio.h"
#include "spectrum.h"
#include "rng.h"
#include "paramset.h"
#include "progressreporter.h"
#include "interaction.h"
#include "sampling.h"
#include "samplers/halton.h"
#include "stats.h"

namespace pbrt {

	void PMIntegrator::Render(const Scene &scene){
		
	}

	Integrator *CreatePMIntegrator(const ParamSet &params,
                                 std::shared_ptr<const Camera> camera){
		int nPhotons = params.FindOneInt("nPhotons", 100000);
        int maxDepth = params.FindOneInt("maxdepth", 5);
	    int K = params.FindOneInt("K", 100);
	    return new PMIntegrator(camera, nPhotons, K, maxDepth);
	}
}