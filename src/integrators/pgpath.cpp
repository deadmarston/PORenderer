/* This file is to implement the "Practical Path Guiding for Efficient Light-Transport Simulation", Thomas Muller, et al., 2017.
 * url: https://tom94.net/
 */


#include "integrators/pgpath.h"
#include "bssrdf.h"
#include "camera.h"
#include "film.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"

//12.7, implement the original path tracer
namespace namespace pbrt {
	PathGuidingIntegrator::PathGuidingIntegrator(int maxDepth, std::shared_ptr<const Camera> camera, 
  							  std::shared_ptr<Sampler> sampler,
  							  const Bounds2i &pixelBounds, Float rrThreshold = 1,
  							  const std::string &lightSampleStrategy = "spatial",
  							  const Float quadThreshold = 0.01f, const Float c = 12000)
		: SamplerIntegrator(camera, sampler, pixelBounds),
		  maxDepth(maxDepth),
		  rrThreshold(rrThreshold),
		  lightSampleStrategy(lightSampleStrategy),
		  quadThreshold(quadThreshold),
		  c(c)
		   {}

  	void Preprocess(const Scene &scene, Sampler &sampler)
  	{
  		lightDistribution = CreateLikghtSampleDistribution(lightSampleStrategy, scene);
  	}

  	Spectrum Li(const RayDifferential &ray, const Scene &scene, 
  					Sampler &sampler, MemoryArena &arena, int depth) const
  	{
  		//variable during the path tracing
  		bool specularBounce = false;
  		Spectrum L(0.0f), beta(1.0f);
  		int bounces;

  		//tracing loop
  		for (bounces = 0; ; bounces++)
  		{
	  		//intersection
	  		SurfaceInteraction isect;
	  		bool foundIntersection = scene.Intesect(ray, &isect);

	  		if (bounces == 0 || specularBounce){
	  			//add emitted light at path vertex
	  			if (foundIntersection){
	  				L += beta * isec.Le(-ray.d);
	  			}else{
	  				for (const auto& light : scene.infiniteLights){
	  					L += beta * light->Le(ray);
	  				}
	  			}
	  		}

	  		if (!foundIntersection || bounces > maxDepth){
	  			//not found intersection or exceed the maximum depth, terminate it
	  			break;
	  		}

	  		isect.ComputeScatteringFunctions(ray, arena, true);
	  		if (!isect.bsdf){
				//ignore surface like participating media
				//simply skip over such surface
	  			ray = isect.Spawn(ray.d);
	  			bounces--;
	  			continue;
	  		}

	  		//direct lighting computation
	  		L += beta * UniformSampleOneLight(isect, scene, arena, sampler);

	  		//indirect lighting computation
	  		//generate the next ray to trace
	  		Vector3f wo = -ray.d, wi;
	  		Float pdf;
	  		BxDFType flags;
	  		Spectrum f = isect.bsdf->(wo, &wi, sampler.Get2D(), &pdf, BSDF_ALL, &flags);

	  		if (f.IsBlack() || pdf == 0.f)
	  			break;

	  		//spawn the next ray and compute the beta value
	  		beta *= f * AbsDot(wi, isect.shading.n) / pdf;
	  		specularBounce = (flags & BSDF_SPECULAR) != 0;
	  		ray = isect.SpawnRay(wi);

	  		//terminate the tracing with russian roulette
	  		if (bounces > 3){
	  			FLoat q = std::max((Float).05, 1 - beta.y());
	  			if (sampler.Get1D() < q))
	  				break;
	  			beta /= 1-q;//update the weight with the russian roulette failure weight 1-q
	  		}
  		}
  		return L;
  	}

	PathGuidingIntegrator *CreatePGPathIntegrator(const ParamSet &params,
												  std::shared_ptr<Sampler> sampler,
												  std::shared_ptr<const Camera> camera)
	{
		int maxDepth = params.FindOneInt("maxdepth", 5);
		int np;
		const int *pb = params.FindInt("pixelBounds", &np);
		Bounds2i pixelBounds = camera->film->GetSampleBounds();

		if (pb){
			if (np != 4)
				Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
			else{
				pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
	            if (pixelBounds.Area() == 0)
	                Error("Degenerate \"pixelbounds\" specified.");
			}
		}
		Float rrThreshold = params.FindOneFloat("rrthreshold", 1.);
		std::string lightStrategy = params.FindOneString("lightsamplestrategy", "spatial");
		return new PathGuidingIntegrator(maxDepth, camera, sampler, pixelBounds,
										 rrThreshold, lightStrategy);
	}

}