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
#include <iostream>
using namespace std;
//12.7, implement the original path tracer
namespace pbrt {
	PathGuidingIntegrator::PathGuidingIntegrator(int maxDepth, std::shared_ptr<const Camera> camera, 
  							  std::shared_ptr<Sampler> sampler,
  							  const Bounds2i &pixelBounds, Float rrThreshold,
  							  const std::string &lightSampleStrategy,
  							  const Float quadThreshold, const Float c)
		: SamplerIntegrator(camera, sampler, pixelBounds),
		  maxDepth(maxDepth),
		  rrThreshold(rrThreshold),
		  lightSampleStrategy(lightSampleStrategy),
		  quadThreshold(quadThreshold),
		  c(c)
		   {}

  	void PathGuidingIntegrator::Preprocess(const Scene &scene, Sampler &sampler)
  	{
  		lightDistribution = CreateLightSampleDistribution(lightSampleStrategy, scene);
  	}

  	void PathGuidingIntegrator::Render(const Scene& scene)
  	{
  		//todo: divide the whole rendering process into several small rendering process
  		//generate one SD-Tree for guiding, and one SD-Tree for storing the light field
  		//consider the parallel programming

  		return;
  	}

  	Spectrum PathGuidingIntegrator::Li(const RayDifferential &r, const Scene &scene, 
  					Sampler &sampler, MemoryArena &arena, int depth) const
  	{
  		//variable during the path tracing
  		bool specularBounce = false;
  		Spectrum L(0.0f), beta(1.0f);
  		int bounces;
  		RayDifferential ray(r);
  		Float etaScale = 1;
  		//tracing loop
  		for (bounces = 0; ; bounces++)
  		{
	  		//intersection
	  		SurfaceInteraction isect;
	  		bool foundIntersection = scene.Intersect(ray, &isect);

	  		if (bounces == 0 || specularBounce){
	  			//add emitted light at path vertex
	  			if (foundIntersection){
	  				L += beta * isect.Le(-ray.d);
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
	  			ray = isect.SpawnRay(ray.d);
	  			bounces--;
	  			continue;
	  		}

	  		//direct lighting computation
	  		//L += beta * UniformSampleOneLight(isect, scene, arena, sampler);

	  		//direct lighting computation with light distribution
	  		const Distribution1D *distrib = lightDistribution->Lookup(isect.p);

	  		//sample illumination from lights to find path contribution
	  		//skip this for perfectly specular BSDFS
	  		if (isect.bsdf->NumComponents(BxDFType(BSDF_ALL & ~BSDF_SPECULAR)) > 0){

	  			Spectrum Ld = beta * UniformSampleOneLight(isect, scene, arena, sampler, false, distrib);
	  			L += Ld;
	  		}

	  		//indirect lighting computation
	  		//generate the next ray to trace
	  		Vector3f wo = -ray.d, wi;
	  		Float pdf;
	  		BxDFType flags;
	  		Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, BSDF_ALL, &flags);

	  		if (f.IsBlack() || pdf == 0.f)
	  			break;

	  		//spawn the next ray and compute the beta value
	  		beta *= f * AbsDot(wi, isect.shading.n) / pdf;
	  		specularBounce = (flags & BSDF_SPECULAR) != 0;
	  		
	  		//account for eta through specular reflection and refraction
	  		if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)){
	  			Float eta = isect.bsdf->eta;

	  			etaScale *= (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
	  		}

	  		ray = isect.SpawnRay(wi);

	  		//todo: deal with the case of bssrdf

	  		Spectrum rrBeta = beta * etaScale;
	  		//terminate the tracing with russian roulette
	  		if (rrBeta.MaxComponentValue() < rrThreshold && bounces > 3){
	  			Float q = std::max((Float).05, 1 - rrBeta.MaxComponentValue());
	  			if (sampler.Get1D() < q)
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
		std::cout << "pgpath integrator created\n";
		return new PathGuidingIntegrator(maxDepth, camera, sampler, pixelBounds,
										 rrThreshold, lightStrategy);
	}

}