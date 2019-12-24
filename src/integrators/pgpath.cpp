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
	Points2i DNode::dirToCanonical(Vector3f dir){
		Points2i can;
		float cisTheta = max(min(d.z, 1), -1);
		float phi = atan(dir.y, dir.x);
		can.x = (dir.x+1)/2.f;
		can.y = phi/(2*M_PI);
		return can;
	}

	Vector3f DNode::canonicalToDir(Points2i can){
		float cosTheta = can.x*2-1;
		float sinTheta = sqrt(1 - cosTheta*cosTheta);
		float phi = can.y*2*M_PI;
		return Vector3f(sinTheta*cos(phi), sinTheta.sin(phi), cosTheta);
	}

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

  		//first, we preprocess the scene
  		this->preprocess(scene, *sampler);

  		//compute nTiles, for parallel programming
  		const int tileSize = 16;
  		Bounds2i sampleBounds = camera->film->GetSampleBounds();
  		Vector2i sampleExtent = sampleBounds.Diagonal();

  		Points2i nTile((sampleExtent.x+tileSize-1)/tileSize, (sampleExtent.y+tileSize-1)/tileSize); 
  		
  		ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
    	{
    		ParallelFor2D([&](Point2i tile) {//solve it parallelly

    			//allocate memory arena for tile
    			MemoryArena arena;

    			//seed
    			int seed = tile.y * nTile.x + tile.x;
    			std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

    			//compute bounds
    			int x0 = sampleBounds.pMin.x + tile.x * tileSize;
    			int y0 = sampleBounds.pMin.y + tile.y * tileSize;
    			int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
    			int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
    			Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
	            LOG(INFO) << "Starting image tile " << tileBounds;

	            // Get _FilmTile_ for tile
	            std::unique_ptr<FilmTile> filmTile =
	                camera->film->GetFilmTile(tileBounds);

    			//lookover pixels to render
    			for (Point2i pixel : tileBounds)
    			{
    				{
	    				ProfilePhase pp(Prof::StartPixel);
	                    tileSampler->StartPixel(pixel);
                	}

                	// Do this check after the StartPixel() call; this keeps
	                // the usage of RNG values from (most) Samplers that use
	                // RNGs consistent, which improves reproducability /
	                // debugging.
	                if (!InsideExclusive(pixel, pixelBounds))
	                    continue;

	                //rendering process
	                do{
	                	//camera sampler
	                	CameraSample cameraSample = tileSampler->GetCameraSample(pixel);

	                	RayDifferential ray;
	                	Float rayWeight = camera->GenerateRayDifferential(cameraSample, &ray);
	                	//todo: figure out why scale
	                	ray.ScaleDifferentials(1 / std::sqrt((Float)tileSampler->samplesPerPixel));
                    	++nCameraRays;

                    	Spectrum L(0.f);
                    	if (rayWeight > 0) L = Li(ray, scene, *tileSampler, arena);

                    	if (L.HasNaNs()){
                    		LOG(ERROR) << StringPrintf(
	                            "Not-a-number radiance value returned "
	                            "for pixel (%d, %d), sample %d. Setting to black.",
	                            pixel.x, pixel.y,
	                            (int)tileSampler->CurrentSampleNumber());
                    		L = Spectrum(0.f);
                    	}else if (L.y() < -1e-5){
                    		LOG(ERROR) << StringPrintf(
	                            "Negative luminance value, %f, returned "
	                            "for pixel (%d, %d), sample %d. Setting to black.",
	                            L.y(), pixel.x, pixel.y,
	                            (int)tileSampler->CurrentSampleNumber());
	                        L = Spectrum(0.f);
                    	} else if (std::isinf(L.y())){
                    		LOG(ERROR) << StringPrintf(
	                            "Infinite luminance value returned "
	                            "for pixel (%d, %d), sample %d. Setting to black.",
	                            pixel.x, pixel.y,
	                            (int)tileSampler->CurrentSampleNumber());
	                        L = Spectrum(0.f);
                    	}
                    	VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
                        	ray << " -> L = " << L;

                        //add ray's contribution to image
                        filmTile->AddSample(cameraSample.pFilm, L, rayWeight);

                        arena.Reset();
	                }while(tileSampler->StartNextSample())
    			}
    			LOG(INFO) << "Finished image tile " << tileBounds;

				camera->file->MergeFilmTile(std::move(filmTile));
				reporter.Update();    		
    		}, nTiles);
    		reporter.Done();
    	}
    	LOG(INFO) << "Rendering finished";

	    // Save final image after rendering
	    camera->film->WriteImage();
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