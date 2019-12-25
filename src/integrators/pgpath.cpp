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
#include "progressreporter.h"
#include <iostream>
#include <vector>
using namespace std;
//12.7 , implement the original path tracer
//12.23, add comments for mis with nee, preparing to finish it asap
namespace pbrt {
	STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);
	void testDNode(){
		std::vector<DNode> nodes;
		DNode dnode;
		nodes.push_back(dnode);
		assert(dnode.sum(0) == 0);
		assert(dnode.sum(1) == 0);
		assert(dnode.sum(2) == 0);
		assert(dnode.sum(3) == 0);
		assert(dnode.isLeaf(0));
		assert(dnode.isLeaf(1));
		assert(dnode.isLeaf(2));
		assert(dnode.isLeaf(3));
		cout << "dnode init pass" << endl;
		//create a child
		DNode child;
		nodes.push_back(child);
		dnode.setChild(0, 1);
		assert(!dnode.isLeaf(0));
		assert(child.isLeaf(0));
		cout << "dnode setchild pass" << endl;
		assert(dnode.childIndex(Point2i(0.25, 0.25)) == 0);
		assert(dnode.childIndex(Point2i(0.75, 0.75)) == 3);
		assert(dnode.childIndex(Point2i(0.25, 0.75)) == 1);
		assert(dnode.childIndex(Point2i(0.75, 0.25)) == 2);
		assert(dnode.childIndex(Point2i(0.0, 0.0)) == 0);
		assert(dnode.childIndex(Point2i(1.0, 1.0)) == 3);
		assert(dnode.childIndex(Point2i(0.5, 0.5)) == 3);
		assert(dnode.childIndex(Point2i(0.0, 0.5)) == 2);
		assert(dnode.childIndex(Point2i(0.5, 0.0)) == 1);
		cout << "dnode childindex pass" << endl;
		dnode.setSum(1, 0.5);
		assert(dnode.sum(1) == 0.5);
		cout << "dnode setsum pass" << endl;
		assert(dnode.depthAt(0.8, 0.8) == 1);
		assert(dnode.depthAt(0.2, 0.2) == 2);
		cout << "dnode depthat pass" << endl;
		Point2i p = Point2i(0.1, 0.1);
		dnode.record(p, 0.3, nodes);
		assert(child.sum(0) == 0.3);
		p = Point2i(0.4, 0.4);
		dnode.record(p, 0.5, nodes);
		assert(child.sum(3) == 0.5);
		cout << "dnode record pass" << endl;
	}

	void testSDTree(){
		cout << "===============test begin===============" << endl;
		testDNode();
		cout << "================test end================" << endl;
	}

	static void addToAtomicFloat(std::atomic<Float>& var, Float val){
		auto current = var.load();
		while(!var.compare_exchange_weak(current, current + val));
	}

	DNode::DNode(){
		for (size_t i = 0; i < m_sum.size(); i++){
			m_sum[i].store(0, std::memory_order_relaxed);
			m_nodes[i] = LEAFINDEX;
		}
	}

	DNode::DNode(const DNode& node){
		for (int i = 0; i < node.m_sum.size(); i++){
			setSum(i, node.sum(i));
			m_nodes[i] = node.child(i);
		}
	}

	void DNode::setSum(int index, Float val){
		m_sum[index].store(val, std::memory_order_relaxed);
	}

	Float DNode::sum(int index) const{
		return m_sum[index].load(std::memory_order_relaxed);
	}

	void DNode::setChild(int index, uint16_t val){
		m_nodes[index] = val;
	}

	uint16_t DNode::child(int index) const{
		return m_nodes[index];
	}

	void DNode::record(Point2i& can, float irradiance, std::vector<DNode>& nodes){
		int index = childIndex(can);
		//didn't add to the sum of current node for efficiency
		if (isLeaf(index)){
			addToAtomicFloat(m_sum[index], irradiance);//make it atomic due to parallel programming
		}else{
			nodes[child(index)].record(can, irradiance, nodes);
		}
	}

	int DNode::depthAt(Point2i& can, const std::vector<DNode>& nodes) const{
		int index = childIndex(can);

		if (isLeaf(index)){
			return 1;
		}else{
			return 1 + nodes[child(index)].depthAt(can, nodes);
		}
	}

	Float DNode::pdf(Point2i& can, const std::vector<DNode>& nodes) const{
		int index = childIndex(can);

		if (sum(index) <= 0){
			return 0;
		}

		Float factor = 4*sum(index) / (sum(0)+sum(1)+sum(2)+sum(3));
		if (isLeaf(index)){
			return factor;
		}else{
			return factor*pdf(can, nodes);
		}
	}

	Point2f DNode::sample(Sampler* sampler, const std::vector<DNode>& nodes) const{
		Float bottomLeft = sum(0);
		Float bottom = bottomLeft + sum(1);
		Float bottom_upperLeft = bottom + sum(2);
		Float total = bottom_upperLeft + sum(3);

		Float sample = sampler->Get1D();
		int index = 0;

		Point2f origin;
		if (sample > bottom/total){
			if (sample > bottomLeft/total){
				index = 1;
				origin.x = 0.5;
			}
		}else{
			if (sample > bottom_upperLeft/total){
				index = 3;
				origin.x = 0.5;
			}else{
				index = 2;
			}
			origin.y = 0.5;
		}
		if (isLeaf(index)){
			return origin+0.5*sampler->Get2D();
		}else{
			return origin+0.5*nodes[child(index)].sample(sampler, nodes);
		}
	}

	uint16_t DNode::childIndex(Point2i& can) const{
		/*
		-------
		|2 |3 |
		-------
		|0 |1 |
		-------
		*/
		assert(can.x >= 0 && can.y >= 0 && can.x <= 1 && can.y <= 1);
		uint16_t index = 0;
		for (int i = 0;  i < 2; i++){
			if (can[i] > 0.5){
				can[i] = (can[i]-0.5)*2;
				index |= (1 << i);
			}else{
				can[i] *= 2;
			}
		}
		return index;
	}

	Point2i DTree::dirToCanonical(Vector3f dir){
		Point2i can;
		float cosTheta = max(min(dir.z, 1.f), -1.f);
		float phi = atan2(dir.y, dir.x);
		can.x = (cosTheta+1)/2.f;
		can.y = phi/(2*M_PI);
		return can;
	}

	Vector3f DTree::canonicalToDir(Point2i can){
		float cosTheta = can.x*2-1;
		float sinTheta = sqrt(1 - cosTheta*cosTheta);
		float phi = can.y*2*M_PI;
		return Vector3f(sinTheta*cos(phi), sinTheta*sin(phi), cosTheta);
	}
	
	Spectrum UniformSampleOneLight_PG(const Interaction &it, const Scene &scene,
					       MemoryArena &arena, Sampler &sampler,
					       bool handleMedia, const Distribution1D *lightDistrib) {
	    ProfilePhase p(Prof::DirectLighting);
	    // Randomly choose a single light to sample, _light_
	    int nLights = int(scene.lights.size());
	    if (nLights == 0) return Spectrum(0.f);
	    int lightNum;
	    Float lightPdf;
	    if (lightDistrib) {
		lightNum = lightDistrib->SampleDiscrete(sampler.Get1D(), &lightPdf);
		if (lightPdf == 0) return Spectrum(0.f);
	    } else {
		lightNum = std::min((int)(sampler.Get1D() * nLights), nLights - 1);
		lightPdf = Float(1) / nLights;
	    }
	    const std::shared_ptr<Light> &light = scene.lights[lightNum];
	    Point2f uLight = sampler.Get2D();
	    Point2f uScattering = sampler.Get2D();
	    return EstimateDirect(it, uScattering, *light, uLight,
				  scene, sampler, arena, handleMedia) / lightPdf;
	}

	Spectrum EstimateDirect_PG(const Interaction &it, const Point2f &uScattering,
				const Light &light, const Point2f &uLight,
				const Scene &scene, Sampler &sampler,
				MemoryArena &arena, bool handleMedia, bool specular) {
	    BxDFType bsdfFlags =
		specular ? BSDF_ALL : BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
	    Spectrum Ld(0.f);
	    // Sample light source with multiple importance sampling
	    Vector3f wi;
	    Float lightPdf = 0, scatteringPdf = 0;
	    VisibilityTester visibility;
	    Spectrum Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &visibility);
	    VLOG(2) << "EstimateDirect uLight:" << uLight << " -> Li: " << Li << ", wi: "
		    << wi << ", pdf: " << lightPdf;
	    if (lightPdf > 0 && !Li.IsBlack()) {
		// Compute BSDF or phase function's value for light sample
		Spectrum f;
		if (it.IsSurfaceInteraction()) {
		    // Evaluate BSDF for light sampling strategy
		    const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
		    f = isect.bsdf->f(isect.wo, wi, bsdfFlags) *
			AbsDot(wi, isect.shading.n);
		    scatteringPdf = isect.bsdf->Pdf(isect.wo, wi, bsdfFlags);
		    VLOG(2) << "  surf f*dot :" << f << ", scatteringPdf: " << scatteringPdf;
		} else {//todo: skip phase function for now
		    // Evaluate phase function for light sampling strategy
		    const MediumInteraction &mi = (const MediumInteraction &)it;
		    Float p = mi.phase->p(mi.wo, wi);
		    f = Spectrum(p);
		    scatteringPdf = p;
		    VLOG(2) << "  medium p: " << p;
		}
		if (!f.IsBlack()) {
		    // Compute effect of visibility for light source sample
		    if (handleMedia) {
			Li *= visibility.Tr(scene, sampler);
			VLOG(2) << "  after Tr, Li: " << Li;
		    } else {
		      if (!visibility.Unoccluded(scene)) {
			VLOG(2) << "  shadow ray blocked";
			Li = Spectrum(0.f);
		      } else
			VLOG(2) << "  shadow ray unoccluded";
		    }

		    // Add light's contribution to reflected radiance
		    if (!Li.IsBlack()) {
			if (IsDeltaLight(light.flags))
			    Ld += f * Li / lightPdf;
			else {
						//todo: consider the sd-tree, and utilize mixed pdf to compute the power heuristic
						//ex: mixed_pdf, d_pdf, b_pdf = pdf(direction)
						//    weight = powerheuristic(mixed_pdf, lightPdf)
						//    compute the Ld
			    Float weight =
				PowerHeuristic(1, lightPdf, 1, scatteringPdf);
			    Ld += f * Li * weight / lightPdf;
			}
		    }
		}
	    }

	    // Sample BSDF with multiple importance sampling
	    if (!IsDeltaLight(light.flags)) {
		Spectrum f;
		bool sampledSpecular = false;
		if (it.IsSurfaceInteraction()) {
		    // Sample scattered direction for surface interactions
		    BxDFType sampledType;
		    const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
				//todo: sample the bsdf and sd-tree with a fraction
				//ex:   direction, bsdf = sampleBSDF&SDTree()
				//      mixed_pdf, d_pdf, b_pdf = pdf(direction)
		    f = isect.bsdf->Sample_f(isect.wo, &wi, uScattering, &scatteringPdf,
					     bsdfFlags, &sampledType);
		    f *= AbsDot(wi, isect.shading.n);
		    sampledSpecular = (sampledType & BSDF_SPECULAR) != 0;
		} else {//todo: skip the medium interactions for now
		    // Sample scattered direction for medium interactions
		    const MediumInteraction &mi = (const MediumInteraction &)it;
		    Float p = mi.phase->Sample_p(mi.wo, &wi, uScattering);
		    f = Spectrum(p);
		    scatteringPdf = p;
		}
		VLOG(2) << "  BSDF / phase sampling f: " << f << ", scatteringPdf: " <<
		    scatteringPdf;
		if (!f.IsBlack() && scatteringPdf > 0) {
		    // Account for light contributions along sampled direction _wi_
		    Float weight = 1;
		    if (!sampledSpecular) {
			lightPdf = light.Pdf_Li(it, wi);
			if (lightPdf == 0) return Ld;
			weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
		    }

		    // Find intersection and compute transmittance
		    SurfaceInteraction lightIsect;
		    Ray ray = it.SpawnRay(wi);
		    Spectrum Tr(1.f);
		    bool foundSurfaceInteraction =
			handleMedia ? scene.IntersectTr(ray, sampler, &lightIsect, &Tr)
				    : scene.Intersect(ray, &lightIsect);

		    // Add light contribution from material sampling
		    Spectrum Li(0.f);
		    if (foundSurfaceInteraction) {
			if (lightIsect.primitive->GetAreaLight() == &light)
			    Li = lightIsect.Le(-wi);
		    } else
			Li = light.Le(ray);
		    if (!Li.IsBlack()) Ld += f * Li * Tr * weight / scatteringPdf;
		}
	    }
	    return Ld;
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
  		//unit test code here
  		// testSDTree();
  		// return;
  		//todo: divide the whole rendering process into several small rendering process
  		//generate one SD-Tree for guiding, and one SD-Tree for storing the light field
  		//consider the parallel programming

  		//first, we preprocess the scene
  		Preprocess(scene, *sampler);

  		//compute nTiles, for parallel programming
  		const int tileSize = 16;
  		Bounds2i sampleBounds = camera->film->GetSampleBounds();
  		Vector2i sampleExtent = sampleBounds.Diagonal();

  		Point2i nTile((sampleExtent.x+tileSize-1)/tileSize, (sampleExtent.y+tileSize-1)/tileSize); 
  		
  		ProgressReporter reporter(nTile.x * nTile.y, "Rendering");
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
                    	//++nCameraRays;

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
	                }while(tileSampler->StartNextSample());
    			};
    			LOG(INFO) << "Finished image tile " << tileBounds;

				camera->film->MergeFilmTile(std::move(filmTile));
				reporter.Update();    		
    		}, nTile);
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
