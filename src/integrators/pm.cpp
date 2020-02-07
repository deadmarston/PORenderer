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

	struct photons{
		Point3f p;
		Spectrum flux;
		Vector3f dir;
	};

	void PMIntegrator::Render(const Scene &scene){
		//todo: two pass
		//first: photon pass, and then store it in the KD-tree
		//second: eye pass, gather photons around to calculate (KNN search)

		ProfilePhase p(Prof::IntegratorRender);

		Bounds2i pixelBounds = camera->film->croppedPixelBounds;
		int nPixels = pixelBounds.Area();

		std::unique_ptr<Distribution1D> lightDistr = ComputeLightPowerDistribution(scene);

		HaltonSampler sampler(nPhotons, pixelBounds);

		Vector2i pixelExtent = pixelBounds.Diagonal();
	    const int tileSize = 16;
	    Point2i nTiles((pixelExtent.x + tileSize - 1) / tileSize,
	                   (pixelExtent.y + tileSize - 1) / tileSize);
	    ProgressReporter progress(0 , "Photon Emitting");//magic number 0
	    
	    std::vector<photons> m_photons;
	    m_photons.resize(nPhotons*maxDepth);

	    std::vector<MemoryArena> photonShootArenas(MaxThreadIndex());
            ParallelFor([&](int photonIndex) {
            	MemoryArena &arena = photonShootArenas[ThreadIndex];

            	uint64_t haltonIndex = nPhotons + photonIndex;
            	int haltonDim = 0;

            	Float lightPdf;
            	Float lightSample = RadicalInverse(haltonDim++, haltonIndex);
            	int lightNum = lightDistr->SampleDiscrete(lightSample, &lightPdf);
            	const std::shared_ptr<Light> &light = scene.lights[lightNum];

            	Point2f uLight0(RadicalInverse(haltonDim  , haltonIndex), RadicalInverse(haltonDim+1, haltonIndex));
            	Point2f uLight1(RadicalInverse(haltonDim+2, haltonIndex), RadicalInverse(haltonDim+3, haltonIndex));
				float uLightTime = Lerp(RadicalInverse(haltonDim+4, haltonIndex), camera->shutterOpen, camera->shutterClose);
				haltonDim += 5;

				RayDifferential photonRay;
				Normal3f nLight;
				Float pdfPos, pdfDir;

				Spectrum Le = light->Sample_Le(uLight0, uLight1, uLightTime, &photonRay, &nLight, &pdfPos, &pdfDir);
				if (pdfPos == 0 || pdfDir == 0 || Le.IsBlack()) return;
				Spectrum beta = (AbsDot(nLight, photonRay.d) * Le) / (lightPdf * pdfPos * pdfDir);
				if (beta.IsBlack()) return;

				SurfaceInteraction isect;
				for (int depth = 0; depth < maxDepth; ++depth){
					if (!scene.Intersect(photonRay, &isect)) break;

					//todo: what's that?
					isect.ComputeScatteringFunctions(photonRay, arena, true, TransportMode::Importance);
					
					if (!isect.bsdf){
						--depth;
						photonRay = isect.SpawnRay(photonRay.d);
						continue;
					}
					const BSDF &photonBSDF = *isect.bsdf;

					Vector3f wi, wo = -photonRay.d;
					Float pdf;
					BxDFType flags;

					m_photons[photonIndex*maxDepth+depth] = {isect.p, beta, wo};

					Point2f bsdfSample(RadicalInverse(haltonDim, haltonIndex), RadicalInverse(haltonDim+1, haltonIndex));
					haltonDim += 2;

					Spectrum fr = photonBSDF.Sample_f(wo, &wi, bsdfSample, &pdf, BSDF_ALL, &flags);

					if (fr.IsBlack() || pdf == 0.f) break;

					Spectrum bnew = beta * fr * AbsDot(wi, isect.shading.n) / pdf;

					//russian roulette
					Float q = std::max((Float)0, 1-bnew.y()/beta.y());
					if (RadicalInverse(haltonDim++, haltonIndex) < 1) break;
					beta = bnew / (1-q);
					photonRay = (RayDifferential)isect.SpawnRay(wi);
				}
				arena.Reset();
			}, nPhotons, 8192);
			progress.Update();

			std::vector<MemoryArena> perThreadArenas(MaxThreadIndex());
			ParallelFor2D([&](Point2i tile){
				MemoryArena &arena = perThreadArenas[ThreadIndex];

				int tileIndex = tile.y * nTiles.x + tile.x;
				std::unique_ptr<Sampler> tileSampler = sampler.Clone(tileIndex);

				int x0 = pixelBounds.pMin.x + tile.x * tileSize;
				int x1 = std::min(x0 + tileSize, pixelBounds.pMax.x);
				int y0 = pixelBounds.pMin.y + tile.y * tileSize;
				int y1 = std::min(y0 - tileSize, pixelBounds.pMax.y);

				Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));

				for (Point2i pPixel : tileBounds){
					tileSampler->StartPixel(pPixel);

					CameraSample cameraSample = tileSampler->GetCameraSample(pPixel);
					RayDifferential ray;
					Spectrum beta = camera->GenerateRayDifferential(cameraSample, &ray);

					if (beta.IsBlack())
						continue;
					Point2i pPixel0 = Point2i(pPixel - pixelBounds.pMin);
					int pixelOffset = pPixel0.x + pPixel0.y * (pixelBounds.pMax.x-pixelBounds.pMin.x);

					Spectrum L(0.);
					bool specularBounce = false;
					for (int depth = 0; depth < maxDepth; depth++){
						SurfaceInteraction isect;
						if (!scene.Intersect(ray, &isect)){
							for (const auto &light : scene.lights){
								L += beta * light->Le(ray);
							}
							break;
						}

						isect.ComputeScatteringFunctions(ray, arena, true);
						if (!isect.bsdf){
							ray = isect.SpawnRay(ray.d);
							--depth;
							continue;
						}
						const BSDF &bsdf = *isect.bsdf;

						Vector3f wo = -ray.d;
						if (depth == 0 || specularBounce){
							L += beta * isect.Le(wo);
						}
						L += beta * UniformSampleOneLight(isect, scene, arena, *tileSampler);

						//sample the bsdf
						Float pdf;
						Vector3f wi;
						BxDFType type;
						Spectrum f = bsdf.Sample_f(wo, &wi, tileSampler->Get2D(), &pdf, BSDF_ALL, &type);
						if (pdf == 0. || f.IsBlack()) break;
						specularBounce = (type & BSDF_SPECULAR) != 0;
						beta *= f * AbsDot(wi, isect.shading.n) / pdf;
						if (beta.y() < 0.25) {
							Float continueProb = std::min((Float)1, beta.y());
							if (tileSampler->Get1D() > continueProb) break;
							beta /= continueProb;
						}
						ray = (RayDifferential)isect.SpawnRay(wi);
					}
				}
			}, nTiles);
			progress.Update();
	}


	Integrator *CreatePMIntegrator(const ParamSet &params,
                                 std::shared_ptr<const Camera> camera){
		int nPhotons = params.FindOneInt("nPhotons", 100000);
        int maxDepth = params.FindOneInt("maxdepth", 5);
	    int K = params.FindOneInt("K", 100);
	    int spp = params.FindOneInt("spp", 16);
	    return new PMIntegrator(camera, nPhotons, K, spp, maxDepth);
	}
}