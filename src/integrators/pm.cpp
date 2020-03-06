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
#include "samplers/random.h"
#include "stats.h"
#include "omp.h"
#include <queue>
#include <mutex>
#include <chrono>
#include <algorithm>
namespace pbrt {

	//todo: implement a balanced kd tree
	Float mean_r = 0.0f;

	void testpq(){
		std::priority_queue<int> pq;
		pq.push(1);
		pq.push(4);
		pq.push(5);
		pq.push(2);
		pq.push(9);
		pq.push(3);
		while(!pq.empty()){
			std::cout << pq.top() << std::endl;
			pq.pop();
		}
	}

	/*int findMedianCounting(std::vector<photon>& photons, int axis, int l, int r){

	}*/

	void exchange(std::vector<photon>& photons, int i1, int i2){
		photon tmp = photons[i1];
		photons[i1] = photons[i2];
		photons[i2] = tmp;
	}

	void kdTree::build(std::vector<photon>& photons){
		m_nodes.resize(photons.size());
		Point3f _mean = Point3f(0.f, 0.f, 0.f);
		Point3f _var = Point3f(0.f, 0.f, 0.f);
		int n = photons.size();
		for (int i = 1; i < photons.size(); i++){
			for (int j = 0; j < 3; j++){
				_mean[j] += photons[i].p[j];
				_var[j] += photons[i].p[j]*photons[i].p[j];
			}
		}
		for (int j = 0; j < 3; j++){
			_mean[j] /= (Float)n;
			_var[j] /= (Float)n;
			_var[j] -= _mean[j]*_mean[j];
		}
		if (_var[0] >= _var[1] && _var[0] >= _var[2]){
			m_axis = 0;
		}else if (_var[1] >= _var[2]){
			m_axis = 1;
		}else{
			m_axis = 2;
		}
		recursiveBuild(photons, photons.begin(), photons.end(), m_axis, 0);
		//get maximun axis
		/*Point3f _min = photons[0].p;
		Point3f _max = photons[0].p;
		for (int i = 0; i < photons.size(); i++){
			for (int j = 0; j < 3; j++){
				_min[j] = std::min(_min[j], photons[i].p[j]);
				_max[j] = std::max(_max[j], photons[i].p[j]);
			}
		}
		Float xx = _max[0]-_min[0];
		Float yy = _max[1]-_min[1];
		Float zz = _max[2]-_max[2];
		int max_axis = 0;
		if (xx >= yy && xx >= zz){
			max_axis = 0;
		}else if (yy >= zz){
			max_axis = 1;
		}else{
			max_axis = 2;
		}
		m_axis = max_axis;
		std::cout << "main axis is " << m_axis << std::endl;
		sort(photons.begin(), photons.end(), compare_func[max_axis]);
		int mid = (photons.size()+1)/2;
		recursiveBuild(photons, 0, mid, (max_axis+1)%3, 1);
		recursiveBuild(photons, mid, photons.size(), (max_axis+1)%3, 1);*/
	}	

	int kdTree::recursiveBuild(std::vector<photon>& photons, std::vector<photon>::iterator l, std::vector<photon>::iterator r, int axis, int depth){
		if (l >= r) { m_depth = std::max(m_depth, depth); return -1;}
		int _axis = axis%3;
		std::size_t len = r-l;
		auto mid = l + len/2;
		std::nth_element(l, mid, r, compare_func[_axis]);
		int count = 0;
		while (mid > l && (mid-1)->p[_axis] == mid->p[_axis]){
			--mid;
			count++;
		}
		int index = currentIndex++;
		int ret1 = recursiveBuild(photons, l, mid, _axis+1, depth+1);
		int ret2 = recursiveBuild(photons, mid+1, r, _axis+1, depth+1);
		//if (ret1 == -1 && ret2 == -1) m_nodes[index].isLeaf = true;
		m_nodes[index].left_idx = ret1;
		m_nodes[index].right_idx = ret2;
		m_nodes[index].idx = mid-photons.begin();
		return index;
	}

	void kdTree::dump(const std::vector<photon>& photons, int index, int layer){
		if (index == 0) {
			std::cout << "first axis is " << m_axis << std::endl;
			std::cout << "max depth of kdtree is " << m_depth << std::endl;
		}
		if (index < 0) return;
		std::cout << photons[m_nodes[index].idx].p << " " << layer << std::endl;
		dump(photons, m_nodes[index].left_idx, layer+1);
		dump(photons, m_nodes[index].right_idx, layer+1);
	}

	void kdTree::dump(const std::vector<photon>& photons, int index, int layer, Point3f& p){
		if (index == 0) {
			std::cout << "first axis is " << m_axis << std::endl;
			std::cout << "max depth of kdtree is " << m_depth << std::endl;
		}
		if (index < 0) return;
		std::cout << photons[m_nodes[index].idx].p << " " << layer << " " << DistanceSquared(photons[m_nodes[index].idx].p, p) << std::endl;
		dump(photons, m_nodes[index].left_idx, layer+1, p);
		dump(photons, m_nodes[index].right_idx, layer+1, p);
	}

	void kdTree::ksearch(const Point3f& p, int K, const std::vector<photon>& photons, std::priority_queue<PhotonIdx>& ids){
		Float radius = m_radius;
		ksearch(p, K, 0, photons, ids, m_axis, radius);
	}

	void kdTree::ksearch(const Point3f& p, int K, int index, const std::vector<photon>& photons, std::priority_queue<PhotonIdx>& ids, int axis, Float& radius){
		if (index < 0) return;
		int phtId = m_nodes[index].idx;

		Float pDis = DistanceSquared(p, photons[phtId].p);
	
		Float _axis = axis%3;
		Float signedDis = p[_axis]-photons[phtId].p[_axis];
		Float distance = signedDis*signedDis;

		//std::cout << photons[phtId].p << " " << pDis << " " << distance << " " << signedDis << std::endl;

		if (signedDis > 0){
			ksearch(p, K, m_nodes[index].right_idx, photons, ids, axis+1, radius);
			if (distance < radius){
				ksearch(p, K, m_nodes[index].left_idx, photons, ids, axis+1, radius);
			}
		}else{
			ksearch(p, K, m_nodes[index].left_idx, photons, ids, axis+1, radius);
			if (distance < radius){
				ksearch(p, K, m_nodes[index].right_idx, photons, ids, axis+1, radius);
			}
		}
		if (pDis < radius){
			ids.push({pDis, phtId});
			if (ids.size() > K){
				ids.pop();
				radius = ids.top().dis;
			}
		}

	}

	Float photonMap::linearSearch(const Point3f& p, int K, std::priority_queue<PhotonIdx>& ids){
		Float ret;
		for (int i = 0; i < m_photons.size(); i++){
			Float dis = DistanceSquared(p, m_photons[i].p);
			if (dis < initialRadius){
				ids.push({dis, i});
				if (ids.size() > K){
					ids.pop();
				}
			}
		}
		ret = ids.empty()?initialRadius:ids.top().dis;
		return ret;
	}

	Spectrum photonMap::gatherKPhotons(const SurfaceInteraction* isect, std::priority_queue<PhotonIdx>& ids){
		Spectrum L(0.);
		//std::cout << ids.size() << std::endl;
		while (!ids.empty()){
			Spectrum f = isect->bsdf->f(isect->wo, m_photons[ids.top().id].wi);
			L += f * m_photons[ids.top().id].flux;
			ids.pop();
		}
		return L;
	}

	/*Spectrum PMIntegrator::KNN(const Interaction &it, std::vector<photon> photonList){
		std::priority_queue<PhotonIdx> pq;
		int minSize = std::min(K, int(photonList.size()));
		int maxSize = std::max(K, int(photonList.size()));
		float radiusSquared = RADIUS_THERSHOLD * RADIUS_THERSHOLD;
		for (int i = 0; i < minSize; i++){
			if (photonList[i].isInited = false){
				continue;
			}
			float disSquared = DistanceSquared(it.p, photonList[i].p);
			if (disSquared > radiusSquared){
				continue;
			}
			pq.push({disSquared, i});
		}
		radiusSquared = pq.top().dis; 
		for (int i = minSize; i < maxSize; i++){
			float disSquared = DistanceSquared(it.p, photonList[i].p);
			if (disSquared > radiusSquared){
				continue;
			}
			pq.pop();
			pq.push({disSquared, i});
			radiusSquared = pq.top().dis;
		}
		Spectrum L(0.);
		const SurfaceInteraction* isect = (SurfaceInteraction*)(&it);
		while(!pq.empty()){
			PhotonIdx tmp = pq.top();
			pq.pop();
			Spectrum f = isect->bsdf->f(isect->wo, photonList[tmp.id].wi);
			L += f * photonList[tmp.id].flux;
		}
		if (radiusSquared == 0) return Spectrum(0.);
		L /= (M_PI * radiusSquared * nPhotons);
		return L;
	}*/

	/*Spectrum PMIntegrator::KNN(const Interaction &it, std::vector<photon> photonList){
		std::priority_queue<PhotonIdx> pq;
		int minSize = std::min(K, int(photonList.size()));
		int maxSize = std::max(K, int(photonList.size()));
		for (int i = 0; i < minSize; i++){
			if (photonList[i].isInited = false){
				continue;
			}
			Float disSquared = DistanceSquared(it.p, photonList[i].p);
			pq.push({disSquared, i});
		}

		for (int i = minSize; i < photonList.size(); i++){
			Float disSquared = DistanceSquared(it.p, photonList[i].p);
			pq.push({disSquared, i});
			pq.pop();
		}
		Spectrum L(0.);
		const SurfaceInteraction* isect = (SurfaceInteraction*)(&it);
		Float radius = std::min(pq.top().dis, initialRadius);
		Float radiusSquared = radius * radius;
		mean_r += radius;
		int count = 0;
		for (int i = 0; i < minSize; i++){
			PhotonIdx tmp = pq.top();
			pq.pop();
			if (tmp.dis > radiusSquared) continue;
			Spectrum f = isect->bsdf->f(isect->wo, photonList[tmp.id].wi);
			L += f * photonList[tmp.id].flux;
			count ++;
		}
		Float inv = 1/(M_PI * radiusSquared * nPhotons);
		L *= inv;
		return L;
	}*/

	Spectrum PMIntegrator::KNN(const Interaction &it){
		std::priority_queue<PhotonIdx> ids;
		Float radius = m_globalMap.ksearch(it.p, K, ids);
		//
		const SurfaceInteraction* isect = (SurfaceInteraction*)(&it);
		Spectrum L = m_globalMap.gatherKPhotons(isect, ids);
		Float inv = 1/(M_PI * radius * nPhotons);
		L *= inv;
		return L;
	}

	void PMIntegrator::Render(const Scene &scene){
		//todo: two pass
		//first: photon pass, and then store it in the KD-tree
		//second: eye pass, gather photons around to calculate (KNN search)
		//testpq();
		//return;
		ProfilePhase p(Prof::IntegratorRender);

		Bounds2i pixelBounds = camera->film->croppedPixelBounds;
		int nPixels = pixelBounds.Area();

		std::unique_ptr<Distribution1D> lightDistr = ComputeLightPowerDistribution(scene);

		HaltonSampler sampler(spp, pixelBounds);
		//RandomSampler sampler(spp);

		Vector2i pixelExtent = pixelBounds.Diagonal();
	    const int tileSize = 16;
	    Point2i nTiles((pixelExtent.x + tileSize - 1) / tileSize,
	                   (pixelExtent.y + tileSize - 1) / tileSize);
	    

	    //std::vector<photon> m_photons;
	    //long pages = sysconf(_SC_PHYS_PAGES);
	    //long page_size = sysconf(_SC_PAGE_SIZE);
	    //std::cout << pages << " " << page_size << " " << pages*page_size << std::endl;
	    //std::cout << m_photons.max_size() << " " << sizeof(photon) << " " << nPhotons*(maxDepth-1) << " " << sizeof(photon)*nPhotons*(maxDepth-1) << std::endl;
	    //m_photons.resize(nPhotons*(maxDepth-1));
	    //std::cout << "allocated success\n";
	    std::timed_mutex mtx;
	    std::vector<MemoryArena> photonShootArenas(MaxThreadIndex());
            ParallelFor([&](int photonIndex) {
            	MemoryArena &arena = photonShootArenas[ThreadIndex];
            	std::vector<photon> tmp_photons;
            	//uint64_t haltonIndex = nPhotons + photonIndex;
            	//int haltonDim = 0;
            	Float lightPdf;
            	RNG rng;
            	rng.SetSequence(photonIndex);
            	//Float lightSample = RadicalInverse(haltonDim++, haltonIndex);
            	Float lightSample = rng.UniformFloat();
            	int lightNum = lightDistr->SampleDiscrete(lightSample, &lightPdf);
            	const std::shared_ptr<Light> &light = scene.lights[lightNum];

            	/*Point2f uLight0(RadicalInverse(haltonDim  , haltonIndex), RadicalInverse(haltonDim+1, haltonIndex));
            	Point2f uLight1(RadicalInverse(haltonDim+2, haltonIndex), RadicalInverse(haltonDim+3, haltonIndex));
				float uLightTime = Lerp(RadicalInverse(haltonDim+4, haltonIndex), camera->shutterOpen, camera->shutterClose);
				haltonDim += 5;*/

            	Point2f uLight0 = {rng.UniformFloat(), rng.UniformFloat()};
            	Point2f uLight1 = {rng.UniformFloat(), rng.UniformFloat()};
            	float uLightTime = rng.UniformFloat();

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

					if (depth > 0){
						tmp_photons.emplace_back(photon{isect.p, beta, wo});
							//m_globalMap.put(photon{isect.p, beta, wo, true}, photonIndex*(maxDepth-1)+depth-1);
					}

					//Point2f bsdfSample(RadicalInverse(haltonDim, haltonIndex), RadicalInverse(haltonDim+1, haltonIndex));
					//haltonDim += 2;
					Point2f bsdfSample = {rng.UniformFloat(), rng.UniformFloat()};

					Spectrum fr = photonBSDF.Sample_f(wo, &wi, bsdfSample, &pdf, BSDF_ALL, &flags);

					if (fr.IsBlack() || pdf == 0.f) break;

					Spectrum bnew = beta * fr * AbsDot(wi, isect.shading.n) / pdf;

					//russian roulette
					Float q = std::max((Float)0, 1-bnew.y()/beta.y());
					//if (RadicalInverse(haltonDim++, haltonIndex) < q) break;
					if (rng.UniformFloat() < q) break;
					beta = bnew / (1-q);
					photonRay = (RayDifferential)isect.SpawnRay(wi);
				}
				while(!mtx.try_lock_for(std::chrono::milliseconds(10))){
				}
				for (int i = 0; i < tmp_photons.size(); i++){
					m_globalMap.put(tmp_photons[i]);
				}
				mtx.unlock();
				arena.Reset();
			}, nPhotons, 8192);
            std::cout << "Photon Pass finished" << std::endl;
            //std::cout << m_photons.size() << std::endl;
            LOG(INFO) << "Photon Pass finished";

            m_globalMap.build();
            /*for test
            std::priority_queue<PhotonIdx>ids;
            Point3f pp = Point3f(0.5f, 0.8f, 0.3f);
            m_globalMap.ksearch(pp, K, ids);
            m_globalMap.dump(pp, ids);
            m_globalMap.dump(pp);
            return;*/
			std::vector<MemoryArena> perThreadArenas(MaxThreadIndex());
			const Float invSqrtSpp = 1.f/spp;
			ProgressReporter reporter(nTiles.x * nTiles.y , "Camera pass");
			ParallelFor2D([&](Point2i tile){
				MemoryArena &arena = perThreadArenas[ThreadIndex];

				int tileIndex = tile.y * nTiles.x + tile.x;
				std::unique_ptr<Sampler> tileSampler = sampler.Clone(tileIndex);

				int x0 = pixelBounds.pMin.x + tile.x * tileSize;
				int x1 = std::min(x0 + tileSize, pixelBounds.pMax.x);
				int y0 = pixelBounds.pMin.y + tile.y * tileSize;
				int y1 = std::min(y0 + tileSize, pixelBounds.pMax.y);

				Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
				std::unique_ptr<FilmTile> filmTile = camera->film->GetFilmTile(tileBounds);

				for (Point2i pPixel : tileBounds){
					tileSampler->StartPixel(pPixel);
					do{
						CameraSample cameraSample = tileSampler->GetCameraSample(pPixel);
						RayDifferential ray;
						Spectrum beta = camera->GenerateRayDifferential(cameraSample, &ray);

						if (beta.IsBlack())
							continue;
						ray.ScaleDifferentials(invSqrtSpp);
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

							bool isDiffuse = bsdf.NumComponents(BxDFType(BSDF_DIFFUSE | BSDF_REFLECTION | BSDF_TRANSMISSION)) > 0;
							bool isGlossy = bsdf.NumComponents(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)) > 0;
							if (isDiffuse || (isGlossy && depth == maxDepth-1)){
								Spectrum densityEstimation = KNN(isect);
								L += beta * densityEstimation;
								break;
							}
							//sample the bsdf 
							if (depth < maxDepth - 1){
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
						filmTile->AddSample(cameraSample.pFilm, L);
						arena.Reset();
					}while(tileSampler->StartNextSample());
				}
				camera->film->MergeFilmTile(std::move(filmTile));
				reporter.Update();
			}, nTiles);
			std::cout << mean_r/(256.0*256.0) << std::endl;
 			reporter.Done();
			LOG(INFO) << "Rendering finished";
			camera->film->WriteImage();
	}


	Integrator *CreatePMIntegrator(const ParamSet &params,
                                 std::shared_ptr<const Camera> camera){
		int nPhotons = params.FindOneInt("nPhotons", 100000);
        int maxDepth = params.FindOneInt("maxdepth", 5);
	    int K = params.FindOneInt("K", 100);
	    int spp = params.FindOneInt("spp", 16);
	    Float initialRadius = params.FindOneFloat("initialRadius", RADIUS_THERSHOLD);
	    return new PMIntegrator(camera, nPhotons, K, spp, maxDepth, initialRadius*initialRadius);
	}
}