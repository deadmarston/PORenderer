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
#include "interaction.h"

#include <iostream>
#include <vector>
#include <queue>

namespace pbrt {

	#define RADIUS_THERSHOLD 1.0f

	struct PhotonIdx{
		Float dis;
		int id;
		friend bool operator > (const struct PhotonIdx& photonIdx1, const struct PhotonIdx& photonIdx2){
			return photonIdx1.dis > photonIdx2.dis;
		}
		friend bool operator < (const struct PhotonIdx& photonIdx1, const struct PhotonIdx& photonIdx2){
			return photonIdx1.dis < photonIdx2.dis;
		}
	};

	struct photon{
		photon(){
		}
		photon(const photon& pho){
			p = pho.p;
			flux = pho.flux;
			wi = pho.wi;
		}
		photon(Point3f& _p, Spectrum& _flux, Vector3f& _wi)
			:p(_p), flux(_flux), wi(_wi)
		{		
		}
		Point3f p;
		Spectrum flux;
		Vector3f wi;
	};

	struct node{
		node(){
			idx = -1;
			left_idx = -1;
			right_idx = -1;
			//isLeaf = false;
		}
		int left_idx;
		int right_idx;
		int idx;
		//int isLeaf;
	};

	inline bool compare_x(const photon& p1, const photon& p2) { return (p1.p[0] < p2.p[0]); }
	inline bool compare_y(const photon& p1, const photon& p2) { return (p1.p[1] < p2.p[1]); }
	inline bool compare_z(const photon& p1, const photon& p2) { return (p1.p[2] < p2.p[2]); }


	//enum kdtype = {balanced=0}
	class kdTree{
	public:
		kdTree(){ 
			m_axis = 0; m_depth = 0;
			currentIndex = 0;
			compare_func[0] = compare_x;
			compare_func[1] = compare_y;
			compare_func[2] = compare_z; 
		};
		void build(std::vector<photon>& photons);
		int recursiveBuild(std::vector<photon>& photons, std::vector<photon>::iterator l, std::vector<photon>::iterator r, int axis, int depth);
		void ksearch(const Point3f& p, int K, int index, const std::vector<photon>& photons, std::priority_queue<PhotonIdx>& id, int axis, Float& radius);
		void ksearch(const Point3f& p, int K, const std::vector<photon>& photons, std::priority_queue<PhotonIdx>& id);
		void dump(const std::vector<photon>& photons, int index, int layer);
		void dump(const std::vector<photon>& photons, int index, int layer, Point3f& p);
		void setSearchRadius(Float square_d){
			m_radius = square_d;
		}
	private:
		int m_depth;
		int m_axis;
		Float m_radius;
		bool (*compare_func[3])(const photon& p1, const photon& p2);
		std::vector<node> m_nodes;
		int currentIndex;
		//kdtype m_type;
	};

	class photonMap{
	private:
		kdTree m_kdtree;
		std::vector<photon> m_photons;
		Float initialRadius;
	public:
		void put(photon phtn){
			//m_photons[index] = phtn;
			m_photons.emplace_back(phtn);
		}
		photonMap(){
		}
		void init(Float _initialRadius, int size){
			initialRadius = _initialRadius;
			//m_photons.resize(size);
		}
		void build(){
			std::cout << m_photons.size() << std::endl;
			m_kdtree.build(m_photons);
			m_kdtree.setSearchRadius(initialRadius);
		}
		void dump(){
			std::cout << "====================\n";
			m_kdtree.dump(m_photons, 0, 0);
			std::cout << std::endl;
		}
		void dump(Point3f& p){
			std::cout << "====================\n";
			m_kdtree.dump(m_photons, 0, 0, p);
			std::cout << std::endl;
		}
		void dump(Point3f& p, std::priority_queue<PhotonIdx>& id){
			std::cout << "====================\n";
			std::cout << "initial radius: " << initialRadius << std::endl;
			std::cout << "center: " << p << std::endl;
			while(!id.empty()){
				std::cout << "ksearch: " << m_photons[id.top().id].p << " " << DistanceSquared(p, m_photons[id.top().id].p) << std::endl;
				id.pop();
			}
		}
		Spectrum gatherKPhotons(const SurfaceInteraction* isect, std::priority_queue<PhotonIdx>& ids);
		Float linearSearch(const Point3f& p, int K, std::priority_queue<PhotonIdx>& id);
		Float ksearch(const Point3f& p, int K, std::priority_queue<PhotonIdx>& id){
			m_kdtree.ksearch(p, K, m_photons, id);
			return id.empty()?initialRadius:id.top().dis;
		}
	};

	// SPPM Declarations
class PMIntegrator : public Integrator {
  public:
    // SPPMIntegrator Public Methods
    PMIntegrator(std::shared_ptr<const Camera> &camera, int nPhotons, int K, int spp, int maxDepth, Float initialRadius)
        : camera(camera),
          maxDepth(maxDepth),
          nPhotons(nPhotons),
          K(K),
          spp(spp),
          initialRadius(initialRadius)
          {
          	m_globalMap.init(initialRadius, nPhotons*(maxDepth-1));
          };
    void Render(const Scene &scene);
    Spectrum KNN(const Interaction &it);
    //Spectrum radiusSearch(const Interaction &it, std::vector<photon> photonList);
    //Spectrum Gaussian(const Interaction &it, std::vector<photon> photonList);
    //Spectrum visualizePMMap(const Interaction &it, std::vector<photon> photonList);
  private:
    // SPPMIntegrator Private Data
    std::shared_ptr<const Camera> camera;
    const int maxDepth;
    const int nPhotons;
    const int K; //KNN search
    const int spp; //spp for aa, motion blur, dof
    const Float initialRadius;
    photonMap m_globalMap;
};

Integrator *CreatePMIntegrator(const ParamSet &params,
                                 std::shared_ptr<const Camera> camera);

}  // namespace pbrt
#endif