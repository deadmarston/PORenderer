/* This file is to implement the "Practical Path Guiding for Efficient Light-Transport Simulation", Thomas Muller, et al., 2017.
 * url: https://tom94.net/
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PGPATH_H
#define PBRT_INTEGRATORS_PGPATH_H

// integrators/path.h*
#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"

namespace pbrt{

//todo: single thread model for SDTree

//spatial binary tree
class SNode{
private:
  uint16_t axis;//change the axis alternative ly
};

class STree{

};
//directional quad tree
#define LEAFINDEX 0//utilize a magic number NODEINDEX to define the index of leaf
class DNode{
public:
  DNode();
  DNode(const DNode& node);
  bool isLeaf(int index) const {return child(index) == LEAFINDEX;}

  //the range of index [0,0] - [1,1]
  //each 0.5*0.5 part is a subnode
  uint16_t childIndex(Point2i& can) const;//get the child index
  void record(Point2i& can, float irradiance, std::vector<DNode>& nodes);
  void setSum(int index, Float val);
  Float sum(int index) const;
  void setChild(int index, uint16_t val);
  uint16_t child(int index) const;
  Float pdf(Point2i& can, const std::vector<DNode>& nodes) const;
  int depthAt(Point2i& can, const std::vector<DNode>& nodes) const;
  Point2f sample(Sampler* sampler, const std::vector<DNode>& nodes) const;
private:
  std::array<std::atomic<Float>, 4> m_sum;//record the irradiance
  std::array<uint16_t, 4> m_nodes;//store the index of children
};

class DTree{
public:
  DTree();
  int getMaxDepth();
  void sample();
  void submit();
  void refine();
  float pdf(Vector3f dir);

  Point2i dirToCanonical(Vector3f dir);//convert the dir vector to canonical 2d
  Vector3f canonicalToDir(Point2i canonical);//convert the canonical into vector
private:
  std::vector<DNode> m_tree;
  int maxDepth;
};

// Path Guiding Path Integrator Declarations
class PathGuidingIntegrator : public SamplerIntegrator{
  public:
  		PathGuidingIntegrator(int maxDepth, std::shared_ptr<const Camera> camera, 
  							  std::shared_ptr<Sampler> sampler,
  							  const Bounds2i &pixelBounds, Float rrThreshold = 1,
  							  const std::string &lightSampleStrategy = "spatial",
  							  const Float quadThreshold = 0.01f, const Float c = 12000);
  		void Preprocess(const Scene &scene, Sampler &sampler);
  		void Render(const Scene &scene);
  		Spectrum Li(const RayDifferential &ray, const Scene &scene, 
  					Sampler &sampler, MemoryArena &arena, int depth=0) const;
  private:
  		//private data for original path intergrator
  		const int maxDepth;//max depth to terminate the tracing
  		const Float rrThreshold;//russian roulette threshold, to terminate the tracing in advance
  		//light sample strategy
  		//uniform, spatial, power
  		//uniform: sample the light sources uniformly
  		//power:   sample the light sources according to their power
  		//spatial: compute the light contributions in regions of the scene and sample from a related distribution
  		const std::string lightSampleStrategy;

  		//light distribution created based on the light sample strategy
  		std::unique_ptr<LightDistribution>lightDistribution;

  		//private data for practical path guiding
  		//todo:
  		Float quadThreshold;//a threshold for dividing the quad
  		Float c;//a constant c trades off convergence of directional quadtrees with spatial resolution of the binary tree
  };

PathGuidingIntegrator *CreatePGPathIntegrator(const ParamSet &params,
											  std::shared_ptr<Sampler> sampler,
											  std::shared_ptr<const Camera> camera);

}//end pbrt

#endif PBRT_INTEGRATORS_PGPATH_H