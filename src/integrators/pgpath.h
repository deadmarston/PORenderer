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

//magic number
enum RecordType{ nearest = 0, filter };

enum BudgetStrategy{ sample = 0, time, naive}; 

//directional quad tree
#define LEAFINDEX 0//utilize a magic number NODEINDEX to define the index of leaf
#define EPSILON 0.00001
#define QUAD_MAX_DEPTH 20

#define VERTEX_MAX_DEPTH 32//according to the original paper
#define LEARNING_MAX_INTERATION 10//max interations for learning and rendering

#define PG_DEBUG

//todo: single thread model for SDTree

class DNode{
public:
  DNode();
  DNode(const DNode& node);
  bool isLeaf(int index) const { return child(index) == LEAFINDEX; }

  //the range of index [0,0] - [1,1]
  //each 0.5*0.5 part is a subnode
  int childIndex(Point2f& can) const;//get the child index
  void record(Point2f& can, Float irradiance, std::vector<DNode>& nodes);
  void setSum(int index, Float val);
  Float sum(int index) const;
  void setChild(int index, uint16_t val);
  uint16_t child(int index) const;
  Float pdf(Point2f& can, const std::vector<DNode>& nodes) const;
  int depthAt(Point2f& can, const std::vector<DNode>& nodes) const;
  Point2f sample(Sampler* sampler, const std::vector<DNode>& nodes) const;
  Float build(std::vector<DNode>& nodes);

  DNode& operator=(const DNode& node);

private:
  std::array<std::atomic<Float>, 4> m_sum;//record the irradiance
  std::array<uint16_t, 4> m_nodes;//store the index of children
};

class DTree{
public:
  DTree();
  int getMaxDepth();
  int depthAt(const Vector3f& dir);
  Vector3f sample(Sampler* sampler);
  void record(const Vector3f& dir, Float irradiance, RecordType type);
  void refine(DTree* sampling, Float quadThreshold);
  void reset();
  Float pdf(const Vector3f& dir);

  Point2f dirToCanonical(const Vector3f& dir);//convert the dir vector to canonical 2d
  Vector3f canonicalToDir(const Point2f& canonical);//convert the canonical into vector

  int numOfChildren() const;
  DNode nodeAtIndex(int index) const;

  void dump() const;
  Float numOfSample() const { return samples.load(std::memory_order_relaxed); };//todo

  DTree& operator=(const DTree& tree);

  std::vector<DNode> m_tree;//maintain a tree to store the index of node, the index of root is 0
  int maxDepth;

  std::atomic<Float> samples;
};

class DTreeWrapper{
public:
    DTreeWrapper();
    DTreeWrapper& operator=(const DTreeWrapper& dwrapper);
    Vector3f sample(Sampler* sampler);
    Float pdf(const Vector3f& dir);
    void record(const Vector3f& dir, Spectrum irradiance, RecordType type);//add recordType

    void dump();
    void rebuild();

    Float getNumOfSamples() const { return building.numOfSample(); };

    void refine(Float quadThreshold);

    DTree sampling;
    DTree building;
};

//spatial binary tree
class SNode{
public:
  SNode();
  SNode(uint16_t _axis);
  SNode(const SNode& node);
  SNode& operator=(const SNode& node);

  //const SNode* acquire(Point3f& pos, std::vector<SNode> nodes) const;
  DTreeWrapper* acquireDTreeWrapper(Point3f& pos, std::vector<SNode>& nodes);
  bool isLeaf(int index) const { return child(index) == LEAFINDEX; };
  uint32_t child(int index) const;
  int childIndex(Point3f& pos) const; 
  void setChild(int i, uint32_t index);
  int depthAt(Point3f& pos, std::vector<SNode>& nodes) const;

  Vector3f sample(Sampler* sampler);
  Float pdf(const Vector3f& dir);
  void record(const Vector3f& dir, Spectrum& irradiance, RecordType type);
  void refine();
  void dump(std::vector<SNode>& nodes);
  bool shallDivide(Float divideThreshold);

  DTreeWrapper wrapper;
  uint16_t axis;//change the axis alternatively
  bool isleaf;
private:
  std::array<uint32_t, 2> m_nodes;
};

class STree{
public:
  STree(Bounds3f bounds);
  int getMaxDepth() const;
  int depthAt(Point3f& pos);
  const Bounds3f& bounds() const{
    return m_bounds;
  }
  DTreeWrapper* acquireDTreeWrapper(Point3f pos);
  void dump();
  void refine();

  void refineDTree();
  void divideSTree();
  void divideSNode(int id);
private:
  void normalize(Point3f& pos) const;

  std::vector<SNode> nodes;
  int maxDepth;
  Bounds3f m_bounds;
  Float m_extent;
};

struct RecordVertex{
  Spectrum radiance;
  DTreeWrapper* dwrapper;
  Vector3f wi;
  RecordType type;

  RecordVertex(DTreeWrapper* _wrapper, Vector3f _wi) 
              : dwrapper(_wrapper), wi(_wi){
                radiance = Spectrum(0.f);
                type = nearest;
              }

  RecordVertex(DTreeWrapper* _wrapper, Spectrum _radiance, Vector3f  _wi)
                : dwrapper(_wrapper),
                  radiance(_radiance),
                  wi(_wi)
                {
                  type = nearest;
                }

  RecordVertex(DTreeWrapper* _wrapper, Spectrum _radiance, Vector3f  _wi, RecordType _type)
                : dwrapper(_wrapper),
                  radiance(_radiance),
                  wi(_wi),
                  type(_type)
                {
                }

  RecordVertex(){
    dwrapper = nullptr;
    radiance = Spectrum(0.f);
    type = nearest;
    wi = Vector3f(0.f, 0.f, 0.f);
  }

  void commit(){
    //todo: submit the irrad
    if (type == nearest){
      dwrapper->record(wi, radiance.Average(), type);
    }
  }

  void record(Spectrum _radiance){
    radiance += _radiance;
  }
};

// Path Guiding Path Integrator Declarations
class PathGuidingIntegrator : public Integrator{
  public:
  		PathGuidingIntegrator(int maxDepth, std::shared_ptr<const Camera> camera, 
                  const int spp, 
  							  const Bounds2i &pixelBounds, Float rrThreshold = 1,
  							  const std::string &lightSampleStrategy = "spatial",
  							  const Float quadThreshold = 0.01f, const Float c = 12000);
  		void Preprocess(const Scene &scene, Sampler &sampler);
  		void Render(const Scene &scene);
  		Spectrum Li(const RayDifferential &ray, const Scene &scene, 
  					Sampler &sampler, MemoryArena &arena, int depth=0) const;
      void figureRenderBudgets();
  private:
      std::shared_ptr<const Camera> camera;
      const Bounds2i pixelBounds;

      STree* m_sdtree;

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

      int spp; 
      std::array<int, LEARNING_MAX_INTERATION> iterations;
      int numOfIterations;

  		//private data for practical path guiding
  		//todo:
  		const Float quadThreshold;//a threshold for dividing the quad
  		const Float c;//a constant c trades off convergence of directional quadtrees with spatial resolution of the binary tree
      
      BudgetStrategy budgetStrategy;
  };

PathGuidingIntegrator *CreatePGPathIntegrator(const ParamSet &params,
											  std::shared_ptr<const Camera> camera);

}//end pbrt

#endif PBRT_INTEGRATORS_PGPATH_H