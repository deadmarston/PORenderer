/* This file is to implement the "Practical Path Guiding for Efficient Light-Transport Simulation", Thomas Muller, et al., 2017.
 * url: https://tom94.net/
 */


#include "integrators/pgpath.h"
#include "samplers/random.h"
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
#include <stack>
using namespace std;
//12.7 , implement the original path tracer
//12.23, add comments for mis with nee, preparing to finish it asap
namespace pbrt {
	STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);

	//debug
  	static int iter = 0;

	bool areSame(Float a, Float b){
	  	return fabs(a-b) < EPSILON;
	}

	bool areFloatSame(Float a, Float b){
		cout << a << " " << b << endl;
		return areSame(a, b);
	}

	bool arePoint2Same(Point2f a, Point2f b){
		cout << a << " " << b << endl;
		return areSame(a.x, b.x) && areSame(a.y, b.y);
	}

	bool areVector3fSame(Vector3f a, Vector3f b){
		cout << a << " " << b << endl;
		return areSame(a.x, b.x) && areSame(a.y, b.y) && areSame(a.z, b.z);
	}

	void testDNode(){
		cout << "===============test dnode===============" << endl;
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
		DNode chi;
		nodes.push_back(chi);
		dnode.setChild(0, 1);
		assert(!dnode.isLeaf(0));
		assert(chi.isLeaf(0));
		assert(dnode.child(0) == 1);
		chi.setSum(0, 0.3);
		cout << "dnode setchild pass" << endl;
		Point2f p1 = Point2f(0.25, 0.25);
		Point2f p2 = Point2f(0.75, 0.75);
		Point2f p3 = Point2f(0.25, 0.75);
		Point2f p4 = Point2f(0.75, 0.25);
		Point2f p5 = Point2f(0.0, 0.0);
		Point2f p6 = Point2f(1.0, 1.0);
		Point2f p7 = Point2f(0.5, 0.5);
		Point2f p8 = Point2f(0.0, 0.5);
		Point2f p9 = Point2f(0.5, 0.0);
		assert(dnode.childIndex(p1) == 0);
		assert(dnode.childIndex(p2) == 3);
		assert(dnode.childIndex(p3) == 2);
		assert(dnode.childIndex(p4) == 1);
		assert(dnode.childIndex(p5) == 0);
		assert(dnode.childIndex(p6) == 3);
		assert(dnode.childIndex(p7) == 0);
		assert(dnode.childIndex(p8) == 0);
		assert(dnode.childIndex(p9) == 0);
		cout << "dnode childindex pass" << endl;
		dnode.setSum(1, 0.5);
		assert(dnode.sum(1) == 0.5);
		cout << "dnode setsum pass" << endl;
		p1 = Point2f(0.8, 0.8);
		p2 = Point2f(0.2, 0.2);
		assert(dnode.depthAt(p1, nodes) == 1);
		assert(dnode.depthAt(p2, nodes) == 2);
		cout << "dnode depthat pass" << endl;
		Point2f p = Point2f(0.1, 0.1);
		dnode.record(p, 0.3, nodes);
		assert(nodes[1].sum(0) == 0.3f);
		p = Point2f(0.4, 0.4);
		dnode.record(p, 0.5, nodes);
		assert(nodes[1].sum(3) == 0.5f);
		cout << "dnode record pass" << endl;

		dnode.setSum(0, 0.0);
		dnode.setSum(1, 0.0);
		dnode.setSum(2, 0.0);
		dnode.setSum(3, 0.0);

		nodes[1].setSum(0, 0.0);
		nodes[1].setSum(1, 0.0);
		nodes[1].setSum(2, 0.0);
		nodes[1].setSum(3, 0.0);
		p1 = Point2f(0.1, 0.1);
		p2 = Point2f(0.3, 0.3);
		p3 = Point2f(0.75, 0.75);
		dnode.record(p1, 0.8, nodes);
		dnode.record(p2, 0.3, nodes);
		dnode.record(p3, 0.4, nodes);
		Float sum = dnode.build(nodes);
		assert(sum == 1.5f);
		assert(dnode.sum(0) == 1.1f);
		cout << "dnode build pass" << endl;
		p1 = Point2f(0.1, 0.1);
		p3 = Point2f(0.75, 0.75);
		assert(areFloatSame(dnode.pdf(p1, nodes), 8.53333));
		assert(areFloatSame(dnode.pdf(p3, nodes), 1.06667));
		cout << "dnode pdf pass" << endl;
		//todo:
		//add testcase for sample function 
		cout << "============test dnode finish============" << endl;
	}

	void testDTree(){
		cout << "===============test dtree===============" << endl;
		DTree tree;
		assert(tree.getMaxDepth() == 1);
		cout << "dtree maxdepth pass" << endl;
		Vector3f v1 = Vector3f(0.0, 0.0, 0.0);
		Vector3f v2 = Vector3f(1.0, 1.0, 1.0);
		Vector3f v3 = Vector3f(-1.0, -1.0, -1.0);
		Vector3f v4 = Vector3f(0.0, 0.0, 5.0);
		Vector3f v5 = Vector3f(0.0, 0.0, -5.0);
		assert(arePoint2Same(tree.dirToCanonical(v1), Point2f(0.5, 0.0)));
		assert(arePoint2Same(tree.dirToCanonical(v2), Point2f(1.0, 0.125)));
		assert(arePoint2Same(tree.dirToCanonical(v3), Point2f(0.0, 0.625)));
		assert(arePoint2Same(tree.dirToCanonical(v4), Point2f(1.0, 0.0)));
		assert(arePoint2Same(tree.dirToCanonical(v5), Point2f(0.0, 0.0)));
		cout << "dtree dirToCanonical pass" << endl;
		Point2f p1 = Point2f(0.0, 0.0);
		Point2f p2 = Point2f(1.0, 1.0);
		Point2f p3 = Point2f(0.5, 0.5);
		assert(areVector3fSame(tree.canonicalToDir(p1), Vector3f(0.0, 0.0, -1.0)));
		assert(areVector3fSame(tree.canonicalToDir(p2), Vector3f(0.0, 0.0, 1.0)));
		assert(areVector3fSame(tree.canonicalToDir(p3), Vector3f(-1.0, 0.0, 0.0)));
		cout << "dtree canonicalToDir pass" << endl;
		//todo:
		//add testcase for sample function
		cout << "============test dtree finish============" << endl;
	}

	void testSNode(){
		cout << "===============test snode===============" << endl;
		SNode snode;
		std::vector<SNode> nodes;
		nodes.push_back(snode);
		assert(snode.isLeaf(0));
		assert(snode.isLeaf(1));
		assert(snode.isLeaf(2));
		assert(snode.isLeaf(3));
		cout << "snode isleaf pass" << endl;

		cout << "============test snode finish============" << endl;
	}	

	void testSTree(){
		cout << "===============test stree===============" << endl;
		cout << "============test stree finish============" << endl;
	}

	void testSDTree(){
		cout << "===============test begin===============" << endl;
		testDNode();
		testDTree();
		testSNode();
		testSTree();
		cout << "================test end================" << endl;
	}

	static void addToAtomicFloat(std::atomic<Float>& var, Float val){
		auto current = var.load();
		while(!var.compare_exchange_weak(current, current + val));
		current = var.load();
	}

	DTreeWrapper::DTreeWrapper(){
	}

	DTreeWrapper& DTreeWrapper::operator=(const DTreeWrapper& dwrapper){
		sampling = dwrapper.sampling;
		building = dwrapper.building;
	}

	Vector3f DTreeWrapper::sample(Sampler* sampler){
		return sampling.sample(sampler);
	}

	Float DTreeWrapper::pdf(const Vector3f& dir){
		return sampling.pdf(dir);
	}

	void DTreeWrapper::record(const Vector3f& dir, Spectrum irradiance, RecordType type = nearest){
		building.record(dir, irradiance.Average(), type);
	}

	void DTreeWrapper::dump(){
		building.dump();
	}

	STree::STree(Bounds3f bounds){
		m_bounds = bounds;
		Vector3f size = m_bounds.pMax - m_bounds.pMin;
		m_extent = max(max(size.x, size.y), size.z);
		m_bounds.pMax = m_bounds.pMin + Vector3f(m_extent, m_extent, m_extent);//todo: make it a cube, needs to figure out why
		maxDepth = 0;
		nodes.emplace_back();
	}

	void STree::refine(){
		//first, divide the stree
		divideSTree();
		//then, for all the dtreewrapper in the leaf node, divide the dtree
	}

	void STree::divideSNode(int id){
		nodes.resize(nodes.size() + 2);//todo: for safety, we need to check the size of nodes

		uint16_t axis = (nodes[id].axis+1)%3;

		for (int j = 0; j < 2; j++){
			int index = nodes.size() - 2 + j;
			nodes[id].setChild(j, index);
			nodes[index].axis = axis;
			nodes[index].wrapper = nodes[id].wrapper;
			nodes[index].isleaf = true;
		}
		nodes[id].wrapper = {};
		nodes[id].isleaf = false;
	}

	void STree::divideSTree(){
		//todo: chech whether the memory usage is over the limit
		std::stack<int> ids;//record the index of nodes
		ids.push(0);

		while (!ids.empty()){//we try to visit all nodes in a bfs way
			int id = ids.top();
			ids.pop();
			if (nodes[id].isleaf){//divide the leaf
				if (nodes[id].shallDivide(1)){//todo: the threshold should be a number controlled by user
					divideSNode(id);
				}
			}else{
				for (int j = 0; j < 2; j++){
					ids.push(nodes[id].child(j));
				}
			}
		}
	}

	void STree::dump(){
		cout << "stree size: " << nodes.size() << endl;
		//dfs to find the leaf
		cout << "============dump sdtree============" << endl;
		for (int i = 0; i < nodes.size(); i++){
			cout << "dump node " << i << endl;
			cout << "isleaf: " << nodes[i].isleaf
				 << " axis: " << nodes[i].axis
				 << endl; 
			if (!nodes[i].isleaf){
				cout << "child1: " << nodes[i].child(0)
					 << " child2: " << nodes[i].child(1)
					 << endl;
			}
 		}

		cout << "============dump sdtree finished============" << endl;
		nodes[0].dump(nodes); 
	}

	int STree::getMaxDepth() const{
		return maxDepth;
	}

	void STree::normalize(Point3f& pos) const{
		Vector3f v = pos-m_bounds.pMin;
		v.x /= m_extent;
		v.y /= m_extent;
		v.z /= m_extent;
		pos.x = v.x;
		pos.y = v.y;
		pos.z = v.z;
	}

	int STree::depthAt(Point3f& pos){
		Point3f _pos = pos;
		normalize(pos);
		return nodes[0].depthAt(_pos, nodes);
	}

	/*const SNode* STree::acquireSNode(Point3f& pos){
		Point3f _pos = pos;
		normalize(pos);
		return nodes[0].acquire(pos, nodes);
	}*/

	DTreeWrapper* STree::acquireDTreeWrapper(Point3f pos){
		normalize(pos);
		return nodes[0].acquireDTreeWrapper(pos, nodes);
	}

	SNode::SNode(){
		axis = 0;
		isleaf = true;
		for (size_t i = 0; i < m_nodes.size(); i++){
			m_nodes[i] = LEAFINDEX;
		}
	}

	SNode::SNode(uint16_t _axis) : axis(_axis){
		isleaf = true;
		for (size_t i = 0; i < m_nodes.size(); i++){
			m_nodes[i] = LEAFINDEX;
		}
	}

	SNode::SNode(const SNode& node){
		axis = node.axis;
		wrapper = node.wrapper;
		isleaf = node.isleaf;
		//current = node.current;
		//previous = node.previous;
		for (size_t i = 0; i < m_nodes.size(); i++){
			m_nodes[i] = node.child(i);
		}
	}

	SNode& SNode::operator=(const SNode& node){
		axis = node.axis;
		wrapper = node.wrapper;
		isleaf = node.isleaf;
		//current = node.current;
		//previous = node.previous;
		for (size_t i = 0; i < m_nodes.size(); i++){
			m_nodes[i] = node.child(i);
		}
	}

	void SNode::setChild(int i, uint32_t index){
		m_nodes[i] = index;
	}

	bool SNode::shallDivide(Float divideThreshold){
		return wrapper.getNumOfSamples() > divideThreshold;
	}

	void SNode::dump(std::vector<SNode>& nodes){
		if (isleaf){
			cout << "dump stree" << endl;
			cout << "dtreewrapper address: " << &wrapper << endl;
			wrapper.dump();
		}else{
			cout << "divide stree into two nodes" << endl;
			cout << "dtreewrapper address: " << &wrapper << endl;
			for (int i = 0; i < 2; i++){
				nodes[child(i)].dump(nodes);
			}
		}
	}

	/*const SNode* SNode::acquire(Point3f& pos, std::vector<SNode> nodes) const{
		int index = childIndex(pos);
		if (isLeaf(index)){
			return this;
		}else{
			return nodes[child(index)].acquire(pos, nodes);
		}
	}*/
	DTreeWrapper* SNode::acquireDTreeWrapper(Point3f& pos, std::vector<SNode>& nodes){
		int index = childIndex(pos);
		if (isleaf){
			return &wrapper;
		}else{
			return nodes[child(index)].acquireDTreeWrapper(pos, nodes);
		}
	}

	uint32_t SNode::child(int index) const{
		return m_nodes[index];
	}

	int SNode::childIndex(Point3f& pos) const{
		/*
		left <---axis 0.5---> right
		*/
		if (pos[axis] > 0.5){
			pos[axis] = (pos[axis]-0.5)*2;
			return 1;
		}else{
			pos[axis] = pos[axis]*2;
			return 0;
		}
	}

	int SNode::depthAt(Point3f& pos, std::vector<SNode>& nodes) const{
		int index = childIndex(pos);

		if (isleaf){
			return 1;
		}else{
			return 1 + nodes[child(index)].depthAt(pos, nodes);
		}
	}

	/*Vector3f SNode::sample(Sampler* sampler){
		return previous.sample(sampler);
	}

	Float SNode::pdf(const Vector3f& dir){
		return previous.pdf(dir)/(4*M_PI);
	}

	void SNode::record(const Vector3f& dir, Spectrum& irradiance, RecordType type = nearest){
		Float average = irradiance.Average();;
		current.record(dir, average, type);
	}*/

	DNode::DNode(){
		for (size_t i = 0; i < m_sum.size(); i++){
			m_sum[i].store(0, std::memory_order_relaxed);
			m_nodes[i] = LEAFINDEX;
		}
	}

	Float DNode::build(std::vector<DNode>& nodes){
		Float _sum = 0.f;
		for (size_t i = 0; i < m_nodes.size(); i++){
			if (isLeaf(i)){
				_sum += sum(i);
			}else{
				Float current = nodes[child(i)].build(nodes);
				setSum(i, current);
				_sum += current;
			}
		}
		return _sum;
	}

	DNode::DNode(const DNode& node){
		for (int i = 0; i < node.m_sum.size(); i++){
			setSum(i, node.sum(i));
			m_nodes[i] = node.child(i);
		}
	}

	DNode& DNode::operator=(const DNode& node){
		for (int i = 0; i < m_sum.size(); i++){
			setSum(i, node.sum(i));
			if (m_nodes[i] > 0) printf("wtf\n");
			m_nodes[i] = node.child(i);
		}
		return *this;
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

	void DNode::record(Point2f& can, Float irradiance, std::vector<DNode>& nodes){
		int index = childIndex(can);
		//didn't add to the sum of current node for efficiency
		if (isLeaf(index)){
			addToAtomicFloat(m_sum[index], irradiance);//make it atomic due to parallel programming
		}else{
			nodes[child(index)].record(can, irradiance, nodes);
		}
	}

	int DNode::depthAt(Point2f& can, const std::vector<DNode>& nodes) const{
		int index = childIndex(can);

		if (isLeaf(index)){
			return 1;
		}else{
			return 1 + nodes[child(index)].depthAt(can, nodes);
		}
	}

	Float DNode::pdf(Point2f& can, const std::vector<DNode>& nodes) const{
		int index = childIndex(can);

		if (sum(index) <= 0){
			return 0;
		}
		Float factor = 4*sum(index) / (sum(0)+sum(1)+sum(2)+sum(3));
		if (isLeaf(index)){
			return factor;
		}else{
			return factor*nodes[child(index)].pdf(can, nodes);
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

	int DNode::childIndex(Point2f& can) const{
		/*
		-------
		|2 |3 |
		-------
		|0 |1 |
		-------
		*/
		assert(can.x >= 0 && can.y >= 0 && can.x <= 1 && can.y <= 1);
		int index = 0;
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

	DTree::DTree(){
		maxDepth = 1;
		m_tree.emplace_back();
		samples.store(0., std::memory_order_relaxed);
	}

	DTree& DTree::operator=(const DTree& dtree){
		m_tree = dtree.m_tree;
		/*m_tree = std::vector<DNode>(dtree.m_tree.size());
		for (int i = 0; i < dtree.m_tree.size(); i++){
			m_tree[i] = dtree.m_tree[i];
		}*/
		maxDepth = dtree.maxDepth;
		samples.store(dtree.numOfSample());
		return *this;
	}

	int DTree::numOfChildren() const{
		return m_tree.size();
	}

	DNode DTree::nodeAtIndex(int index) const{
		return m_tree[index];
	}

	int DTree::getMaxDepth(){
		return maxDepth;
	}

	int DTree::depthAt(const Vector3f& dir){
		Point2f can = dirToCanonical(dir);
		return m_tree[0].depthAt(can, m_tree);
	}

	Vector3f DTree::sample(Sampler* sampler){
		Point2f can = m_tree[0].sample(sampler, m_tree);
		return canonicalToDir(can);
	}

	void DTree::record(const Vector3f& dir, Float irradiance, RecordType type = nearest){
		Point2f can = dirToCanonical(dir);
		addToAtomicFloat(samples, 1.0);
		//todo: improvement for filter recording
		if (type == nearest){
			m_tree[0].record(can, irradiance, m_tree);
		}
	}

	void DTree::refine(){

	}

	void DTree::dump() const{
		for (int i = 0; i < m_tree.size(); i++){
			cout << "dnode:" << i << " address: " << &m_tree[i] << endl;
			for (int j = 0; j < 4; j++){
				cout << "child: " << j 
					 << " isLeaf: " << m_tree[i].isLeaf(j)
					 << " sum: " << m_tree[i].sum(j)
					 << " index: " << m_tree[i].child(j)
					 << endl;
			}
		}
	}

	Float DTree::pdf(const Vector3f& dir){
		Point2f can = dirToCanonical(dir);
		return m_tree[0].pdf(can, m_tree);
	}

	Point2f DTree::dirToCanonical(const Vector3f& dir){
		/*Point2f can;
		float cosTheta = max(min(dir.z, 1.f), -1.f);
		float phi = atan2(dir.y, dir.x);
		while(phi < 0) phi += 2*M_PI;
		can.x = (cosTheta+1)/2.f;
		can.y = phi/(2*M_PI);
		return can;*/
		Point2f can;
		float cosTheta = max(min(dir.y, 1.f), -1.f);
		float phi = atan2(dir.z, dir.x);
		while(phi < 0) phi += 2*M_PI;

		can.x = (cosTheta+1)/2.f;
		can.y = phi/(2*M_PI);

		return can;
	}

	Vector3f DTree::canonicalToDir(const Point2f& can){
		assert(can.x >= 0 && can.y >= 0 && 1 >= can.x &&  1>= can.y);
		/*float cosTheta = can.x*2-1;
		float sinTheta = sqrt(1 - cosTheta*cosTheta);
		float phi = can.y*2*M_PI;
		return Vector3f(sinTheta*cos(phi), sinTheta*sin(phi), cosTheta);*/
		float cosTheta = can.x*2-1;
		float sinTheta = sqrt(1-cosTheta*cosTheta);
		float phi = can.y*2*M_PI;
		return Vector3f(sinTheta*cos(phi), sinTheta*sin(phi), cosTheta);
	}
	
	Spectrum EstimateDirect_PG(const Interaction &it, const Point2f &uScattering,
				const Light &light, const Point2f &uLight,
				const Scene &scene, Sampler &sampler,
				MemoryArena &arena, DTreeWrapper* dwrapper, bool handleMedia = false, bool specular = false) {
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
				    }else{
						VLOG(2) << "  shadow ray unoccluded";
				    }
			    }

			    // Add light's contribution to reflected radiance
			    if (!Li.IsBlack()) {
					if (IsDeltaLight(light.flags)){
					    Ld += f * Li / lightPdf;
					    RecordVertex v(dwrapper, f * Li/ lightPdf, wi);
					    v.commit();
					}
					else{
						//todo: consider the sd-tree, and utilize mixed pdf to compute the power heuristic
						//ex: mixed_pdf, d_pdf, b_pdf = pdf(direction)
						//    weight = powerheuristic(mixed_pdf, lightPdf)
						//    compute the Ld
					    Float weight =
						PowerHeuristic(1, lightPdf, 1, scatteringPdf);
					    Ld += f * Li * weight / lightPdf;
					    RecordVertex v(dwrapper, f * Li * weight / lightPdf, wi);
					    v.commit();
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
			    if (!Li.IsBlack()){
			    	Ld += f * Li * Tr * weight / scatteringPdf;
			    	dwrapper->record(wi, Li * Tr * weight / scatteringPdf);
			    }
			}
	    }
	    return Ld;
	}

	Spectrum UniformSampleOneLight_PG(const Interaction &it, const Scene &scene,
					       MemoryArena &arena, Sampler &sampler,
					       bool handleMedia, DTreeWrapper* dwrapper, const Distribution1D *lightDistrib = nullptr) {
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
	    return EstimateDirect_PG(it, uScattering, *light, uLight,
                          scene, sampler, arena, dwrapper, handleMedia) / lightPdf;
	}
	
	PathGuidingIntegrator::PathGuidingIntegrator(int maxDepth, std::shared_ptr<const Camera> camera, 
  							  const int spp,
  							  const Bounds2i &pixelBounds, Float rrThreshold,
  							  const std::string &lightSampleStrategy,
  							  const Float quadThreshold, const Float c)
		: camera(camera),
		  pixelBounds(pixelBounds),
		  spp(spp),
		  maxDepth(maxDepth),
		  rrThreshold(rrThreshold),
		  lightSampleStrategy(lightSampleStrategy),
		  quadThreshold(quadThreshold),
		  c(c)
		   {
		   		budgetStrategy = naive;
		   		figureRenderBudgets();
		   }

    void PathGuidingIntegrator::figureRenderBudgets(){
    	switch (budgetStrategy){//todo: finish the configure of rendering budget 
    		case naive: default: {
    			int budget = 1;
    			int iter = 0;
    			int remain = spp;
    			do{
    				iterations[iter] = budget;
    				remain -= budget;
    				budget *= 2;
    			} while (iter++ < LEARNING_MAX_INTERATION && remain > 0);
    			numOfIterations = iter;
    			if (numOfIterations == LEARNING_MAX_INTERATION){
    				iterations[numOfIterations-1] += remain;
    			}
    		}	
    	}
    }

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
  		m_sdtree = new STree(scene.WorldBound());//build the initial sdtree

  		for (iter = 0; iter < numOfIterations; iter++){ //progressive rendering
  			cout << "the " << iter+1 << "th rendering with budget " <<  iterations[iter] << endl;
  			std::shared_ptr<Sampler> sampler = std::shared_ptr<Sampler>(CreateRandomSampler(iterations[iter]));

	  		//first, we preprocess the scene
	  		Preprocess(scene, *sampler);

	  		//compute nTiles, for parallel programming
	  		const int tileSize = 16;
	  		Bounds2i sampleBounds = camera->film->GetSampleBounds();
	  		Vector2i sampleExtent = sampleBounds.Diagonal();

	  		Point2i nTile((sampleExtent.x+tileSize-1)/tileSize, (sampleExtent.y+tileSize-1)/tileSize); 

	  		m_sdtree->dump();

	  		ProgressReporter reporter(nTile.x * nTile.y, "Rendering");
	    	{
	    		ParallelFor2D([&](Point2i tile) {//solve it parallelly

	    			//allocate memory arena for tile
	    			MemoryArena arena;

	    			//seed
	    			int seed = (tile.y * nTile.x + tile.x)*iterations[iter];
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

	                        if (iter == numOfIterations-1){
		                        //add ray's contribution to image
		                        filmTile->AddSample(cameraSample.pFilm, L, rayWeight);
		                    }

	                        arena.Reset();
		                }while(tileSampler->StartNextSample());
	    			};
	    			LOG(INFO) << "Finished image tile " << tileBounds;

					camera->film->MergeFilmTile(std::move(filmTile));
					reporter.Update();    		
	    		}, nTile);
	    		//todo: collect the sd tree, refine it
	    		//print the information of sdtree
	    		//
	    		reporter.Done();
	    	};
	    	m_sdtree->dump();
	    	cout << "begin refining" << endl;
	    	m_sdtree->refine();
	    	cout << "refining finished\n";
	    	m_sdtree->dump();


	    	LOG(INFO) << "Rendering finished";
	    }
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

  		//todo: add vertex to record the irradiance, and the submit it to the dtree finally
  		//the struct of vertex requires stored irradiance, dtreewrapper
  		
  		int nVertex = 0;

  		std::array<RecordVertex, VERTEX_MAX_DEPTH> vertex;

  		auto recordRadianceForVertex = [&](Spectrum& radiance){
  			for (int i = 0; i < nVertex; i++){
  				vertex[i].record(radiance);
  			}
  		};

  		//tracing loop
  		for (bounces = 0; ; bounces++)
  		{
	  		//intersection
	  		SurfaceInteraction isect;
	  		bool foundIntersection = scene.Intersect(ray, &isect);

	  		if (bounces == 0 || specularBounce){//todo: need to figure out
	  			//add emitted light at path vertex
	  			if (foundIntersection){
	  				L += beta * isect.Le(-ray.d);
	  			}else{
	  				for (const auto& light : scene.infiniteLights){
	  					L += beta * light->Le(ray);
	  				}
	  			}
	  		}

	  		if (!foundIntersection || bounces >= maxDepth){
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

	  		//acquire the dtree from the sdtree
	  		DTreeWrapper* dwrapper = m_sdtree->acquireDTreeWrapper(isect.p);
	  		//direct lighting computation
	  		//L += beta * UniformSampleOneLight(isect, scene, arena, sampler);

	  		//direct lighting computation with light distribution
	  		const Distribution1D *distrib = lightDistribution->Lookup(isect.p);
	  		//sample illumination from lights to find path contribution
	  		//skip this for perfectly specular BSDFS
	  		if (isect.bsdf->NumComponents(BxDFType(BSDF_ALL & ~BSDF_SPECULAR)) > 0){
	  			//todo: sample the sdtree and compute the contribution, after sampling the direct light, we record
	  			//it back to the sdtree
	  			Spectrum Ld = beta * UniformSampleOneLight_PG(isect, scene, arena, sampler, false, dwrapper, distrib);
	  			L += Ld;
	  			recordRadianceForVertex(Ld);
	  		}
	  		//indirect lighting computation
	  		//generate the next ray to trace

	  		//todo: sample the sdtree and compute the new pdf
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

	  		//generate a new vertex by this wi and dwrapper
	  		vertex[nVertex++] = RecordVertex(dwrapper, wi);

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

  		for (int i = 0; i < nVertex; i++){
  			vertex[i].commit();
  		}

  		return L;
  	}

	PathGuidingIntegrator *CreatePGPathIntegrator(const ParamSet &params,
												  std::shared_ptr<const Camera> camera)
	{
		int maxDepth = params.FindOneInt("maxdepth", 5);
		int spp = params.FindOneInt("spp", 1023);
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
		return new PathGuidingIntegrator(maxDepth, camera, spp, pixelBounds,
										 rrThreshold, lightStrategy);
	}

}
