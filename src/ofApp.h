//
//  RayCaster - Set of simple classes to create a camera/view setup for our Ray Tracer HW Project
//
//  I've included these classes as a mini-framework for our introductory ray tracer.
//  You are free to modify/change.   
//
//  These classes provide a simple render camera which can can return a ray starting from
//  it's position to a (u, v) coordinate on the view plane.
//
//  The view plane is where we can locate our photorealistic image we are rendering.
//  The field-of-view of the camera by moving it closer/further 
//  from the view plane.  The viewplane can be also resized.  When ray tracing an image, the aspect
//  ratio of the view plane should the be same as your image. So for example, the current view plane
//  default size is ( 6.0 width by 4.0 height ).   A 1200x800 pixel image would have the same
//  aspect ratio.
//
//  This is not a complete ray tracer - just a set of skelton classes to start.  The current
//  base scene object only stores a value for the diffuse/specular color of the object (defaut is gray).
//  at some point, we will want to replace this with a Material class that contains these (and other 
//  parameters)
//  
//  (c) Kevin M. Smith  - 24 September 2018
//
#pragma once

#include "ofMain.h"
#include <glm/gtx/intersect.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include <fstream>

enum Shape { SphereShape, PlaneShape, TorusShape, BoxShape, InfiniteRepShape, MeshShape };
//  General Purpose Ray class 
//
class Face {
public:
	glm::vec3 vert_ind;
	glm::vec3 tex;
	glm::vec3 normal;
	Face(glm::vec3 v, glm::vec3 t, glm::vec3 n) {
		vert_ind = v;
		tex = t;
		normal = n;
	}
	Face() {}
};

class Node {
public:
	glm::vec3 bbStart;
	glm::vec3 dim;
	std::vector<Node> children;
	std::vector<int> triangles;
	Node() {}
	Node(glm::vec3 bbStart, glm::vec3 dim) {
		this->bbStart = bbStart;
		this->dim = dim;
	}
	bool contains(glm::vec3 point) {
		return point.x >= bbStart.x && point.y >= bbStart.y && point.z >= bbStart.z && 
				point.x <= bbStart.x+dim.x && point.y <= bbStart.y + dim.y && point.z <= bbStart.z + dim.z;
	}
};

static std::vector<glm::vec3> bbParts = { glm::vec3(0.0f,0.0f,0.0f),glm::vec3(0.5f,0.0f,0.0f), glm::vec3(0.0f,0.5f,0.0f), glm::vec3(0.0f,0.0f,0.5f), glm::vec3(0.5f,0.5f,0.0f), glm::vec3(0.0f,0.5f,0.5f), glm::vec3(0.5f,0.0f,0.5f), glm::vec3(0.5f,0.5f,0.5f)};

class Ray {
public:
	Ray(glm::vec3 p, glm::vec3 d) { this->p = p; this->d = d; }
	void draw(float t) { ofDrawLine(p, p + t * d); }

	glm::vec3 evalPoint(float t) {
		return (p + t * d);
	}

	glm::vec3 p, d;
};

//  Base class for any renderable object in the scene
//
class SceneObject {
public:
	virtual void draw() {}    // pure virtual funcs - must be overloaded
	virtual bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal, bool debug) { //cout << "SceneObject::intersect" << endl; 
		return false;
	}
	virtual float sdf(glm::vec3& point, bool debug) { //cout << "SceneObject::intersect" << endl; 
		return 0;
	}
	virtual glm::vec3 getNormal(glm::vec3& point, bool debug) { //cout << "SceneObject::intersect" << endl; 
		return glm::vec3();
	}
	// any data common to all scene objects goes here
	glm::vec3 position = glm::vec3(0, 0, 0);
	Shape shape;
	// material properties (we will ultimately replace this with a Material class - TBD)
	//
	ofColor diffuseColor = ofColor::grey;    // default colors - can be changed.
	ofColor specularColor = ofColor::lightGray;
	glm::vec3 rotation;
	glm::mat4 rotationMatrix;
	glm::mat4 inverseRotationMatrix;
};

//  General purpose sphere  (assume parametric)
//
class Sphere : public SceneObject {
public:
	Sphere(glm::vec3 p, float r, ofColor diffuse = ofColor::lightGray) { position = p; radius = r; diffuseColor = diffuse; shape = SphereShape; }
	Sphere() {}
	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal, bool debug) {
		bool result = glm::intersectRaySphere(ray.p, ray.d, position, radius, point, normal);
		if (debug)
			std::cout << "[Sphere::intersect] "<<ray.p<<" "<<ray.d<<" | sphere center- "<< position<<" radius- "<<radius<<", intersection at- "<<point<<", normal- "<< normal<<", result- "<< result << std::endl;
		return result;
	}
	float sdf(glm::vec3& point, bool debug) {
		return glm::length(point-position)-radius;
	}
	glm::vec3 getNormal(glm::vec3& point, bool debug) { //cout << "SceneObject::intersect" << endl; 
		return point-position;
	}
	void draw() {
		ofDrawSphere(position, radius);
	}

	float radius = 1.0;
};

//  Mesh class (will complete later- this will be a refinement of Mesh from Project 1)
//
class Mesh : public SceneObject {
public:
	std::vector<glm::vec3> v;
	std::vector<glm::vec3> vt;
	std::vector<glm::vec3> vn;
	std::vector<Face> f;
	glm::vec3 scale;
	glm::vec3 position;
	std::vector<std::vector<int>> v_normals_list;
	std::vector<glm::vec3> new_vn;
	Node root;
	int depth = 0;
	int totalDepth = 0;
	int totalSize = 0;

	Mesh(std::string filename, glm::vec3 p, glm::vec3 size, ofColor diffuse = ofColor::lightGray, float w = 20, float h = 20) {
		diffuseColor = diffuse;
		shape = MeshShape;
		scale = size;
		position = p;
		std::ifstream infile(filename);
		std::string line;
		glm::vec3 sum(0.0f, 0.0f, 0.0f);
		root.bbStart = glm::vec3(0.0f, 0.0f, 0.0f);
		root.dim = glm::vec3(0.0f, 0.0f, 0.0f);
		int triangleNo = 0;
		while (std::getline(infile, line))
		{
			std::istringstream iss(line);
			std::string lineType;
			iss >> lineType;
			if (lineType.compare("v")==0) {
				glm::vec3 vec;
				iss >> vec.x;
				iss >> vec.y;
				iss >> vec.z;
				sum = sum + vec;
				v.push_back(glm::diagonal3x3(scale) * vec + position);
				root.bbStart = glm::min(root.bbStart, v[v.size() - 1]);
				root.dim = glm::max(root.dim, v[v.size() - 1]-root.bbStart);
				std::vector<int> new_vec;
				v_normals_list.push_back(new_vec);
			}
			else if (lineType.compare("vt") == 0) {
				glm::vec3 tex;//texture coord
				iss >> tex.x;
				iss >> tex.y;
				iss >> tex.z;
				vt.push_back(tex);
			}
			else if (lineType.compare("vn") == 0) {
				glm::vec3 vecNormal;
				iss >> vecNormal.x;
				iss >> vecNormal.y;
				iss >> vecNormal.z;
				//std::cout << vecNormal << endl;
				vn.push_back(vecNormal);
			}
			else if (lineType.compare("f") == 0) {
				std::replace(line.begin(), line.end(), '/', ' ');
				//std::cout << line << endl;
				//vertex_index/texture_index/normal_index 
				std::istringstream iss1(line);
				iss1 >> lineType;
				glm::vec3 vert;
				glm::vec3 tex;
				glm::vec3 normal;
				int t;
				iss1 >> t;
				vert.x = t;
				iss1 >> tex.x;
				iss1 >> normal.x;
				iss1 >> t;
				vert.y = t;
				iss1 >> tex.y;
				iss1 >> normal.y;
				iss1 >> vert.z;
				iss1 >> tex.z;
				iss1 >> normal.z;
				vert -= 1;
				normal -= 1;
				tex -= 1;
				v_normals_list[vert.x].push_back(normal.x);
				v_normals_list[vert.y].push_back(normal.x);
				v_normals_list[vert.z].push_back(normal.x);
				root.triangles.push_back(triangleNo);
				triangleNo++;
				Face face(vert, tex, normal);
				f.push_back(face);
			}
		}

		for (int i = 0; i < v.size(); i++) {
			glm::vec3 sum1(0.0f,0.0f,0.0f);
			for (int j = 0; j < v_normals_list[i].size(); j++) {
				sum1 += vn[v_normals_list[i][j]];
			}
			sum1 /= v_normals_list[i].size();
			new_vn.push_back(sum1);
		}
		int level = 0;
		generateOctree(&root, level);

		infile.close();
		std::cout << "\nmesh size:\n\tvertices: "<<v.size()<<"\n\ttriangles: "<<f.size() << "\n\tvertex normals: " << vn.size() << "\n\tvert avg: "<<sum<<"\n";
	}

	void generateOctree(Node* temp, int level) {
		//cout << "\nstart: " << temp->bbStart << ", dim: " << temp->dim << "\n"; 
		//split into 8
		std::vector<Node> children;
		glm::vec3 center = temp->bbStart + temp->dim/2.0f;
		for (int i = 0; i < bbParts.size(); i++) {
			Node newNode(temp->bbStart + temp->dim*bbParts[i]-glm::vec3(0.1f,0.1f,0.1f), temp->dim / 2 + glm::vec3(0.1f, 0.1f, 0.1f));
			for (int t = 0; t < temp->triangles.size(); t++) {
				int triIndex = temp->triangles[t];
				if (newNode.contains(v[f[triIndex].vert_ind.r]) || newNode.contains(v[f[triIndex].vert_ind.g]) || newNode.contains(v[f[triIndex].vert_ind.b])) {
					newNode.triangles.push_back(triIndex);
				}
			}
			if (newNode.triangles.size() > 0) {
				std::cout << "Node level " << level << ":\n\tstart: " << newNode.bbStart << "\n\tend: " << newNode.bbStart + newNode.dim << "\n\n";
				std::cout << "\tcount: " << newNode.triangles.size() << "\n";
				children.push_back(newNode);
			}
		}
		//recursively check and call generateOctree
		if (children.size() == 0) return;
		level++;
		temp->children.insert(temp->children.end(), children.begin(), children.end());
		for (int i = 0; i < temp->children.size(); i++) {
			if(temp->children[i].triangles.size()>50)
				generateOctree(&(temp->children[i]), level);
		}

	}

	void setRotation(glm::vec3 rot) {
		rotation = rot;
		inverseRotationMatrix = glm::rotate(glm::mat4(1.0f), -rotation.z, glm::vec3(0, 0, 1));
		inverseRotationMatrix *= glm::rotate(inverseRotationMatrix, -rotation.y, glm::vec3(0, 1, 0));
		inverseRotationMatrix *= glm::rotate(inverseRotationMatrix, -rotation.x, glm::vec3(1, 0, 0));
		rotationMatrix = glm::rotate(glm::mat4(1.0f), rotation.x, glm::vec3(1, 0, 0));
		rotationMatrix *= glm::rotate(rotationMatrix, rotation.y, glm::vec3(0, 1, 0));
		rotationMatrix *= glm::rotate(rotationMatrix, rotation.z, glm::vec3(0, 0, 1));
	}

	bool intersectRayBox(Node node, const Ray& ray) {
		//https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
		bool intersects = false;
		std::vector<float> t_values;
		//std::cout << "[intersectRayBox] ray dir: " << ray.d << "\n";
		float t0x = (node.bbStart.x - ray.p.x) / ray.d.x;
		if(t0x>0){
			intersects |= node.contains(ray.p + (t0x+0.01f) * ray.d);
		}
		float t0y = (node.bbStart.y - ray.p.y) / ray.d.y;
		if (t0y > 0) {
			intersects |= node.contains(ray.p + (t0y + 0.01f) * ray.d);
		}
		float t0z = (node.bbStart.z - ray.p.z) / ray.d.z;
		if (t0z > 0) {
			intersects |= node.contains(ray.p + (t0z + 0.01f) * ray.d);
		}
		float t1x = (node.bbStart.x + node.dim.x - ray.p.x) / ray.d.x;
		if (t1x > 0) {
			intersects |= node.contains(ray.p + (t1x + 0.01f) * ray.d);
		}
		float t1y = (node.bbStart.y + node.dim.y - ray.p.y) / ray.d.y;
		if (t1y > 0) {
			intersects |= node.contains(ray.p + (t1y + 0.01f) * ray.d);
		}
		float t1z = (node.bbStart.z + node.dim.z - ray.p.z) / ray.d.z;
		if (t1z > 0) {
			intersects |= node.contains(ray.p + (t1z + 0.01f) * ray.d);
		}
		//std::cout << "[intersectRayBox] end\n";
		return intersects;
	}

	std::vector<Node*>  parseTree(const Ray& ray, glm::vec3& point, glm::vec3& normal, Node* node) {
		std::vector<Node*> intersectingBoxes;
		depth++;
		//std::cout << "\n[mesh.parseTree]"<<ray.d<<"\n";
		for (int i = 0; i < node->children.size(); i++) {
			//std::cout << "\n[mesh.parseTree] child: " << i << "\n";
			if (intersectRayBox(node->children[i], ray)) {
				if (node->children[i].children.size() == 0)
					intersectingBoxes.push_back(&(node->children[i]));
				else {
					std::vector<Node*> temp = parseTree(ray, point, normal, &(node->children[i]));
					intersectingBoxes.insert(intersectingBoxes.end(), temp.begin(), temp.end());
				}
			}
		}
		return intersectingBoxes;
	}

	Mesh() { }
	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal, bool debug) {
		//https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution
		int closestTriangle = -1;
		float minDist = 99999.0f;
		std::vector<Node*> list;
		if (debug) {
			std::cout << "\n[mesh.intersect]\n";
		}
		//if (!intersectRayBox(root, ray)) {
		//	return false;
		//}
		depth = 0;
		//list = parseTree(ray, point, normal, &root);
		//if (list.size()==0) return false;
		//std::cout << " [" << list.size() << ", "<<depth<<"] ";
//		totalDepth += depth;
	//	for (int b = 0; b < list.size(); b++) {
	//		totalSize += list[b]->triangles.size();
	//		for (int i = 0; i < list[b]->triangles.size(); i++) {
			for (int i = 0; i < f.size(); i++) {
			//	int triIndex = list[b]->triangles[i];
				int triIndex = i;
				glm::vec2 baryPos(0.0f, 0.0f);
				float dist;
				//std::cout << f[triIndex].vert_ind << std::endl;
				if (glm::intersectRayTriangle(ray.p, ray.d, v[f[triIndex].vert_ind.x], v[f[triIndex].vert_ind.y], v[f[triIndex].vert_ind.z], baryPos, dist)) {
					if (baryPos.x <= 1 && baryPos.y <= 1 && dist < minDist) {
						normal = vn[f[triIndex].normal.x];
						point = ray.p + dist * ray.d;
						minDist = dist;
						closestTriangle = triIndex;
					}
				}

//			}
		}
		if (false && closestTriangle >= 0) {
			/*std::vector<int> neighbors;
			for (int i = 0; i < f.size(); i++) {
				bool matchingEdgeX = (f[i].vert_ind.x == f[closestTriangle].vert_ind.x || f[i].vert_ind.x == f[closestTriangle].vert_ind.y || f[i].vert_ind.x == f[closestTriangle].vert_ind.z);
				bool matchingEdgeY = (f[i].vert_ind.y == f[closestTriangle].vert_ind.x || f[i].vert_ind.y == f[closestTriangle].vert_ind.y || f[i].vert_ind.y == f[closestTriangle].vert_ind.z);
				bool matchingEdgeZ = (f[i].vert_ind.z == f[closestTriangle].vert_ind.x || f[i].vert_ind.z == f[closestTriangle].vert_ind.y || f[i].vert_ind.z == f[closestTriangle].vert_ind.z);
				if (i != closestTriangle && ((matchingEdgeX&& matchingEdgeY)|| (matchingEdgeY && matchingEdgeZ)|| (matchingEdgeX && matchingEdgeZ))) {
					neighbors.push_back(i);

				}
			}
			glm::vec3 closestAB = glm::closestPointOnLine(point, v[f[closestTriangle].vert_ind.r], v[f[closestTriangle].vert_ind.g]- v[f[closestTriangle].vert_ind.r]);
			glm::vec3 closestBC = glm::closestPointOnLine(point, v[f[closestTriangle].vert_ind.g], v[f[closestTriangle].vert_ind.b] - v[f[closestTriangle].vert_ind.g]);
			glm::vec3 closestAC = glm::closestPointOnLine(point, v[f[closestTriangle].vert_ind.r], v[f[closestTriangle].vert_ind.b] - v[f[closestTriangle].vert_ind.r]);
			normal = normal + glm::distance(point, v[f[closestTriangle].vert_ind.b]) / glm::distance(closestAB, point) * vn[f[neighbors[0]].normal.x]
				+ glm::distance(point, v[f[closestTriangle].vert_ind.b]) / glm::distance(closestBC, point) * vn[f[neighbors[0]].normal.x]
				+ glm::distance(point, v[f[closestTriangle].vert_ind.b]) / glm::distance(closestAC, point) * vn[f[neighbors[0]].normal.x];*/
			//Approach from https://www.flipcode.com/archives/Interpolating_Normals_For_Ray-Tracing.shtml
			int p1, p2, p3;
			if (glm::distance(point, v[f[closestTriangle].vert_ind.x]) > glm::distance(point, v[f[closestTriangle].vert_ind.y]) &&
				glm::distance(point, v[f[closestTriangle].vert_ind.x]) > glm::distance(point, v[f[closestTriangle].vert_ind.z])) {
				p1=f[closestTriangle].vert_ind.x;
				p2=f[closestTriangle].vert_ind.y;
				p3= f[closestTriangle].vert_ind.z;
			}else if (glm::distance(point, v[f[closestTriangle].vert_ind.y]) > glm::distance(point, v[f[closestTriangle].vert_ind.x]) &&
					  glm::distance(point, v[f[closestTriangle].vert_ind.y]) > glm::distance(point, v[f[closestTriangle].vert_ind.z])) {
				p1 = f[closestTriangle].vert_ind.y;
				p2 = f[closestTriangle].vert_ind.x;
				p3 = f[closestTriangle].vert_ind.z;
			}
			else {
				p1 = f[closestTriangle].vert_ind.z;
				p2 = f[closestTriangle].vert_ind.x;
				p3 = f[closestTriangle].vert_ind.y;
			}
			if (debug) {
				std::cout << "\np1: " << v[p1] << "\np2: " << v[p2] << "\np3: " << v[p3] << "\n";
				std::cout << "\nn1: " << new_vn[p1] << "\nn2: " << new_vn[p2] << "\nn3: " << new_vn[p3] << "\n";
			}
			float pointToEdgeDist;
			glm::intersectRayPlane(v[p1], point - v[p1], v[p3], glm::cross(v[p3] - v[p2], vn[f[closestTriangle].normal.x]), pointToEdgeDist);
			glm::vec3 intersectionPoint = v[p1] + pointToEdgeDist * (point - v[p1]);

			float dist2q = glm::distance(v[p2], intersectionPoint);
			float dist23 = glm::distance(v[p2], v[p3]);
			normal = new_vn[p2] + (new_vn[p3] - new_vn[p2]) * dist2q / dist23;
			if (debug) {
				std::cout << "\ndist2q: " << dist2q << "\n";
				std::cout << "dist23: " << dist23 << "\n";
				std::cout << "3->p intersectionPoint: " << intersectionPoint << "\n";
				std::cout << "\nq normal: " << normal << "\n";
			}
			float distpq = glm::distance(point, intersectionPoint);
			float dist1q = glm::distance(v[p1], intersectionPoint);
			normal = normal + (new_vn[p1] - normal) * distpq / dist1q;
			if (debug) {
				std::cout << "\nq normal: " << normal << "\n";
			}
			//normal=vn[f[closestTriangle].normal.x];
		}
		return closestTriangle!=-1;
	}
	void draw() { }
};


//  General purpose plane 
//
class Plane : public SceneObject {
public:
	Plane(glm::vec3 p, glm::vec3 n, ofColor diffuse = ofColor::darkOliveGreen, float w = 20, float h = 20) {
		position = p; normal = n;
		width = w;
		height = h;
		diffuseColor = diffuse;
		plane.rotateDeg(90, 1, 0, 0);
		shape = PlaneShape;
		//std::cout << "plane init\n";
	}
	Plane() { }
	glm::vec3 normal = glm::vec3(0, 1, 0);
	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal, bool debug);
	
	float sdf(glm::vec3& point, bool debug) {
		return point.y - position.y;
	}
	
	glm::vec3 getNormal(glm::vec3& point, bool debug) { //cout << "SceneObject::intersect" << endl; 
		return normal;
	}
	void draw() {
		plane.setPosition(position);
		plane.setWidth(width);
		plane.setHeight(height);
		plane.setResolution(4, 4);
		plane.drawWireframe();
	}
	ofPlanePrimitive plane;
	float width = 20;
	float height = 20;
};

// view plane for render camera
// 
class  ViewPlane : public Plane {
public:
	ViewPlane(glm::vec2 p0, glm::vec2 p1) { min = p0; max = p1; }

	ViewPlane() {                         // create reasonable defaults (6x4 aspect)
		min = glm::vec2(-3, -2);
		max = glm::vec2(3, 2);
		position = glm::vec3(0, 0, 5);
		normal = glm::vec3(0, 0, 1);      // viewplane currently limited to Z axis orientation
	}

	void setSize(glm::vec2 min, glm::vec2 max) { this->min = min; this->max = max; }
	float getAspect() { return width() / height(); }

	glm::vec3 toWorld(float u, float v);   //   (u, v) --> (x, y, z) [ world space ]

	void draw() {
		ofDrawRectangle(glm::vec3(min.x, min.y, position.z), width(), height());
	}


	float width() {
		return (max.x - min.x);
	}
	float height() {
		return (max.y - min.y);
	}

	// some convenience methods for returning the corners
	//
	glm::vec2 topLeft() { return glm::vec2(min.x, max.y); }
	glm::vec2 topRight() { return max; }
	glm::vec2 bottomLeft() { return min; }
	glm::vec2 bottomRight() { return glm::vec2(max.x, min.y); }

	//  To define an infinite plane, we just need a point and normal.
	//  The ViewPlane is a finite plane so we need to define the boundaries.
	//  We will define this in terms of min, max  in 2D.  
	//  (in local 2D space of the plane)
	//  ultimately, will want to locate the ViewPlane with RenderCam anywhere
	//  in the scene, so it is easier to define the View rectangle in a local'
	//  coordinate system.
	//
	glm::vec2 min, max;
};

class Torus : public SceneObject {
public:
	Torus(glm::vec3 p, float r, float r2, ofColor diffuse = ofColor::lightGray) { position = p; radius = r; radius2 = r2; diffuseColor = diffuse; shape = TorusShape; }
	Torus() {}

	void setRotation(glm::vec3 rot) {
		rotation = rot;
		inverseRotationMatrix =		glm::rotate(glm::mat4(1.0f),		-rotation.z,	glm::vec3(0, 0, 1));
		inverseRotationMatrix *=	glm::rotate(inverseRotationMatrix,	-rotation.y,	glm::vec3(0, 1, 0));
		inverseRotationMatrix *=	glm::rotate(inverseRotationMatrix,	-rotation.x,	glm::vec3(1, 0, 0));
		rotationMatrix = glm::rotate(glm::mat4(1.0f), rotation.x, glm::vec3(1, 0, 0));
		rotationMatrix *= glm::rotate(rotationMatrix, rotation.y, glm::vec3(0, 1, 0));
		rotationMatrix *= glm::rotate(rotationMatrix, rotation.z, glm::vec3(0, 0, 1));
	}

	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal, bool debug) {
		//glm::vec2 q = glm::vec2(glm::length(p.xz) - t.x, p.y);
		return false;// length(q) - t.y;
	}
	
	float sdf(glm::vec3& point, bool debug) {
		//do inverse transform on point
		glm::vec3 temp(point);
		temp -= position;
		temp = inverseRotationMatrix * glm::vec4(temp,1.0f);
		//point = rot * point;
		//find distance between 'point' and radius, then check if it is within radius1 distance
		glm::vec2 q = glm::vec2(glm::length(glm::vec2(temp.x, temp.z)) - radius, temp.y);

		//point = rotationMatrix * glm::vec4(point, 1.0f);
		//point = rot * point;
		//point += position;
		return glm::length(q) - radius2;
	}
	
	glm::vec3 getNormal(glm::vec3& point, bool debug) { //cout << "SceneObject::intersect" << endl; 
		return point - position;
	}
	
	void draw() {
		//ofDrawTorus(position, radius);
	}

	float radius = 2.0;
	float radius2 = 1.0;
};

class Box : public SceneObject {
public:
	Box(glm::vec3 p, glm::vec3 size, ofColor diffuse = ofColor::lightGray) { position = p; s=size; diffuseColor = diffuse; shape = BoxShape; }
	Box() {}

	void setRotation(glm::vec3 rot) {
		rotation = rot;
		inverseRotationMatrix = glm::rotate(glm::mat4(1.0f), -rotation.z, glm::vec3(0, 0, 1));
		inverseRotationMatrix *= glm::rotate(inverseRotationMatrix, -rotation.y, glm::vec3(0, 1, 0));
		inverseRotationMatrix *= glm::rotate(inverseRotationMatrix, -rotation.x, glm::vec3(1, 0, 0));
		rotationMatrix = glm::rotate(glm::mat4(1.0f), rotation.x, glm::vec3(1, 0, 0));
		rotationMatrix *= glm::rotate(rotationMatrix, rotation.y, glm::vec3(0, 1, 0));
		rotationMatrix *= glm::rotate(rotationMatrix, rotation.z, glm::vec3(0, 0, 1));
	}

	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal, bool debug) {
		//glm::vec2 q = glm::vec2(glm::length(p.xz) - t.x, p.y);
		return false;// length(q) - t.y;
	}

	float sdf(glm::vec3& point, bool debug) {
		//do inverse transform on point
		glm::vec3 temp(point);
		temp -= position;
		temp = inverseRotationMatrix * glm::vec4(temp, 1.0f);
		//point = rot * point;
		//find distance between 'point' and radius, then check if it is within radius1 distance
		glm::vec3 q = abs(temp) - s;
		return glm::length(glm::max(q, 0.0f)) + glm::min(glm::max(q.x, max(q.y, q.z)), 0.0f);
	}

	glm::vec3 getNormal(glm::vec3& point, bool debug) { //cout << "SceneObject::intersect" << endl; 
		return point - position;
	}

	void draw() {
		//ofDrawTorus(position, radius);
	}

	glm::vec3 s;
};


class BooleanShape : public SceneObject {
public:
	BooleanShape(SceneObject* object1, Shape a, SceneObject* object2, Shape b, int op, ofColor diffuse = ofColor::lightGray) { A = object1; aShape = a; B = object2; bShape = b; operation = op; diffuseColor = diffuse; shape = BoxShape; }
	BooleanShape() {}

	void setRotation(glm::vec3 rot) {
		rotation = rot;
		inverseRotationMatrix = glm::rotate(glm::mat4(1.0f), -rotation.z, glm::vec3(0, 0, 1));
		inverseRotationMatrix *= glm::rotate(inverseRotationMatrix, -rotation.y, glm::vec3(0, 1, 0));
		inverseRotationMatrix *= glm::rotate(inverseRotationMatrix, -rotation.x, glm::vec3(1, 0, 0));
		rotationMatrix = glm::rotate(glm::mat4(1.0f), rotation.x, glm::vec3(1, 0, 0));
		rotationMatrix *= glm::rotate(rotationMatrix, rotation.y, glm::vec3(0, 1, 0));
		rotationMatrix *= glm::rotate(rotationMatrix, rotation.z, glm::vec3(0, 0, 1));
	}

	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal, bool debug) {
		//glm::vec2 q = glm::vec2(glm::length(p.xz) - t.x, p.y);
		return false;// length(q) - t.y;
	}

	float sdf(glm::vec3& point, bool debug) {
		//do inverse transform on point
		float sdf1 = ((Box*)A)->sdf(point, debug);
		float sdf2 = ((Sphere*)B)->sdf(point, debug);
		//if (operation == 0)//difference
			return glm::max(sdf1, sdf2);
		//return glm::length(glm::max(q, 0.0f)) + glm::min(glm::max(q.x, max(q.y, q.z)), 0.0f);
	}

	glm::vec3 getNormal(glm::vec3& point, bool debug) { //cout << "SceneObject::intersect" << endl; 
		return point - position;
	}

	void draw() {
		//ofDrawTorus(position, radius);
	}

	glm::vec3 s;
	SceneObject* A;
	Shape aShape;
	SceneObject* B;
	Shape bShape;
	int operation;
};

class InfiniteRep : public SceneObject {
public:
	InfiniteRep(SceneObject* objectToRepeat, glm::vec3 gap, ofColor diffuse = ofColor::darkOliveGreen) {
		diffuseColor = diffuse;
		shape = InfiniteRepShape;
		objToRepeat = objectToRepeat;
		c = gap;
		//std::cout << "plane init\n";
	}
	InfiniteRep() { }
	glm::vec3 normal = glm::vec3(0, 1, 0);
	//bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal, bool debug);

	float sdf(glm::vec3& point, bool debug) {
		glm::vec3 q = glm::mod(point + 0.5 * c, c) - 0.5 * c;
		return objToRepeat->sdf(q, false);
	}

	glm::vec3 getNormal(glm::vec3& point, bool debug) { //cout << "SceneObject::intersect" << endl; 
		return normal;
	}
	void draw() {
	}
	glm::vec3 c;
	SceneObject* objToRepeat;
};


//  render camera  - currently must be z axis aligned (we will improve this in project 4)
//
class RenderCam : public SceneObject {
public:
	RenderCam() {
		position = glm::vec3(0, 0, 10);
		aim = glm::vec3(0, 0, -1);
	}
	Ray getRay(float u, float v);
	void draw() { ofDrawBox(position, 1.0); };
	void drawFrustum();

	glm::vec3 aim;
	ViewPlane view;          // The camera viewplane, this is the view that we will render 
};



class ofApp : public ofBaseApp {

public:
	void setup();
	void update();
	void draw();

	void keyPressed(int key);
	void keyReleased(int key);
	void mouseMoved(int x, int y);
	void mouseDragged(int x, int y, int button);
	void mousePressed(int x, int y, int button);
	void mouseReleased(int x, int y, int button);
	void mouseEntered(int x, int y);
	void mouseExited(int x, int y);
	void windowResized(int w, int h);
	void dragEvent(ofDragInfo dragInfo);
	void gotMessage(ofMessage msg);
	void rayTrace();
	void drawGrid();
	void drawAxis(glm::vec3 position);




	bool bHide = true;
	bool bShowImage = false;

	ofEasyCam  mainCam;
	ofCamera sideCam;
	ofCamera previewCam;
	ofCamera* theCam;    // set to current camera either mainCam or sideCam

	// set up one render camera to render image throughn
	//
	RenderCam renderCam;
	ofImage image;

	vector<SceneObject*> scene;

	int imageWidth = 600;
	int imageHeight = 400;
};


void renderRow(float h, int width, int height, std::vector<SceneObject*> scene, ofImage* img);
void renderRow(float h, float w, int width, int height, SceneObject** scene, int size, glm::vec3** img);