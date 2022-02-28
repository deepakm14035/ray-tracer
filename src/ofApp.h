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
enum Shape { SphereShape, PlaneShape, TorusShape, BoxShape };
//  General Purpose Ray class 
//
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
	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal, bool debug) { return false; }
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