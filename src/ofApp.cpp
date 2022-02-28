#include "ofApp.h"

bool isDebug = false;
glm::vec3 camera = glm::vec3(0, 8, 20);
std::vector<glm::vec3> lightSource = { glm::vec3(-8, 18, 0) , glm::vec3(8, 18, 0) };

float viewportWidth = 2.0f;
float viewportHeight = 2.0f;

glm::vec3 cameraMultiplier = { 0.02f, 0.02f, 1.0f };
float cameraToScreenDist = 1.5f;
bool renderImage = false;

int MAX_RAY_STEPS = 100;
float DIST_THRESHOLD = 0.01f;

std::vector<SceneObject*> generateScene() {
	std::vector<SceneObject*> scene;
	Plane* floor	= new Plane	(glm::vec3(0, -5, 0), glm::vec3(0, 1, 0));
	//Sphere* sphere	= new Sphere(glm::vec3(0, 5, 0), 2.0f);
	//Sphere* sphere1 = new Sphere(glm::vec3(0, 5, -10), 4.0f, ofColor::blueViolet);
	//Sphere* sphere2 = new Sphere(glm::vec3(0, 5, -20), 6.0f, ofColor::fireBrick);
	Torus* torus1 = new Torus(glm::vec3(-8, 0, -10), 6.0f, 1.5f, ofColor::fireBrick);
	Torus* torus2 = new Torus(glm::vec3(8, 0, -10), 6.0f, 1.5f, ofColor::fireBrick);
	torus1->setRotation(glm::vec3(0, 0, 30*3.1428f / 180.0f));
	scene.push_back(floor);
	//scene.push_back(sphere);
	//scene.push_back(sphere1);
	//scene.push_back(sphere2);
	scene.push_back(torus1);
	//scene.push_back(torus2);
	return scene;
}

// Intersect Ray with Plane  (wrapper on glm::intersect*
//
bool Plane::intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal, bool debug) {
	float dist;
	bool hit = glm::intersectRayPlane(ray.p, ray.d, position, this->normal, dist);
	normal.x = this->normal.x;
	normal.y = this->normal.y;
	normal.z = this->normal.z;
	//std::cout << "Plane::intersect" << std::endl;
	if (hit) {
		Ray r = ray;
		point = r.evalPoint(dist);
	}

	//custom ray plane intersection
	//float aDotB = glm::dot(ray.d, this->normal);
	//float cosTheta = aDotB / (glm::l2Norm(ray.d) * glm::l2Norm(this->normal));
	//std::cout << "intersecting - " << cosTheta << std::endl;
	//if () 
	return hit;
	//return (glm::abs(cosTheta) > 0.001f);
}


// Convert (u, v) to (x, y, z) 
// We assume u,v is in [0, 1]
//
glm::vec3 ViewPlane::toWorld(float u, float v) {
	float w = width();
	float h = height();
	return (glm::vec3((u * w) + min.x, (v * h) + min.y, position.z));
}

// Get a ray from the current camera position to the (u, v) position on
// the ViewPlane
//
Ray RenderCam::getRay(float u, float v) {
	glm::vec3 pointOnPlane = view.toWorld(u, v);
	return(Ray(position, glm::normalize(pointOnPlane - position)));
}




//--------------------------------------------------------------
void ofApp::setup() {

}

//--------------------------------------------------------------
void ofApp::update() {

}


float lightMultiplier(glm::vec3 position, glm::vec3 lightSource, glm::vec3 normal) {
	//return glm::dot(normal, lightSource - position) / (glm::l2Norm(normal), glm::l2Norm(lightSource - position));
	return  glm::clamp( glm::dot(glm::normalize(normal), glm::normalize(lightSource-position)), 0.0f, 1.0f);// (glm::l2Norm(normal), glm::l2Norm(lightSource-position));
}

float handleLighting(glm::vec3 intersectionPoint, std::vector<SceneObject*> scene, glm::vec3 normalDirection, int objNo) {
	float result = 0.0f;
	glm::vec3 intersectionPoint1;
	for (int i = 0; i < lightSource.size(); i++) {
		float value = lightMultiplier(intersectionPoint, lightSource[i], normalDirection);
		Ray ray1(intersectionPoint, glm::normalize(lightSource[i] - intersectionPoint));
		for (int objNo1 = 0; objNo1 < scene.size(); objNo1++) {
			if (objNo1 != objNo && scene[objNo1]->intersect(ray1, intersectionPoint1, glm::normalize(normalDirection), false)) {
				value = 0.0f;
				break;
			}
		}
		result += value;
	}
	
	return glm::clamp( result, 0.0f, 1.0f);
}

//using raymarching
float handleLighting1(glm::vec3 intersectionPoint, std::vector<SceneObject*> scene, glm::vec3 normalDirection, int objNo) {
	float result = 0.0f;
	
	for (int i = 0; i < lightSource.size(); i++) {
		float value = lightMultiplier(intersectionPoint, lightSource[i], normalDirection);
		Ray ray1(intersectionPoint, glm::normalize(lightSource[i] - intersectionPoint));
		int iterations = 0;
		glm::vec3 intersectionPoint1;
		glm::vec3 nextStep;
		int closestObjectDistance = 99999;
		nextStep = ray1.p + ray1.d;
		while (iterations < MAX_RAY_STEPS && closestObjectDistance>DIST_THRESHOLD) {
			closestObjectDistance = 999999;
			for (int objNo1 = 0; objNo1 < scene.size(); objNo1++) {
				//if (isDebug) std::cout << scene[objNo]->sdf(nextStep, false) << ", " << ray1.d << std::endl;
				if (objNo!=objNo1 && scene[objNo1]->sdf(nextStep, false) < closestObjectDistance) {
					closestObjectDistance = scene[objNo1]->sdf(nextStep, false);
				}
			}
			nextStep = nextStep + ray1.d * closestObjectDistance;
			iterations++;
		}
		if(iterations== MAX_RAY_STEPS)//no other object between this object and light source
			result += value;
	}

	return glm::clamp(result, 0.0f, 1.0f);
}

void rayTracing(std::vector<SceneObject*> scene, Ray ray, glm::vec3 pixelPosition, ofImage& img) {
	glm::vec3 intersectionPoint;
	glm::vec3 normalDirection;
	bool setColor = false;
	glm::vec3 minIntersectionPoint = glm::vec3(99999, 99999, 99999);
	glm::vec3 minNormalDirection = glm::vec3(99999, 99999, 99999);
	int minObjectNo = -1;
	for (int objNo = 0; objNo < scene.size(); objNo++) {
		if (scene[objNo]->intersect(ray, intersectionPoint, normalDirection, false) && glm::l2Norm(camera, intersectionPoint) < 200.0f) {
			if (glm::l2Norm(minIntersectionPoint, camera) > glm::l2Norm(camera, intersectionPoint)) {
				minIntersectionPoint = intersectionPoint;
				minNormalDirection = normalDirection;
				//if (isDebug) myfile << (int)(w + width / 2.0f)<<", "<< height - (int)(h + height / 2.0f) << "\n";
				//std::cout << (int)(w + width / 2.0f) << ", " << height - (int)(h + height / 2.0f) << std::endl;

				minObjectNo = objNo;
				setColor = true;
			}
		}
	}

	if (minObjectNo != -1) {
		scene[minObjectNo]->intersect(ray, intersectionPoint, normalDirection, false);
		float shadedColor = handleLighting(intersectionPoint, scene, normalDirection, minObjectNo);
		ofColor c = scene[minObjectNo]->diffuseColor;
		c.r *= shadedColor;
		c.g *= shadedColor;
		c.b *= shadedColor;

		img.setColor(pixelPosition.x, pixelPosition.y, c);
	}
}

bool march(std::vector<SceneObject*> scene, Ray ray, int* minObjectNo, glm::vec3* nextStep) {
	*nextStep = ray.p + ray.d;
	float closestObjectDistance = 999999;
	*(minObjectNo) = -1;
	int iterations = 0;
	while (iterations < MAX_RAY_STEPS && closestObjectDistance>DIST_THRESHOLD) {
		closestObjectDistance = 999999;
		for (int objNo = 0; objNo < scene.size(); objNo++) {
			if (isDebug) std::cout << scene[objNo]->sdf(*nextStep, false) << ", " << ray.d << std::endl;
			if (scene[objNo]->sdf(*nextStep, false) < closestObjectDistance) {
				closestObjectDistance = scene[objNo]->sdf(*nextStep, false);
				*minObjectNo = objNo;
			}
		}
		*nextStep = *nextStep + ray.d * closestObjectDistance;
		iterations++;
	}
	if (isDebug) std::cout << "iterations-" << iterations << ", closestPoint-" << nextStep << std::endl;
	return iterations < MAX_RAY_STEPS;
}


void rayMarching(std::vector<SceneObject*> scene, Ray ray, glm::vec3 pixelPosition, ofImage& img, bool isDebug) {
	int* minObjectNo=(int*)calloc(1,sizeof(int));
	glm::vec3* nextStep = new glm::vec3();
	if (march(scene, ray, minObjectNo, nextStep)) {
		//glm::vec3 intersectionPoint;
		//glm::vec3 normalDirection;
		//scene[minObjectNo]->intersect(ray, intersectionPoint, normalDirection, false);
		float shadedColor = handleLighting1(*nextStep, scene, scene[*minObjectNo]->getNormal(*nextStep, false), *minObjectNo);
		ofColor c = scene[*minObjectNo]->diffuseColor;
		c.r *= shadedColor;
		c.g *= shadedColor;
		c.b *= shadedColor;

		img.setColor(pixelPosition.x, pixelPosition.y, c);
	}
}


void drawScene() {
	//ofDrawRectangle(64, 64, 64, 64);
	int width = ofGetWindowWidth() - 2;
	int height = ofGetWindowHeight() - 2;


	std::vector<SceneObject*> scene = generateScene();
	glm::vec3 screenOffset = glm::vec3();
	ofImage img;
	img.allocate(ofGetWindowWidth(), ofGetWindowHeight(), OF_IMAGE_COLOR);
	//std::cout << ofGetWindowWidth() << ", " << ofGetWindowHeight() << "\n";
	//ofPixels& pixels = img.getPixels();
	img.setColor(ofColor::black);
	
	ofstream myfile;
	if (isDebug)	myfile.open("C:\\Users\\Deepak\\OneDrive\\Desktop\\a.txt");

	//myfile << "Writing this to a file.\n";
	for (float h = -height / 2.0f; h < height / 2.0f; h++) {
		for (float w = -width / 2.0f; w < width / 2.0f; w++) {
			float u = w / width;
			float v = h / height;
			glm::vec3 pixelPos = glm::vec3(u * viewportWidth, v * viewportHeight, -cameraToScreenDist) + camera;
			glm::vec3 rayDirection = pixelPos - camera;
			Ray ray(camera, glm::normalize(rayDirection));
			//if (isDebug) myfile << "[" << h << " " << w << "]ray start- " << camera<<" direction- "<<rayDirection << std::endl;// "\t" << rayDirection << std::endl;
			//rayTracing(scene, ray, glm::vec3((int)(w + width / 2.0f), height - (int)(h + height / 2.0f),0), img);
			rayMarching(scene, ray, glm::vec3((int)(w + width / 2.0f), height - (int)(h + height / 2.0f), 0), img, false);
		}
		std::cout << h <<" ";
	}
	if (isDebug) myfile.close();
	//ofSetColor(255);
	img.update();
	img.draw(0, 0);
}

void drawbox() {
	ofImage img;
	img.allocate(ofGetWidth(), ofGetHeight(), OF_IMAGE_COLOR_ALPHA);
	//img.getPixels().setColor(x, y, ofColor(r, g, b, a));//place your pixels
	img.update();
	for (int y = 0; y < 200; y++) {
		for (int x = 0; x < 200; x++) {
			img.getPixels().setColor(x, y, ofColor::yellow);
		}
	}
	ofSetColor(255);// if you dont set the color to white the image will be tinted by whichever color was set before.
	img.draw(0, 0);
}

//--------------------------------------------------------------
void ofApp::draw() {
	if (renderImage) {
		drawScene();
		ofSetBackgroundAuto(false);
		renderImage = false;
	}
	//drawbox();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button) {
	renderImage = true;
	//drawScene();
	int width = ofGetWindowWidth() - 2;
	int height = ofGetWindowHeight() - 2;

	float w = x - width / 2.0f;
	float h = height / 2.0f - y;
	float u = w / width;
	float v = h / height;
	ofImage img;
	glm::vec3 pixelPos = glm::vec3(u * viewportWidth, v * viewportHeight, -cameraToScreenDist) + camera;
	glm::vec3 rayDirection = pixelPos - camera;
	Ray ray(camera, glm::normalize(rayDirection));
	std::cout << rayDirection<<"\n";
	std::vector<SceneObject*> scene = generateScene();
	//rayMarching(scene, ray, glm::vec3((int)(w + width / 2.0f), height - (int)(h + height / 2.0f), 0), img, true);

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h) {

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg) {

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo) {

}