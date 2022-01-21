#include "ofApp.h"

bool isDebug = true;
glm::vec3 camera = glm::vec3(0, 5, 5);
glm::vec3 screenCenter = glm::vec3(0, 5, 10);
glm::vec3 lightSource = glm::vec3(5, 15, 10);
float cameraMultiplier = 0.02f;
float cameraToScreenDist = 0.8f;


// Intersect Ray with Plane  (wrapper on glm::intersect*
//
bool Plane::intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal, bool debug) {
	float dist;
	bool hit = glm::intersectRayPlane(ray.p, ray.d, position, this->normal, dist);
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

std::vector<SceneObject*> generateScene() {
	std::vector<SceneObject*> scene;
	//Plane* floor = new Plane(glm::vec3(0, -5, 0), glm::vec3(0, 1, 0));
	Plane* floor = new Plane(glm::vec3(0, 0, 0), glm::vec3(0, 0, 1));
	Sphere* sphere = new Sphere(glm::vec3(0, 5, 0), 3.0f);
	//Ray ray(glm::vec3(4, 3, 2), glm::vec3(4, 3, 2));
	//floor.intersect(ray, glm::vec3(4, 3, 2), glm::vec3(4, 3, 2));
	//scene.push_back(floor);
	scene.push_back(sphere);
	return scene;
}

float lightMultiplier(glm::vec3 position, glm::vec3 lightSource, glm::vec3 normal) {
	//return glm::dot(normal, lightSource - position) / (glm::l2Norm(normal), glm::l2Norm(lightSource - position));
	return  glm::clamp( glm::dot(glm::normalize(normal), glm::normalize(lightSource-position)), 0.0f, 1.0f);// (glm::l2Norm(normal), glm::l2Norm(lightSource-position));
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
	glm::vec3 intersectionPoint;
	glm::vec3 normalDirection;
	ofstream myfile;
	if (isDebug)	myfile.open("C:\\Users\\Deepak\\OneDrive\\Desktop\\a.txt");
	Plane* floor = new Plane(glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
	Sphere* sphere = new Sphere(glm::vec3(0, 5, 0), 4.0f);
	//myfile << "Writing this to a file.\n";
	for (float h = -height / 2.0f; h < height / 2.0f; h++) {
		for (float w = -width / 2.0f; w < width / 2.0f; w++) {
			glm::vec3 pixelPos = glm::vec3(w * cameraMultiplier, h * cameraMultiplier, -cameraToScreenDist) + camera;
			glm::vec3 rayDirection = pixelPos - camera;
			Ray ray(camera, glm::normalize(rayDirection));
			//if (isDebug) myfile << "[" << h << " " << w << "]ray start- " << camera<<" direction- "<<rayDirection << std::endl;// "\t" << rayDirection << std::endl;
			bool setColor = false;
			glm::vec3 minIntersectionPoint = glm::vec3(99999, 99999, 99999);
			glm::vec3 minNormalDirection = glm::vec3(99999, 99999, 99999);
			int objectNo = 0;
			//for (int objNo = 0; objNo < scene.size(); objNo++) {
				//if (scene[objNo]->shape == PlaneShape) {
					//std::cout << "plane\n";
			if (floor->intersect(ray, intersectionPoint, normalDirection, false) && glm::l2Norm(camera, intersectionPoint) < 200.0f) {
				if (glm::l2Norm(minIntersectionPoint, camera) > glm::l2Norm(camera, intersectionPoint)) {
					minIntersectionPoint = intersectionPoint;
					minNormalDirection = normalDirection;
					//if (isDebug) myfile << (int)(w + width / 2.0f)<<", "<< height - (int)(h + height / 2.0f) << "\n";
					//std::cout << (int)(w + width / 2.0f) << ", " << height - (int)(h + height / 2.0f) << std::endl;

					float shadedColor = lightMultiplier(intersectionPoint, lightSource, floor->normal);
					ofColor c = floor->diffuseColor;
					c.r *= shadedColor;
					c.g *= shadedColor;
					c.b *= shadedColor;
					Ray ray1(intersectionPoint, glm::normalize(lightSource - intersectionPoint));
					if (sphere->intersect(ray1, intersectionPoint, glm::normalize(normalDirection), false)) {
						c.r *= 0;
						c.g *= 0;
						c.b *= 0;
					}

					img.setColor((int)(w + width / 2.0f), height - (int)(h + height / 2.0f), c);
					objectNo = 0;
					setColor = true;
				}
			}
			//}
			//if (sphere->shape == SphereShape) {
			if (sphere->intersect(ray, intersectionPoint, normalDirection, false) && glm::l2Norm(camera, intersectionPoint) < 200.0f) {
				if (isDebug) myfile << glm::l2Norm(minIntersectionPoint, camera) << ", " << glm::l2Norm(camera, intersectionPoint) << "\n";
				if (glm::l2Norm(minIntersectionPoint, camera) > glm::l2Norm(camera, intersectionPoint)) {
					minIntersectionPoint = intersectionPoint;
					minNormalDirection = normalDirection;
					if (isDebug) myfile << (int)(w + width / 2.0f) << ", " << height - (int)(h + height / 2.0f) << "\n";
					ofColor c = sphere->diffuseColor;
					float shadedColor = lightMultiplier(intersectionPoint, lightSource, normalDirection);
					c.r *= shadedColor;
					c.g *= shadedColor;
					c.b *= shadedColor;
					img.setColor((int)(w + width / 2.0f), height - (int)(h + height / 2.0f), c);
					objectNo = 1;
					setColor = true;
				}
			}
			//}
			//if (isDebug) myfile << "[" << "abc" << "] " << glm::l2Norm(camera, intersectionPoint) << "\n";
		//}
			if (setColor) {
				//if (isDebug) myfile << "[" << "abc" << "] " << "intersectionPoint- " << intersectionPoint <<" "<< glm::l2Norm(camera, intersectionPoint)<< "\n";
			}
			//std::cout << objectNo;
		}
		//std::cout << endl;
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
	drawScene();
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
	int width = ofGetWindowWidth() - 2;
	int height = ofGetWindowHeight() - 2;

	int w = x - width / 2.0f;
	int h = height / 2.0f - y;

	Sphere* sphere = new Sphere(glm::vec3(0, 5, 0), 4.0f);
	Plane* floor = new Plane(glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));

	glm::vec3 pixelPos = glm::vec3(w * cameraMultiplier, h * cameraMultiplier, -cameraToScreenDist) + camera;
	glm::vec3 rayDirection = pixelPos - camera;
	glm::vec3 intersectionPoint;
	glm::vec3 normalDirection;
	Ray ray(camera, glm::normalize(rayDirection));
	std::cout << "\n[Debug] " << w << ", " << h << "\n";
	if (sphere->intersect(ray, intersectionPoint, normalDirection, false) && glm::l2Norm(camera, intersectionPoint) < 200.0f) {
		ofColor c = sphere->diffuseColor;
		float shadedColor = lightMultiplier(intersectionPoint, lightSource, normalDirection);
		c.r *= shadedColor;
		c.g *= shadedColor;
		c.b *= shadedColor;
		std::cout << "[Debug] shaded color-"<<shadedColor<<", color-"<<c << "\n";
		std::cout << "[Debug] light vector- " << glm::normalize(lightSource-intersectionPoint) << ", normalDirection- " << glm::normalize(normalDirection) << "\n";

	}

	if (floor->intersect(ray, intersectionPoint, normalDirection, false) && glm::l2Norm(camera, intersectionPoint) < 200.0f) {
		ofColor c = sphere->diffuseColor;
		float shadedColor = lightMultiplier(intersectionPoint, lightSource, floor->normal);
		c.r *= shadedColor;
		c.g *= shadedColor;
		c.b *= shadedColor;
		std::cout << "[Debug plane] shaded color-" << shadedColor << "\n";
		std::cout << "[Debug plane] intersectionPoint- " << intersectionPoint << ", normalDirection- " << floor->normal << "\n";

		Ray ray1(intersectionPoint, glm::normalize( lightSource - intersectionPoint));
		std::cout << "[Debug plane] ray to light-  start-" << ray1.p << ", dir- " << ray1.d<< "\n";
		if (sphere->intersect(ray1, intersectionPoint, normalDirection, true)) {
			std::cout << "[Debug plane] after intersection- " << intersectionPoint << ", normalDirection- " << floor->normal << "\n";
			c.r *= 0;
			c.g *= 0;
			c.b *= 0;
		}



	}
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