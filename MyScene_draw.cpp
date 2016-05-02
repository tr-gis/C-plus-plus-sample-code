#include "../cse452.h"
#include "MyScene.h"

void MyScene::resize(int w, int h) {
    // resize the film plane to the specified width/height
    camera.setWidthHeight(w, h);
}

/// Note: your camera and screen clear, etc, will be set up by
/// SceneviewUI.cpp *before* this gets called
void MyScene::draw() {
	// render the scene using OpenGL
	if (!isLoaded) // Don't draw if loadSceneFile hasn't been called yet
		return;


	// Turn off all lights
	for (int i = 0; i < 7; i++)
		glDisable(GL_LIGHT0 + i);
	//print_master_subgraph();

	//  .. and reset
	// TODO: draw the rest of the scene here
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, &ambientLight[0]);
	for (unsigned int i = 0; i < lights.size(); i++) {
		lights[i].setOpenGLLight(GL_LIGHT0 + i);
	}
	if (matersubgraph_list.find("root") != matersubgraph_list.end())
	{
		//std::cout << "Root found";
		glEnable(GL_COLOR_MATERIAL);
		drawTree(matersubgraph_list["root"]);
		glDisable(GL_COLOR_MATERIAL);
	}
    
}
/*#################################################Codes that I added ###########################################################*/
void MyScene::drawTree(Tree*t)
{
	//std::cout << "drawTree called";
	for (int i = 0; i < t->transblocks.size(); i++)
	{
		//t->Node*->obj*/Tree*/Vector3/Matrix4
		//std::cout << "before if statement \t";
		if (t->transblocks[i]->hasobject){
			//std::cout << "Object name" << t->transblocks[i]->obj->obj_name;

			if (t->transblocks[i]->obj->obj_name == "cube")
			{
				
				
				glPushMatrix();
				set_object_attributes(t->transblocks[i]->obj);
				for (int j = 0; j < t->transblocks[i]->matrices.size(); j++){
					glMultMatrixd(&(t->transblocks[i]->matrices[j])(0, 0));
				}
				drawcube();
				//std::cout << "cube is drawn";
				glPopMatrix();

			}
			else if (t->transblocks[i]->obj->obj_name == "cone")
			{
				
				glPushMatrix();
				set_object_attributes(t->transblocks[i]->obj);
				for (int j = 0; j < t->transblocks[i]->matrices.size(); j++){
					glMultMatrixd(&(t->transblocks[i]->matrices[j])(0, 0));
				}
				drawcone();
				//std::cout << "cone is drawn";
				glPopMatrix();
			}
			else if (t->transblocks[i]->obj->obj_name == "cylinder")
			{
				
				glPushMatrix();
				set_object_attributes(t->transblocks[i]->obj);
				for (int j = 0; j < t->transblocks[i]->matrices.size(); j++){
					glMultMatrixd(&(t->transblocks[i]->matrices[j])(0, 0));
				}
				drawcylinder();
				//std::cout << "cylinder is drawn";
				glPopMatrix();
			}
			else if (t->transblocks[i]->obj->obj_name == "sphere")
			{
				
				glPushMatrix();
				set_object_attributes(t->transblocks[i]->obj);
				for (int j = 0; j < t->transblocks[i]->matrices.size(); j++){
					glMultMatrixd(&(t->transblocks[i]->matrices[j])(0, 0));
				}
				drawsphere();
				//std::cout << "sphere is drawn";
				glPopMatrix();
			}
		}
		else {
			//std::cout << "object is a subgraph";
			glPushMatrix();
			for (int j = 0; j < t->transblocks[i]->matrices.size(); j++){
				glMultMatrixd(&(t->transblocks[i]->matrices[j])(0, 0));
			}
			//std::cout << "subgraph name is " << t->transblocks[i]->subgraph_name;
			//matersubgraph_list[t->transblocks[i]->subgraph_name];
			glEnable(GL_COLOR_MATERIAL);
			drawTree(matersubgraph_list[t->transblocks[i]->subgraph_name]);
			glDisable(GL_COLOR_MATERIAL);
			//std::cout << "after drawsubgraph";
			glPopMatrix();

		}
	}
	
}
void MyScene::set_object_attributes(Object* O){
	
	if (O->diffuse_bool){
		GLfloat diffuse[] = { O->diffuse.r, O->diffuse.g, O->diffuse.b, 1.0 };
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	}

	if (O->ambient_bool){
		GLfloat ambient[] = { O->ambient.r, O->ambient.g, O->ambient.b, 1.0 };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	}

	if (O->specular_bool){
		GLfloat specular[] = { O->specular.r, O->specular.g, O->specular.b, 1.0 };
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	}

	if (O->emit_bool){
		GLfloat emition[] = { O->emit.r, O->emit.g, O->emit.b, 1.0 };
		glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emition);
	}

	if (O->shine_bool){
		const GLfloat s = GLfloat(O->shine);
		//glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, s);
	}
}


//********************************************drawcodes**************************************************
void MyScene::drawcube(){
	int X, Y, Z;
	X = 0; Y = 0; Z = 0;
	float Ax, Ay, Az;

	/*####Face 1 #######*/
	Cube.attr.nx = 1; Cube.attr.ny = 0; Cube.attr.nz = 0;
	Cube.attr.qx = 0.5; Cube.attr.qy = 0; Cube.attr.qz = 0;
	Cube.QNs.push_back(Cube.attr);
	/*####Face 2 #######*/
	Cube.attr.nx = -1; Cube.attr.ny = 0; Cube.attr.nz = 0;
	Cube.attr.qx = -0.5; Cube.attr.qy = 0; Cube.attr.qz = 0;
	Cube.QNs.push_back(Cube.attr);
	/*####Face 3 #######*/
	Cube.attr.nx = 0; Cube.attr.ny = 1; Cube.attr.nz = 0;
	Cube.attr.qx = 0; Cube.attr.qy = 0.5; Cube.attr.qz = 0;
	Cube.QNs.push_back(Cube.attr);
	/*####Face 4 #######*/
	Cube.attr.nx = 0; Cube.attr.ny = -1; Cube.attr.nz = 0;
	Cube.attr.qx = 0; Cube.attr.qy = -0.5; Cube.attr.qz = 0;
	Cube.QNs.push_back(Cube.attr);
	/*####Face 5 #######*/
	Cube.attr.nx = 0; Cube.attr.ny = 0; Cube.attr.nz = 1;
	Cube.attr.qx = 0; Cube.attr.qy = 0; Cube.attr.qz = 0.5;
	Cube.QNs.push_back(Cube.attr);
	/*####Face 6 #######*/
	Cube.attr.nx = 0; Cube.attr.ny = 0; Cube.attr.nz = -1;
	Cube.attr.qx = 0; Cube.attr.qy = 0; Cube.attr.qz = -0.5;
	Cube.QNs.push_back(Cube.attr);

	Cube.normal.x = 1; Cube.normal.y = 0; Cube.normal.z = 0;
	while (X < 2)//X=0 means x=0.5 face,X=1 means x=-0.5 face
	{
		float dz = 1.0 / 20;
		float dy = pow(-1, X + 1)*dz;
		Ax = 0.5 - X; Ay = 0.5 - X; Az = -0.5;



		for (int i = 1; i <= 20; i++)
		{
			for (int j = 1; j <= 20; j++)
			{
				//vertices for left bottom triangle
				Xs.push_back(Ax); Ys.push_back(Ay); Zs.push_back(Az);//left bottom
				Xs.push_back(Ax); Ys.push_back(Ay); Zs.push_back(Az + (dz));//left top
				Xs.push_back(Ax); Ys.push_back(Ay + (dy)); Zs.push_back(Az);//right bottom
				Xs.push_back(Ax); Ys.push_back(Ay + (dy)); Zs.push_back(Az + (dz));//right top
				Cube.drawquad(Xs, Ys, Zs);
				Xs.clear(); Ys.clear(); Zs.clear();
				//move starting point 
				Az += dz;
			}
			Ay += dy; Az = -0.5;
		}

		++X;
		Cube.normal.x = -1; Cube.normal.y = 0; Cube.normal.z = 0;
	}


	Cube.normal.x = 0; Cube.normal.y = 1; Cube.normal.z = 0;

	while (Y < 2)//Y=0 means y=0.5 face,y=1 means y=-0.5 face
	{
		float dz = 1.0 / 20;
		float dx = pow(-1, Y)*dz;


		Ax = -0.5 + Y; Ay = 0.5 - Y; Az = -0.5;

		for (int i = 1; i <= 20; i++)
		{
			for (int j = 1; j <= 20; j++)
			{

				Xs.push_back(Ax); Ys.push_back(Ay); Zs.push_back(Az);//left bottom vertex
				Xs.push_back(Ax); Ys.push_back(Ay); Zs.push_back(Az + (dz));//left top vertex
				Xs.push_back(Ax + dx); Ys.push_back(Ay); Zs.push_back(Az);//right bottom vertex
				Xs.push_back(Ax + dx); Ys.push_back(Ay); Zs.push_back(Az + (dz));//right top vertex
				Cube.drawquad(Xs, Ys, Zs);
				Xs.clear(); Ys.clear(); Zs.clear();
				Az += dz;
			}

			Ax += dx; Az = -0.5;
		}

		++Y;
		Cube.normal.x = 0; Cube.normal.y = -1; Cube.normal.z = 0;


	}


	Cube.normal.x = 0; Cube.normal.y = 0; Cube.normal.z = 1;

	while (Z < 2)//Z=0 means z=0.5 face,Z=1 means z=-0.5 face
	{
		float dx = 1.0 / 20;
		float dy = pow(-1, Z)*dx;

		Ax = -0.5; Ay = -0.5 + Z; Az = 0.5 - Z;


		for (int i = 1; i <= 20; i++)
		{
			for (int j = 1; j <= 20; j++)
			{

				Xs.push_back(Ax); Ys.push_back(Ay); Zs.push_back(Az);
				Xs.push_back(Ax + dx); Ys.push_back(Ay); Zs.push_back(Az);
				Xs.push_back(Ax); Ys.push_back(Ay + dy); Zs.push_back(Az);
				Xs.push_back(Ax + dx); Ys.push_back(Ay + dy); Zs.push_back(Az);
				Cube.drawquad(Xs, Ys, Zs);
				Xs.clear(); Ys.clear(); Zs.clear();
				Ax += dx;
			}
			Ay += dy; Ax = -0.5;
		}

		++Z;
		Cube.normal.x = 0; Cube.normal.y = 0; Cube.normal.z = -1;


	}

}
void MyScene::drawcylinder(){
	int N = 60;
	float dy = 1.0 / 10;//this will be used to move point down the cylinder
	Xs.clear(); Ys.clear(); Zs.clear();
	//top cap

	for (int i = 0; i < N; i++){	//i is the i'th triangle of a total of N-1 triangles on each cap
		if (i < (N - 1)){

			Xs.push_back(0); Ys.push_back(0.5); Zs.push_back(0);//center of cap
			Xs.push_back(0.5*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(0.5); Zs.push_back(0.5*sin((i + 1) * 2 * 3.14 / N));
			Xs.push_back(0.5*cos(i * 2 * 3.14 / N)); Ys.push_back(0.5); Zs.push_back(0.5*sin(i * 2 * 3.14 / N));
			Cylinder.normal.x = 0; Cylinder.normal.y = 1; Cylinder.normal.z = 0;
			Cylinder.drawtriangle(Xs, Ys, Zs);
			Xs.clear(); Ys.clear(); Zs.clear();
			for (int j = 0; j < 10; j++)
			{
				Xs.push_back(0.5*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(0.5 - (dy*(j + 1))); Zs.push_back(0.5*sin((i + 1) * 2 * 3.14 / N));
				Xs.push_back(0.5*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(0.5 - (dy*j)); Zs.push_back(0.5*sin((i + 1) * 2 * 3.14 / N));
				Xs.push_back(0.5*cos(i * 2 * 3.14 / N)); Ys.push_back(0.5 - (dy*(1 + j))); Zs.push_back(0.5*sin(i * 2 * 3.14 / N));
				Xs.push_back(0.5*cos(i * 2 * 3.14 / N)); Ys.push_back(0.5 - (dy*j)); Zs.push_back(0.5*sin(i * 2 * 3.14 / N));
				Cylinder.normal.x = Xs[0]; Cylinder.normal.y = 0; Cylinder.normal.z = Zs[0];
				Cylinder.drawbarrel(Xs, Ys, Zs);
				Xs.clear(); Ys.clear(); Zs.clear();
			}

		}
		else
		{
			Xs.push_back(0); Ys.push_back(0.5); Zs.push_back(0);//center of cap
			Xs.push_back(0.5*cos(0)); Ys.push_back(0.5); Zs.push_back(0.5*sin(0));
			Xs.push_back(0.5*cos(i * 2 * 3.14 / N)); Ys.push_back(0.5); Zs.push_back(0.5*sin(i * 2 * 3.14 / N));
			Cylinder.normal.x = 0; Cylinder.normal.y = 1; Cylinder.normal.z = 0;
			Cylinder.drawtriangle(Xs, Ys, Zs);
			Xs.clear(); Ys.clear(); Zs.clear();
			for (int j = 0; j < 10; j++)
			{
				Xs.push_back(0.5*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(0.5 - (dy*(j + 1))); Zs.push_back(0.5*sin((i + 1) * 2 * 3.14 / N));
				Xs.push_back(0.5*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(0.5 - (dy*j)); Zs.push_back(0.5*sin((i + 1) * 2 * 3.14 / N));
				Xs.push_back(0.5*cos(i * 2 * 3.14 / N)); Ys.push_back(0.5 - (dy*(1 + j))); Zs.push_back(0.5*sin(i * 2 * 3.14 / N));
				Xs.push_back(0.5*cos(i * 2 * 3.14 / N)); Ys.push_back(0.5 - (dy*j)); Zs.push_back(0.5*sin(i * 2 * 3.14 / N));
				Cylinder.normal.x = Xs[0]; Cylinder.normal.y = 0; Cylinder.normal.z = Zs[0];
				Cylinder.drawbarrel(Xs, Ys, Zs);
				Xs.clear(); Ys.clear(); Zs.clear();
			}
		}
	}
	Cylinder.normal.x = 0; Cylinder.normal.y = -1; Cylinder.normal.z = 0;
	for (int i = 0; (i) < N; i++){	//i is the i'th triangle of a total of N-1 triangles on each cap
		if (i < (N - 1)){

			Xs.push_back(0); Ys.push_back(-0.5); Zs.push_back(0);//center of cap
			Xs.push_back(0.5*cos(i * 2 * 3.14 / N)); Ys.push_back(-0.5); Zs.push_back(0.5*sin(i * 2 * 3.14 / N));
			Xs.push_back(0.5*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(-0.5); Zs.push_back(0.5*sin((i + 1) * 2 * 3.14 / N));
			Cylinder.drawtriangle(Xs, Ys, Zs);
			Xs.clear(); Ys.clear(); Zs.clear();
		}
		else
		{
			Xs.push_back(0); Ys.push_back(-0.5); Zs.push_back(0);//center of cap
			Xs.push_back(0.5*cos(i * 2 * 3.14 / N)); Ys.push_back(-0.5); Zs.push_back(0.5*sin(i * 2 * 3.14 / N));
			Xs.push_back(0.5*cos(0)); Ys.push_back(-0.5); Zs.push_back(0.5*sin(0));
			Cylinder.drawtriangle(Xs, Ys, Zs);
			Xs.clear(); Ys.clear(); Zs.clear();
		}
	}
}
void MyScene::drawcone(){
	int N = 40;
	float radius = 0.5;
	float dy = 1.0 / 10;//this will be used to move up through the cone
	float dr = 0.5 / 10;
	Xs.clear(); Ys.clear(); Zs.clear();
	for (int i = 0; i < N; i++){	//i is the (i+1)'th triangle of a total of N-1 triangles on each cap
		if (i < (N - 1)){
			Xs.push_back(0); Ys.push_back(-0.5); Zs.push_back(0);//center of cap
			Xs.push_back(0.5*cos(i * 2 * 3.14 / N)); Ys.push_back(-0.5); Zs.push_back(0.5*sin(i * 2 * 3.14 / N));
			Xs.push_back(0.5*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(-0.5); Zs.push_back(0.5*sin((i + 1) * 2 * 3.14 / N));
			Cone.normal.x = 0; Cone.normal.y = -1; Cone.normal.z = 0;
			Cone.drawtriangle(Xs, Ys, Zs);
			Xs.clear(); Ys.clear(); Zs.clear();

			for (int j = 0; j < 10; j++)//shapesUI->getTessel2() 
			{
				Xs.push_back((0.5 - (dr*j))*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(-0.5 + (dy*j)); Zs.push_back((0.5 - (dr*j))*sin((i + 1) * 2 * 3.14 / N));
				Xs.push_back((0.5 - (dr*(j + 1)))*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(-0.5 + (dy*(j + 1))); Zs.push_back((0.5 - (dr*(j + 1)))*sin((i + 1) * 2 * 3.14 / N));
				Xs.push_back((0.5 - (dr*j))*cos(i * 2 * 3.14 / N)); Ys.push_back(-0.5 + (dy*j)); Zs.push_back((0.5 - (dr*j))*sin(i * 2 * 3.14 / N));
				Xs.push_back((0.5 - (dr*(j + 1)))*cos(i * 2 * 3.14 / N)); Ys.push_back(-0.5 + (dy*(j + 1))); Zs.push_back((0.5 - (dr*(j + 1)))*sin(i * 2 * 3.14 / N));
				Cone.normal.x = Xs[0]; Cone.normal.y = 0.5*sqrt(pow(Xs[0], 2) + pow(Zs[0], 2)); Cone.normal.z = Zs[0];
				Cone.drawbarrel(Xs, Ys, Zs);
				Xs.clear(); Ys.clear(); Zs.clear();

			}


		}
		else{
			Xs.push_back(0); Ys.push_back(-0.5); Zs.push_back(0);//center of cap
			Xs.push_back(0.5*cos(i * 2 * 3.14 / N)); Ys.push_back(-0.5); Zs.push_back(0.5*sin(i * 2 * 3.14 / N));
			Xs.push_back(0.5*cos(0)); Ys.push_back(-0.5); Zs.push_back(0.5*sin(0));
			Cone.normal.x = 0; Cone.normal.y = -1; Cone.normal.z = 0;
			Cone.drawtriangle(Xs, Ys, Zs);
			Xs.clear(); Ys.clear(); Zs.clear();
			for (int j = 0; j < 10; j++)//
			{
				Xs.push_back((0.5 - (dr*j))*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(-0.5 + (dy*j)); Zs.push_back((0.5 - (dr*j))*sin((i + 1) * 2 * 3.14 / N));
				Xs.push_back((0.5 - (dr*(j + 1)))*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(-0.5 + (dy*(j + 1))); Zs.push_back((0.5 - (dr*(j + 1)))*sin((i + 1) * 2 * 3.14 / N));
				Xs.push_back((0.5 - (dr*j))*cos(i * 2 * 3.14 / N)); Ys.push_back(-0.5 + (dy*j)); Zs.push_back((0.5 - (dr*j))*sin(i * 2 * 3.14 / N));
				Xs.push_back((0.5 - (dr*(j + 1)))*cos(i * 2 * 3.14 / N)); Ys.push_back(-0.5 + (dy*(j + 1))); Zs.push_back((0.5 - (dr*(j + 1)))*sin(i * 2 * 3.14 / N));
				Cone.normal.x = Xs[0]; Cone.normal.y = 0.5*sqrt(pow(Xs[0], 2) + pow(Zs[0], 2)); Cone.normal.z = Zs[0];
				Cone.drawbarrel(Xs, Ys, Zs);
				Xs.clear(); Ys.clear(); Zs.clear();

			}

		}

	}

}
void MyScene::drawsphere(){
	float dtheta = 5;
	float r = 0.5;
	float dphi = 0;
	dphi = 360 / (100);

	float d2r = 3.14 / 180;//convertion factor from degree to radians
	for (int theta = -90; theta <= 90 - dtheta; theta += dtheta) {//
		for (int phi = 0; phi <= 360 - dphi; phi += dphi) {//
			Xs.push_back(r*cos((theta + dtheta)*d2r) * cos(phi*d2r)); Ys.push_back(r*cos((theta + dtheta)*d2r) * sin(phi*d2r)); Zs.push_back(r*sin((theta + dtheta)*d2r));
			Xs.push_back(r*cos((theta + dtheta)*d2r) * cos((phi + dphi)*d2r)); Ys.push_back(r*cos((theta + dtheta)*d2r) * sin((phi + dphi)*d2r)); Zs.push_back(r*sin((theta + dtheta)*d2r));
			Xs.push_back(r*cos(theta*d2r) * cos(phi*d2r)); Ys.push_back(r*cos(theta*d2r) * sin(phi*d2r)); Zs.push_back(r*sin(theta*d2r));
			S.normal.x = -Xs[1]; S.normal.y = -Ys[1]; S.normal.z = -Zs[1];
			S.drawtriangle(Xs, Ys, Zs);
			if (theta > -90 && theta < 90) {
				Xs.push_back(r*cos(theta*d2r) * cos((phi + dphi)*d2r)); Ys.push_back(r*cos(theta*d2r) * sin((phi + dphi)*d2r)); Zs.push_back(r*sin(theta*d2r));
				S.drawquad(Xs, Ys, Zs);

			}
			Xs.clear(); Ys.clear(); Zs.clear();
		}
	}


}
