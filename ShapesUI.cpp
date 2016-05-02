#include "../cse452.h"
#include "ShapesUI.h"
#include "ShapesInterface.h"
#include "../Color.h"
#include <FL/Fl.H>
#include <FL/gl.h>
#include <GL/glu.h>
#include <Vector3.h>
#include <deque>
#include "../intersection/HitRecord.h"


ShapesUI::ShapesUI() {
    width = height = 0;

    // ToDo: initialize your variables here
}

ShapesUI::~ShapesUI() {
    // ToDo: delete your variables here
}

void ShapesUI::resize(int w, int h) {
    width = w;
    height = h;
}

void ShapesUI::draw() {
	// Sets up the viewport and background color
	setup3DDrawing(Color(0, 0, 0), width, height, true);

	// Changes the way triangles are drawn
	switch (shapesUI->getDisplayType()) {
	case DISPLAY_WIREFRAME: {
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT, GL_LINE);
		glColor3f(1.0f, 1.0f, 1.0f);
	} break;
	case DISPLAY_FLAT_SHADING: {
		glEnable(GL_LIGHTING);
		glPolygonMode(GL_FRONT, GL_FILL);
		glColor3f(1.0f, 1.0f, 1.0f);
		glShadeModel(GL_FLAT);
	} break;
	case DISPLAY_SMOOTH_SHADING: {
		glEnable(GL_LIGHTING);
		glPolygonMode(GL_FRONT, GL_FILL);
		glColor3f(1.0f, 1.0f, 1.0f);
		glShadeModel(GL_SMOOTH);
	} break;
	default: break;
	}

	// Setup the camera
	gluLookAt(3.5 * cos(shapesUI->getYRot()) * cos(shapesUI->getXRot()),
		3.5 * sin(shapesUI->getYRot()),
		3.5 * cos(shapesUI->getYRot()) * sin(shapesUI->getXRot()), 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

	
/*#################################################Codes that I added ###########################################################*/
	Shape Cube;
	Cyl Cylinder;
	Cone Cone;
	if (shapesUI->getShapeType() == 3)//when the shape to be drawn is cube
	{

		int X, Y, Z;
		X = 0; Y = 0; Z = 0;
		float Ax, Ay, Az;

		Cube.normal.x = 1; Cube.normal.y = 0; Cube.normal.z = 0;
		
		while (X < 2)//X=0 means x=0.5 face,X=1 means x=-0.5 face
		{
			float dz = 1.0 / shapesUI->getTessel1();
			float dy = pow(-1, X + 1)*dz;
			Ax = 0.5 - X; Ay = 0.5 - X; Az = -0.5;
			
			for (int i = 1; i <= shapesUI->getTessel1(); i++)
			{
				for (int j = 1; j <= shapesUI->getTessel1(); j++)
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
		Cube.attr.nx = Cube.normal.x; Cube.attr.ny = Cube.normal.y; Cube.attr.nz = Cube.normal.z;
		while (Y < 2)//Y=0 means y=0.5 face,y=1 means y=-0.5 face
		{
			float dz = 1.0 / shapesUI->getTessel1();
			float dx = pow(-1, Y)*dz;


			Ax = -0.5 + Y; Ay = 0.5 - Y; Az = -0.5;
			Cube.attr.qx = Ax; Cube.attr.qy = Ay; Cube.attr.qz = Az;

			for (int i = 1; i <= shapesUI->getTessel1(); i++)
			{
				for (int j = 1; j <= shapesUI->getTessel1(); j++)
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
			Cube.QNs.push_back(Cube.attr);
			Cube.normal.x = 0; Cube.normal.y = -1; Cube.normal.z = 0;
			Cube.attr.nx = Cube.normal.x; Cube.attr.ny = Cube.normal.y; Cube.attr.nz = Cube.normal.z;
			++Y;
		}
		Cube.QNs.push_back(Cube.attr);

		Cube.normal.x = 0; Cube.normal.y = 0; Cube.normal.z = 1;
		while (Z < 2)//Z=0 means z=0.5 face,Z=1 means z=-0.5 face
		{
			float dx = 1.0 / shapesUI->getTessel1();
			float dy = pow(-1, Z)*dx;

			Ax = -0.5; Ay = -0.5 + Z; Az = 0.5 - Z;

			for (int i = 1; i <= shapesUI->getTessel1(); i++)
			{
				for (int j = 1; j <= shapesUI->getTessel1(); j++)
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
			Cube.normal.x = 0; Cube.normal.y = 0; Cube.normal.z = -1;
			++Z;
		}


	}

	if (shapesUI->getShapeType() == 2)//the shape is cylinder
	{

		int N = shapesUI->getTessel1();
		float dy = 1.0 / shapesUI->getTessel2();//this will be used to move point down the cylinder
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
				for (int j = 0; j < shapesUI->getTessel2(); j++)
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
				for (int j = 0; j < shapesUI->getTessel2(); j++)
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

		//bottom cap
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

	if (shapesUI->getShapeType() == 1)//Cone is the shape chosen
	{
		int N = shapesUI->getTessel1();
		float radius = 0.5;
		float dy = 1.0 / shapesUI->getTessel2();//this will be used to move up the cone
		float dr = 0.5 / shapesUI->getTessel2();
		Xs.clear(); Ys.clear(); Zs.clear();
		for (int i = 0; i < N; i++){	//i is the (i+1)'th triangle of a total of N-1 triangles on each cap
			if (i < (N - 1)){

				Xs.push_back(0); Ys.push_back(-0.5); Zs.push_back(0);//center of cap
				Xs.push_back(0.5*cos(i * 2 * 3.14 / N)); Ys.push_back(-0.5); Zs.push_back(0.5*sin(i * 2 * 3.14 / N));
				Xs.push_back(0.5*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(-0.5); Zs.push_back(0.5*sin((i + 1) * 2 * 3.14 / N));
				Cone.normal.x = 0; Cone.normal.y = -1; Cone.normal.z = 0;
				Cone.drawtriangle(Xs, Ys, Zs);
				Xs.clear(); Ys.clear(); Zs.clear();

				for (int j = 0; j < shapesUI->getTessel2(); j++)//shapesUI->getTessel2() 
				{
					Xs.push_back((0.5 - (dr*j))*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(-0.5 + (dy*j)); Zs.push_back((0.5 - (dr*j))*sin((i + 1) * 2 * 3.14 / N));
					Xs.push_back((0.5 - (dr*(j + 1)))*cos((i + 1) * 2 * 3.14 / N)); Ys.push_back(-0.5 + (dy*(j + 1))); Zs.push_back((0.5 - (dr*(j + 1)))*sin((i + 1) * 2 * 3.14 / N));
					Xs.push_back((0.5 - (dr*j))*cos(i * 2 * 3.14 / N)); Ys.push_back(-0.5 + (dy*j)); Zs.push_back((0.5 - (dr*j))*sin(i * 2 * 3.14 / N));
					Xs.push_back((0.5 - (dr*(j + 1)))*cos(i * 2 * 3.14 / N)); Ys.push_back(-0.5 + (dy*(j + 1))); Zs.push_back((0.5 - (dr*(j + 1)))*sin(i * 2 * 3.14 / N));
					Cone.normal.x = Xs[0]; Cone.normal.y = 0.5*sqrt(pow(Xs[0],2)+pow(Zs[0],2)); Cone.normal.z = Zs[0];
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
				for (int j = 0; j < shapesUI->getTessel2(); j++)//
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
	Batman bat;
	if (shapesUI->getShapeType() == 4)//batman time !!! 
	{
		std::vector<float> b_Xs;
		std::vector<float> b_Ys;
		std::vector<float> b_Zs;
		int N = 14;
		Xs.clear(); Ys.clear(); Zs.clear();
		//top cap

		for (int i = 0; i < N; i++){	//i is the i'th triangle of a total of N-1 triangles on each cap
			Xs.push_back(0.5*cos(i * 2 * 3.14 / N)); Ys.push_back(1.0); Zs.push_back(0.5*sin(i * 2 * 3.14 / N));
			}
		//left wing tip
		b_Xs.push_back(Xs[5]); b_Ys.push_back(0.5); b_Zs.push_back(Zs[5]);
		b_Xs.push_back(Xs[5]); b_Ys.push_back(0.9); b_Zs.push_back(Zs[5]);
		b_Xs.push_back(Xs[6]); b_Ys.push_back(0.7); b_Zs.push_back(Zs[6]);
		bat.normal.x = b_Xs[0]; bat.normal.y=0; bat.normal.z = b_Zs[0];
		bat.drawtriangle(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();
		//left wing upper portion
		b_Xs.push_back(Xs[5]); b_Ys.push_back(0.5); b_Zs.push_back(Zs[5]);
		b_Xs.push_back(Xs[5]); b_Ys.push_back(0.9); b_Zs.push_back(Zs[5]);
		b_Xs.push_back(Xs[3]); b_Ys.push_back(0.5); b_Zs.push_back(Zs[3]);
		b_Xs.push_back(Xs[3]); b_Ys.push_back(0.9); b_Zs.push_back(Zs[3]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawbarrel(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();
		//left wing lower portion
		b_Xs.push_back(Xs[4]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[4]);
		b_Xs.push_back(Xs[4]); b_Ys.push_back(0.5); b_Zs.push_back(Zs[4]);
		b_Xs.push_back(Xs[3]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[3]);
		b_Xs.push_back(Xs[3]); b_Ys.push_back(0.5); b_Zs.push_back(Zs[3]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawbarrel(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();
		//left wing joint
		b_Xs.push_back(Xs[3]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[3]);
		b_Xs.push_back(Xs[3]); b_Ys.push_back(0.6); b_Zs.push_back(Zs[3]);
		b_Xs.push_back(Xs[2]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[2]);
		b_Xs.push_back(Xs[2]); b_Ys.push_back(0.6); b_Zs.push_back(Zs[2]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawbarrel(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();
		//main body left half
		b_Xs.push_back(Xs[2]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[2]);
		b_Xs.push_back(Xs[2]); b_Ys.push_back(0.7); b_Zs.push_back(Zs[2]);
		b_Xs.push_back(Xs[0]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[0]);
		b_Xs.push_back(Xs[0]); b_Ys.push_back(0.7); b_Zs.push_back(Zs[0]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawbarrel(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();
		//main body right half
		b_Xs.push_back(Xs[0]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[0]);
		b_Xs.push_back(Xs[0]); b_Ys.push_back(0.7); b_Zs.push_back(Zs[0]);
		b_Xs.push_back(Xs[12]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[12]);
		b_Xs.push_back(Xs[12]); b_Ys.push_back(0.7); b_Zs.push_back(Zs[12]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawbarrel(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();

		//tail left half
		b_Xs.push_back(Xs[0]); b_Ys.push_back(0); b_Zs.push_back(Zs[0]);
		b_Xs.push_back(Xs[0]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[0]);
		b_Xs.push_back(Xs[2]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[2]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawtriangle(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();
		//tail right half
		b_Xs.push_back(Xs[0]); b_Ys.push_back(0); b_Zs.push_back(Zs[0]);
		b_Xs.push_back(Xs[12]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[12]);
		b_Xs.push_back(Xs[0]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[0]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawtriangle(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();
		//left ear
		b_Xs.push_back(Xs[2]); b_Ys.push_back(0.7); b_Zs.push_back(Zs[2]);
		b_Xs.push_back(Xs[0]); b_Ys.push_back(0.7); b_Zs.push_back(Zs[0]);
		b_Xs.push_back(Xs[1]); b_Ys.push_back(0.9); b_Zs.push_back(Zs[1]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawtriangle(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();
		//right ear
		b_Xs.push_back(Xs[0]); b_Ys.push_back(0.7); b_Zs.push_back(Zs[0]);
		b_Xs.push_back(Xs[12]); b_Ys.push_back(0.7); b_Zs.push_back(Zs[12]);
		b_Xs.push_back(Xs[13]); b_Ys.push_back(0.9); b_Zs.push_back(Zs[13]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawtriangle(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();
		//right wing joint
		b_Xs.push_back(Xs[12]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[12]);
		b_Xs.push_back(Xs[12]); b_Ys.push_back(0.6); b_Zs.push_back(Zs[12]);
		b_Xs.push_back(Xs[11]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[11]);
		b_Xs.push_back(Xs[11]); b_Ys.push_back(0.6); b_Zs.push_back(Zs[11]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawbarrel(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();
		//right wing lower portion
		b_Xs.push_back(Xs[11]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[11]);
		b_Xs.push_back(Xs[11]); b_Ys.push_back(0.5); b_Zs.push_back(Zs[11]);
		b_Xs.push_back(Xs[10]); b_Ys.push_back(0.3); b_Zs.push_back(Zs[10]);
		b_Xs.push_back(Xs[10]); b_Ys.push_back(0.5); b_Zs.push_back(Zs[10]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawbarrel(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();
		//right wing upper portion
		b_Xs.push_back(Xs[11]); b_Ys.push_back(0.5); b_Zs.push_back(Zs[11]);
		b_Xs.push_back(Xs[11]); b_Ys.push_back(0.9); b_Zs.push_back(Zs[11]);
		b_Xs.push_back(Xs[9]); b_Ys.push_back(0.5); b_Zs.push_back(Zs[9]);
		b_Xs.push_back(Xs[9]); b_Ys.push_back(0.9); b_Zs.push_back(Zs[9]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawbarrel(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();
		//right wing tip
		b_Xs.push_back(Xs[9]); b_Ys.push_back(0.9); b_Zs.push_back(Zs[9]);
		b_Xs.push_back(Xs[9]); b_Ys.push_back(0.5); b_Zs.push_back(Zs[9]);
		b_Xs.push_back(Xs[8]); b_Ys.push_back(0.7); b_Zs.push_back(Zs[8]);
		bat.normal.x = b_Xs[0]; bat.normal.y = 0; bat.normal.z = b_Zs[0];
		bat.drawtriangle(b_Xs, b_Ys, b_Zs);
		b_Xs.clear(); b_Ys.clear(); b_Zs.clear();



	}
	//*********************************************************************************************
	Sphere S;
	if (shapesUI->getShapeType() == 0)//the object selected is sphere
	{
		float dtheta = 5;
		float r = 0.5;
		float dphi=0;
		if (shapesUI->getTessel1()!=0)
			dphi = 360 / (shapesUI->getTessel1());
		else
			dphi = 360 / (10);

		float d2r = 3.14 / 180;//convertion factor from degree to radians
		for (int theta = -90; theta <= 90 -dtheta; theta += dtheta) {//
			for (int phi = 0; phi <= 360 - dphi; phi += dphi) {//
				Xs.push_back(r*cos((theta + dtheta)*d2r) * cos(phi*d2r)); Ys.push_back(r*cos((theta + dtheta)*d2r) * sin(phi*d2r)); Zs.push_back(r*sin((theta + dtheta)*d2r));
				Xs.push_back(r*cos((theta + dtheta)*d2r) * cos((phi + dphi)*d2r)); Ys.push_back(r*cos((theta + dtheta)*d2r) * sin((phi + dphi)*d2r)); Zs.push_back(r*sin((theta + dtheta)*d2r));
				Xs.push_back(r*cos(theta*d2r) * cos(phi*d2r)); Ys.push_back(r*cos(theta*d2r) * sin(phi*d2r)); Zs.push_back(r*sin(theta*d2r));
				S.normal.x=-Xs[1]; S.normal.y=-Ys[1]; S.normal.z=-Zs[1];
				S.drawtriangle(Xs, Ys, Zs);
				if (theta > -90 && theta < 90) {
					Xs.push_back(r*cos(theta*d2r) * cos((phi + dphi)*d2r)); Ys.push_back(r*cos(theta*d2r) * sin((phi + dphi)*d2r)); Zs.push_back(r*sin(theta*d2r));
					S.drawquad(Xs, Ys, Zs);
					
									}
				Xs.clear(); Ys.clear(); Zs.clear();
				}
			}
		
		}
	
	endDrawing();
}


int ShapesUI::handle(int event) {
    return 0;
}

void ShapesUI::changedShape()
{
    // ToDo: Switch shapes

    
    RedrawWindow();
}

void ShapesUI::changedTessel( ) {
    // ToDo: tessellate your shape here

    
    RedrawWindow();
}

/*INCOMPLETE ISOCAHEDRON TESSELLATION OF SPHERE
std::deque<Triangle> triangle;
double theta = 26.56505117707799 * 3.14 / 180.0;
double stheta = sin(theta);
double ctheta = 0.5;
int tri_ver[20][3] = { { 0, 1, 2 },{ 0, 2, 3 }, { 0, 3, 4 }, { 0, 4, 5 }, { 0, 5, 1 },
{ 11, 7, 6 }, { 11, 8, 7 }, { 11, 9, 8 }, { 11, 10, 9 }, { 11, 6, 10 },
{ 1, 10, 6 }, { 6, 2, 1 }, { 2, 6, 7 }, { 7, 3, 2 }, { 3, 7, 8 }, { 8, 3, 4 }, { 4, 8, 9 }, { 9, 4, 5 }, { 5, 9, 10 }, {10,5,1} };//

S.icosaVertices.push_back(Point3(0.0, 0.0, 0.5));//top vertex
//upper pentagon coordinates
double phi = 0;
for (int i = 1; i < 6; ++i) {
S.icosaVertices.push_back(Point3(ctheta * cos(phi), ctheta * sin(phi), stheta));
phi += (2.0 * 3.14 / 5.0);
}

// the lower pentagon
phi = 3.14 / 5.0;
for (int i = 1; i < 6; ++i) {
S.icosaVertices.push_back(Point3(ctheta * cos(phi), ctheta * sin(phi), -stheta));
phi += 2.0 * 3.14 / 5.0;
}
S.icosaVertices.push_back(Point3(0.0, 0.0, -0.5));//bottom vertex
for (int i = 0; i < 20; i++)//decides the number of triagnles to be drawn
{
Triangle tri(S.icosaVertices[tri_ver[i][0]], S.icosaVertices[tri_ver[i][1]], S.icosaVertices[tri_ver[i][2]]);
//tri.print();
triangle.push_back(tri);
S.drawtriangle(tri);
}
*/
/*
for (int i = 0; i < S.icosaVertices.size(); i++)
{
std::cout << "Point " << i << "\t";
S.icosaVertices[i].print();
std::cout << std::endl;
}
system("PAUSE");
*/
