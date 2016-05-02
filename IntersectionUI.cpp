#include "../cse452.h"
#include "IntersectionInterface.h"
#include "HitRecord.h"
#include "../shapes/ShapesUI.h"
#include <FL/Fl.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Value_Input.H>
#include <FL/gl.h>
#include <GL/glu.h>
#include <fstream>

IntersectionUI::IntersectionUI() {
    width = height = 0;
}

IntersectionUI::~IntersectionUI() {
}

void IntersectionUI::resize(int w, int h) {
    width = w;
    height = h;
}

void IntersectionUI::draw() {
    setup3DDrawing( Color(1,1,1), width, height, true );

    glMatrixMode(GL_MODELVIEW);
    gluLookAt( 3.5 * cos( intersectionUI->getYRot() ) * cos( intersectionUI->getXRot() ), 
               3.5 * sin( intersectionUI->getYRot() ), 
               3.5 * cos( intersectionUI->getYRot() ) * sin( intersectionUI->getXRot() ), 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);


    glDisable(GL_LIGHTING);
    if ( intersectionUI->m_bGrid->value() ) {
        glBegin(GL_LINES);
        // draw grid
        glColor3f(0.0f, 0.0f, 0.0f);
        for (int i = 0; i <= 10; i++) {
            float s = -2.0f + i / 2.5f;
            glVertex3f(s, 0.0f, -2.0f);
            glVertex3f(s, 0.0f,  2.0f);
            glVertex3f(-2.0f, 0.0f, s);
            glVertex3f( 2.0f, 0.0f, s);
        }
        glEnd();
        // draw (X,Y,Z) axes
        glLineWidth(3.0f);
        glBegin(GL_LINES);
        glColor3f(1.0f, 0.0f, 0.0f);
        glVertex3f(-2.0f, 0.0f, 0.0f);
        glVertex3f( 2.0f, 0.0f, 0.0f);
        glColor3f(0.0f, 1.0f, 0.0f);
        glVertex3f(0.0f, -2.0f, 0.0f);
        glVertex3f(0.0f,  2.0f, 0.0f);
        glColor3f(0.0f, 0.0f, 1.0f);
        glVertex3f(0.0f, 0.0f, -2.0f);
        glVertex3f(0.0f, 0.0f,  2.0f);
        glEnd();
        glLineWidth(1.0f);
    }

    // compute ray origin from parameters
    Point3 pAt( intersectionUI->m_dXAt->value(), intersectionUI->m_dYAt->value(), intersectionUI->m_dZAt->value() );

    // compute ray direction from parameters
    Vector3 dir;
    dir[0] = cos(intersectionUI->getPhi()) * cos(intersectionUI->getTheta());
    dir[1] = sin(intersectionUI->getPhi());
    dir[2] = cos(intersectionUI->getPhi()) * sin(intersectionUI->getTheta());

    const Point3 pE1 = pAt - dir * 2.0;
    const Point3 pE2 = pAt + dir * 2.0;
    
    if (intersectionUI->m_bRay->value()) {
        glPointSize(6.0f);
        glLineWidth(3.0f);
        glColor3f(0.5f, 0.5f, 0.0f);
        glBegin(GL_POINTS);
        glVertex3dv( &pE1[0]);
        glEnd();
        glColor3f(0.5f, 0.0f, 0.5f);
        glBegin(GL_LINES);
        glVertex3dv( &pE1[0]);
        glVertex3dv( &pE2[0]);
        glEnd();
        glLineWidth(1.0f);
        glPointSize(1.0f);
    }
    if (intersectionUI->m_bRayShadow->value()) {
        glPointSize(6.0f);
        glLineWidth(2.0f);
        glColor3f(0.1f, 0.1f, 0.0f);
        glBegin(GL_POINTS);
        glVertex3d( pE1[0], 0.0, pE1[2] );
        glEnd();
        glColor3f(0.0f, 0.0f, 0.0f);
        glEnable(GL_LINE_STIPPLE);
        glLineStipple(1, 0xF0F0);
        glBegin(GL_LINES);
        glVertex3d( pE1[0], 0.0, pE1[2] );
        glVertex3d( pE2[0], 0.0, pE2[2] );
        glVertex3d( pE1[0], 0.0, pE1[2] );
        glVertex3dv( &pE1[0]);
        glVertex3d( pE2[0], 0.0, pE2[2] );
        glVertex3dv( &pE2[0]);
        glEnd();
        glLineWidth(1.0f);
        glPointSize(1.0f);
        glDisable(GL_LINE_STIPPLE);
    }
    glEnable(GL_LIGHTING);

    glEnable( GL_BLEND );
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    const float colPlane[4] = {0.5, 0.5, 0.75, 0.5};
    const float colObj[4] = {0.5, 0.25, 0.25, 0.5};
    glColor4fv( colPlane );
    glMaterialfv(GL_FRONT, GL_DIFFUSE , colPlane);

    glBegin( GL_POLYGON );
    glVertex3f( -2.0, 0.0, -2.0 );
    glVertex3f(  2.0, 0.0, -2.0 );
    glVertex3f(  2.0, 0.0,  2.0 );
    glVertex3f( -2.0, 0.0,  2.0 );
    glEnd();

    glBegin( GL_POLYGON );
    glVertex3f( -2.0, 0.0001f, -2.0 );
    glVertex3f( -2.0, 0.0001f,  2.0 );
    glVertex3f(  2.0, 0.0001f,  2.0 );
    glVertex3f(  2.0, 0.0001f, -2.0 );
    glEnd();

    
    glColor4fv( colObj );
    glMaterialfv(GL_FRONT, GL_DIFFUSE , colObj);
    glMaterialfv(GL_FRONT, GL_SPECULAR, colObj);

    // ToDo: draw your shape here and perform the intersection
    // then call drawHits so you can see where the ray has hit the shape
    // the origin is in variable 'p' and direction in variable 'dir'

	
	if (current_shape == 0){
		drawsphere();
		drawHits(S.intersect(pE1, dir));
		
	}
	if (current_shape == 2){
		drawcone();
		drawHits(Cone.intersect(pE1, dir));
	}
	if (current_shape == 1)
	{
		drawcylinder();
		drawHits(Cylinder.intersect(pE1, dir));
	}
	if (current_shape == 3)
	{
		drawcube();
		drawHits(Cube.intersect(pE1, dir));
	}
	

    //Call HitRecord hr = intersect(pE1, dir);
    //drawHits(hr);

    endDrawing();
	
}

void IntersectionUI::drawHits(HitRecord& hr) {
    double t, u, v;
    Point3 p;
    Vector3 n;
    glDisable(GL_LIGHTING);
    while (hr.getFirstHit(t, u, v, p, n)) {
        glPointSize(8.0f);
        glLineWidth(6.0f);
        glColor3f(1.0f, 0.0f, 0.0f);
        glBegin(GL_POINTS);
        glVertex3d(p[0], p[1], p[2]);
        glEnd();
        glColor3f(0.2f, 0.2f, 0.2f);
        glBegin(GL_LINES);
        glVertex3d(p[0], p[1], p[2]);
        glVertex3d(p[0] + 0.5 * n[0], p[1] + 0.5 * n[1], p[2] + 0.5 * n[2]);
        glEnd();
        glLineWidth(1.0f);
        glPointSize(1.0f);
        hr.removeFirstHit();
    }
}

void IntersectionUI::changeShape( ShapesUI::ShapeType type )
{
    // ToDo: Change which shape
	
    switch ( type ) {
    case ShapesUI::SHAPE_SPHERE :
		//std::cout << "Current shape is sphere \n";
		current_shape = 0;
		break;
    case ShapesUI::SHAPE_CYLINDER : 
		//std::cout << "Current shape is cylinder\n";
		current_shape = 1;
		break;
    case ShapesUI::SHAPE_CONE :
		//std::cout << "Current shape is cone\n";
		current_shape = 2;
		break;
    case ShapesUI::SHAPE_CUBE : 
		//std::cout << "Current shape is cube\n";
		current_shape = 3;
		break;
    }
}


int IntersectionUI::handle(int event) {
    return 0;
}


void IntersectionUI::writeTest() {
    // creates a deterministic sequence of ray positions and directions
    // and writes the resulting intersections to a file
    // you must add the proper intersect calls for this file to be generated
    
    double invBase[5] = {1.0 / 2.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 11.0};
    double values[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    std::ofstream file("../intersections.txt");
    file.precision(4);

    const int seed = static_cast<int>(intersectionUI->m_iSeed->value());
    // generate a halton sequence to pick position/ray combinations
    // skip the first 'seed' values
    for (int i = 0; i < seed; i++) {
        for (int j = 0; j < 5; j++) {
            double r = 1.0 - values[j] - 1e-10;
            if (invBase[j] < r)
                values[j] += invBase[j];
            else {
                double hh;
                double h = invBase[j];
                do {
                    hh = h;
                    h *= invBase[j];
                } while (h >= r);
                values[j] += ((hh + h) - 1.0);
            }
        }
    }
    for (int i = seed; i < (seed + 1638); i++) {
        for (int j = 0; j < 5; j++) {
            double r = 1.0 - values[j] - 1e-10;
            if (invBase[j] < r)
                values[j] += invBase[j];
            else {
                double hh;
                double h = invBase[j];
                do {
                    hh = h;
                    h *= invBase[j];
                } while (h >= r);
                values[j] += ((hh + h) - 1.0);
            }
        }
        // create the ray from the five random values
        // compute ray origin
        Point3 p;
        p[0] = values[4] * sin(values[0] * M_PI) * cos(values[1] * 2.0 * M_PI);
        p[1] = values[4] * sin(values[0] * M_PI) * sin(values[1] * 2.0 * M_PI);
        p[2] = values[4] * cos(values[0] * M_PI);
        // compute ray direction
        Vector3 dir;
        dir[0] = sin(values[2] * M_PI) * cos(values[3] * 2.0 * M_PI);
        dir[1] = sin(values[2] * M_PI) * sin(values[3] * 2.0 * M_PI);
        dir[2] = cos(values[2] * M_PI);
        
        HitRecord cubeHr, cylinderHr, coneHr, sphereHr;
        // ToDo: intersect with your shapes here and store the result
        // in the appropriate hit record
		
        cubeHr=Cube.intersect(p, dir);
		cylinderHr=Cylinder.intersect(p, dir);
        coneHr = Cone.intersect(p, dir);
        sphereHr = S.intersect(p, dir);

        // write out
        file << i << " Cube     " << cubeHr     << std::endl;
        file << i << " Cylinder " << cylinderHr << std::endl;
        file << i << " Cone     " << coneHr     << std::endl;
        file << i << " Sphere   " << sphereHr   << std::endl;
    }
    file.close();
}
/*###################################  Draw Routines that I added #############################################*/
void IntersectionUI::drawsphere(){
	
	float dtheta = 5;
	float r = 0.5;
	float dphi = 0;
	dphi = 360 /(100);

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
void IntersectionUI::drawcone()
{
	
	int N = 20;
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

void IntersectionUI::drawcylinder(){
	
	int N = 20;
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
void IntersectionUI::drawcube(){
	
	int X, Y, Z;
	X = 0; Y = 0; Z = 0;
	float Ax, Ay, Az;

	/*####Face 1 #######*/
	Cube.attr.nx = 1; Cube.attr.ny = 0; Cube.attr.nz = 0;
	Cube.attr.qx = 0.5; Cube.attr.qy =0; Cube.attr.qz = 0;
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


