#ifndef _INTERSECTION_UI_H_
#define _INTERSECTION_UI_H_

#include "../cse452.h"
#include "../shapes/ShapesUI.h"
#include "HitRecord.h"
#include "../UIInterface.h"
#include <FL/Fl_Window.H>

class IntersectionInterface;
class IntersectionUI : public UIInterface {
public:
    IntersectionUI();
    ~IntersectionUI();

	//vectors that I have declared
	std::vector<float> Xs;
	std::vector<float> Ys;
	std::vector<float> Zs;
	//shapes that I declared
	Sphere S; Cone Cone; Cyl Cylinder; Shape Cube;

    // Inherited from userInterface
    void resize(int width, int height);
    void draw();
    int handle(int event);
	//changes I made
	void drawsphere();   //0 
	void drawcylinder(); //1
	void drawcone();     //2
	void drawcube();     //3

	int current_shape;

    // Link to the intersection widget
    void setUI( const IntersectionInterface *in_ui ) { intersectionUI = in_ui; }
    void changeShape( ShapesUI::ShapeType type );

    void writeTest();

private:
    const IntersectionInterface *intersectionUI;
    int width, height;
	

    void drawHits(HitRecord& hr);

    // declare your variables here
};

#endif /* _INTERSECTION_UI_H_ */
