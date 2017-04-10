/*
  CSC418 - RayTracer code - Winter 2017 - Assignment 3&4
  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.
  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c
  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.
  You only need to modify or add code in sections
  clearly marked "TO DO"
*/

#include "utils.h"
#include "float.h"
#include "omp.h"

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
struct colourRGB background;  
int MAX_DEPTH;

void buildScene(void)
{
 // Sets up all objects in the more complex scene for Assignment 4 ( a pool table scene)

 struct object3D *o;
 struct pointLS *l;
 struct point3D p; 

  //foreground texture
  const char pooltable[] = "pool_table_green.ppm";
  const char wood[] = "wood.ppm";
  const char border[] = "wood_border.ppm";
  
  //pool ball texture filenames  
  const char ball1[] = "1_full.ppm";
  const char ball3[] = "3_full.ppm";
  const char ball4[] = "4_full.ppm";
  const char ball5[] = "5_full.ppm";
  const char ball6[] = "6_full.ppm";
  const char ball8[] = "8_full.ppm";
  const char ball9[] = "9_stripe.ppm";
  const char ball10[] = "10_stripe.ppm";
  const char ball13[] = "13_stripe.ppm";
  const char ball15[] = "15_stripe.ppm";
 
  //environment mapping texture
  const char wall[] = "woodwall.ppm";
  const char ceiling[] = "ceiling.ppm";
 
 
   // This is the rotation offset
  double offset = PI/6.0;
  
  // This is the room size
  double wallsize = 300;
  double baseheight = 0;
 
  //Creating Room
 
  //~ // Left wall
  o=newPlane(.05,.75,0,1,.55,.55,.75,1,1,0); 
  Scale(o,wallsize,wallsize,1);      
  RotateY(o,-PI/2);
  Translate(o,-wallsize,0, 0);
  RotateY(o,offset);
  Translate(o,0,-3 + baseheight ,10);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, wall);
  insertObject(o,&object_list);    

  //~ // Right wall
  o=newPlane(.05,.75,0,1,.55,.55,.75,1,1,0);  
  Scale(o,wallsize,wallsize,1);     
  RotateY(o,PI/2);
  Translate(o,wallsize,0, 0);
  RotateY(o,offset);
  Translate(o,0,-3 + baseheight ,10);
  invert(&o->T[0][0],&o->Tinv[0][0]);   
  loadTexture(o, wall);
  insertObject(o,&object_list);  

  // Back Wall
  o=newPlane(.05,.75,0,1,.55,.55,.75,1,1,0);  
  Scale(o,wallsize,wallsize,1);      
  Translate(o,0, 0, wallsize);
  RotateY(o,offset);
  Translate(o,0,-3 + baseheight ,10);
  invert(&o->T[0][0],&o->Tinv[0][0]);   
  loadTexture(o, wall);
  insertObject(o,&object_list);     

  // Front Wall
  o=newPlane(.05,.75,0,1,.55,.55,.75,1,1,0);  
  Scale(o,wallsize,wallsize,1);     
  Translate(o,0, 0, -wallsize);
  RotateY(o,offset);
  Translate(o,0,-3 + baseheight ,10);
  invert(&o->T[0][0],&o->Tinv[0][0]);   
  loadTexture(o, wall);
  insertObject(o,&object_list);   
 
  //~ // Ceiling
  o=newPlane(.05,.75,0,0.55,.55,.8,.75,1,1,0);  
  Scale(o,wallsize,wallsize,1);    
  RotateX(o,-PI/2);
  Translate(o,0,wallsize, 0);
  RotateY(o,offset);
  invert(&o->T[0][0],&o->Tinv[0][0]);   
  loadTexture(o, ceiling);
  insertObject(o,&object_list);     
 
  // Creating Table 
  
  // Felt Table Top
  o=newPlane(.05,.75,.25,0.55,.55,.8,.75,1,1,0); 
  Scale(o,12,12,1);       
  RotateX(o,PI/2);
  RotateY(o,offset);
  Translate(o,0,-3 + baseheight ,10);
  invert(&o->T[0][0],&o->Tinv[0][0]);  
  loadTexture(o, pooltable);
  insertObject(o,&object_list);   
 
 
  //back table border
  o=newCylinder(.05,.95,.15,.35,1,.25,.25,1,1,3);
  Scale(o,.3,.3,12); 
  RotateX(o,PI/2);
  RotateZ(o,-PI/2); 
  Translate(o, 0, 0, 12);
  RotateY(o,offset);
  Translate(o,0,-3 + baseheight ,10);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, border);
  insertObject(o,&object_list);
 
  //left table border
  o=newCylinder(.05,.95,.15,.35,1,.25,.25,1,1,3);
  Scale(o,.3,.3,12); 
  Translate(o, -12, 0, 0);
  RotateY(o,offset );
  Translate(o,0,-3 + baseheight ,10);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, border);
  insertObject(o,&object_list);
 
  //right table border
  o=newCylinder(.05,.95,.15,.35,1,.25,.25,1,1,3);
  Scale(o,.3,.3,12); 
  Translate(o, 12, 0, 0);
  RotateY(o,offset );
  Translate(o,0,-3 + baseheight ,10);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, border);
  insertObject(o,&object_list);
 
  //Creating billiard balls

  // 13 orange
  o=newSphere(.05,.95,.25,1,1,.25,.25,1,1,10);
  Scale(o,0.7,0.7,0.7); 
  RotateX(o,PI/2);
  Translate(o,-3.45,-2.25 + baseheight ,3.5);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, ball13);
  insertObject(o,&object_list);
 
  // // white 
  o=newSphere(.05,.95,.25,1,.95,.95,.95,1,1,10);
  Scale(o,0.7,0.7,0.7); 
  RotateX(o,PI/2);
  RotateY(o,PI);
  Translate(o, 0.45, -2.25 + baseheight, 3.5);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  insertObject(o,&object_list);
 
  //8 black
  o=newSphere(.05,.95,.25,1,1,.25,.25,1,1,10);
  Scale(o,0.7,0.7,0.7); 
  RotateX(o,PI/2);
  RotateY(o,PI);
  Translate(o, 1.45, -2.25 + baseheight, 5);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, ball8);
  insertObject(o,&object_list);
 
  //4 purple
  o=newSphere(.05,.95,.25,1,1,.25,.25,1,1,10);
  Scale(o,0.7,0.7,0.7); 
  RotateX(o,PI/2);
  RotateY(o,PI);
  Translate(o, 2.45, -2.25 + baseheight, 2);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, ball4);
  insertObject(o,&object_list);
 
  //10 blue
  o=newSphere(.05,.95,.25,1,1,.25,.25,1,1,10);
  Scale(o,0.7,0.7,0.7); 
  RotateX(o,PI/2);
  RotateY(o,PI);
  Translate(o, 3.45, -2.25 + baseheight , 4.5);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, ball10);
  insertObject(o,&object_list);

  // Creating the stick 
  o=newCylinder(.05,.95,.15,.35,1,.25,.25,1,1,3);
  Scale(o,.3,.3,15); 
  RotateY(o,-3*PI/4);
  RotateZ(o,-PI/15);
  Translate(o,-10.5,0,-8);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, wood);
  insertObject(o,&object_list);

 
 addAreaLight(3, 3, 0, -1, 0, 0, 10, 3, 20, 20, 1, 1, 0.9, &object_list, &light_list);

}
void buildScene1(void)
{
 // Sets up all objects in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 // Light sources must be defined, positioned, and their colour defined.
 // All objects must be inserted in the object_list. All light sources
 // must be inserted in the light_list.
 //
 // To create hierarchical objects:
 //   Copy the transform matrix from the parent node to the child, and
 //   apply any required transformations afterwards.
 //
 // NOTE: After setting up the transformations for each object, don't
 //       forget to set up the inverse transform matrix!

 struct object3D *o;
 struct pointLS *l;
 struct point3D p;

 // Simple scene for Assignment 3:
 // Insert a couple of objects. A plane and two spheres
 // with some transformations.

 // Let's add a plane
 // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)
 o=newPlane(1,0,0,1,.55,.8,.75,1,1,2);	
 Scale(o,6,6,1);			
 RotateZ(o,PI/1.20);
 RotateX(o,PI/2.25);
 Translate(o,0,-3,10);
 invert(&o->T[0][0],&o->Tinv[0][0]);		
 insertObject(o,&object_list);		

 // Let's add a couple spheres
 o=newSphere(1,0,0,1,1,.25,.25,1,1,6);
 Scale(o,.75,.5,1.5);
 RotateY(o,PI/2);
 Translate(o,-1.45,1.1,3.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(1,0,0,1,.75,.95,.55,1,1,6);
 Scale(o,.5,2.0,1.0);
 RotateZ(o,PI/1.5);
 Translate(o,1.75,1.25,5.0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 // Insert a single point light source.
 p.px=0;
 p.py=15.5;
 p.pz=-5.5;
 p.pw=1;
 l=newPLS(&p,.95,.95,.95);
 insertPLS(l,&light_list);

}

void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
 // This function implements the shading model as described in lecture. It takes
 // - A pointer to the first object intersected by the ray (to get the colour properties)
 // - The coordinates of the intersection point (in world coordinates)
 // - The normal at the point
 // - The ray (needed to determine the reflection direction to use for the global component, as well as for
 //   the Phong specular component)
 // - The current racursion depth
 // - The (a,b) texture coordinates (meaningless unless texture is enabled)
 //
 // Returns:
 // - The colour for this ray (using the col pointer)
 //
 
 struct colourRGB tmp_col, refl_col;  // Accumulator for colour components

 double R,G,B;      // Colour for the object in R G and B

 // This will hold the colour as we process all the components of
 // the Phong illumination model
 tmp_col.R=0;
 tmp_col.G=0;
 tmp_col.B=0;
 
 if (obj->texImg==NULL)   // Not textured, use object colour
 {
  R=obj->col.R;
  G=obj->col.G;
  B=obj->col.B;
 }
 else
 {
  // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
  // for the object. Note that we will use textures also for Photon Mapping.
  obj->textureMap(obj->texImg,a,b,&R,&G,&B);
 }


double shinyness = obj->shinyness;
double ra = obj->alb.ra;
double rd = obj->alb.rd;
double rs = obj->alb.rs;
double rg = obj->alb.rg;
double count = 0;
double totalAmb;
double disToLight, disHit;
 
double *lambda = (double *) malloc (sizeof (double));
struct object3D *objHit;
struct point3D *pHit = newPoint(0, 0, 0);
struct point3D *nHit = newPoint(0, 0, 0);


struct point3D *r = newPoint(dot(&ray->d , n) * n->px, dot(&ray->d , n) * n->py, dot(&ray->d , n) * n->pz);
struct point3D *u, *v ;
r->px = -r->px;
r->py = -r->py;
r->pz = -r->pz;
addVectors(r, r);
addVectors(&ray->d, r);
normalize(r);
struct ray3D *reflectedRay = newRay(p, r);

struct pointLS *currLight = light_list;
struct point3D *direction = newPoint(0, 0, 0);

    // counting number of lights to make scene 
	while (currLight != NULL){
		  count = count + 1;
		  currLight = currLight->next;
	}
	
	u = cross(r, n);
	v = cross(r, u);
	currLight = light_list;
	
  // go through all the light sources and figure out if a shadow is cast then add the correct light contribution. 
  while (currLight != NULL)
  {    

		direction->px = currLight->p0.px;
		direction->py = currLight->p0.py;
		direction->pz = currLight->p0.pz;

		subVectors(p, direction);
		double intensityAmbient = 1;
		disToLight =  length(direction);
		normalize(direction);
		
		struct ray3D *shadowRay = newRay(p, direction);
	  
		double intensityDiffuse = 1;
		double intensitySpecular = 1;

		findFirstHit(shadowRay, lambda, obj, &objHit, pHit, nHit, &a, &b);
		double ambient = ra * intensityAmbient;
		totalAmb += ambient;
		 
		// check if valid hit 
		if(*lambda == DBL_MAX)
		{
			 double specular = rs * pow(max(0, dot(&(shadowRay->d) , r)), shinyness) * intensitySpecular;
			 double diffuse = rd * max(0, dot(n, &shadowRay->d)) * intensityDiffuse;
		   
			 tmp_col.R += (diffuse + specular) * R * currLight->col.R;
			 tmp_col.G += (diffuse + specular) * G * currLight->col.G;
			 tmp_col.B += (diffuse + specular) * B * currLight->col.B;
		}
		else 
		{
			subVectors(p, pHit);
			double disHit =  length(pHit);
				
			if (disHit > disToLight)
			{
				double specular = rs * pow(max(0, dot(&(shadowRay->d) , r)), shinyness) * intensitySpecular;
				double diffuse = rd * max(0, dot(n, &shadowRay->d)) * intensityDiffuse;

				tmp_col.R += (diffuse + specular) * R * currLight->col.R;
				tmp_col.G += (diffuse + specular) * G * currLight->col.G;
				tmp_col.B += (diffuse + specular) * B * currLight->col.B;
		   
			}
		}
		
		currLight = currLight->next;
		free(shadowRay);
  }    
  
  if (depth > 0)
    {
	  // create reflections recursively 	
      rayTrace(reflectedRay, depth-1, &refl_col, obj);

       tmp_col.R +=  rs*refl_col.R;
       tmp_col.G +=  rs*refl_col.G;
       tmp_col.B +=  rs*refl_col.B;

    }
	
  //  divide by count to adjust for scenario if too many lights are added. Ensure scene is not overexposed
  tmp_col.R += R*totalAmb/count;
  tmp_col.G += G*totalAmb/count;
  tmp_col.B += B*totalAmb/count;
  
  col->R = min(tmp_col.R * rg, 1);
  col->G = min(tmp_col.G * rg, 1);
  col->B = min(tmp_col.B * rg, 1);
  
  free(u);
  free(v);
  free(reflectedRay);
  free(lambda);
  free(nHit);
  free(pHit);
 

}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Find the closest intersection between the ray and any objects in the scene.
 // It returns:
 //   - The lambda at the intersection (or < 0 if no intersection)
 //   - The pointer to the object at the intersection (so we can evaluate the colour in the shading function)
 //   - The location of the intersection point (in p)
 //   - The normal at the intersection point (in n)
 //
 // Os is the 'source' object for the ray we are processing, can be NULL, and is used to ensure we don't 
 // return a self-intersection due to numerical errors for recursive raytrace calls.
 //

 /////////////////////////////////////////////////////////////
 // TO DO: Implement this function. See the notes for
 // reference of what to do in here
 /////////////////////////////////////////////////////////////
 
 struct object3D* curr_object = object_list;
 double closest, curr_distance, curr_lambda;
 struct point3D *closest_p, *closest_n;
 struct point3D *curr_p, *curr_n;
 struct point3D distance_test;
 double a_temp, b_temp;
 double closest_a, closest_b;

 curr_p = newPoint(0, 0, 0);
 curr_n = newPoint(0, 0, 0); 
 closest_p = newPoint(0, 0, 0);
 closest_n = newPoint(0, 0, 0);

 closest = DBL_MAX;
 while (curr_object != NULL)
 {
   if (curr_object != Os)
   {
    curr_object->intersect(curr_object, ray, &curr_lambda, curr_p, curr_n, &a_temp, &b_temp);

    if (curr_lambda < DBL_MAX )
    { 
      distance_test = ray->p0;
      subVectors(curr_p, &distance_test);
      curr_distance = length(&distance_test);
      if (curr_distance < closest)
      {

      *obj = curr_object;
      closest = curr_distance;
	  closest_a = a_temp;
	  closest_b = b_temp;
      memcpy(closest_n,curr_n,sizeof(point3D));
      memcpy(closest_p,curr_p,sizeof(point3D));
      }
    }
  }
  curr_object = curr_object->next;
 }
 
 *lambda = closest;
 memcpy(n,closest_n,sizeof(point3D));
 memcpy(p,closest_p,sizeof(point3D));  
 memcpy(a,&closest_a,sizeof(double));
 memcpy(b,&closest_b,sizeof(double));    

 free(closest_n);
 free(closest_p);
 free(curr_n);
 free(curr_p);

}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // Ray-Tracing function. It finds the closest intersection between
 // the ray and any scene objects, calls the shading function to
 // determine the colour at this intersection, and returns the
 // colour.
 //
 // Os is needed for recursive calls to ensure that findFirstHit will
 // not simply return a self-intersection due to numerical
 // errors. For the top level call, Os should be NULL. And thereafter
 // it will correspond to the object from which the recursive
 // ray originates.
 //

 double *lambda;   // Lambda at intersection
 double a,b;    // Texture coordinates
 struct object3D *obj;  // Pointer to object at intersection
 struct point3D *p;  // Intersection point
 struct point3D *n;  // Normal at intersection
 struct colourRGB I;  // Colour returned by shading function

 
 lambda =  (double *) malloc (sizeof (double));

 p = newPoint(0, 0, 0);
 n = newPoint(0, 0, 0);
 findFirstHit(ray, lambda, Os, &obj, p, n, &a, &b);

 if(*lambda < DBL_MAX){

  I.R = obj->col.R;
  I.G = obj->col.G;
  I.B = obj->col.B;

   rtShade(obj, p, n, ray, depth, a, b, &I);
 }
 else
 {
   I.R = background.R;
   I.G = background.G;
   I.B = background.B;
 }
 memcpy(col, &I, sizeof(colourRGB));
 free(n);
 free(p);
 free(lambda);

}

 int main(int argc, char *argv[])
 {
  // Main function for the raytracer. Parses input parameters,
  // sets up the initial blank image, and calls the functions
  // that set up the scene and do the raytracing.
  struct image *im;  // Will hold the raytraced image
  struct view *cam;  // Camera and view for this scene
  int sx;    // Size of the raytraced image
  int antialiasing;  // Flag to determine whether antialiaing is enabled or disabled
  char output_name[1024];  // Name of the output file for the raytraced .ppm image
  struct point3D e;    // Camera view parameters 'e', 'g', and 'up'
  struct point3D g;
  struct point3D up;
  double du, dv;     // Increase along u and v directions for pixel coordinates
  struct point3D pc,d;   // Point structures to keep the coordinates of a pixel and
         // the direction or a ray
  struct ray3D *ray;   // Structure to keep the ray from e to a pixel
  struct colourRGB areaCol;    // Return colour for raytraced pixels
  struct colourRGB col;
  int i,j, ii, jj;     // Counters for pixel coordinates
  unsigned char *rgbIm;
  if (argc<5)
  {
   fprintf(stderr,"RayTracer: Can not parse input parameters\n");
   fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
   fprintf(stderr,"   size = Image size (both along x and y)\n");
   fprintf(stderr,"   rec_depth = Recursion depth\n");
   fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
   fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
   exit(0);
  }
  sx=atoi(argv[1]);
  MAX_DEPTH=atoi(argv[2]);
  if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
  strcpy(&output_name[0],argv[4]);

  fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
  fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
  if (!antialiasing) fprintf(stderr,"Antialising is off\n");
  else fprintf(stderr,"Antialising is on\n");
  fprintf(stderr,"Output file name: %s\n",output_name);

  object_list=NULL;
  light_list=NULL;

  // Allocate memory for the new image
  im=newImage(sx, sx);
  if (!im)
  {
   fprintf(stderr,"Unable to allocate memory for raytraced image\n");
   exit(0);
  }
  else rgbIm=(unsigned char *)im->rgbdata;

  ///////////////////////////////////////////////////
  // TO DO: You will need to implement several of the
  //        functions below. For Assignment 3, you can use
  //        the simple scene already provided. But
  //        for Assignment 4 you need to create your own
  //        *interesting* scene.
  ///////////////////////////////////////////////////
  buildScene();    // Create a scene. This defines all the
       // objects in the world of the raytracer

  //////////////////////////////////////////
  // TO DO: For Assignment 3 you can use the setup
  //        already provided here. For Assignment 4
  //        you may want to move the camera
  //        and change the view parameters
  //        to suit your scene.
  //////////////////////////////////////////

  // Mind the homogeneous coordinate w of all vectors below. DO NOT
  // forget to set it to 1, or you'll get junk out of the
  // geometric transformations later on.

  // Camera center is at (0,0,-1)
  e.px=0;
  e.py=0;
  e.pz=-3;
  e.pw=1;

  // To define the gaze vector, we choose a point 'pc' in the scene that
  // the camera is looking at, and do the vector subtraction pc-e.
  // Here we set up the camera to be looking at the origin, so g=(0,0,0)-(0,0,-1)
  g.px=0;
  g.py=0;
  g.pz=1;
  g.pw=0;

  // Define the 'up' vector to be the Y axis
  up.px=0;
  up.py=1;
  up.pz=0;
  up.pw=0;

  // Set up view with given the above vectors, a 4x4 window,
  // and a focal length of -1 (why? where is the image plane?)
  // Note that the top-left corner of the window is at (-2, 2)
  // in camera coordinates.  
cam=setupView(&e, &g, &up, -3, -1, 1, 2);

  if (cam==NULL)
  {
   fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
   cleanup(object_list,light_list);
   deleteImage(im);
   exit(0);
 }

  // Set up background colour here
  background.R=0.2;
  background.G=0.2;
  background.B=0.2;
  areaCol.R = 0;
  areaCol.G = 0;
  areaCol.B = 0;

  // Do the raytracing
  //////////////////////////////////////////////////////
  // TO DO: You will need code here to do the raytracing
  //        for each pixel in the image. Refer to the
  //        lecture notes, in particular, to the
  //        raytracing pseudocode, for details on what
  //        to do here. Make sure you undersand the
  //        overall procedure of raytracing for a single
  //        pixel.
  //////////////////////////////////////////////////////
  du=cam->wsize/(sx-1);    // du and dv. In the notes in terms of wl and wr, wt and wb,
  dv=-cam->wsize/(sx-1);   // here we use wl, wt, and wsize. du=dv since the image is
         // and dv is negative since y increases downward in pixel
         // coordinates and upward in camera coordinates.

  fprintf(stderr,"View parameters:\n");
  fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
  fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
  printmatrix(cam->C2W);
  fprintf(stderr,"World to camera conversion matrix\n");
  printmatrix(cam->W2C);
  fprintf(stderr,"\n");

   fprintf(stderr,"Rendering row: ");
   //since the origin at 0, can easily find direction without subtracting origin from image point.

   
   double superSampleScale = pow(antialiasing + 1, atoi(argv[3])); // assume power of two 
   fprintf(stderr,"%d", superSampleScale);
   omp_set_dynamic(0);     // Explicitly disable dynamic teams
   omp_set_num_threads(8);
   
   // setup for loop to be paralized in only the outer loop and explicitly. State which variables are shared between threads and which need to be independent here:
   #pragma omp parallel for collapse(2) shared(rgbIm) private(pc) private(ray) private(d) private(i) private(col) private(areaCol)  firstprivate(MAX_DEPTH) num_threads(8) schedule(static)
	
  for (j=0;j<sx;j++)   // For each of the pixels in the image
  {
	
   for (i=0;i<sx;i++)
   {
     col.R = 0;
     col.G = 0;
     col.B = 0;
     areaCol.R = 0;
     areaCol.G = 0;
     areaCol.B = 0;
	 
	 // antialising to create more rays then average the results
     for (int jj = 0;jj<superSampleScale; jj++)
     {
       for  (int ii = 0; ii<superSampleScale; ii++)
       {
        pc.px = 0;
        pc.py = 0;
        pc.pz = 0;
        pc.pw = 1;
        d.px = (-sx/2 + (i + 0.5) + ( -superSampleScale/2 + (ii + 0.5))/superSampleScale)*du;
        d.py = (-sx/2 + (j + 0.5) + ( -superSampleScale/2 + (jj + 0.5))/superSampleScale)*dv;
        d.pz = -1;
        d.pw = 0;

        normalize(&d);
      
        matVecMult(cam->C2W, &d);
        matVecMult(cam->C2W, &pc);
       
        ray = newRay(&pc, &d);
        rayTrace(ray, MAX_DEPTH, &col, NULL);
        free(ray);
        areaCol.R += col.R;
        areaCol.G += col.G;
        areaCol.B += col.B;
       }
     }
     areaCol.R = areaCol.R/(superSampleScale*superSampleScale);
     areaCol.G = areaCol.G/(superSampleScale*superSampleScale);
     areaCol.B = areaCol.B/(superSampleScale*superSampleScale);
	 
	 // assign colors to output
     *(rgbIm+3*(sx*j+i)) = (unsigned char)(areaCol.R *255);
     *(rgbIm+(3*(sx*j+i)+1)) = (unsigned char)(areaCol.G*255);
     *(rgbIm+(3*(sx*j+i)+2)) = (unsigned char)(areaCol.B*255);
	 
     areaCol.R = 0;
     areaCol.G = 0;
     areaCol.B = 0;
   } 

  } 

  imageOutput(im,output_name);

  // Exit section. Clean up and return.
  cleanup(object_list,light_list);   // Object and light lists
  deleteImage(im);       // Rendered image
  free(cam);         // camera view
  exit(0);
 }
