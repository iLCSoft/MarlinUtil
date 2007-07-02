/* "C" event display.
 * Client side elements definitions.
 *
 * Alexey Zhelezov, DESY/ITEP, 2005 */
#include <string.h>

#include <ced_cli.h>
#include <ced.h>


/*
 * Hit element
 */

static unsigned HIT_ID=0;

void ced_hit(float x,float y,float z,unsigned type,unsigned size,unsigned color){
  CED_Hit *h=(CED_Hit *)ced_add(HIT_ID);
  if(!h)
    return;
  h->p.x=x;
  h->p.y=y;
  h->p.z=z;
  h->type=type;
  h->size=size;
  h->color=color;
}

/*
 * Line element
 */

static unsigned LINE_ID=0;

void ced_line(float x0,float y0,float z0,
	      float x1,float y1,float z1,
	      unsigned type, unsigned width,unsigned color){
  CED_Line *l=(CED_Line *)ced_add(LINE_ID);
  if(!l)
    return;
  l->p0.x=x0;
  l->p0.y=y0;
  l->p0.z=z0;
  l->p1.x=x1;
  l->p1.y=y1;
  l->p1.z=z1;
  l->type=type;
  l->width=width;
  l->color=color;
}

/*
 * GeoCylinder
 */
static unsigned GEOC_ID=0;

void ced_geocylinder(float d,unsigned sides,float rotate,float z,float shift,
		     unsigned color){
  CED_GeoCylinder *c=(CED_GeoCylinder *)ced_add(GEOC_ID);
  if(!c)
    return;
  c->d=d;
  c->sides=sides;
  c->rotate=rotate;
  c->z=z;
  c->shift=shift;
  c->color=color;
}
 
void ced_geocylinders(unsigned n,CED_GeoCylinder *all){
  CED_GeoCylinder *c;
  unsigned i;
  for(i=0;i<n;i++){
    c=(CED_GeoCylinder *)ced_add(GEOC_ID);
    if(!c)
      return;
    memcpy(c,all+i,sizeof(CED_GeoCylinder));
  }
}


static unsigned GEOB_ID=0;

void ced_geobox(double * sizes, double * center, unsigned int color ) {
  int iDim;
  CED_GeoBox * box = (CED_GeoBox*) ced_add(GEOB_ID);
  if ( ! box ) return;
  for ( iDim = 0; iDim < 3; iDim ++ ) {
    box->sizes[iDim]   = sizes[iDim];
    box->center[iDim]  = center[iDim];
  }
  box->color   = color;

}

void ced_geoboxes(unsigned int nBox, CED_GeoBox * allBoxes ) {
  
  CED_GeoBox * box;
  unsigned int iBox;
  for ( iBox = 0; iBox < nBox ; iBox++ ) {
    box = (CED_GeoBox *) ced_add(GEOB_ID);
    if ( ! box ) return;
    memcpy( box, allBoxes + iBox, sizeof(CED_GeoBox) );
  }
}

void ced_register_elements(void){
  GEOC_ID  =ced_register_element(sizeof(CED_GeoCylinder),0);
  LINE_ID  =ced_register_element(sizeof(CED_Line),0);
  HIT_ID   =ced_register_element(sizeof(CED_Hit),0);
  GEOB_ID  =ced_register_element(sizeof(CED_GeoBox), 0);
}

