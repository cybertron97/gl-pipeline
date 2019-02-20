#include "driver_state.h"
#include <cstring>
using namespace std;

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    state.image_color = new pixel[width * height];
   for (int i = 0; i < width * height; ++i)
   {
      state.image_color[i] = make_pixel(0, 0, 0);
   }

    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
  //  std::cout<<"TODO: implement rendering."<<std::endl;

switch (type) {

  case render_type:: triangle:{

//  for(int v = 0; v < state.num_vertices; v +=3)
//   {
//     //data_geometry const* arr = new arr[3];
//
// //     some_type const* const_some_array = some_array;
// // f(&const_some_array);
//    //data_geometry** arr = new data_geometry*[3];
//


data_geometry arr[3];
data_vertex ver[3];
const data_geometry *g[3];
for (int i =0 ; i<3 ; ++i)
{
  for (int j = 0 ; j< state.num_vertices*state.floats_per_vertex; j += 3*state.floats_per_vertex)
  {
    ver[i].data= &state.vertex_data[(state.floats_per_vertex * i) + j];
  }
  state.vertex_shader(ver[i], arr[i], state.uniform_data);
  g[i] = &arr[i];
}
rasterize_triangle(state, g);
}

break;




  case render_type::indexed: break;
  case render_type::strip: break ;
  case render_type::fan: break;

	//	break;
    default: std::cout<<"This is hard! Send help bros!"<<std::endl;
}}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
  //  std::cout<<"TODO: implement rasterization"<<std::endl;
  float ax = ((state.image_width/2) * in[0]->gl_Position[0]);
  float ay = ((state.image_height/2) * in[0]->gl_Position[1]);

  float bx = ((state.image_width/2) * in[1]->gl_Position[0]);
  float by = ((state.image_height/2) * in[1]->gl_Position[1]);

  float cx = ((state.image_width/2) * in[2]->gl_Position[0]);
  float cy = ((state.image_height/2) * in[2]->gl_Position[1]);

  float i = ((state.image_width/2) - 0.5);
  float j = ((state.image_height/2) - 0.5);


  float ai = ax + i;
  float aj = ay  + j;

  float bi =  bx + i ;
  float bj = by + j;

  float ci =  cx + i;
  float cj = cy + j;

  float area_Triangle = 0.5 * (((bi * cj) - (ci * bj)) - ((ai * cj) - (ci * aj)) + ((ai * bj) - (bi * aj)));


   for (int i = min(min(ai,bi),ci); i < max(max(ai,bi),ci); i++)
   {

     for (int j=min(min(aj,bj),cj); j<max(max(aj,bj),cj) ; j++)
 {
       float alpha = 0.5 * (((bi * cj)-(ci * bj))-((i * cj) - (ci * j))+((i * bj)-(bi * j)));
       float beta = 0.5 * (((i * cj)-(ci * j))-((ai * cj) - (ci * aj))+((ai * j)-(i * aj)));
       float gamma = 0.5 * (((bi * j)-(i * bj))-((ai * j)-(i * aj))+((ai * bj) - (bi * aj)));


   alpha = alpha/area_Triangle;
   beta = beta/area_Triangle;
   gamma= gamma/area_Triangle;

   if (alpha >= 0 && beta>= 0 && gamma <= 1)
   {

     int index = (state.image_width  * j) + i;
     state.image_color[index] = make_pixel(255,255,255);
   }
 }
   }
}
