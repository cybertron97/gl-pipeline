#include "driver_state.h"
#include <cstring>
 #include <vector>
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
    state.image_depth= new float[width * height];
    state.image_color = new pixel[width * height];
   for (int i = 0; i < width * height; ++i)
   {
      state.image_color[i] = make_pixel(0, 0, 0);
      state.image_depth[i]=1.0;
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

  case render_type:: triangle:
  {

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
      for (int j =0 ; j<state.num_vertices*state.floats_per_vertex; j += 3*state.floats_per_vertex)
        {
          for (int i = 0 ; i< 3; i++)
          {
            ver[i].data= &state.vertex_data[(state.floats_per_vertex * i) + j];
            arr[i].data=ver[i].data;
            state.vertex_shader(ver[i], arr[i], state.uniform_data);
            g[i] = &arr[i];
          }
          // rasterize_triangle(state, g);
          clip_triangle(state,g,0);
        }

   }
   break;
   case render_type::indexed:
   {
     data_geometry arr[3];
     data_vertex ver[3];
     const data_geometry *g[3];
     for (int j =0 ; j<3*state.num_triangles; j += 3)
        {
          for (int i = 0 ; i< 3; i++)
            {
              ver[i].data= &state.vertex_data[state.index_data[i+j]*state.floats_per_vertex];
              arr[i].data=ver[i].data;
              state.vertex_shader(ver[i], arr[i], state.uniform_data);
              g[i] = &arr[i];
            }
            // rasterize_triangle(state, g);
           clip_triangle(state,g,0);
          }

    }
    break;
   case render_type::strip:
   {
     data_geometry arr[3];
     data_vertex ver[3];
     const data_geometry *g[3];
      for (int j =0 ; j<(state.num_vertices-2); j++)
        {
          for (int i = 0 ; i< 3; i++)
           {
             ver[i].data= &state.vertex_data[(state.floats_per_vertex) * (i + j)];
             arr[i].data=ver[i].data;
             state.vertex_shader(ver[i], arr[i], state.uniform_data);
             g[i] = &arr[i];
           }
           // rasterize_triangle(state, g);
            clip_triangle(state,g,0);
        }
   }
  break;
    case render_type::fan:
    {
      data_geometry arr[3];
      data_vertex ver[3];
      const data_geometry *g[3];
      for (int j =0 ; j<state.num_vertices; j ++)
        {
          for (int i = 0 ; i< 3; i++)
            {
              // ver[i].data= &state.vertex_data [(state.floats_per_vertex * i) + j];
              int a = ((i+j)*state.floats_per_vertex);
              if (i==0)
                {
                  a = 0;
                }
              ver[i].data= state.vertex_data + a;
              arr[i].data=ver[i].data;
              state.vertex_shader(ver[i], arr[i], state.uniform_data);
              g[i] = &arr[i];

            }
            // rasterize_triangle(state, g);
            clip_triangle(state,g,0);
        }
      }break;

	//	break;
    default: std::cout<<"This is hard! Send help bros!"<<std::endl;
  }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
  data_geometry arr[3];
  data_geometry o[3];//lol i was running out of alphabets
  data_vertex ver[3];
  const data_geometry *g[3];
    for (int i = 0;i<3;i++)
      {
        // ver[i].data = in[i]-> data;
        // arr[i].data=ver[i].data;
        arr[i].gl_Position= in[i]->gl_Position;
        // if(face==0)
        //   {
        //    // state.vertex_shader(ver[i], arr[i], state.uniform_data);
        //   }
        //g[i] = &arr[i];
      }
      //in[0]->gl_Position[0];
     // auto a = in[0]->gl_Position[3];
     // auto b = in[1]->gl_Position[3];
     // auto c = in[2]->gl_Position[3];
      // vec4 a = {arr[0].gl_Position[0], arr[0].gl_Position[1],arr[0].gl_Position[2], arr[0].gl_Position[3]};
      // vec4 b = {arr[1].gl_Position[0], arr[1].gl_Position[1],arr[1].gl_Position[2], arr[1].gl_Position[3]};
      // vec4 c = {arr[2].gl_Position[0], arr[2].gl_Position[1],arr[2].gl_Position[2], arr[2].gl_Position[3]};
      // // vec3 p,n = new vec3(x,y,z);
  int  position =0, value =0;

    if(face==6)
    {

        rasterize_triangle(state, in);
        return;
    }
      else if (face==0 )
    {
      position = 0;
      value =1;
    }
      else if (face ==1)// set the position values
    {
        position =0;
        value =-1;
    }
    else if (face == 2)
    {
      position =1;
      value =1;
    }
    else if (face == 3 )
    {
      position =1;
      value =-1;
    }
    else if (face == 4 )
    {
      position =2;
      value =1;
    }
    else if (face == 5 )
    {
      position =2;
       value =-1;
    }
    //going to try piazza
vec3 tri = {0,0,0};
int count =0;
for (int i = 0; i < 3; i++)
{
    if (value == -1)
    {
      if(in[i]->gl_Position[position] >= -in[i]->gl_Position[3])
        {
        count++;
        tri[i] = 1;
        }
     }
    else if(in[i]->gl_Position[position] <=  in[i]->gl_Position[3])
     {
      count++;
      tri[i] = 1;
      }
  }

if (count ==3 )
{
//everything inside
  clip_triangle(state, in ,face+1);
  return;
}
if (count ==0 )
{
  //everything is outside, ..... so no cares
  return;
}

const data_geometry* a =0;
const data_geometry* b =0;
const data_geometry* c = 0;
for (int i =0; i<3; i++)
{
  if ((count ==1 && tri[i]==1 ) || (count ==2 && tri[i]==0))
   {
//setting the triangles vertices;
     a = in[i];
     b = in[(i+1)%3];
     c = in[(i+2)%3];
   }
}
float position_A = a->gl_Position[3];
float position_B = b->gl_Position[3];
float position_C = c->gl_Position[3];
float beta = ((value * position_B) - b->gl_Position[position]) / ((a->gl_Position[position] - (value * position_A)) + ((value * position_B) - b->gl_Position[position]));
float gamma = ((value * position_C) - c->gl_Position[position]) / ((a->gl_Position[position] - (value * position_A)) + ((value * position_C) - c->gl_Position[position]));

data_geometry* ba = new data_geometry();
data_geometry* ca = new data_geometry();
 float BA_data[state.floats_per_vertex];
ba->data= BA_data;
float CA_data[state.floats_per_vertex];
ca->data= CA_data;
for (int j =0; j<state.floats_per_vertex;  j++)
{
  switch(state.interp_rules[j])
  {
    case interp_type::flat:
    {
      ba->data[j]=a->data[j];
      ca->data[j]=a->data[j];
    }break;
    case interp_type::smooth:
    {
      ba->data[j]= beta* a->data[j]+ ((1-beta)*b->data[j]);
      ca->data[j]= gamma* a->data[j]+((1-gamma)*c->data[j]);

    }break;
    case interp_type::noperspective:
    {
      float  k = (beta*position_A)+((1-beta)*position_B);
      float alpha_prime = (beta * position_A)/k;
      ba->data[j]= (alpha_prime * a->data[j])+((1-alpha_prime )*b->data[j]);
      float k_gamma = (gamma *position_A)+((1-gamma)*position_C);
      float gamma_prime = (gamma*position_A)/k_gamma;
      ca->data[j]= (gamma_prime*a->data[j])+((1-gamma_prime)*c->data[j]);
    }break;
    default : std::cout<<"GOD!HELP ME!!!!!!"<<std::endl;

  }

}
for (int i = 0;i<4; i++)
{
  ba->gl_Position[i]= (beta * a->gl_Position[i]) + ((1-beta)*b->gl_Position[i]);
  ca->gl_Position[i]= (gamma * a->gl_Position[i]) + ((1-gamma)*c->gl_Position[i]);

}


if (count ==1)
{
  g[0]=a;
  g[1]=ba;
  g[2]=ca;
  clip_triangle(state,g,face+1);

  return;
}
if (count ==2 )
{
  g[0]=ba;
  g[1]=b;
  g[2]=c;
  clip_triangle(state,g,face+1);
  g[0]=c;
  g[1]=ca;
  g[2]=ba;
  clip_triangle(state,g,face+1);
  return;

}



//plane's points to be {0,0, 1}
//clipping on one side
// if (tri[0]==0&&tri[1]==0&&tri[2]==1)
// {
//   float  alpha = (((value * a)-in[0]->gl_Position[index])/(in[2]->gl_Position[index]-(value*c)+(value * a)-in[0]->gl_Position[index]));
//   float beta = (((value *b)-in[1]->gl_Position[index])/(in[2]->gl_Position[index]-(value *c)+(value * b)-in[1]->gl_Position[index]));
// 	vec4 Wa=  (alpha*in[2]->gl_Position)+(1-alpha)*in[0]->gl_Position;
//   vec4 Wb = (beta*in[2]->gl_Position)+(1-beta)*in[1]->gl_Position;
//   o[0].gl_Position = Wa;
//   o[1].gl_Position = Wb;
//   o[2].gl_Position = in[2]->gl_Position;
//   g[0]= &o[0];
//   g[1]=&o[1];
//   g[2]=&o[2];
//   clip_triangle(state,g,face+1);
// 	return;
//
// }
// else if (tri[0]==1&& tri[1]==0&& tri[2]==1)
// {
//   float alpha = (((value * b) - in[1]->gl_Position[index])/(in[0]->gl_Position[index] - (value * a) + (value * b) - in[1]->gl_Position[index]));
//   float gamma = (((value * b) - in[1]->gl_Position[index])/(in[2]->gl_Position[index] - (value * c) + (value * b) - in[1]->gl_Position[index]));
//   vec4 Wa=  (alpha*in[0]->gl_Position)+(1-alpha)*in[1]->gl_Position;
//   vec4 Wc = (gamma*in[2]->gl_Position)+(1-gamma)*in[1]->gl_Position;
//   o[0].gl_Position = in[0]->gl_Position;
//   o[1].gl_Position = Wa;
//   o[2].gl_Position = Wc;
//   g[0]= &o[0];
//   g[1]=&o[1];
//   g[2]=&o[2];
//   clip_triangle(state,g,face+1);
//   o[0].gl_Position = in[0]->gl_Position;
//   o[1].gl_Position = Wc;
//   o[2].gl_Position = in[2]->gl_Position;
//   g[0]= &o[0];
//   g[1]=&o[1];
//   g[2]=&o[2];
//   clip_triangle(state,g,face+1);
//   return;
// }
// else if (tri[0]==1&&tri[1]==1&&tri[2]==0)
// {
//   float  alpha = (((value * c)-in[2]->gl_Position[index])/(in[0]->gl_Position[index]-(value*a)+(value * c)-in[2]->gl_Position[index]));
//   float beta = ((value * a)-in[1]->gl_Position[index])/(in[2]->gl_Position[index]-(value *c)+(value * b)-in[1]->gl_Position[index]);
// 	vec4 Wa=  (alpha*in[0]->gl_Position)+(1-alpha)*in[2]->gl_Position;
//   vec4 Wb = (beta*in[2]->gl_Position)+(1-beta)*in[1]->gl_Position;
//   o[0].gl_Position = in[0]->gl_Position;
//   o[1].gl_Position = Wa;
//   o[2].gl_Position = Wb;
//   g[0]= &o[0];
//   g[1]=&o[1];
//   g[2]=&o[2];
// clip_triangle(state,g,face+1);
// o[0].gl_Position = in[0]->gl_Position;
// o[1].gl_Position = in[1]->gl_Position;
// o[2].gl_Position = Wb;
// g[0]= &o[0];
// g[1]=&o[1];
// g[2]=&o[2];
// clip_triangle(state,g,face+1);
// return;
// }
// else if(tri[0]==0 && tri[1]==1 && tri[2]==0)
// {
//     float alpha = (((value * a) - in[0]->gl_Position[index])/(in[1]->gl_Position[index] - (value * b) + (value * a) - in[0]->gl_Position[index]));
// 		float gamma = (((value * c) - in[2]->gl_Position[index])/(in[1]->gl_Position[index] - (value * b) + (value * c) - in[2]->gl_Position[index]));
//     vec4 Wa=  (alpha*in[1]->gl_Position)+(1-alpha)*in[0]->gl_Position;
//     vec4 Wc = (gamma*in[1]->gl_Position)+(1-gamma)*in[2]->gl_Position;
//     o[0].gl_Position = Wa;
//     o[1].gl_Position = in[1]->gl_Position;
//     o[2].gl_Position = Wc;
//     g[0]= &o[0];
//     g[1]=&o[1];
//     g[2]=&o[2];
//   clip_triangle(state,g,face+1);
//   return;
// }
// else if (tri[0]==0 && tri[1]==1&& tri[2]==1)
// {
//   float beta = (((value * a) - in[0]->gl_Position[index])/(in[1]->gl_Position[index] - (value * b) + (value * a) - in[0]->gl_Position[index]));
//   float gamma  = (((value * c) - in[2]->gl_Position[index])/(in[0]->gl_Position[index] - (value * a) + (value * c) - in[2]->gl_Position[index]));
//   vec4 Wb=  (beta*in[1]->gl_Position)+(1-beta)*in[0]->gl_Position;
//   vec4 Wc = (gamma*in[0]->gl_Position)+(1-gamma)*in[2]->gl_Position;
//   o[0].gl_Position = Wc;
//   o[1].gl_Position = in[1]->gl_Position;
//   o[2].gl_Position = Wb;
//   g[0]= &o[0];
//   g[1]=&o[1];
//   g[2]=&o[2];
//   clip_triangle(state,g,face+1);
//   o[0].gl_Position = Wc;
//   o[1].gl_Position = in[1]->gl_Position;
//   o[2].gl_Position = in[2]->gl_Position;
//   g[0]= &o[0];
//   g[1]=&o[1];
//   g[2]=&o[2];
//   clip_triangle(state,g,face+1);
//   return;
// }
// else if (tri[0]==1&&tri[1]==0&&tri[2]==0)
// {
//   float beta = (((value * b) - in[1]->gl_Position[index])/(in[0]->gl_Position[index] - (value * a) + (value * b) - in[1]->gl_Position[index]));
// 	float gamma = (((value * c) - in[2]->gl_Position[index])/(in[0]->gl_Position[index] - (value * a) + (value  * c) - in[2]->gl_Position[index]));
//   vec4 Wb=  (beta*in[0]->gl_Position)+(1-beta)*in[1]->gl_Position;
//   vec4 Wc = (gamma*in[0]->gl_Position)+(1-gamma)*in[2]->gl_Position;
//   o[0].gl_Position = in[0]->gl_Position;
//   o[1].gl_Position = Wb;
//   o[2].gl_Position = Wc;
//   g[0]= &o[0];
//   g[1]=&o[1];
//   g[2]=&o[2];
//   clip_triangle(state,g,face+1);
//   return;

// }


  //  std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
  //  clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
  auto a = in[0]->gl_Position/in[0]->gl_Position[3];
  auto b = in[1]->gl_Position/in[1]->gl_Position[3];
  auto c = in[2]->gl_Position/in[2]->gl_Position[3];
  //  std::cout<<"TODO: implement rasterization"<<std::endl;
  float ax = ((state.image_width/2) * a[0]);
  float ay = ((state.image_height/2) * a[1]);

  float bx = ((state.image_width/2) * b[0]);
  float by = ((state.image_height/2) * b[1]);

  float cx = ((state.image_width/2) * c[0]);
  float cy = ((state.image_height/2) * c[1]);

  float i = ((state.image_width/2) - 0.5);
  float j = ((state.image_height/2) - 0.5);


  float ai = ax + i;
  float aj = ay  + j;

  float bi =  bx + i ;
  float bj = by + j;

  float ci =  cx + i;
  float cj = cy + j;

  float area_Triangle = 0.5 * (((bi * cj) - (ci * bj)) - ((ai * cj) - (ci * aj)) + ((ai * bj) - (bi * aj)));

      int minx = min(min(ai, bi),ci);
      int maxx =  max(max(ai,bi),ci);
      int miny = min(min(aj,bj),cj);
      int maxy = max(max(aj,bj),cj);

  if(minx < 0)
       minx = 0;
  if(miny < 0)
       miny = 0;
  if(maxx > state.image_width)
       maxx = state.image_width - 1;
  if(maxy > state.image_height)
       maxy = state.image_height - 1;


   for (int i = minx; i <= maxx; i++)
   {

     for (int j=miny; j<=maxy ; j++)
      {
       float alpha = 0.5 * (((bi * cj)-(ci * bj))-((i * cj) - (ci * j))+((i * bj)-(bi * j)));
       float beta = 0.5 * (((i * cj)-(ci * j))-((ai * cj) - (ci * aj))+((ai * j)-(i * aj)));
       float gamma = 0.5 * (((bi * j)-(i * bj))-((ai * j)-(i * aj))+((ai * bj) - (bi * aj)));


       alpha = alpha/area_Triangle;
       beta = beta/area_Triangle;
       gamma= gamma/area_Triangle;

       if (alpha >= -0.001 && beta>= -0.001 && gamma >= -0.001)
            {
                data_fragment frag;
                float frag_data[state.floats_per_vertex];
                       //float no_perspective [MAX_FLOATS_PER_VERTEX];
                frag.data = frag_data;

                  for(int j = 0; j < state.floats_per_vertex; j++)
                      {
                         switch(state.interp_rules[j])
                         {
                           case interp_type::flat:
                                {
                                  frag.data[j]= in[0]->data[j];
                                }
                           break;
                           case interp_type::smooth:
                                {
                                  float  k = (alpha/in[0]->gl_Position[3])+(beta/in[1]->gl_Position[3])
                                          +(gamma/in[2]->gl_Position[3]);
                              //alpha, beta , gamma prime values
                                  float alpha1 = alpha/ (k * in[0]->gl_Position[3]);
                                  float beta1 = beta/(k * in[1]->gl_Position[3]);
                                  float gamma1  = gamma/(k * in[2]->gl_Position[3]);
                                  frag.data[j]= alpha1 * in[0]->data[j] + beta1*in[1]->data[j]+ gamma1*in[2]->data[j];
                                }
                           break;
                           case interp_type::noperspective:
                                  {
                                  // no_perspective[j]=(alpha * in[0]->data[j]) + (beta * in[1]->data[j])+(gamma * in[2]->data[j]);
                                  // frag.data =no_perspective;
                                    frag.data[j] =(alpha * in[0]->data[j]) + (beta * in[1]->data[j])+(gamma * in[2]->data[j]);

                                  }
                                  break;
                           case interp_type::invalid:
                           break;
                         }
                       }

                        //call fragment shader
                       data_output dout;
                       state.fragment_shader(frag, dout,state.uniform_data);

                       int red = 255 * dout.output_color[0];
                       int green = 255 * dout.output_color[1];
                       int blue = 255 * dout.output_color[2];

                         float z = alpha*a[2]+beta*b[2]+gamma*c[2];
                         int index = (state.image_width * j) + i;

                      if ( z<=state.image_depth[index] )
                            {
                              state.image_depth[index] = z;
                              state.image_color[index] = make_pixel(red,green,blue);
                            }


              // int index = state.image_width  * j + i;
              // state.image_color[index] = make_pixel(255,255,255);

             }
       }
   }
}
