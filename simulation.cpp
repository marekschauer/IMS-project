#include <iostream>
#include <vector>
#include "bitmap_image.hpp"

#define IMAGE_WIDTH   (700)
#define IMAGE_HEIGHT  (700)



typedef struct material_struct {
 float structure_combustibility;
 float spread_ability_time;
 float extinction_time;
} material;

class Cell {

public:
  float spread_wind_dir;
  material material_properties;
  float ratio_of_area;
  float spread_probability;
  unsigned int state;
  unsigned int x;
  unsigned int y;
  // Cell();
  // ~Cell();
  
};

class Automata {
  // unsigned int water_utilization_capacity;
  // unsigned int nozzle_volume;
  // unsigned int spraying_time;
  // unsigned int nozzle_range;
  // unsigned int arrive_time;
  std::vector<std::vector<Cell*> > next_grid;



public:
  std::vector<std::vector<Cell*> > actual_grid;
  /**
   * Rows of CA lattice
   */
  unsigned int rows;
    /**
     * Columns of CA lattice
     */
  unsigned int cols;
    /**
     * Overall simulation time in minutes
     */
  unsigned int simulation_time;
    /**
     * Minutes to be passed by every step of simulation
     */
  unsigned int time_step;
    /**
     * Actual time of simulation (in minutes)
     */
  unsigned int actual_time;
    /**
     * Wind velocity in meters per second
     */
  unsigned int wind_velocity;
    /**
     * Number representing wind direction 
     */
  unsigned int wind_direction;
  unsigned int alpha;
  unsigned int beta;
  Automata (unsigned int, unsigned int);
  void process_fire();
  // int area () {return (width*height);}
};

/**********************************************
**********************************************
**********************************************
***************IMAGE PROCESSING START*********
**********************************************
**********************************************
***/
double x_coord(double x) {
 return x-(IMAGE_WIDTH/2);
}

double y_coord(double y) {

 return (IMAGE_HEIGHT/2)-y;
}

// TODO - do nazvu suboru pridavat cas simulacie
void print_ca(Automata* ca) {
  cartesian_canvas canvas(700,700);

  canvas.image().clear(255);

  canvas.pen_width(1);

  ::srand(0xA5A5A5A5);
  float cell_width = IMAGE_WIDTH/ca->cols;
  for (unsigned int i = 0; i < ca->rows; i++) {
    for (unsigned int j = 0; j < ca->cols; j++) {
      int state = ca->actual_grid[i][j]->state;
      switch(state) {
        case 0:
          canvas.pen_color(255,255,255);
          break; //optional
        case 1  :
          canvas.pen_color(200,200,200);
          break; //optional
        case 2  :
          canvas.pen_color(255, 191, 0);
          break; //optional
        case 3  :
          canvas.pen_color(255, 0, 0);
          break; //optional
      
       // you can have any number of case statements.
        default : //Optional
          canvas.pen_color(0,0,0);
      }
      std::cout << i << std::endl;
      canvas.fill_rectangle(x_coord((IMAGE_WIDTH/ca->cols)*j),y_coord((IMAGE_HEIGHT/ca->rows)*i),x_coord((IMAGE_WIDTH/ca->cols)*j+IMAGE_WIDTH/ca->cols),y_coord((IMAGE_HEIGHT/ca->rows)*i+(IMAGE_HEIGHT/ca->rows)));
      // canvas.fill_rectangle(x_coord((IMAGE_WIDTH/ca->cols)*i),y_coord((IMAGE_HEIGHT/ca->rows)*j),x_coord((IMAGE_WIDTH/ca->cols)*i+IMAGE_WIDTH/ca->cols),y_coord((IMAGE_HEIGHT/ca->rows)*j+(IMAGE_HEIGHT/ca->rows)));
      canvas.pen_color(0,0,0);
      canvas.rectangle(x_coord((IMAGE_WIDTH/ca->cols)*j),y_coord((IMAGE_HEIGHT/ca->rows)*i),x_coord((IMAGE_WIDTH/ca->cols)*j+IMAGE_WIDTH/ca->cols),y_coord((IMAGE_HEIGHT/ca->rows)*i+(IMAGE_HEIGHT/ca->rows)));
     }
   }
 canvas.image().save_image("moj_obrazok2.bmp");

}

/**********************************************
**********************************************
**********************************************
***************IMAGE PROCESSING END***********
**********************************************
**********************************************
***/


/*
 *void Box::setLength( double len ) {
   length = len;
}
*/

Automata::Automata(unsigned int r, unsigned int c) {
  rows = r;
  cols = c;
  actual_time = 0;
  time_step = 3;
  std::vector<Cell*> tmp_vec;

  for (unsigned int i = 0; i < r; ++i) {
    for (unsigned int j = 0; j < c; ++j) {
      Cell* to_be_added = new Cell;
      to_be_added->x = i;
      to_be_added->y = j;
      to_be_added->state = 1; // TODO - put it to define
      tmp_vec.push_back(to_be_added);
    }
    actual_grid.push_back(tmp_vec);
    tmp_vec.clear();
  }

  // for (int i = 0; i < r; ++i) {
  //   for (int j = 0; j < c; ++j) {
  //     std::cout << "|" << actual_grid[i][j]->x << "," << actual_grid[i][j]->y << "|";
  //   }
  //   std::cout << std::endl;
  // }
}

void Automata::process_fire() {
  actual_time += time_step;
  // this vector holds the spreading cells
  std::vector<Cell*> spreading_cells;
  // filling spreading_cells vector
  for (unsigned int i = 0; i < rows; ++i) {
    for (unsigned int j = 0; j < cols; ++j) {
      if (actual_grid[i][j]->state == 3) { // TODO - put 3 to define
        spreading_cells.push_back(actual_grid[i][j]);
      }
    }
  }

  
}

int main () {
  Automata* ca = new Automata(15,15);

  ca->actual_grid[9][12]->state = 3; // TODO - put 3 to define
  ca->actual_grid[0][0]->state = 3; // TODO - put 3 to define
  ca->actual_grid[0][12]->state = 3; // TODO - put 3 to define
  ca->time_step = 5;
  ca->simulation_time = 100;
  ca->wind_velocity = 1;
  // wind direction - TODO
  ca->wind_direction = 1;
  ca->alpha = 1;
  ca->beta = 1;


  print_ca(ca);
  ca->process_fire();



  return 0;
}
