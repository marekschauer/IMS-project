/**
 * TODOs:
 *   - pridat podporu Wij pre juhovychod, juhozapad, severovychod a severozapad
 *   - vyriesit otacanie mapky pre Wij nejak automaticky
 */

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <random>
#include "bitmap_image.hpp"
#include "json.hpp"

#define IMAGE_WIDTH   (1500)
#define IMAGE_HEIGHT  (1500)

typedef struct material_struct {
 float structure_combustibility;
 // cas, kedy sa bunka rozhori tak, aby vedela sirit ohen (cas, za ktory bunka zmeni svoj stav z 2 do 3)
 float spread_ability_time;
 // cas do vyhorenia (cas, za ktory bunka zmeni svoj stav z 3 do 4)
 float extinction_time;
} material;

typedef enum {
 north,
 south,
 east,
 west,
 northeast,
 northwest,
 southeast,
 southwest
} WindDirection;


typedef struct automata_parameters_struct {
  unsigned int cell_size;
  unsigned int rows;
  unsigned int cols;
  unsigned int simulation_time;
  unsigned int time_step;
  unsigned int actual_time;
  float wind_velocity;
  WindDirection wind_direction;
  unsigned int alpha;
  unsigned int beta;
} automata_parameters;

class Cell {

public:
  /**
   * Wij for equation of Fij
   */
  float spread_wind_dir;
  material material_properties;
  float ratio_of_area;
  /**
   * Fij
   */
  float spread_probability;
  unsigned int state;
  unsigned int x;
  unsigned int y;
  unsigned int burning_time;
  Cell* copy();
  // Cell();
  // ~Cell();
  
};

class Automata {
  // unsigned int water_utilization_capacity;
  // unsigned int nozzle_volume;
  // unsigned int spraying_time;
  // unsigned int nozzle_range;
  // unsigned int arrive_time;



public:
  std::vector<std::vector<Cell*> > actual_grid;
  std::vector<std::vector<Cell*> > next_grid;
  
  /**
   * Cell size in meters
   */
  unsigned int cell_size;
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
  float wind_velocity;
    /**
     * Number representing wind direction 
     */
  WindDirection wind_direction;
  unsigned int alpha;
  unsigned int beta;
  Automata (automata_parameters);
  void process_fire(int, int);
  float getSpreadingProbability(int, int);
  std::vector<std::vector<Cell*>> getNeighbourhood(unsigned int, unsigned int);
  std::vector<Cell*> getSpreadingCells();
  std::vector<Cell*> getBurningCells();
  bool isInLattice(unsigned int, unsigned int);
  void simulate();
  void actualize_grid();
};

/**********************************************
**********************************************
**********************************************
***************IMAGE PROCESSING START*********
**********************************************
**********************************************
***/
double x_coord(int width, double x) {
 return x-(width/2);
}

double y_coord(int height, double y) {

 return (height/2)-y;
}

// TODO - do nazvu suboru pridavat cas simulacie
void print_ca(Automata* ca, int time) {
  int image_width = ca->cols*5;
  int image_height = ca->rows*5;

  cartesian_canvas canvas(image_width,image_height);

  canvas.image().clear(255);

  canvas.pen_width(1);

  ::srand(0xA5A5A5A5);

  for (unsigned int j = 0; j < ca->cols; j++) {
    for (unsigned int i = 0; i < ca->rows; i++) {
      int state = ca->actual_grid[j][i]->state;
      float compustibility = ca->actual_grid[j][i]->material_properties.structure_combustibility;
      if (compustibility == 0.6) {
        canvas.pen_color(45,24,13);
      } else if (compustibility == 1.0) {
        canvas.pen_color(35,35,0);
      } else if (compustibility == 0) {
        canvas.pen_color(255,255,255);
      }
      switch(state) {
        // case 0:
        // canvas.pen_color(255,255,255);
          // break; //optional
          // case 1  :
          // canvas.pen_color(200,200,200);
          // break; //optional
          case 2  :
          canvas.pen_color(255, 191, 0);
          break; //optional
          case 3  :
          canvas.pen_color(255, 0, 0);
          break; //optional
          case 4  :
          canvas.pen_color(0,0,0);
        }
      //  // you can have any number of case statements.
      //   default : //Optional
      //   canvas.pen_color(0,0,0);
      // }

      // int shade = 255*ca->actual_grid[j][i]->spread_wind_dir;
      // if (ca->actual_grid[j][i]->spread_wind_dir == 42) {
      //   canvas.pen_color(255,255,255);
      // } else {
      //   canvas.pen_color(shade, 0, 0);
      // }
      // std::cout << i << std::endl;
      canvas.fill_rectangle(x_coord(image_width,(image_width/ca->cols)*j),y_coord(image_height,(image_height/ca->rows)*i),x_coord(image_width,(image_width/ca->cols)*j+image_width/ca->cols),y_coord(image_height,(image_height/ca->rows)*i+(image_height/ca->rows)));
      // canvas.fill_rectangle(x_coord(image_width,(image_width/ca->cols)*i),y_coord(image_height,(image_height/ca->rows)*j),x_coord(image_width,(image_width/ca->cols)*i+image_width/ca->cols),y_coord(image_height,(image_height/ca->rows)*j+(image_height/ca->rows)));
      canvas.pen_color(0,0,0);
      canvas.rectangle(x_coord(image_width,(image_width/ca->cols)*j),y_coord(image_height,(image_height/ca->rows)*i),x_coord(image_width,(image_width/ca->cols)*j+image_width/ca->cols),y_coord(image_height,(image_height/ca->rows)*i+(image_height/ca->rows)));
    }
  }
  std::string timeStr = std::to_string(time);
  canvas.image().save_image("moj_obrazok-"+timeStr+".bmp");

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

double randnum (double a, double b)
{
  static std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution (a,b);
  return distribution(generator);
}

/**
 * @brief      Constructor for copying cells
 *
 * @param      pointer  Pointer to the cell to be copied
 */
Cell* Cell::copy() {
  Cell* to_be_returned = new Cell();

  to_be_returned->x = this->x;
  to_be_returned->y = this->y;
  to_be_returned->spread_wind_dir = this->spread_wind_dir;
  to_be_returned->material_properties = this->material_properties;
  to_be_returned->ratio_of_area = this->ratio_of_area;
  to_be_returned->spread_probability = this->spread_probability;
  to_be_returned->state = this->state;
  to_be_returned->burning_time = this->burning_time;

  return to_be_returned;
}


Automata::Automata(automata_parameters params) {
  cell_size = params.cell_size;
  rows = params.rows;
  cols = params.cols;
  simulation_time = params.simulation_time;
  time_step = params.time_step;
  actual_time = 0;
  wind_velocity = params.wind_velocity;
  wind_direction = params.wind_direction;
  alpha = params.alpha;
  beta = params.beta;


  // ---------------ORIG algoritmus----------------------
  std::vector<Cell*> tmp_vec;
  // filling actual grid
  for (unsigned int x = 0; x < cols; ++x) {
    for (unsigned int y = 0; y < rows; ++y) {
      material tmpMat;
      tmpMat.structure_combustibility = 0;
      tmpMat.spread_ability_time = 0;
      tmpMat.extinction_time = 0;

      Cell* to_be_added = new Cell;
      to_be_added->material_properties = tmpMat;
      to_be_added->x = x;
      to_be_added->y = y;
      to_be_added->state = 1; // TODO - put it to define
      to_be_added->spread_wind_dir = 42; // TODO - delete this line
      to_be_added->ratio_of_area = 1;
      tmp_vec.push_back(to_be_added);
    }
    actual_grid.push_back(tmp_vec);
    tmp_vec.clear();
  }
  // ----------------------------------------------------

  bitmap_image image("input-1.bmp");

  if (!image) {
    printf("Error - Failed to open: input.bmp\n");
    exit(1);
  }

  const unsigned int height = image.height();
  const unsigned int width  = image.width();

  rgb_t colour;

  for (std::size_t y = 0; y < height; y+=3) {
    for (std::size_t x = 0; x < width; x+=3) {
      float woodRatio = 0.0;
      float fireproofWoodRatio = 0.0;
      bool fire = false;
      for (int yi = 0; yi < 3; ++yi) {
        for (int xi = 0; xi < 3; ++xi) {
          rgb_t colour;
          image.get_pixel(x+xi, y+yi, colour);
          // std::cout << "Citam z [" << x+xi << "," << y+yi << "], je tu farba: " << colour.red << "," << colour.green << "," << colour.blue << std::endl;
          printf("Citam z [%d,%d], je tu farba %d,%d,%d\n", x+xi, y+y, colour.red, colour.green, colour.blue);
          if (colour.red == 104 && colour.green == 55 && colour.blue == 29) {
            // wood
            woodRatio += 1.0/9;
            std::cout << woodRatio << std::endl;
          } else if (colour.red == 89 && colour.green == 89 && colour.blue == 0) {
            // fireproofWood
            fireproofWoodRatio += 1.0/9;
          } else if (colour.red == 255 && colour.green == 0 && colour.blue == 0) {
            fire = true;
          }
        }
      }


      if (fire) {
        actual_grid[x/3][y/3]->state = 3;
      }
      if (woodRatio != 0.0 || fireproofWoodRatio != 0.0) {
        actual_grid[x/3][y/3]->material_properties.structure_combustibility = woodRatio > fireproofWoodRatio ? 1.0 : 0.6;
        actual_grid[x/3][y/3]->material_properties.spread_ability_time = 2;
        actual_grid[x/3][y/3]->material_properties.extinction_time = 10;

        actual_grid[x/3][y/3]->state = 1; // TODO - put it to define
        actual_grid[x/3][y/3]->ratio_of_area = woodRatio > fireproofWoodRatio ? woodRatio : fireproofWoodRatio;
      }
    }
  }

  // filling next grid
  for (unsigned int x = 0; x < cols; ++x) {
    for (unsigned int y = 0; y < rows; ++y) {
      Cell* to_be_added_to_new = actual_grid[x][y]->copy();
      tmp_vec.push_back(to_be_added_to_new);
    }
    next_grid.push_back(tmp_vec);
    tmp_vec.clear();
  }

}


  // for (int i = 0; i < r; ++i) {
  //   for (int j = 0; j < c; ++j) {
  //     std::cout << "|" << actual_grid[i][j]->x << "," << actual_grid[i][j]->y << "|";
  //   }
  //   std::cout << std::endl;
  // }
// }

void Automata::actualize_grid() {
  for (unsigned int x = 0; x < cols; ++x) {
    for (unsigned int y = 0; y < rows; ++y) {
      delete actual_grid[x][y];
      actual_grid[x][y] = next_grid[x][y]->copy();
    }
  }
}

bool Automata::isInLattice(unsigned int x, unsigned int y) {
  // I'm not verifying, if y >= 0 and x >= 0, because unsigned int is always >= 0
  return (x < cols && y < rows);
}

float Automata::getSpreadingProbability(int x, int y) {
  if (!isInLattice(x,y)) {
    return 0;
  }
  Cell* cell = actual_grid[x][y];

  float ptckl = 1.0;
  float t1 = cell->material_properties.spread_ability_time;
  float t2 = cell->material_properties.extinction_time;
  float t = cell->burning_time;

  if (((t2-t1)/5)+t1 <= t && t <= t2) {
    ptckl = (5*(-t+t2))/(4*(t2-t1));
  } else if (t1 <= t && t <= ((t2-t1)/5)+t1) {
    ptckl = (4.0*t)/(t2-t1)+(0.2*t2-4.2*t1)/(t2-t1);
  }

  // material compustibility
  float sij = cell->material_properties.structure_combustibility;
  // ratio of area occupied by building
  float pij = cell->ratio_of_area;
  // relates to spreading speed and direction
  float wij = cell->spread_wind_dir;

  return alpha*sij*pij*pow(wij,beta)*ptckl;
  // return alpha*cell->material_properties.structure_combustibility*cell->ratio_of_area*cell->spread_wind_dir;
}

/**
 * @brief      Changes probabilities of fire spreading of each cell in neighbourhood of cell specified by x and y
 *
 * @param[in]  x     { parameter_description }
 * @param[in]  y     { parameter_description }
 */
void Automata::process_fire(int x, int y) {

  //JUST PLAYING
  std::vector<std::vector<float>> v;
  if (wind_velocity == 0) {
    // the same for all wind directions
    v = {{0.3,0.3,0.3},{0.3,0.3},{0.3,0.3,0.3}};
  } else if (wind_velocity <= 1) {
      // the same for all wind directions
    v = {{0.3,0.3,0.3},{0.3,0.5,0.5,0.5,0.3},{0.3,0.5,0.5,0.3},{0.3,0.5,0.5,0.5,0.3},{0.3,0.3,0.3}};
  } else if (wind_velocity < 5) {
    if (wind_direction == west || wind_direction == east) {
      // this is from west
      v = {{0.1,0.3,0.3},{0.1,0.2,0.5,0.8,0.4},{0.1,0.2,0.8,0.4},{0.1,0.2,0.5,0.8,0.4},{0.1,0.3,0.3}};
    } else if (wind_direction == north || wind_direction == south) {
      // this is from north
      v = {{0.1,0.1,0.1},{0.1,0.2,0.2,0.2,0.1},{0.3,0.5,0.5,0.3},{0.3,0.8,0.8,0.8,0.3},{0.4,0.4,0.4}};
    }
  } else if (wind_velocity < 8) {
    if (wind_direction == west || wind_direction == east) {
      // this is from west
      v = {{0.05,0.2,0.3},{0.1,0.1,0.3,0.6,0.5},{0.05,0.1,0.2,0.5,0.8,0.6,0.4},{0.05,0.1,0.2,0.8,0.6,0.4},{0.05,0.1,0.2,0.5,0.8,0.6,0.4},{0.1,0.1,0.3,0.6,0.5},{0.05,0.2,0.3}};
    } else if (wind_direction == north || wind_direction == south) {
      // this is from north
      v = {{0.05,0.05,0.05},{0.1,0.1,0.1,0.1,0.1},{0.05,0.1,0.2,0.2,0.2,0.1,0.05},{0.2,0.3,0.5,0.5,0.3,0.2},{0.3,0.6,0.8,0.8,0.8,0.6,0.3},{0.5,0.6,0.6,0.6,0.5},{0.4,0.4,0.4}};
    }
  } else {
    if (wind_direction == west || wind_direction == east) {
      //this is from west
      v = {{0.05,0.05,0.2},{0.05,0.1,0.3,0.4,0.5},{0.05,0.05,0.2,0.3,1,0.6,0.3},{0.05,0.05,0.2,1,0.6,0.3},{0.05,0.05,0.2,0.3,1,0.6,0.3},{0.05,0.1,0.3,0.4,0.5},{0.05,0.05,0.2}};
    } else if (wind_direction == north || wind_direction == south) {
      // this is from north
      v = {{0.05,0.05,0.05},{0.05,0.05,0.05,0.05,0.05},{0.05,0.1,0.2,0.2,0.2,0.1,0.05},{0.05,0.3,0.3,0.3,0.3,0.05},{0.2,0.4,1,1,1,0.4,0.2},{0.5,0.6,0.6,0.6,0.5},{0.3,0.3,0.3}};
    }
  }

  std::vector<std::vector<Cell*>> neighbourhood;

  neighbourhood = getNeighbourhood(x,y);

  // changing the Wij of each cell in neigbourhood
  // row and col has meaning like row and column of neighbourhood
  int row_index, col_index;
  for (unsigned int row = 0; row < neighbourhood.size(); ++row) {
    for (unsigned int col = 0; col < neighbourhood[row].size(); ++col) {
      // std::cout << "[" << neighbourhood[row][col]->x << "," << neighbourhood[row][col]->y << "]";
      if (wind_direction == west || wind_direction == north) {
        row_index = row; col_index = col;
      } else {
        row_index = neighbourhood.size()-row-1;
        col_index = neighbourhood[row].size()-col-1;
      }
      // assigning spread_wind_dir property
      neighbourhood[row][col]->spread_wind_dir  = v[row_index][col_index];
    }
      // std::cout << std::endl << "-----------------------" << std::endl;
  }

  //JUST PLAYING

  // for (int row = 0; row < neighbourhood.size(); ++row) {
  //   for (int col = 0; col < neighbourhood[row].size(); ++col) {
  //     Cell* processed_cell = neighbourhood[row][col];
  //     if (processed_cell->state == 1) {
  //       float ran = 0.5;
  //       if (getSpreadingProbability(processed_cell->x,processed_cell->y) > ran) {
  //         // bunka na suradniciach bunky neighbourhood[row][col] bude v dalsom kroku horiet
  //         // getCellFromNextGrid(processed_cell->x, processed_cell->y)->state = 2;
  //         // getCellFromNextGrid(processed_cell->x, processed_cell->y)->burning_time = 0; // 0 alebo time_step?
  //       }
  //     }
  //   }
  // }

}

std::vector<std::vector<Cell*>> Automata::getNeighbourhood(unsigned int x_given, unsigned int y_given) {
  std::vector<std::vector<Cell*>> to_be_returned;
  std::vector<Cell*> tmp_vec;
  std::vector<int> neighbourhood_map;
  
  if (wind_velocity == 0) {
    neighbourhood_map = {1,1,1};
  } else if (wind_velocity < 5) {
    neighbourhood_map = {1,2,2,2,1};
  } else {
    neighbourhood_map = {1,2,3,3,3,2,1};
  }

  int i = 0;
  for (unsigned int y = y_given - neighbourhood_map.size()/2; y <= y_given + neighbourhood_map.size()/2; ++y) {
    for (unsigned int x = x_given-neighbourhood_map[i]; x <= x_given+neighbourhood_map[i]; ++x) {
      if (!(x == x_given && y == y_given)) {
        if (isInLattice(x,y)) {
          tmp_vec.push_back(actual_grid[x][y]);
        }
      }

    }
    to_be_returned.push_back(tmp_vec);
    tmp_vec.clear();
    i++;
  }

  return to_be_returned;
}

/**
 * @brief      Gets all the cells, that have state = 2 or state = 3
 *
 * @return     The 1 dimensional vector of spreading cells.
 */
std::vector<Cell*> Automata::getBurningCells() {
    // this vector holds the spreading cells (spreading cells have state = 3)
  std::vector<Cell*> burning_cells;
    // filling burning_cells vector
  for (unsigned int i = 0; i < cols; ++i) {
    for (unsigned int j = 0; j < rows; ++j) {
      if (isInLattice(i,j)) {
                if (actual_grid[i][j]->state == 3 || actual_grid[i][j]->state == 2) { // TODO - put 3 to define
                  burning_cells.push_back(actual_grid[i][j]);
                }
              }
            }
          }
          return burning_cells;
        }

/**
 * @brief      Gets all the cells, that have state = 3
 *
 * @return     The 1 dimensional vector of spreading cells.
 */
        std::vector<Cell*> Automata::getSpreadingCells() {
  // this vector holds the spreading cells (spreading cells have state = 3)
          std::vector<Cell*> spreading_cells;
  // filling spreading_cells vector
          for (unsigned int i = 0; i < cols; ++i) {
            for (unsigned int j = 0; j < rows; ++j) {
              if (isInLattice(i,j)) {
        if (actual_grid[i][j]->state == 3) { // TODO - put 3 to define
          spreading_cells.push_back(actual_grid[i][j]);
        }
      }
    }
  }
  return spreading_cells;
}

Cell* getCell(unsigned int x, unsigned int y, std::vector<std::vector<Cell*>> grid) {
  if (x < grid.size() && y < grid[0].size()) {
    return grid[x][y];
  }
  return NULL;
}

void Automata::simulate() {
  for (actual_time = 0; actual_time < simulation_time; actual_time += time_step) {

    // Updating cell states if they are burning too long for their states
    std::vector<Cell*> burning_cells = getBurningCells();
    for (unsigned int i = 0; i < burning_cells.size(); ++i) {
      Cell* tmpCell = burning_cells[i];
      unsigned int x = tmpCell->x;
      unsigned int y = tmpCell->y;

      if (tmpCell->state == 2 || tmpCell->state == 3) {
        next_grid[x][y]->burning_time += time_step;
      }
      if ((tmpCell->state == 2 && tmpCell->burning_time >= tmpCell->material_properties.spread_ability_time) || (tmpCell->state == 3 && tmpCell->burning_time >= tmpCell->material_properties.extinction_time)) {
        next_grid[x][y]->state++;
      }
    }

    // std::vector<Cell*> spreading_cells = getSpreadingCells();
    // for each spreading cell
    for (unsigned int i = 0; i < burning_cells.size(); i++) {
      if (burning_cells[i]->state == 3) {
      // getting the neighbourhood of spreading cell
        Cell* spreading_cell = burning_cells[i];
        std::vector<std::vector<Cell*>> neighbourhood_of_spreading = getNeighbourhood(spreading_cell->x, spreading_cell->y);
      // set Wij parameters for equation of Fij
        process_fire(spreading_cell->x,spreading_cell->y);
      // for each cell from neighbourhood of spreading cell
        for (unsigned int j = 0; j < neighbourhood_of_spreading.size(); j++) {
          for (unsigned int k = 0; k < neighbourhood_of_spreading[j].size(); k++) {
            Cell* tmpCell = getCell(neighbourhood_of_spreading[j][k]->x,neighbourhood_of_spreading[j][k]->y,actual_grid);
            if (tmpCell->state == 1) {
              float spread_prob = getSpreadingProbability(tmpCell->x,tmpCell->y);
              float rnd = randnum(0.0,1.0);
              std::cout << "rnd... " << rnd << std::endl;
              if (spread_prob > rnd) {
                next_grid[tmpCell->x][tmpCell->y]->state = 2;
                next_grid[tmpCell->x][tmpCell->y]->burning_time = 0;
              }
            }
          }
        }
      }
    }
    burning_cells.clear();


    print_ca(this, actual_time);

    actualize_grid();

    /* code */
    /* code */
    /* code */
    // burning_time kazdej horiacej bunky += time_step
    // dealokuj a vymaz to, co je v actual_grid
    // 
    // actual_grid = copy(next_grid);
    // 
    // obsah next_grid bude ten isty, budem na nom stavat dalej
  }
}

int main () {
  automata_parameters params;
  params.cell_size = 3;
  params.rows = 100;
  params.cols = 100;
  params.time_step = 1;
  params.simulation_time = 30;
  params.wind_velocity = 15;
  params.wind_direction = south;
  params.alpha = 1;
  params.beta = 1;
  params.time_step = 1;
  params.simulation_time = 30;
  params.wind_velocity = 15;
  params.wind_direction = south;
  params.alpha = 1;
  params.beta = 1;

  Automata* ca = new Automata(params);
  ca->simulate();

  print_ca(ca, 999);

  // getCell(25,45,ca->actual_grid)->state = 3;
  // getCell(25,45,ca->next_grid)->state = 3;
  // getCell(15,15,ca->actual_grid)->state = 3;
  // getCell(15,15,ca->next_grid)->state = 3;

  // ca->simulate();

  // std::cout << ca->actual_grid[25][25]->state << std::endl;
  // std::cout << ca->next_grid[25][25]->state << std::endl;



  // //////just testing start///////////
  // material tmpMat;
  // tmpMat.structure_combustibility = 1;
  // tmpMat.spread_ability_time = 2;
  // tmpMat.extinction_time = 10;

  // Cell* tmpCell = new Cell;
  // tmpCell->material_properties = tmpMat;
  // tmpCell->x = 5;
  // tmpCell->y = 6;
  // tmpCell->state = 1; // TODO - put it to define
  // tmpCell->spread_wind_dir = 42; // TODO - delete this line
  // tmpCell->ratio_of_area = 1;

  // std::cout << tmpCell << std::endl;
  // Cell* anotherCell = tmpCell->copy();
  // anotherCell->x = 1598;
  // std::cout << anotherCell->x << std::endl;
  // std::cout << tmpCell->x << std::endl;
  //////just testing end/////////////

  // std::vector<std::vector<Cell*>> tmp = ca->getNeighbourhood(5,20);

  // std::cout << tmp.size() << "x" << tmp[0].size() << std::endl;
  // for (int i = 0; i < tmp.size(); ++i) {
    // for (int j = 0; j < tmp[i].size(); ++j) {
  //     std::cout << "|" << actual_grid[i][j]->x << "," << actual_grid[i][j]->y << "|";
      // std::cout << "|" << tmp[i][j]->x << "," << tmp[i][j]->y;
    // }
    // std::cout << std::endl;
  // }

  // std::cout << "================================" << std::endl;

  // for (int i = 0; i < ca->actual_grid.size(); ++i) {
  //   for (int j = 0; j < ca->actual_grid[i].size(); ++j) {
  //     std::cout << "[" << ca->actual_grid[i][j]->x << "," << ca->actual_grid[i][j]->y << "]" << std::endl;
  //   }
  //   std::cout << std::endl;
  // }


  // for (unsigned int i = 0; i < tmp.size(); ++i) {
  //   for (unsigned int j = 0; j < tmp[i].size(); ++j) {
  //     tmp[i][j]->state = 2;
  //   }
  // }

  // ca->process_fire();



  return 0;
}
