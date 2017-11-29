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

#define IMAGE_WIDTH   (700)
#define IMAGE_HEIGHT  (700)

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
  Automata (unsigned int, unsigned int);
  void process_fire(int, int);
  float getSpreadingProbability(int, int);
  std::vector<std::vector<Cell*>> getNeighbourhood(int, int);
  std::vector<Cell*> getSpreadingCells();
  bool isInLattice(int, int);
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
double x_coord(double x) {
 return x-(IMAGE_WIDTH/2);
}

double y_coord(double y) {

 return (IMAGE_HEIGHT/2)-y;
}

// TODO - do nazvu suboru pridavat cas simulacie
void print_ca(Automata* ca, int time) {
  cartesian_canvas canvas(700,700);

  canvas.image().clear(255);

  canvas.pen_width(1);

  ::srand(0xA5A5A5A5);
  float cell_width = IMAGE_WIDTH/ca->cols;
  std::cout << "cols: " << ca->cols << std::endl;
  std::cout << "rows: " << ca->rows << std::endl;
  for (unsigned int j = 0; j < ca->cols; j++) {
    for (unsigned int i = 0; i < ca->rows; i++) {
      int state = ca->actual_grid[j][i]->state;
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
      // int shade = 255*ca->actual_grid[j][i]->spread_wind_dir;
      // if (ca->actual_grid[j][i]->spread_wind_dir == 42) {
      //   canvas.pen_color(255,255,255);
      // } else {
      //   canvas.pen_color(shade, 0, 0);
      // }
      // std::cout << i << std::endl;
      canvas.fill_rectangle(x_coord((IMAGE_WIDTH/ca->cols)*j),y_coord((IMAGE_HEIGHT/ca->rows)*i),x_coord((IMAGE_WIDTH/ca->cols)*j+IMAGE_WIDTH/ca->cols),y_coord((IMAGE_HEIGHT/ca->rows)*i+(IMAGE_HEIGHT/ca->rows)));
      // canvas.fill_rectangle(x_coord((IMAGE_WIDTH/ca->cols)*i),y_coord((IMAGE_HEIGHT/ca->rows)*j),x_coord((IMAGE_WIDTH/ca->cols)*i+IMAGE_WIDTH/ca->cols),y_coord((IMAGE_HEIGHT/ca->rows)*j+(IMAGE_HEIGHT/ca->rows)));
      canvas.pen_color(0,0,0);
      canvas.rectangle(x_coord((IMAGE_WIDTH/ca->cols)*j),y_coord((IMAGE_HEIGHT/ca->rows)*i),x_coord((IMAGE_WIDTH/ca->cols)*j+IMAGE_WIDTH/ca->cols),y_coord((IMAGE_HEIGHT/ca->rows)*i+(IMAGE_HEIGHT/ca->rows)));
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


Automata::Automata(unsigned int c, unsigned int r) {
  rows = r;
  cols = c;
  actual_time = 0;
  time_step = 3;
  std::vector<Cell*> tmp_vec;

  // filling actual grid
  for (unsigned int x = 0; x < c; ++x) {
    for (unsigned int y = 0; y < r; ++y) {
      material tmpMat;
      tmpMat.structure_combustibility = 1;
      tmpMat.spread_ability_time = 2;
      tmpMat.extinction_time = 10;

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

  // filling next grid
  for (unsigned int x = 0; x < c; ++x) {
    for (unsigned int y = 0; y < r; ++y) {
      Cell* to_be_added_to_new = actual_grid[x][y]->copy();
      tmp_vec.push_back(to_be_added_to_new);
    }
    next_grid.push_back(tmp_vec);
    tmp_vec.clear();
  }

  // for (int i = 0; i < r; ++i) {
  //   for (int j = 0; j < c; ++j) {
  //     std::cout << "|" << actual_grid[i][j]->x << "," << actual_grid[i][j]->y << "|";
  //   }
  //   std::cout << std::endl;
  // }
}

void Automata::actualize_grid() {
  for (unsigned int x = 0; x < cols; ++x) {
    for (unsigned int y = 0; y < rows; ++y) {
      delete actual_grid[x][y];
      actual_grid[x][y] = next_grid[x][y]->copy();
    }
  }
}

bool Automata::isInLattice(int x, int y) {
  return (x >= 0 && x < cols &&  y >= 0 && y < rows);
}

float Automata::getSpreadingProbability(int x, int y) {
  if (!isInLattice(x,y)) {
    return 0;
  }
  Cell* cell = actual_grid[x][y];
  // std::cout << "========================" << std::endl;
  // std::cout << "SDSD:" << cell->spread_wind_dir << std::endl;
  // std::cout << "========================" << std::endl;
  // TODO equation (2) - add Wij, beta and p(tckl)
  return alpha*cell->material_properties.structure_combustibility*cell->ratio_of_area*cell->spread_wind_dir;
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
  for (int row = 0; row < neighbourhood.size(); ++row) {
    for (int col = 0; col < neighbourhood[row].size(); ++col) {
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

std::vector<std::vector<Cell*>> Automata::getNeighbourhood(int x_given, int y_given) {
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
  for (int y = y_given - neighbourhood_map.size()/2; y <= y_given + neighbourhood_map.size()/2; ++y) {
    for (int x = x_given-neighbourhood_map[i]; x <= x_given+neighbourhood_map[i]; ++x) {
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

Cell* getCell(int x, int y, std::vector<std::vector<Cell*>> grid) {
  if (x >= 0 && x < grid.size() && y >= 0 && y < grid[0].size()) {
    return grid[x][y];
  }
  return NULL;
}

void Automata::simulate() {
  for (actual_time = 0; actual_time < simulation_time; actual_time += time_step) {

    // Updating cell states if they are burning too long for their states
    for (unsigned int x = 0; x < cols; ++x) {
      for (unsigned int y = 0; y < rows; ++y) {
        Cell* tmpCell = actual_grid[x][y];
        if (tmpCell->state == 2 || tmpCell->state == 3) {
          next_grid[x][y]->burning_time += time_step;
        }
        if ((tmpCell->state == 2 && tmpCell->burning_time >= tmpCell->material_properties.spread_ability_time) || (tmpCell->state == 3 && tmpCell->burning_time >= tmpCell->material_properties.extinction_time)) {
          next_grid[x][y]->state++;
        }
      }
    }

    std::vector<Cell*> spreading_cells = getSpreadingCells();
    // for each spreading cell
    for (int i = 0; i < spreading_cells.size(); i++) {
      // getting the neighbourhood of spreading cell
      Cell* spreading_cell = spreading_cells[i];
      std::vector<std::vector<Cell*>> neighbourhood_of_spreading = getNeighbourhood(spreading_cell->x, spreading_cell->y);
      // set Wij parameters for equation of Fij
      process_fire(spreading_cell->x,spreading_cell->y);
      // for each cell from neighbourhood of spreading cell
      for (int j = 0; j < neighbourhood_of_spreading.size(); j++) {
        for (int k = 0; k < neighbourhood_of_spreading[j].size(); k++) {
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
    spreading_cells.clear();


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
  Automata* ca = new Automata(50,50);

  ca->time_step = 1;
  ca->simulation_time = 30;
  ca->wind_velocity = 0;
  // wind direction - TODO
  ca->wind_direction = west;
  ca->alpha = 1;
  ca->beta = 1;
  
  getCell(25,45,ca->actual_grid)->state = 3;
  getCell(25,45,ca->next_grid)->state = 3;
  getCell(15,15,ca->actual_grid)->state = 3;
  getCell(15,15,ca->next_grid)->state = 3;

  ca->simulate();

  std::cout << ca->actual_grid[25][25]->state << std::endl;
  std::cout << ca->next_grid[25][25]->state << std::endl;



  //////just testing start///////////
  material tmpMat;
  tmpMat.structure_combustibility = 1;
  tmpMat.spread_ability_time = 2;
  tmpMat.extinction_time = 10;

  Cell* tmpCell = new Cell;
  tmpCell->material_properties = tmpMat;
  tmpCell->x = 5;
  tmpCell->y = 6;
  tmpCell->state = 1; // TODO - put it to define
  tmpCell->spread_wind_dir = 42; // TODO - delete this line
  tmpCell->ratio_of_area = 1;

  std::cout << tmpCell << std::endl;
  Cell* anotherCell = tmpCell->copy();
  anotherCell->x = 1598;
  std::cout << anotherCell->x << std::endl;
  std::cout << tmpCell->x << std::endl;
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
