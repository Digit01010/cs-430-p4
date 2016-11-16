#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265

int line = 1;

// Holds an rgb triple of a pixel
typedef struct Pixel {
  unsigned char red, green, blue;
} Pixel;

// Holds information about the header of a ppm file
typedef struct Header {
   unsigned char magicNumber;
   unsigned int width, height, maxColor;
} Header;

// Plymorphism in C

typedef struct {
  int kind; // 0 = camera, 1 = sphere, 2 = plane, 3 = light
  double color[3]; // diffuse color
  double specular_color[3];
  double position[3];
  union {
    struct {
      double width;
      double height;
    } camera;
    struct {
      double radius;
    } sphere;
    struct {
      double normal[3];
    } plane;
    struct {
      double radial_a0;
      double radial_a1;
      double radial_a2;
      double theta;
      double angular_a0;
      double direction[3];
    } light;
  };
} Object;

double sphere_intersection(double*, double*, double*, double);
double plane_intersection(double*, double*, double*, double*);
void writeP3(Pixel *, Header, FILE *);
int next_c(FILE*);
void expect_c(FILE*, int);
void skip_ws(FILE*);
char* next_string(FILE*);
double next_number(FILE*);
double* next_vector(FILE*);
Object** read_scene(char*);

static inline double sqr(double v) {
  return v*v;
}

static inline double v3_len(double x, double y, double z) {
  return sqrt(sqr(x) + sqr(y) + sqr(z));
}

static inline double v3_dot(double* a, double* b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double clamp(double n, double min, double max) {
  if (n < min) n = min;
  if (n > max) n = max;
  return n;
}

static inline void normalize(double* v) {
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

int main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(stderr, "Error: Incorrect number of arguments.\n");
    printf("Usage: raycast width height input.json output.ppm\n");
    return(1);
  }
  
  // Get dimensions
  int N = atoi(argv[1]);
  int M = atoi(argv[2]);
  
  if (N <= 0 || M <= 0) {
    fprintf(stderr, "Error: Invalid dimensions.\n");
    exit(1);
  }
  
  // Read the scene file

  
  Object** cameras = malloc(sizeof(Object*)*129);
  Object** objects = malloc(sizeof(Object*)*129);
  Object** lights = malloc(sizeof(Object*)*129);
  
  Object** json_objects = read_scene(argv[3]);
  int camcnt = 0;
  int objcnt = 0;
  int lgtcnt = 0;
  for (int i = 0; json_objects[i] != NULL; i++) {
    if (json_objects[i]->kind == 0) {
      cameras[camcnt] = json_objects[i];
      camcnt++;
    } else if (json_objects[i]->kind == 1 || json_objects[i]->kind == 2) {
      objects[objcnt] = json_objects[i];
      objcnt++;
    } else if (json_objects[i]->kind == 3) {
      lights[lgtcnt] = json_objects[i];
      lgtcnt++;
    }
  }
  cameras[camcnt] = NULL;
  objects[objcnt] = NULL;
  lights[lgtcnt] = NULL;

  // Find the camera and get the height and width  
  if (camcnt == 0) {
    fprintf(stderr, "Error: No camera found.\n");
    exit(1);
  }
  
  if (camcnt > 1) {
    fprintf(stderr, "Error: Multiple cameras not supported.\n");
    exit(1);
  }
  
  double cx = 0;
  double cy = 0;
  
  double h = cameras[0]->camera.height;
  double w = cameras[0]->camera.width;
   
  // Initialize pixel buffer
  Pixel *buffer = malloc(sizeof(Pixel) * N * M);
  
  double pixheight = h / M;
  double pixwidth = w / N;
  for (int y = M; y > 0; y -= 1) { // Going through y from greatest to least
                                   // for convenient pixel output
    for (int x = 0; x < N; x += 1) {
      // Raycast
      
      double Ro[3] = {0, 0, 0};
      // Rd = normalize(P - Ro)
      double Rd[3] = {
        cx - (w/2) + pixwidth * (x + 0.5),
        cy - (h/2) + pixheight * (y + 0.5),
        1
      };
      normalize(Rd);

      double best_t = INFINITY;
      int best_i = 0;
      // Check intersections
      for (int i=0; objects[i] != NULL; i += 1) {
        double t = 0;
        // Call correct intersection function
        switch(objects[i]->kind) {
        case 0:
          t = -1;
          break;
        case 1:
          t = sphere_intersection(Ro, Rd,
                                    objects[i]->position,
                                    objects[i]->sphere.radius);
          break;
        case 2:
          t = plane_intersection(Ro, Rd,
                                    objects[i]->position,
                                    objects[i]->plane.normal);
          break;
        default:
          fprintf(stderr, "Error: Programmer forgot to implement an intersection.");
          exit(1);
        }
        if (t > 0 && t < best_t) {
          best_t = t;
          best_i = i;
        }
      }
      
      double* color = malloc(sizeof(double)*3);
      color[0] = 0; // ambient_color[0];
      color[1] = 0; // ambient_color[1];
      color[2] = 0; // ambient_color[2];
      
      if (best_t > 0 && best_t != INFINITY) {
        for (int j=0; lights[j] != NULL; j++) {
        // Shadow test
          double Ron[3] = {0, 0, 0};
          double Rdn[3] = {0, 0, 0};
          Ron[0] = best_t * Rd[0] + Ro[0];
          Ron[1] = best_t * Rd[1] + Ro[1];
          Ron[2] = best_t * Rd[2] + Ro[2];
          Rdn[0] = lights[j]->position[0] - Ron[0];
          Rdn[1] = lights[j]->position[1] - Ron[1];
          Rdn[2] = lights[j]->position[2] - Ron[2];
          double lgtdist = v3_len(Rdn[0], Rdn[1], Rdn[2]);
          normalize(Rdn);
          int closest_i = -1;
          double closest_t = INFINITY;
          //closest_shadow_object = ...;
          for (int i=0; objects[i] != NULL; i += 1) {
            double t = 0;
            if (i == best_i) continue; // Skip own object
            // Call correct intersection function
            switch(objects[i]->kind) {
            case 0:
              t = -1;
              break;
            case 1:
              t = sphere_intersection(Ron, Rdn,
                                        objects[i]->position,
                                        objects[i]->sphere.radius);
              break;
            case 2:
              t = plane_intersection(Ron, Rdn,
                                        objects[i]->position,
                                        objects[i]->plane.normal);
              break;
            default:
              fprintf(stderr, "Error: Programmer forgot to implement an intersection in light test.");
              exit(1);
            }
            if (t > 0 && t < lgtdist) {
              closest_t = t;
              closest_i = i;
            }
          }
          
          if (closest_i == -1) {
           // N, L, R, V
           
            double N[3];
            switch (objects[best_i]->kind) {
              case 1: // sphere
                N[0] = Ron[0] - objects[best_i]->position[0];
                N[1] = Ron[1] - objects[best_i]->position[1];
                N[2] = Ron[2] - objects[best_i]->position[2];
                break;
              case 2: // plane
                N[0] = objects[best_i]->plane.normal[0];
                N[1] = objects[best_i]->plane.normal[1];
                N[2] = objects[best_i]->plane.normal[2];
                break;
              default:
                fprintf(stderr, "Error: Programmer forgot to implement object normal %d.\n", line);
                exit(1);
                break;
            }
            
            
            normalize(N);
            double* L = Rdn; // light_position - Ron;
            normalize(L);
                     
            double fang;
            double theta = lights[j]->light.theta;
            if (theta == 0) { // Not a spotlight
              fang = 1;
            } else {
              double nRdn[3];
              nRdn[0] = -Rdn[0];
              nRdn[1] = -Rdn[1];
              nRdn[2] = -Rdn[2];
              normalize(nRdn);
              normalize(lights[j]->light.direction);
              double cos_a = v3_dot(lights[j]->light.direction, nRdn);
              if (cos_a < cos(theta*PI/180.0)) {
                fang = 0;
              }
              else {
                fang = pow(cos_a, lights[j]->light.angular_a0);
              }
            }
            
            
            double a_0 = lights[j]->light.radial_a0;
            double a_1 = lights[j]->light.radial_a0;
            double a_2 = lights[j]->light.radial_a0;
            double frad = 1/(a_0 + a_1 * lgtdist + a_2 * sqr(lgtdist));
            
            double NdotL = v3_dot(N, L);
            
            double diffuse = 0;
            if (NdotL > 0) {
              diffuse = NdotL;
            }
            
            double V[3];
            V[0] = -Rd[0];
            V[1] = -Rd[1];
            V[2] = -Rd[2];
            normalize(V);
            
            double R[3]; // Reflection of L
            R[0] = L[0] + 2*(NdotL*N[0] - L[0]);
            R[1] = L[1] + 2*(NdotL*N[1] - L[1]);
            R[2] = L[2] + 2*(NdotL*N[2] - L[2]);
            normalize(R);
            
            double specular = 0;
            double VdotR = v3_dot(V, R);
            if (VdotR > 0) {
              specular = pow(VdotR, 20);
            }
            
            double* dc = objects[best_i]->color;
            double* sc = objects[best_i]->specular_color;
            double* lc = lights[j]->color;
            
            color[0] += frad * fang * lc[0]*(diffuse*dc[0] + specular*sc[0]);
            color[1] += frad * fang * lc[1]*(diffuse*dc[1] + specular*sc[1]);
            color[2] += frad * fang * lc[2]*(diffuse*dc[2] + specular*sc[2]);
          }
        }
        
      }
      // Note: Going through y in reverse, so adjust index accordingly
      int p = (M - y)*N + x; // Index of buffer
      buffer[p].red = (int) (255.0 *  clamp(color[0], 0.0, 1.0));
      buffer[p].green = (int) (255.0 * clamp(color[1], 0.0, 1.0));
      buffer[p].blue = (int) (255.0 * clamp(color[2], 0.0, 1.0));
      
    }
  }
  
  FILE* output = fopen(argv[4], "w");
  
  if (output == NULL) {
    fprintf(stderr, "Error: Could not write to file \"%s\"\n", argv[4]);
    exit(1);
  }
  
  Header outHeader;
  outHeader.magicNumber = 3;
  outHeader.maxColor = 255;
  outHeader.width = N;
  outHeader.height = M;
  
  writeP3(buffer, outHeader, output);
  
  return 0;
}

double sphere_intersection(double* Ro, double* Rd,
                             double* C, double r) {
  double a = (sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]));
  double b = 2 * (Rd[0] * (Ro[0] - C[0]) + Rd[1] * (Ro[1] - C[1]) + Rd[2] * (Ro[2] - C[2]));  
  double c = sqr(Ro[0] - C[0]) + sqr(Ro[1] - C[1]) + sqr(Ro[2] - C[2]) - sqr(r);

  double det = sqr(b) - 4 * a * c;
  if (det < 0) return -1;

  det = sqrt(det);
  
  double t0 = (-b - det) / (2*a);
  if (t0 > 0) return t0;

  double t1 = (-b + det) / (2*a);
  if (t1 > 0) return t1;

  return -1;
}

double plane_intersection(double* Ro, double* Rd,
                             double* P, double* n) {
  double t = -(n[0]*(Ro[0]-P[0]) + n[1]*(Ro[1]-P[1]) + n[2]*(Ro[2]-P[2])) / (n[0]*Rd[0] + n[1]*Rd[1] + n[2]*Rd[2]);
  return t;
}


// Writes P3 data
void writeP3(Pixel *buffer, Header h, FILE *fh) {
  // Write the header
  fprintf(fh, "P%d\n%d %d\n%d\n", h.magicNumber, h.width, h.height, h.maxColor);
  // Write the ascii data
  for (int i = 0; i < h.width * h.height; i++) {
     fprintf(fh, "%d\n%d\n%d\n", buffer[i].red, buffer[i].green, buffer[i].blue);
  }
}

// next_c() wraps the getc() function and provides error checking and line
// number maintenance
int next_c(FILE* json) {
  int c = fgetc(json);
#ifdef DEBUG
  printf("next_c: '%c'\n", c);
#endif
  if (c == '\n') {
    line += 1;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
    exit(1);
  }
  return c;
}


// expect_c() checks that the next character is d.  If it is not it emits
// an error.
void expect_c(FILE* json, int d) {
  int c = next_c(json);
  if (c == d) return;
  fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
  exit(1);    
}


// skip_ws() skips white space in the file.
void skip_ws(FILE* json) {
  int c = next_c(json);
  while (isspace(c)) {
    c = next_c(json);
  }
  ungetc(c, json);
}


// next_string() gets the next string from the file handle and emits an error
// if a string can not be obtained.
char* next_string(FILE* json) {
  char buffer[129];
  int c = next_c(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line %d.\n", line);
    exit(1);
  }  
  c = next_c(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
      exit(1);      
    }
    if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
      exit(1);      
    }
    if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
      exit(1);
    }
    buffer[i] = c;
    i += 1;
    c = next_c(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}

double next_number(FILE* json) {
  double value;
  int cnt = fscanf(json, "%lf", &value);
  if (cnt != 1) {
    fprintf(stderr, "Error: Could not read number on line %d.\n", line);
    exit(1);
  }
  return value;
}

double* next_vector(FILE* json) {
  double* v = malloc(3*sizeof(double));
  expect_c(json, '[');
  skip_ws(json);
  v[0] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[1] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[2] = next_number(json);
  skip_ws(json);
  expect_c(json, ']');
  return v;
}


Object** read_scene(char* filename) {
  int c;
  
  Object** objects;
  objects = malloc(sizeof(Object*)*129);
  
  FILE* json = fopen(filename, "r");

  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", filename);
    exit(1);
  }
  
  skip_ws(json);
  
  // Find the beginning of the list
  expect_c(json, '[');

  skip_ws(json);

  // Find the objects
  int objcnt = 0;
  while (1) {
    c = fgetc(json);
    if (c == ']') {
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      exit(1);
    }
    if (c == '{') {
      skip_ws(json);
    
      // Parse the object
      char* key = next_string(json);
      if (strcmp(key, "type") != 0) {
        fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
        exit(1);
      }

      skip_ws(json);

      expect_c(json, ':');

      skip_ws(json);

      char* value = next_string(json);
      if (strcmp(value, "camera") == 0) {     
        objects[objcnt] = malloc(sizeof(Object));
        objects[objcnt]->kind = 0;
      } else if (strcmp(value, "sphere") == 0) {
        objects[objcnt] = malloc(sizeof(Object));
        objects[objcnt]->kind = 1;
      } else if (strcmp(value, "plane") == 0) {
        objects[objcnt] = malloc(sizeof(Object));
        objects[objcnt]->kind = 2;
      } else if (strcmp(value, "light") == 0) {
        objects[objcnt] = malloc(sizeof(Object));
        objects[objcnt]->kind = 3;
        objects[objcnt]->light.theta = 0; // Must be zero if a point light
      } else {
        fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
        exit(1);
      }
      
      // Defaults
      objects[objcnt]->color[0] = 0;
      objects[objcnt]->color[1] = 0;
      objects[objcnt]->color[2] = 0;

      objects[objcnt]->specular_color[0] = 1;
      objects[objcnt]->specular_color[1] = 1;
      objects[objcnt]->specular_color[2] = 1;

      skip_ws(json);

      int valcnt = 0;
      while (1) {
        // , }
        c = next_c(json);
        if (c == '}') {
          // stop parsing this object
          switch (objects[objcnt]->kind) {
          case 0:
            if (valcnt != 2 ) {
              fprintf(stderr, "Error: Bad value count.");
              exit(1);
            }
            break;
          case 1:
            if (valcnt < 3 || valcnt > 4) {
              fprintf(stderr, "Error: Bad value count.");
              exit(1);
            }
            break;
          case 2:
            if (valcnt < 3 || valcnt > 4) {
              fprintf(stderr, "Error: Bad value count.");
              exit(1);
            }
            break;
          case 3:
            if (valcnt < 6 || valcnt > 9) {
              fprintf(stderr, "Error: Bad value count.");
              exit(1);
            }
            break;
          default:
            fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
            exit(1);
            break;
          }
          
          
          objcnt += 1;
          if (objcnt > 128) {
            fprintf(stderr, "Error: 128 object count exceeded.");
            exit(1);
          }
          break;
        } else if (c == ',') {
          // read another field
          skip_ws(json);
          char* key = next_string(json);
          skip_ws(json);
          expect_c(json, ':');
          skip_ws(json);
          
          
          
          if (strcmp(key, "width") == 0) {
            double value = next_number(json);
            switch (objects[objcnt]->kind) {
            case 0:
              objects[objcnt]->camera.width = value;
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          } else if (strcmp(key, "height") == 0) {
            double value = next_number(json);
            switch (objects[objcnt]->kind) {
            case 0:
              objects[objcnt]->camera.height = value;
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          } else if (strcmp(key, "radius") == 0) {
            double value = next_number(json);
            switch (objects[objcnt]->kind) {
            case 1:
              objects[objcnt]->sphere.radius = value;
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          } else if (strcmp(key, "diffuse_color") == 0 || strcmp(key, "color") == 0) {
            double* value = next_vector(json);
            switch (objects[objcnt]->kind) {
            case 0:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            default:
              objects[objcnt]->color[0] = value[0];
              objects[objcnt]->color[1] = value[1];
              objects[objcnt]->color[2] = value[2];
              break;
            }
          } else if (strcmp(key, "specular_color") == 0) {
            double* value = next_vector(json);
            switch (objects[objcnt]->kind) {
            case 1:
              objects[objcnt]->specular_color[0] = value[0];
              objects[objcnt]->specular_color[1] = value[1];
              objects[objcnt]->specular_color[2] = value[2];
              break;           
            case 2:            
              objects[objcnt]->specular_color[0] = value[0];
              objects[objcnt]->specular_color[1] = value[1];
              objects[objcnt]->specular_color[2] = value[2];
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          } else if (strcmp(key, "position") == 0){
            double* value = next_vector(json);
            switch (objects[objcnt]->kind) {
            case 1:
              objects[objcnt]->position[0] = value[0];
              objects[objcnt]->position[1] = value[1];
              objects[objcnt]->position[2] = value[2];
              break;
            case 2:
              objects[objcnt]->position[0] = value[0];
              objects[objcnt]->position[1] = value[1];
              objects[objcnt]->position[2] = value[2];
              break;
            case 3:
              objects[objcnt]->position[0] = value[0];
              objects[objcnt]->position[1] = value[1];
              objects[objcnt]->position[2] = value[2];
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          } else if (strcmp(key, "normal") == 0) {
            double* value = next_vector(json);
            switch (objects[objcnt]->kind) {
            case 2:
              objects[objcnt]->plane.normal[0] = value[0];
              objects[objcnt]->plane.normal[1] = value[1];
              objects[objcnt]->plane.normal[2] = value[2];
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          } 
          else if (strcmp(key, "radial-a0") == 0) {
            double value = next_number(json);
            switch (objects[objcnt]->kind) {
            case 3:
              objects[objcnt]->light.radial_a0 = value;
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          } else if (strcmp(key, "radial-a1") == 0) {
            double value = next_number(json);
            switch (objects[objcnt]->kind) {
            case 3:
              objects[objcnt]->light.radial_a1 = value;
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          }
          else if (strcmp(key, "radial-a2") == 0) {
            double value = next_number(json);
            switch (objects[objcnt]->kind) {
            case 3:
              objects[objcnt]->light.radial_a2 = value;
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          }
          else if (strcmp(key, "theta") == 0) {
            double value = next_number(json);
            switch (objects[objcnt]->kind) {
            case 3:
              objects[objcnt]->light.theta = value;
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          }
          else if (strcmp(key, "angular-a0") == 0) {
            double value = next_number(json);
            switch (objects[objcnt]->kind) {
            case 3:
              objects[objcnt]->light.radial_a0 = value;
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          }
          else if (strcmp(key, "direction") == 0) {
            double* value = next_vector(json);
            switch (objects[objcnt]->kind) {
            case 3:
              objects[objcnt]->light.direction[0] = value[0];
              objects[objcnt]->light.direction[1] = value[1];
              objects[objcnt]->light.direction[2] = value[2];
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          }
          else {
            fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n",
                    key, line);
            //char* value = next_string(json);
          }
          valcnt += 1;
          skip_ws(json);
        } else {
          fprintf(stderr, "Error: Unexpected value on line %d\n", line);
          exit(1);
        }
      }
      skip_ws(json);
      c = next_c(json);
      if (c == ',') {
        // noop
        skip_ws(json);
      } else if (c == ']') {
        objects[objcnt] = NULL;
        fclose(json);
        return objects;
      } else {
        fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
        exit(1);
      }
    }

  }

  
  return NULL;
}

