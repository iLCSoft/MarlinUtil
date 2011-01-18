#ifndef __ELEMENTS_H
#define __ELEMENTS_H

typedef struct {
  int         a_number;  // atomic number of element
  const char     *name;
  double        z_by_a;  // Z/A for each element
  double     exit_enrg;  // [MeV]
  double       density;  // [g/cm^3]
} BooElement;

extern BooElement BooElements[];

#endif /* __ELEMENTS_H */
