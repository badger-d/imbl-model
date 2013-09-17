#ifndef StackParameterisation_H
#define StackParameterisation_H 1

#include <vector>
#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4RotationMatrix.hh"

using namespace std;


class G4VPhysicalVolume;
class G4VTouchable;
class G4Box;
class G4Material;
class G4String;

// Dummy declarations to get rid of warnings ...
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;

class StackParameterisation : public G4VPVParameterisation
{
public:

  StackParameterisation(G4int num_layers,                         // Number of layers.
						std::vector<G4Material*> &mat_store,      // Vector of materials.
						std::vector<G4ThreeVector> &pos_store,    // Vector of positions.
						std::vector<G4ThreeVector> &shape_store); // Vector of shapes (x, y, z).

  virtual ~StackParameterisation();

  // Compute the location of the current layer.
  virtual void ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const;

  // Compute the shape of the current layer.
  virtual void ComputeDimensions (G4Box& box, const G4int copyNo, const G4VPhysicalVolume* physVol) const;

  // Compute the material of the current layer.
  virtual G4Material* ComputeMaterial(const G4int copyNo, G4VPhysicalVolume * physVol,  const G4VTouchable *parentTouch );

private:  // Dummy declarations to get rid of warnings ...

  // void ComputeDimensions (G4Box&, const G4int, const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Tubs&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Trd&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Trap&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Cons&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Sphere&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Orb&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Torus&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Para&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Hype&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Polycone&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Polyhedra&,const G4int,const G4VPhysicalVolume*) const {}

private:

  G4double _num_layers;
  std::vector<G4Material*> _mat_store;
  std::vector<G4ThreeVector> _pos_store;
  std::vector<G4ThreeVector> _shape_store;
};


#endif


