#include "StackParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Material.hh"

StackParameterisation::StackParameterisation(G4int num_layers,                        // Number of layers.
											 std::vector<G4Material*> &mat_store,     // Vector of materials.
											 std::vector<G4ThreeVector> &pos_store,   // Vector of positions.
											 std::vector<G4ThreeVector> &shape_store) // Vector of shapes (x, y, z).
{
	_num_layers = num_layers;
	_mat_store = mat_store;
	_pos_store = pos_store;
	_shape_store = shape_store;
}

StackParameterisation::~StackParameterisation()
{}

void StackParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
	G4ThreeVector pos = _pos_store[copyNo];
	physVol->SetTranslation(pos);
    physVol->SetRotation(0);

}

void StackParameterisation::ComputeDimensions(G4Box& box, const G4int copyNo, const G4VPhysicalVolume*) const
{
	// Set the shape of the current layer.
	box.SetXHalfLength(_shape_store[copyNo].getX() * 0.5);
	box.SetYHalfLength(_shape_store[copyNo].getY() * 0.5);
	box.SetZHalfLength(_shape_store[copyNo].getZ() * 0.5);

}

G4Material* StackParameterisation::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume*, const G4VTouchable* parentTouch)
{
  // Protection for initialization and vis at idle state
  if(parentTouch==0) return _mat_store[copyNo];

  // Otherwise, return the material for the current layer.
  return _mat_store[copyNo];
}
