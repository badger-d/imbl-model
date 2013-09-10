//filename : VisManager
//----------------------------------------------------------------------------

#ifndef VisManager_h
#define VisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

//----------------------------------------------------------------------------
class VisManager: public G4VisManager {

public:

  VisManager ();

private:

  void RegisterGraphicsSystems ();

};

//----------------------------------------------------------------------------
#endif

#endif
