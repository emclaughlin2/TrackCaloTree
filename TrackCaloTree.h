// Track Calorimeter Calibration TreeMaker

#ifndef __TRACKCALOTREE_H__
#define __TRACKCALOTREE_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include <vector>
#include <map>
#include "TTree.h"
#include "TFile.h"
#include "TH2.h"

#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/PHCompositeNode.h>

#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerContainer.h>

#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

class PHCompositeNode;
class RawTowerContainer;
class PHG4TruthInfoContainer;
class RawTowerGeomContainer;
class SvtxTrackMap;
class SvtxTrackEval;
class SvtxEvalStack;
class PHHepMCGenEventMap;
class PHHepMCGenEvent;

class TrackCaloTree: public SubsysReco {
 public:
  
  TrackCaloTree(const std::string &name = "TrackCaloTree.root");
  
  virtual ~TrackCaloTree() {}
  
  int Init(PHCompositeNode *);
  
  int reset_tree_vars();

  int GetNodes(PHCompositeNode *);
  
  int process_event(PHCompositeNode *);
  
  int End(PHCompositeNode *);
  
 private:
  std::string _foutname; // outfile name
  std::string _hcal_raw_name;
  std::string _hcal_sim_name;
  std::string _hcal_calib_name;
  std::string _hcalIN_raw_name;
  std::string _hcalIN_sim_name;
  std::string _hcalIN_calib_name;
  std::string _EMcal_raw_name;
  std::string _EMcal_sim_name;
  std::string _EMcal_calib_name;
  std::string _track_name;

  int _startingSector, _numSectors; 
  RawTowerContainer* _towersSimOH;
  RawTowerContainer* _towersRawOH;
  RawTowerContainer* _towersCalibOH;
  RawTowerGeomContainer* _geomOH;
  RawTowerGeomContainer* _geomIH;
  RawTowerGeomContainer* _geomEM;

  RawTowerContainer* _towersRawIH;
  RawTowerContainer* _towersRawEM;
  RawTowerContainer* _towersSimIH;
  RawTowerContainer* _towersSimEM;
  RawTowerContainer* _towersCalibIH;
  RawTowerContainer* _towersCalibEM;
  SvtxEvalStack* m_svtxEvalStack;
  SvtxTrackMap *trackmap;
  PHG4TruthInfoContainer *truthinfo;
  PHHepMCGenEventMap *hepmceventmap;

  TTree* _tree;
  TFile* _file;
  
  int _ievent, _b_event;

  static const bool _debug = false;
  static const bool include_tracks = false;
  static const bool track_debug = false && include_tracks;
  static const int nTowers = 10000;
  static const int nTracks = 10000;
  static const int nTruth = 15000;
  static const int nHEP = 10000;
  
  // outer Hcal
  
  int _b_tower_sim_n;
  float _b_tower_sim_E[nTowers];
  float _b_tower_sim_eta[nTowers];
  float _b_tower_sim_phi[nTowers];
  int   _b_tower_sim_ieta[nTowers];
  int   _b_tower_sim_iphi[nTowers];
  
  int _b_tower_raw_n;
  float _b_tower_raw_E[nTowers];
  float _b_tower_raw_eta[nTowers];
  float _b_tower_raw_phi[nTowers];
  
  int _b_tower_calib_n;
  float _b_tower_calib_E[nTowers];
  float _b_tower_calib_eta[nTowers];
  float _b_tower_calib_phi[nTowers];
  int _b_tower_calib_ieta  [nTowers];
  int _b_tower_calib_iphi  [nTowers];
  
  // inner HCal
  
  int   _b_hcalIN_sim_n = 0;
  float _b_hcalIN_sim_E[nTowers];
  float _b_hcalIN_sim_eta[nTowers];
  float _b_hcalIN_sim_phi[nTowers];
  int   _b_hcalIN_sim_ieta[nTowers];
  int   _b_hcalIN_sim_iphi[nTowers];

  int   _b_hcalIN_raw_n = 0;
  float _b_hcalIN_raw_E[nTowers];
  float _b_hcalIN_raw_eta[nTowers];
  float _b_hcalIN_raw_phi[nTowers];
  
  int   _b_hcalIN_calib_n = 0;
  float _b_hcalIN_calib_E[nTowers];
  float _b_hcalIN_calib_eta[nTowers];
  float _b_hcalIN_calib_phi[nTowers];
  int   _b_hcalIN_calib_ieta[nTowers];
  int   _b_hcalIN_calib_iphi[nTowers];

  // EMCal
  
  int   _b_EMcal_sim_n = 0;
  float _b_EMcal_sim_E[nTowers];
  float _b_EMcal_sim_eta[nTowers];
  float _b_EMcal_sim_phi[nTowers];
  int   _b_EMcal_sim_ieta[nTowers];
  int   _b_EMcal_sim_iphi[nTowers];

  int   _b_EMcal_raw_n = 0;
  float _b_EMcal_raw_E[nTowers];
  float _b_EMcal_raw_eta[nTowers];
  float _b_EMcal_raw_phi[nTowers];
  
  int   _b_EMcal_calib_n = 0;
  float _b_EMcal_calib_E[nTowers];
  float _b_EMcal_calib_eta[nTowers];
  float _b_EMcal_calib_phi[nTowers];
  int   _b_EMcal_calib_ieta[nTowers];
  int   _b_EMcal_calib_iphi[nTowers];

  // Truth Tracks
  
  int _b_n_truth = 0;
  float _b_truthenergy[nTruth];
  float _b_trutheta[nTruth];
  float _b_truthphi[nTruth];
  float _b_truthpx[nTruth];
  float _b_truthpy[nTruth];
  float _b_truthpz[nTruth];
  float _b_truthpt[nTruth];
  float _b_truthp[nTruth];
  int _b_truthpid[nTruth];
  
  // Track variables
  int _b_n_tracks = 0;
  float _b_tr_px[nTracks];
  float _b_tr_py[nTracks];
  float _b_tr_pz[nTracks];
  float _b_tr_p[nTracks];
  float _b_tr_pt[nTracks];
  float _b_tr_phi[nTracks];
  float _b_tr_eta[nTracks];
  float _b_charge[nTracks];

  /// HEPMC Tree variables
  int m_num = 0;
  float m_truthenergy[nHEP];
  float m_trutheta[nHEP];
  float m_truthphi[nHEP];
  float m_truthpx[nHEP];
  float m_truthpy[nHEP];
  float m_truthpz[nHEP];
  float m_truthpt[nHEP];
  float m_truthp[nHEP];
  int m_truthpid[nHEP];

};

#endif
