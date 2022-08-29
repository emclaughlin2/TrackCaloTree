// Track Calorimeter Calibration TreeMaker
#include <fun4all/Fun4AllBase.h>
#include <iostream>
#include "TrackCaloTree.h"

// General F4A includes 
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4VtxPoint.h>

// Calorimeter includes
#include <calobase/RawTowerv2.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

// Tracking includes
#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

// HEPMC truth includes 
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

// ROOT includes 
#include "TLorentzVector.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH2.h"

TrackCaloTree::TrackCaloTree(const std::string &name) : SubsysReco("TRACK_CALO_TREE")
{
  _foutname = name;
  _hcal_sim_name = "TOWER_SIM_HCALOUT";
  _hcal_raw_name = "TOWER_RAW_HCALOUT";
  _hcal_calib_name = "TOWER_CALIB_HCALOUT";
  _hcalIN_raw_name = "TOWER_RAW_HCALIN";
  _hcalIN_calib_name = "TOWER_CALIB_HCALIN";
  _hcalIN_sim_name = "TOWER_SIM_HCALIN";
  _EMcal_sim_name = "TOWER_SIM_CEMC";
  _EMcal_raw_name = "TOWER_RAW_CEMC";
  _EMcal_calib_name = "TOWER_CALIB_CEMC";
  _track_name = "SvtxTrackMap";

}

int TrackCaloTree::reset_tree_vars() {

  _b_event = -99; 
  _b_tower_sim_n = -99;
  _b_tower_raw_n = -99;
  _b_tower_calib_n = -99;
  _b_EMcal_sim_n = -99;
  _b_EMcal_raw_n = -99;
  _b_EMcal_calib_n = -99;
  _b_hcalIN_sim_n = -99;
  _b_hcalIN_raw_n = -99;
  _b_hcalIN_calib_n = -99;
  m_num = -99;

  for (int i = 0; i <nTowers; i++){
    
    _b_tower_sim_E[i] = -99;
    _b_tower_sim_eta[i] = -99;
    _b_tower_sim_phi[i] = -99;
    _b_tower_sim_ieta[i] = -99;
    _b_tower_sim_iphi[i] = -99;
    
    _b_tower_raw_E[i] = -99;
    _b_tower_raw_eta[i] = -99;
    _b_tower_raw_phi[i] = -99;
    
    _b_EMcal_sim_E[i] = -99;
    _b_EMcal_sim_eta[i] = -99;
    _b_EMcal_sim_phi[i] = -99;
    _b_EMcal_sim_iphi[i] = -99;
    _b_EMcal_sim_ieta[i] = -99;

    _b_EMcal_raw_E[i] = -99;
    _b_EMcal_raw_eta[i] = -99;
    _b_EMcal_raw_phi[i] = -99;

    _b_hcalIN_sim_E[i] = -99;
    _b_hcalIN_sim_eta[i] = -99;
    _b_hcalIN_sim_phi[i] = -99;
    _b_hcalIN_sim_iphi[i] = -99;
    _b_hcalIN_sim_ieta[i] = -99;

    _b_hcalIN_raw_E[i] = -99;
    _b_hcalIN_raw_eta[i] = -99;
    _b_hcalIN_raw_phi[i] = -99;
    
    _b_tower_calib_E[i] = -99;
    _b_tower_calib_eta[i] = -99;
    _b_tower_calib_phi[i] = -99;
    _b_tower_calib_ieta[i] = -99;
    _b_tower_calib_iphi[i] = -99;

    _b_EMcal_calib_E[i] = -99;
    _b_EMcal_calib_eta[i] = -99;
    _b_EMcal_calib_phi[i] = -99;
    _b_EMcal_calib_iphi[i] = -99;
    _b_EMcal_calib_ieta[i] = -99;

    _b_hcalIN_calib_E[i] = -99;
    _b_hcalIN_calib_eta[i] = -99;
    _b_hcalIN_calib_phi[i] = -99;
    _b_hcalIN_calib_iphi[i] = -99;
    _b_hcalIN_calib_ieta[i] = -99;
  }
  
  _b_n_truth = -99;
  for (int i = 0; i < nTruth; i++) {
    _b_truthenergy[i] = -99;
    _b_trutheta[i] = -99;
    _b_truthphi[i] = -99;
    _b_truthpx[i] = -99;
    _b_truthpy[i] = -99;
    _b_truthpz[i] = -99;
    _b_truthpt[i] = -99;
    _b_truthp[i] = -99;
    _b_truthpid[i] = -99;
  }
  

  for (int i = 0; i < nHEP; i++) {
    m_truthenergy[i] = -99;
    m_trutheta[i] = -99;
    m_truthphi[i] = -99;
    m_truthpx[i] = -99;
    m_truthpy[i] = -99;
    m_truthpz[i] = -99;
    m_truthpt[i] = -99;
    m_truthp[i] = -99;
    m_truthpid[i] = -99;
  }

  if (include_tracks) {
    _b_n_tracks = -99;
    for (int i = 0; i < nTracks; i++) {
      _b_tr_py[i] = -99;
      _b_tr_px[i] = -99;
      _b_tr_pz[i] = -99;
      _b_tr_p[i] = -99;
      _b_tr_pt[i] = -99;
      _b_tr_phi[i] = -99;
      _b_tr_eta[i] = -99;
      _b_charge[i] = -99;
    }
  }

  return 1;
}

int TrackCaloTree::GetNodes(PHCompositeNode *topNode) {

  if (_debug) std::cout<<"GettingNodes..."<<std::endl;
  _geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if(!_geomOH) std::cout<<"No TOWERGeOM_HCALOUT"<<std::endl;

  _geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  if(!_geomIH) std::cout<<"No TOWERGeOM_HCALIN"<<std::endl;

  _geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  if(!_geomEM) std::cout<<"No TOWERGeOM_CEMC"<<std::endl;
  
  _towersSimOH = findNode::getClass<RawTowerContainer>(topNode, _hcal_sim_name);
  if (!_towersSimOH) std::cout<<"No TOWER_SIM_HCALOUT Node"<<std::endl;

  _towersRawOH = findNode::getClass<RawTowerContainer>(topNode, _hcal_raw_name);
  if (!_towersRawOH) std::cout<<"No TOWER_RAW_HCALOUT Node"<<std::endl;
  
  _towersCalibOH = findNode::getClass<RawTowerContainer>(topNode, _hcal_calib_name);
  if (!_towersCalibOH) std::cout<<"No TOWER_CALIB_HCALOUT Node"<<std::endl;
  
  _towersRawIH = findNode::getClass<RawTowerContainer>(topNode, _hcalIN_raw_name);
  if (!_towersRawIH) std::cout<<"No TOWER_RAW_HCALIN Node"<<std::endl;

  _towersRawEM = findNode::getClass<RawTowerContainer>(topNode, _EMcal_raw_name);
  if (!_towersRawEM) std::cout<<"No TOWER_RAW_CEMC Node"<<std::endl;
  
  _towersCalibIH = findNode::getClass<RawTowerContainer>(topNode, _hcalIN_calib_name);
  if (!_towersCalibIH) std::cout<<"No TOWER_CALIB_HCALIN Node"<<std::endl;

  _towersCalibEM = findNode::getClass<RawTowerContainer>(topNode, _EMcal_calib_name);
  if (!_towersCalibEM) std::cout<<"No TOWER_CALIB_CEMC Node"<<std::endl;
  
  _towersSimIH = findNode::getClass<RawTowerContainer>(topNode, _hcalIN_sim_name);
  if (!_towersSimIH) std::cout<<"No TOWER_SIM_HCALIN Node"<<std::endl;

  _towersSimEM = findNode::getClass<RawTowerContainer>(topNode, _EMcal_sim_name);
  if (!_towersSimEM) std::cout<<"No TOWER_SIM_CEMC Node"<<std::endl;
  
  if (include_tracks) {
    trackmap = findNode::getClass<SvtxTrackMap>(topNode, _track_name);
    if (!trackmap) std::cout <<"SvtxTrackMap node is missing, can't collect tracks"<< std::endl;
    if (track_debug && trackmap->size() == 0) std::cout << "SvtxTrackMap has zero entries" << std::endl;
  }

  truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthinfo) std::cout << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"<< std::endl;

  hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!hepmceventmap) std::cout <<"HEPMC event map node is missing, can't collected HEPMC truth particles"<< std::endl;

  return 1;
}

int TrackCaloTree::Init(PHCompositeNode *topNode) {
  _ievent = 0;
  _b_event = -1;

  if (_debug) std::cout<<"Initiating..."<<std::endl;

  reset_tree_vars();

  _file = new TFile(_foutname.c_str(), "RECREATE");
  
  _tree = new TTree("T", "keep on giving tree");
  
  _tree->Branch("EMcal_calib_n",   &_b_EMcal_calib_n,     "EMcal_calib_n/I");
  _tree->Branch("EMcal_calib_E",    _b_EMcal_calib_E,     "EMcal_calib_E[EMcal_calib_n]/F");
  _tree->Branch("EMcal_calib_eta",  _b_EMcal_calib_eta,   "EMcal_calib_eta[EMcal_calib_n]/F");
  _tree->Branch("EMcal_calib_phi",  _b_EMcal_calib_phi,   "EMcal_calib_phi[EMcal_calib_n]/F");
  _tree->Branch("EMcal_calib_iphi",  _b_EMcal_calib_iphi,   "EMcal_calib_iphi[EMcal_calib_n]/I");
  _tree->Branch("EMcal_calib_ieta",  _b_EMcal_calib_ieta,   "EMcal_calib_ieta[EMcal_calib_n]/I");

  _tree->Branch("hcalIN_calib_n",   &_b_hcalIN_calib_n,     "hcalIN_calib_n/I");
  _tree->Branch("hcalIN_calib_E",    _b_hcalIN_calib_E,     "hcalIN_calib_E[hcalIN_calib_n]/F");
  _tree->Branch("hcalIN_calib_eta",  _b_hcalIN_calib_eta,   "hcalIN_calib_eta[hcalIN_calib_n]/F");
  _tree->Branch("hcalIN_calib_phi",  _b_hcalIN_calib_phi,   "hcalIN_calib_phi[hcalIN_calib_n]/F");
  _tree->Branch("hcalIN_calib_iphi",  _b_hcalIN_calib_iphi,   "hcalIN_calib_iphi[hcalIN_calib_n]/I");
  _tree->Branch("hcalIN_calib_ieta",  _b_hcalIN_calib_ieta,   "hcalIN_calib_ieta[hcalIN_calib_n]/I");

  _tree->Branch("tower_calib_n",&_b_tower_calib_n, "tower_calib_n/I");
  _tree->Branch("tower_calib_E",_b_tower_calib_E, "tower_calib_E[tower_calib_n]/F");
  _tree->Branch("tower_calib_eta",_b_tower_calib_eta, "tower_calib_eta[tower_calib_n]/F");
  _tree->Branch("tower_calib_phi",_b_tower_calib_phi, "tower_calib_phi[tower_calib_n]/F");
  _tree->Branch("tower_calib_ieta",_b_tower_calib_ieta, "tower_calib_ieta[tower_calib_n]/I");
  _tree->Branch("tower_calib_iphi",_b_tower_calib_iphi, "tower_calib_iphi[tower_calib_n]/I");
  
  if (include_tracks) {
  _tree->Branch("n_tracks",&_b_n_tracks,"n_tracks/I");
  _tree->Branch("tr_py",_b_tr_py,"tr_py[n_tracks]/F");
  _tree->Branch("tr_px",_b_tr_px,"tr_px[n_tracks]/F");
  _tree->Branch("tr_pz",_b_tr_pz,"tr_pz[n_tracks]/F");
  _tree->Branch("tr_p",_b_tr_p,"tr_p[n_tracks]/F");
  _tree->Branch("tr_pt",_b_tr_pt,"tr_pt[n_tracks]/F");
  _tree->Branch("tr_phi",_b_tr_phi,"tr_phi[n_tracks]/F");
  _tree->Branch("tr_eta",_b_tr_eta,"tr_eta[n_tracks]/F");
  _tree->Branch("charge",_b_charge,"charge[n_tracks]/F");
  }
  
  _tree->Branch("n_truth",&_b_n_truth,"n_truth/I");
  _tree->Branch("truthenergy",_b_truthenergy,"truthenergy[n_truth]/F");
  _tree->Branch("trutheta",_b_trutheta,"trutheta[n_truth]/F");
  _tree->Branch("truthphi",_b_truthphi,"truthphi[n_truth]/F");
  _tree->Branch("truthpx",_b_truthpx,"truthpx[n_truth]/F");
  _tree->Branch("truthpy",_b_truthpy,"truthpy[n_truth]/F");
  _tree->Branch("truthpz",_b_truthpz,"truthpz[n_truth]/F");
  _tree->Branch("truthpt",_b_truthpt,"truthpt[n_truth]/F");
  _tree->Branch("truthp",_b_truthp,"truthp[n_truth]/F");
  _tree->Branch("truthpid",_b_truthpid,"truthpid[n_truth]/I");
  
  _tree->Branch("m_num",&m_num,"m_num/I");
  _tree->Branch("m_truthenergy",m_truthenergy,"m_truthenergy[m_num]/F");
  _tree->Branch("m_trutheta",m_trutheta,"m_trutheta[m_num]/F");
  _tree->Branch("m_truthphi",m_truthphi,"m_truthphi[m_num]/F");
  _tree->Branch("m_truthp",m_truthp,"m_truthp[m_num]/F");
  _tree->Branch("m_truthpt",m_truthpt,"m_truthpt[m_num]/F");
  _tree->Branch("m_truthpx",m_truthpx,"m_truthpx[m_num]/F");
  _tree->Branch("m_truthpy",m_truthpy,"m_truthpy[m_num]/F");
  _tree->Branch("m_truthpz",m_truthpz,"m_truthpz[m_num]/F");
  _tree->Branch("m_truthpid",m_truthpid,"m_truthpid[m_num]/F");

  _tree->Branch("EMcal_sim_n",   &_b_EMcal_sim_n,     "EMcal_sim_n/I");
  _tree->Branch("EMcal_sim_E",    _b_EMcal_sim_E,     "EMcal_sim_E[EMcal_sim_n]/F");
  _tree->Branch("EMcal_sim_eta",  _b_EMcal_sim_eta,   "EMcal_sim_eta[EMcal_sim_n]/F");
  _tree->Branch("EMcal_sim_phi",  _b_EMcal_sim_phi,   "EMcal_sim_phi[EMcal_sim_n]/F");
  _tree->Branch("EMcal_sim_iphi",  _b_EMcal_sim_iphi,   "EMcal_sim_iphi[EMcal_sim_n]/I");
  _tree->Branch("EMcal_sim_ieta",  _b_EMcal_sim_ieta,   "EMcal_sim_ieta[EMcal_sim_n]/I");

  _tree->Branch("hcalIN_sim_n",   &_b_hcalIN_sim_n,     "hcalIN_sim_n/I");
  _tree->Branch("hcalIN_sim_E",    _b_hcalIN_sim_E,     "hcalIN_sim_E[hcalIN_sim_n]/F");
  _tree->Branch("hcalIN_sim_eta",  _b_hcalIN_sim_eta,   "hcalIN_sim_eta[hcalIN_sim_n]/F");
  _tree->Branch("hcalIN_sim_phi",  _b_hcalIN_sim_phi,   "hcalIN_sim_phi[hcalIN_sim_n]/F");
  _tree->Branch("hcalIN_sim_iphi",  _b_hcalIN_sim_iphi,   "hcalIN_sim_iphi[hcalIN_sim_n]/I");
  _tree->Branch("hcalIN_sim_ieta",  _b_hcalIN_sim_ieta,   "hcalIN_sim_ieta[hcalIN_sim_n]/I");

  _tree->Branch("tower_sim_n",&_b_tower_sim_n, "tower_sim_n/I");
  _tree->Branch("tower_sim_E",_b_tower_sim_E, "tower_sim_E[tower_sim_n]/F");
  _tree->Branch("tower_sim_eta",_b_tower_sim_eta, "tower_sim_eta[tower_sim_n]/F");
  _tree->Branch("tower_sim_phi",_b_tower_sim_phi, "tower_sim_phi[tower_sim_n]/F");
  _tree->Branch("tower_sim_ieta",_b_tower_sim_ieta, "tower_sim_ieta[tower_sim_n]/I");
  _tree->Branch("tower_sim_iphi",_b_tower_sim_iphi, "tower_sim_iphi[tower_sim_n]/I");
  /*
  _tree->Branch("EMcal_raw_n",   &_b_EMcal_raw_n,     "EMcal_raw_n/I");
  _tree->Branch("EMcal_raw_E",    _b_EMcal_raw_E,     "EMcal_raw_E[EMcal_raw_n]/F");
  _tree->Branch("EMcal_raw_eta",  _b_EMcal_raw_eta,   "EMcal_raw_eta[EMcal_raw_n]/F");
  _tree->Branch("EMcal_raw_phi",  _b_EMcal_raw_phi,   "EMcal_raw_phi[EMcal_raw_n]/F");

  _tree->Branch("hcalIN_raw_n",   &_b_hcalIN_raw_n,     "hcalIN_raw_n/I");
  _tree->Branch("hcalIN_raw_E",    _b_hcalIN_raw_E,     "hcalIN_raw_E[hcalIN_raw_n]/F");
  _tree->Branch("hcalIN_raw_eta",  _b_hcalIN_raw_eta,   "hcalIN_raw_eta[hcalIN_raw_n]/F");
  _tree->Branch("hcalIN_raw_phi",  _b_hcalIN_raw_phi,   "hcalIN_raw_phi[hcalIN_raw_n]/F");

  _tree->Branch("tower_raw_n",&_b_tower_raw_n, "tower_raw_n/I");
  _tree->Branch("tower_raw_E",_b_tower_raw_E, "tower_raw_E[tower_raw_n]/F");
  _tree->Branch("tower_raw_eta",_b_tower_raw_eta, "tower_raw_eta[tower_raw_n]/F");
  _tree->Branch("tower_raw_phi",_b_tower_raw_phi, "tower_raw_phi[tower_raw_n]/F");
  */
  return 0;
}

int TrackCaloTree::process_event(PHCompositeNode *topNode) {

  GetNodes(topNode);
  
  _b_event = _ievent;
  if (_ievent %5000==0) std::cout<<"Event: "<<_ievent<<std::endl;


  if (_debug) std::cout<<"Processing Event: "<< _b_event<<std::endl;
  if (_debug) std::cout << "hello";
  
  //////////////////////////////
  // OUTER HCAL
  //////////////////////////////
  
  _b_tower_raw_n = 0;
  RawTowerContainer::ConstRange begin_end_raw = _towersRawOH->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_raw.first; rtiter != begin_end_raw.second; ++rtiter) {
    if (_debug) std::cout << "looking at wother " << _b_tower_raw_n << std::endl;
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = _geomOH->get_tower_geometry(tower->get_key());
    if (_debug) std::cout << "looking at wother geometry " << _b_tower_raw_n << std::endl;
      _b_tower_raw_E[ _b_tower_raw_n ]   =tower->get_energy();
      _b_tower_raw_eta[ _b_tower_raw_n ] =tower_geom->get_eta();
      _b_tower_raw_phi[ _b_tower_raw_n ] =tower_geom->get_phi();
      _b_tower_raw_n++;
      if (_debug) std::cout << "going to tower " << _b_tower_raw_n << std::endl; 
  }

  if (_debug) std::cout<<"Got raw n: "<< _b_tower_raw_n<<std::endl;
  _b_tower_sim_n=0;
  RawTowerContainer::ConstRange begin_end_sim = _towersSimOH->getTowers();
  if (_debug) std::cout<<"Got the iterator"<<std::endl;

  for (RawTowerContainer::ConstIterator rtiter = begin_end_sim.first; rtiter != begin_end_sim.second; ++rtiter) {
    if ((_b_tower_sim_n%10==0)&&(_debug)) std::cout<<"At sim tower: "<< _b_tower_sim_n<<std::endl;

    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = _geomOH->get_tower_geometry(tower->get_key());
      _b_tower_sim_E   [ _b_tower_sim_n ] = tower->get_energy();
      _b_tower_sim_eta [ _b_tower_sim_n ] = tower_geom->get_eta();
      _b_tower_sim_phi [ _b_tower_sim_n ] = tower_geom->get_phi();
      _b_tower_sim_ieta[ _b_tower_sim_n ] = tower_geom->get_bineta();
      _b_tower_sim_iphi[ _b_tower_sim_n ] = tower_geom->get_binphi();
      _b_tower_sim_n++;
    if(_b_tower_sim_n >= nTowers){
      std::cout << __FILE__ << " ERROR: _b_tower_sim_n has hit cap of " << nTowers << "!!!" << std::endl;
    }
    
  }
  
  if (_debug) std::cout<<"Got sims n: "<< _b_tower_sim_n<<std::endl;
  _b_tower_calib_n = 0;
  
  RawTowerContainer::ConstRange begin_end_calib = _towersCalibOH->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_calib.first; rtiter != begin_end_calib.second; ++rtiter) {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = _geomOH->get_tower_geometry(tower->get_key());
      _b_tower_calib_E[ _b_tower_calib_n ] = tower->get_energy();
      _b_tower_calib_eta[ _b_tower_calib_n ] = tower_geom->get_eta();
      _b_tower_calib_phi[ _b_tower_calib_n ] = tower_geom->get_phi();
      _b_tower_calib_ieta [ _b_tower_calib_n ] = tower_geom->get_bineta();
      _b_tower_calib_iphi [ _b_tower_calib_n ] = tower_geom->get_binphi();
      _b_tower_calib_n++;

    if(_b_tower_calib_n >= nTowers){
      std::cout << __FILE__ << " ERROR: _b_tower_calib_n has hit cap of " << nTowers << "!!!" << std::endl;
    }
    
  }
  if (_debug) std::cout<<"Got calib n: "<< _b_tower_calib_n<<std::endl;


  //////////////////////
  // INNER HCAL
  /////////////////////
  
  _b_hcalIN_sim_n = 0;
  RawTowerContainer::ConstRange begin_end_simIN = _towersSimIH->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_simIN.first; rtiter != begin_end_simIN.second; ++rtiter) {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = _geomIH->get_tower_geometry(tower->get_key());
      _b_hcalIN_sim_E    [ _b_hcalIN_sim_n ] = tower->get_energy();
      _b_hcalIN_sim_eta  [ _b_hcalIN_sim_n ] = tower_geom->get_eta();
      _b_hcalIN_sim_phi  [ _b_hcalIN_sim_n ] = tower_geom->get_phi();
      _b_hcalIN_sim_ieta [ _b_hcalIN_sim_n ] = tower_geom->get_bineta();
      _b_hcalIN_sim_iphi [ _b_hcalIN_sim_n ] = tower_geom->get_binphi();
      _b_hcalIN_sim_n++;

  }

  _b_hcalIN_raw_n = 0;
  RawTowerContainer::ConstRange begin_end_rawIN = _towersRawIH->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_rawIN.first; rtiter != begin_end_rawIN.second; ++rtiter) {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = _geomIH->get_tower_geometry(tower->get_key());
      _b_hcalIN_raw_E    [ _b_hcalIN_raw_n ] = tower->get_energy();
      _b_hcalIN_raw_eta  [ _b_hcalIN_raw_n ] = tower_geom->get_eta();
      _b_hcalIN_raw_phi  [ _b_hcalIN_raw_n ] = tower_geom->get_phi();
      _b_hcalIN_raw_n++;

  }
  
  _b_hcalIN_calib_n = 0;
  RawTowerContainer::ConstRange begin_end_calibIN = _towersCalibIH->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_calibIN.first; rtiter != begin_end_calibIN.second; ++rtiter) {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = _geomIH->get_tower_geometry(tower->get_key());
      _b_hcalIN_calib_E    [ _b_hcalIN_calib_n ] = tower->get_energy();
      _b_hcalIN_calib_eta  [ _b_hcalIN_calib_n ] = tower_geom->get_eta();
      _b_hcalIN_calib_phi  [ _b_hcalIN_calib_n ] = tower_geom->get_phi();
      _b_hcalIN_calib_ieta [ _b_hcalIN_calib_n ] = tower_geom->get_bineta();
      _b_hcalIN_calib_iphi [ _b_hcalIN_calib_n ] = tower_geom->get_binphi();
      _b_hcalIN_calib_n++;

  }

  //////////////////////
  // EM CAL
  //////////////////////
  
  _b_EMcal_sim_n = 0;
  RawTowerContainer::ConstRange begin_end_simEM = _towersSimEM->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_simEM.first; rtiter != begin_end_simEM.second; ++rtiter) {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = _geomEM->get_tower_geometry(tower->get_key());
      _b_EMcal_sim_E    [ _b_EMcal_sim_n ] = tower->get_energy();
      _b_EMcal_sim_eta  [ _b_EMcal_sim_n ] = tower_geom->get_eta();
      _b_EMcal_sim_phi  [ _b_EMcal_sim_n ] = tower_geom->get_phi();
      _b_EMcal_sim_ieta [ _b_EMcal_sim_n ] = tower_geom->get_bineta();
      _b_EMcal_sim_iphi [ _b_EMcal_sim_n ] = tower_geom->get_binphi();
      _b_EMcal_sim_n++;

  }

  _b_EMcal_raw_n = 0;
  RawTowerContainer::ConstRange begin_end_rawEM = _towersRawEM->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_rawEM.first; rtiter != begin_end_rawEM.second; ++rtiter) {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = _geomEM->get_tower_geometry(tower->get_key());
      _b_EMcal_raw_E    [ _b_EMcal_raw_n ] = tower->get_energy();
      _b_EMcal_raw_eta  [ _b_EMcal_raw_n ] = tower_geom->get_eta();
      _b_EMcal_raw_phi  [ _b_EMcal_raw_n ] = tower_geom->get_phi();
      _b_EMcal_raw_n++;

  }
  
  _b_EMcal_calib_n = 0;
  RawTowerContainer::ConstRange begin_end_calibEM = _towersCalibEM->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_calibEM.first; rtiter != begin_end_calibEM.second; ++rtiter) {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = _geomEM->get_tower_geometry(tower->get_key());
      _b_EMcal_calib_E    [ _b_EMcal_calib_n ] = tower->get_energy();
      _b_EMcal_calib_eta  [ _b_EMcal_calib_n ] = tower_geom->get_eta();
      _b_EMcal_calib_phi  [ _b_EMcal_calib_n ] = tower_geom->get_phi();
      _b_EMcal_calib_ieta [ _b_EMcal_calib_n ] = tower_geom->get_bineta();
      _b_EMcal_calib_iphi [ _b_EMcal_calib_n ] = tower_geom->get_binphi();
      _b_EMcal_calib_n++;

  }

  ///////////////////////////
  // RECONSTRUCTED TRACKS
  //////////////////////////
  if (include_tracks) {
    if(track_debug) {
      std::cout << "getting tracks from map" << std::endl;
      std::cout << "track map size is " << trackmap->size() << std::endl;
      std::cout << "eta values from nodetree" << std::endl;
    }

    _b_n_tracks = 0;
    for (SvtxTrackMap::Iter iter = trackmap->begin();
      iter != trackmap->end(); ++iter) {
      SvtxTrack *track = iter->second;

      if (track_debug) std::cout << track->get_eta() << " ";
      _b_tr_px[_b_n_tracks] = track->get_px();
      _b_tr_py[_b_n_tracks] = track->get_py();
      _b_tr_pz[_b_n_tracks] = track->get_pz();
      _b_tr_p[_b_n_tracks] = sqrt(_b_tr_px[_b_n_tracks] * _b_tr_px[_b_n_tracks]
                                + _b_tr_py[_b_n_tracks] * _b_tr_py[_b_n_tracks] 
                                + _b_tr_pz[_b_n_tracks] * _b_tr_pz[_b_n_tracks]);
      _b_tr_pt[_b_n_tracks] = sqrt(_b_tr_px[_b_n_tracks] * _b_tr_px[_b_n_tracks] 
                                + _b_tr_py[_b_n_tracks] * _b_tr_py[_b_n_tracks]);
      _b_tr_phi[_b_n_tracks] = track->get_phi();
      _b_tr_eta[_b_n_tracks] = track->get_eta();
      _b_charge[_b_n_tracks] = track->get_charge();
      _b_n_tracks++;
      if (_b_n_tracks >= 10000) break;
    }
    if (track_debug) { 
      std::cout << std::endl;
      std::cout << "eta values saved to tree variables" << std::endl;
      for (int it = 0; it < _b_n_tracks; it++) {
        std::cout << _b_tr_eta[it] << " ";
      }
      std::cout << std::endl;
    }
  }
  
  ////////////////////
  // TRUTH INFO
  ////////////////////
  
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  //PHG4TruthInfoContainer::Range range = truthinfo->GetParticleRange();
  /// Loop over the G4 truth (stable) particles
  _b_n_truth = 0;
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second; ++iter) {
    /// Get this truth particle
    const PHG4Particle *truth = iter->second;
    if (!truthinfo->is_primary(truth)) continue;

    /// Get this particles momentum, etc.
    _b_truthpx[_b_n_truth] = truth->get_px();
    _b_truthpy[_b_n_truth] = truth->get_py();
    _b_truthpz[_b_n_truth] = truth->get_pz();
    _b_truthp[_b_n_truth] = sqrt(_b_truthpx[_b_n_truth] * _b_truthpx[_b_n_truth]
                           + _b_truthpy[_b_n_truth] * _b_truthpy[_b_n_truth]
                           + _b_truthpz[_b_n_truth] * _b_truthpz[_b_n_truth]);
    _b_truthpt[_b_n_truth] = sqrt(_b_truthpx[_b_n_truth] * _b_truthpx[_b_n_truth]
                            + _b_truthpy[_b_n_truth] * _b_truthpy[_b_n_truth]);
    _b_truthenergy[_b_n_truth] = truth->get_e();
    _b_truthpt[_b_n_truth] = sqrt(_b_truthpx[_b_n_truth] * _b_truthpx[_b_n_truth] 
                            + _b_truthpy[_b_n_truth] * _b_truthpy[_b_n_truth]);
    _b_truthphi[_b_n_truth] = atan(_b_truthpy[_b_n_truth] / _b_truthpx[_b_n_truth]);
    _b_trutheta[_b_n_truth] = atanh(_b_truthpz[_b_n_truth] / _b_truthenergy[_b_n_truth]);

    /// Check for nansx
    if (_b_trutheta[_b_n_truth] != _b_trutheta[_b_n_truth])
      _b_trutheta[_b_n_truth] = -99;
    _b_truthpid[_b_n_truth] = truth->get_pid();
    //_b_truthpid[_b_n_truth] = truth->get_barcode();
    _b_n_truth++;
    if (_b_n_truth >= 15000) break;
  }
  

  /// You can iterate over the number of events in a hepmc event
  /// for pile up events where you have multiple hard scatterings per bunch crossing
  for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin();
       eventIter != hepmceventmap->end();
       ++eventIter)
  {
    /// Get the event
    PHHepMCGenEvent *hepmcevent = eventIter->second;
    std::cout << "embedding id for HEPMC event" << hepmcevent->get_embedding_id() << std::endl;
    // To fill TTree, require that the event be the primary event (embedding_id > 0)
    if (hepmcevent && hepmcevent->get_embedding_id() == 0)
    {
      /// Get the event characteristics, inherited from HepMC classes
      HepMC::GenEvent *truthevent = hepmcevent->getEvent();
      if (!truthevent)
      {
        std::cout << PHWHERE
             << "no evt pointer under phhepmvgeneventmap found "
             << std::endl;
      }

      if (Verbosity() > 2)
      {
        std::cout << " Iterating over an event" << std::endl;
      }
      /// Loop over all the truth particles and get their information
      m_num = 0;
      for (HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin();
           iter != truthevent->particles_end();
           ++iter)
      {
        /// Get each pythia particle characteristics
        m_truthenergy[m_num] = (*iter)->momentum().e();
        m_truthpid[m_num] = (*iter)->pdg_id();

        m_trutheta[m_num] = (*iter)->momentum().pseudoRapidity();
        m_truthphi[m_num] = (*iter)->momentum().phi();
        m_truthpx[m_num] = (*iter)->momentum().px();
        m_truthpy[m_num] = (*iter)->momentum().py();
        m_truthpz[m_num] = (*iter)->momentum().pz();
        m_truthpt[m_num] = sqrt(m_truthpx[m_num] * m_truthpx[m_num] 
                            + m_truthpy[m_num] * m_truthpy[m_num]);
        /// Fill the truth tree
        m_num++;
        if (m_num >= 10000) break;
      }
      break;
    } 
  }

  _file->cd();
  if (track_debug) { 
    std::cout << std::endl;
    std::cout << "eta values in tree variables just prior to tree fill" << std::endl;
    for (int it = 0; it < _b_n_tracks; it++) {
      std::cout << _b_tr_eta[it] << " ";
    }
    std::cout << std::endl;
  }
  _tree->Fill();
  if (track_debug) { 
    std::cout << std::endl;
    std::cout << "eta values in tree variables just after to tree fill" << std::endl;
    for (int it = 0; it < _b_n_tracks; it++) {
      std::cout << _b_tr_eta[it] << " ";
    }
    std::cout << std::endl;
  }
  _ievent++;
  
  return 0;
}

int TrackCaloTree::End(PHCompositeNode *topNode)
{
  if (_debug) std::cout<<"Writing File"<<std::endl;
  _file->cd();
  if (track_debug) { 
    std::cout << std::endl;
    std::cout << "eta values in tree variables just prior to writing to file" << std::endl;
    for (int it = 0; it < _b_n_tracks; it++) {
      std::cout << _b_tr_eta[it] << " ";
    }
    std::cout << std::endl;
  }
  _file->Write();
  if (track_debug) { 
    std::cout << std::endl;
    std::cout << "eta values in tree variables just after writing to file" << std::endl;
    for (int it = 0; it < _b_n_tracks; it++) {
      std::cout << _b_tr_eta[it] << " ";
    }
    std::cout << std::endl;
  }
  _file->Close();

  return 0;
}
