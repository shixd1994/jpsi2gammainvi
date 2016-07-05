// -*- C++ -*-
//
//
// Description: psi(2S)-> J/psi pi+ pi-,J/psi -> rho pi 
//
// Original Author:  SHI Xiaodong <wherenpc@mail.ustc.edu.cn>
//         Created:  [2016-06-22 Wed 09:12] 
//         Inspired by Shi Xin's Psip2diPiRhoPi code 
//

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/LoadFactoryEntries.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"

#include "EventModel/EventHeader.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/Helix.h"
#include "VertexFit/WTrackParameter.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/KalmanKinematicFit.h"

#include "ParticleID/ParticleID.h"
#include "McTruth/McParticle.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McDecayModeSvc/IMcDecayModeSvc.h"


#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <vector>

class Psip2diPiRhoPi: public Algorithm {
  
public:
  Psip2diPiRhoPi(const std::string&, ISvcLocator*);
  ~Psip2diPiRhoPi(); 
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:
  // declare some cuts for charged tracks
	double m_vz0cut;
	double m_vr0cut;
	double m_cha_costheta_cut;
	double m_prob_pion_min;
	double m_min_emctime;
	double m_max_emctime;
	double m_costheta_barrel_max;
	double m_costheta_endcap_min;
	double m_costheta_endcap_max;
	double m_energy_barrel_min;
	double m_energy_endcap_min;
	double m_gammaCosCut;
	double m_photon_iso_angle_min;
	double m_chisq4C_cut;
	double m_pi0mass_cut;
	double m_Jpsimass_cut;

  //  MC truth info
	int m_numParticle;
	const int kMaxId; 
	int m_pdgid[500],m_motherindex[500];

	//const ecms_4p
	HepLorentzVector ECMS_4P;

  // output file
  std::string m_output_filename;
  //bool m_isMonteCarlo; 
  TFile* m_fout; 

	McDecayModeSvc* m_svc;
  // define Histograms
  TH1F* h_evtflw; 
  
  // define Trees
  TTree* m_tree;

  // common info 
  int m_run;
  int m_event;

  // charged tracks
  int m_ncharged;
  int m_nptrk;
  int m_nmtrk;
  
  // neutral tracks
	int m_ig1,m_ig2;
  int m_ngam;
	vector<double> *m_shower_energy;
	vector<double> *m_shower_phi;
	vector<double> *m_shower_theta;
	vector<double> *m_shower_e3;
	vector<double> *m_shower_e5;
	vector<int> *m_shower_status;
	vector<int> *m_shower_numHits;
	vector<double> *m_shower_secondMoment;
	vector<double> *m_shower_latMoment;
	vector<double> *m_shower_a42Moment;
	vector<double> *m_shower_a20Moment;
	vector<double> *m_shower_dE;

	// pion0 info
	double m_pi0_m;

	// Jpsi info
	double m_Jpsi_m;



  // functions		
	int psi2s_sample(int runNo);
	int psi2s_status(int runNo);
  void book_histogram();
  void book_tree(); 
	bool runStatusCheck();
  bool buildPsip2diPiRhoPi();
  bool saveGenInfo(SmartDataPtr<Event::McParticleCol>); 
  int selectChargedTracks(SmartDataPtr<EvtRecEvent>,
			  SmartDataPtr<EvtRecTrackCol>,
			  std::vector<int> &,
			  std::vector<int> &); 
  int PionPid(SmartDataPtr<EvtRecTrackCol>,
			      std::vector<int> &,
			      std::vector<int> &);
  void calcTrackPID(EvtRecTrackIterator,
		    double& ,
		    double&);
  int selectNeutralTracks(SmartDataPtr<EvtRecEvent>,
			  SmartDataPtr<EvtRecTrackCol>,
				std::vector<int> &);
  void clearNeutralTracks();
  void storeNeutralTracks(SmartDataPtr<EvtRecEvent>,
			  SmartDataPtr<EvtRecTrackCol>,
				std::vector<int> &);
  bool passVertexSelection(CLHEP::Hep3Vector,
			   RecMdcKalTrack* ); 
  CLHEP::Hep3Vector getOrigin();
  bool FourC(SmartDataPtr<EvtRecEvent>,
			  SmartDataPtr<EvtRecTrackCol>,
				std::vector<int> &,
				std::vector<int> &,
				std::vector<int> &,
			  int&,
			  int&,
				std::vector<HepLorentzVector> &,
				std::vector<HepLorentzVector> &,
				std::vector<HepLorentzVector> &,
				double&); 
  bool selectPion0(SmartDataPtr<EvtRecEvent>,
			  SmartDataPtr<EvtRecTrackCol>,
				std::vector<HepLorentzVector> &); 
  bool selectJpsi(SmartDataPtr<EvtRecEvent>,
			  SmartDataPtr<EvtRecTrackCol>,
				std::vector<HepLorentzVector> &,
				std::vector<HepLorentzVector> &); 
}; 


//
// module declare
//

DECLARE_ALGORITHM_FACTORY( Psip2diPiRhoPi )
DECLARE_FACTORY_ENTRIES( Psip2diPiRhoPiAlg ) {
  DECLARE_ALGORITHM(Psip2diPiRhoPi);
}

LOAD_FACTORY_ENTRIES( Psip2diPiRhoPiAlg )

//
// constants
//

const double PION_MASS = 0.139570;
const double PION0_MASS = 0.134977;
const double JPSI_MASS = 3.096916;


//
// member functions
//
  
Psip2diPiRhoPi::Psip2diPiRhoPi(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator),
	m_shower_energy(0),
	m_shower_phi(0),
	m_shower_theta(0),
	m_shower_e3(0),
	m_shower_e5(0),
	m_shower_status(0),
	m_shower_numHits(0),
	m_shower_secondMoment(0),
	m_shower_latMoment(0),
	m_shower_a42Moment(0),
	m_shower_a20Moment(0),
	m_shower_dE(0),
	kMaxId(500)
	//m_pdgid(0),
	//m_motherindex(0)
	{			
  declareProperty("OutputFileName", m_output_filename);
  //declareProperty("IsMonteCarlo", m_isMonteCarlo);
  declareProperty("Vr0cut", m_vr0cut=1.0);
  declareProperty("Vz0cut", m_vz0cut=10.0);
  declareProperty("ChaCosthetaCut", m_cha_costheta_cut=0.93);
  declareProperty("MinEstCut", m_min_emctime=0.0);
  declareProperty("MaxEstCut", m_max_emctime=14.0);
  declareProperty("GammaCosCut",  m_gammaCosCut= 0.93); 
  declareProperty("CosthetaBarrelMax", m_costheta_barrel_max=0.8);
  declareProperty("CosthetaEndcapMin", m_costheta_endcap_min=0.86);
  declareProperty("CosthetaEndcapMax", m_costheta_endcap_max=0.92);
  declareProperty("EnergyBarrelMin", m_energy_barrel_min=0.025); 
  declareProperty("EnergyEndcapMin", m_energy_endcap_min=0.050); 
  declareProperty("PhotonIsoAngleMin", m_photon_iso_angle_min=20);
  declareProperty("ProbPionMin", m_prob_pion_min=0.001);
  declareProperty("chisq4C", m_chisq4C_cut=200);
  declareProperty("pi0mass", m_pi0mass_cut=0.5);
  declareProperty("Jpsimass", m_Jpsimass_cut=0.2);
}


StatusCode Psip2diPiRhoPi::initialize(){
	ECMS_4P.set(3.686*sin(0.011),0.,0.,3.686);     
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << ">>>>>>> in initialize()" << endmsg;

  m_fout = new TFile(m_output_filename.c_str(), "RECREATE");
  m_fout->cd(); 

  book_histogram(); 
  book_tree(); 

  log << MSG::INFO << "successfully return from initialize()" <<endmsg;
  return StatusCode::SUCCESS;
}


StatusCode Psip2diPiRhoPi::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

  h_evtflw->Fill(0); // raw 
  SmartDataPtr<Event::EventHeader>eventHeader(eventSvc(),"/Event/EventHeader");
  if(!eventHeader) return StatusCode::FAILURE;

  m_run = eventHeader->runNumber();
  m_event = eventHeader->eventNumber();

	if(!runStatusCheck())  return StatusCode::SUCCESS;
  
  if(buildPsip2diPiRhoPi() == true) m_tree->Fill();

  return StatusCode::SUCCESS; 
}

StatusCode Psip2diPiRhoPi::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endmsg;

  m_fout->cd();
  m_tree->Write();
  h_evtflw->Write();
  m_fout->Close();
  
  return StatusCode::SUCCESS;
}

Psip2diPiRhoPi::~Psip2diPiRhoPi() {
}


void Psip2diPiRhoPi::book_histogram() {

  h_evtflw = new TH1F("hevtflw", "eventflow", 7, 0, 7);
  if (!h_evtflw) return;
  h_evtflw->GetXaxis()->SetBinLabel(1, "raw");
  h_evtflw->GetXaxis()->SetBinLabel(2, "N_{Good}=2&&Charge=0");
  h_evtflw->GetXaxis()->SetBinLabel(3, "Pid for pion");
  h_evtflw->GetXaxis()->SetBinLabel(4, "N_{#gamma}>=2");
  h_evtflw->GetXaxis()->SetBinLabel(5, "4C");
  h_evtflw->GetXaxis()->SetBinLabel(6, "Pion mass window"); 
  h_evtflw->GetXaxis()->SetBinLabel(7, "Jpsi mass window");
}


void Psip2diPiRhoPi::book_tree() {

  m_tree=new TTree("tree", "Psip2diPiRhoPi");
  if (!m_tree) return; 

  //commom info
  m_tree->Branch("run",&m_run,"run/I");
  m_tree->Branch("event",&m_event,"event/I");
	  
  //netual tracks
  m_tree->Branch("ig1",&m_ig1,"ig1/I");
  m_tree->Branch("ig2",&m_ig2,"ig2/I");
  m_tree->Branch("ngam",&m_ngam,"ngam/I");
	m_tree->Branch("sh_e", &m_shower_energy);
	m_tree->Branch("sh_phi", &m_shower_phi);
	m_tree->Branch("sh_the", &m_shower_theta);
	m_tree->Branch("sh_e3", &m_shower_e3);
	m_tree->Branch("sh_e5", &m_shower_e5);
	m_tree->Branch("sh_status", &m_shower_status);
	m_tree->Branch("sh_numHits", &m_shower_numHits);
	m_tree->Branch("sh_iiMoment", &m_shower_secondMoment);
	m_tree->Branch("sh_latMoment", &m_shower_latMoment);
	m_tree->Branch("sh_a42Moment", &m_shower_a42Moment);
	m_tree->Branch("sh_a20Moment", &m_shower_a20Moment);
	m_tree->Branch("sh_dE", &m_shower_dE);

	// pion0 info
  m_tree->Branch("pi0m",&m_pi0_m,"pi0m/D");
	
	// Jpsi info
  m_tree->Branch("Jpsim",&m_Jpsi_m,"Jpsim/D");
	
	//MC info 
	//if(!m_isMonteCarlo) return;
  m_tree->Branch("indexmc",&m_numParticle,"indexmc/I");
	m_tree->Branch("pdgid", m_pdgid, "pdgid[indexmc]/I");
	m_tree->Branch("motheridx", m_motherindex, "motheridx[indexmc]/I");
}

bool Psip2diPiRhoPi::runStatusCheck(){
	int runNo = m_run;
	if(runNo<0) runNo=0-runNo;
	int statuS=-1;
	statuS=psi2s_status(runNo);
	if(statuS==0) return false;
	else if(statuS!=3) cout<<"error!!"<<endl;
	return true;
}

bool Psip2diPiRhoPi::buildPsip2diPiRhoPi() {
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
 	if (m_run < 0) {
		if(!saveGenInfo(mcParticleCol)) return false;
		}

  SmartDataPtr<EvtRecEvent>evtRecEvent(eventSvc(),"/Event/EvtRec/EvtRecEvent");
  if(!evtRecEvent) return false;

  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  if(!evtRecTrkCol) return false;

	std::vector<int> iPGood,iMGood;
  selectChargedTracks(evtRecEvent, evtRecTrkCol, iPGood, iMGood);

	if (m_ncharged != 4) return false;		
	if (m_nptrk != 2) return false;
	if (m_nmtrk != 2) return false;
	h_evtflw->Fill(1); // N_{Good} = 4 && N_{P} = 2 && N_{M} = 2

	if (PionPid(evtRecTrkCol, iPGood, iMGood) !=4) return false;
	h_evtflw->Fill(2); // N_{Pi} = 4

  std::vector<int> iGam;
  iGam.clear();
  selectNeutralTracks(evtRecEvent, evtRecTrkCol, iGam);
  if (m_ngam < 2) return false;
  h_evtflw->Fill(3); // N_{#gamma} >= 2
	clearNeutralTracks();
  storeNeutralTracks(evtRecEvent, evtRecTrkCol, iGam);
	m_ig1=-1;
	m_ig2=-1;
	double chisq4C;
	std::vector<HepLorentzVector> pip4p(2),pim4p(2),gam4p(2);
	pip4p.clear();
	pim4p.clear();
	gam4p.clear();
	if(FourC(evtRecEvent, evtRecTrkCol,iGam ,iPGood, iMGood, 
		m_ig1, m_ig2, pip4p, pim4p, gam4p, chisq4C) == false) return false;
  h_evtflw->Fill(4); // 4C
	if(selectPion0(evtRecEvent, evtRecTrkCol ,gam4p) == false) return false;
  h_evtflw->Fill(5); //	Pion mass window cut 
	if(selectJpsi(evtRecEvent, evtRecTrkCol, pip4p, pim4p) == false) return false;
  h_evtflw->Fill(6); //	Jpsi mass window cut 
	return true;

}

bool Psip2diPiRhoPi::saveGenInfo(SmartDataPtr<Event::McParticleCol> mcParticleCol) {
  MsgStream log(msgSvc(), name());
	if (!mcParticleCol){
		std::cout << "Could not retrieve McParticelCol" << std::endl;
		return false;
	}
  IMcDecayModeSvc* i_svc;
  StatusCode sc_DecayModeSvc = service("McDecayModeSvc", i_svc);
  if ( sc_DecayModeSvc.isFailure() ){
      log << MSG::FATAL << "Could not load McDecayModeSvc!" << endreq;
      return false;
  }
  m_svc = dynamic_cast<McDecayModeSvc*>(i_svc);

	std::vector<int> pdgid, motherindex;
	pdgid.clear();
	motherindex.clear();
  Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
  for (; iter_mc != mcParticleCol->end(); iter_mc++){
    if ((*iter_mc)->primaryParticle()) continue;
  //  if (!(*iter_mc)->decayFromGenerator()) continue;

		if ((*iter_mc)->particleProperty()==100443){
			int mode = m_svc->extract(*iter_mc, pdgid, motherindex);
			m_numParticle = pdgid.size();
			for(int i = 0; i < m_numParticle; i++){
				m_pdgid[i]=(pdgid[i]);
				m_motherindex[i]=(motherindex[i]);
			}
		}
  } 
	return true;
}

CLHEP::Hep3Vector Psip2diPiRhoPi::getOrigin() {
  CLHEP::Hep3Vector xorigin(0,0,0);
  IVertexDbSvc*  vtxsvc;
  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  if(vtxsvc->isVertexValid()){
    double *dbv = vtxsvc->PrimaryVertex(); 
    xorigin.setX(dbv[0]);
    xorigin.setY(dbv[1]);
    xorigin.setZ(dbv[2]);
  }
  return xorigin; 
}

bool Psip2diPiRhoPi::passVertexSelection(CLHEP::Hep3Vector xorigin,
				    RecMdcKalTrack* mdcTrk ) {
  HepVector a = mdcTrk->helix();
  HepSymMatrix Ea = mdcTrk->err();
  HepPoint3D point0(0.,0.,0.);
  VFHelix helixip(point0,a,Ea);
  HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
  helixip.pivot(IP);
  HepVector vecipa = helixip.a();
  
  double vz0 = vecipa[3];
  double vr0 = vecipa[0];
  
  if(fabs(vz0) >= m_vz0cut) return false;
  if(fabs(vr0) >= m_vr0cut) return false;
  
  return true;
}

int Psip2diPiRhoPi::selectChargedTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
				   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
				   std::vector<int> & iPGood,
				   std::vector<int> & iMGood) {

  CLHEP::Hep3Vector xorigin = getOrigin(); 

  std::vector<int> iGood;
  iGood.clear();
  iPGood.clear();
  iMGood.clear();
  
  // loop through charged tracks 
  for(int i = 0; i < evtRecEvent->totalCharged(); i++){
    // get mdcTrk 
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;

    // Good Kalman Track 
    if(!(*itTrk)->isMdcKalTrackValid()) continue;

    if(!(*itTrk)->isMdcTrackValid()) continue;
    RecMdcKalTrack* mdcTrk = (*itTrk)->mdcKalTrack();

    // Good Vertex 
    if (!passVertexSelection(xorigin, mdcTrk) ) continue; 

    // Polar angle cut
    if(fabs(cos(mdcTrk->theta())) > m_cha_costheta_cut) continue;

    iGood.push_back((*itTrk)->trackId());
    if(mdcTrk->charge()>0) iPGood.push_back(i);
    if(mdcTrk->charge()<0) iMGood.push_back(i);

  } // end charged tracks

  m_ncharged = iGood.size();
  m_nptrk = iPGood.size();
  m_nmtrk = iMGood.size(); 
  
  return iGood.size(); 
}

int Psip2diPiRhoPi::PionPid(SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
										std::vector<int> & iPGood,
										std::vector<int> & iMGood){
  int npi= 0;
  bool evtflw_filled = false;
  
  for(unsigned int i1 = 0; i1 < iPGood.size(); i1++) {
    EvtRecTrackIterator itTrk_p = evtRecTrkCol->begin() + iPGood[i1];
    RecMdcTrack* mdcTrk_p = (*itTrk_p)->mdcTrack();
    if (mdcTrk_p->charge() < 0) continue; // only positive charged tracks

    // track PID
    double prob_pip, prob_kp;

    calcTrackPID(itTrk_p, prob_pip, prob_kp);  

    if(! (prob_pip > prob_kp &&
			prob_pip > m_prob_pion_min) ) continue;

    npi++;
  } 
  for(unsigned int i2 = 0; i2 < iMGood.size(); i2++) {
    EvtRecTrackIterator itTrk_m = evtRecTrkCol->begin() + iMGood[i2];
    RecMdcTrack* mdcTrk_m = (*itTrk_m)->mdcTrack();
    if (mdcTrk_m->charge() > 0) continue; // only negative charged tracks

    // track PID
		double prob_pim, prob_km; 
    calcTrackPID(itTrk_m, prob_pim, prob_km);
    // printf(">>> %f, %f, %f, %f \n", prob_pip, prob_kp, prob_pim, prob_km);

    if(! (prob_pim > prob_km &&
			prob_pim > m_prob_pion_min) ) continue;

    npi++;
  }

  return npi; 
}
			
void Psip2diPiRhoPi::calcTrackPID(EvtRecTrackIterator itTrk_p,
			     double& prob_pip,
			     double& prob_kp) {
  prob_pip = 999.; 
  prob_kp = 999.; 
  ParticleID * pidp = ParticleID::instance();
  pidp->init();
  pidp->setMethod(pidp->methodProbability());
  pidp->setChiMinCut(4);
  pidp->setRecTrack(*itTrk_p);
  // use PID sub-system
  pidp->usePidSys(pidp->useDedx() | pidp->useTof1() | pidp->useTof2());
  pidp->identify(pidp->onlyPionKaonProton());
  pidp->calculate();
  if(pidp->IsPidInfoValid()) {
    prob_pip = pidp->probPion();
    prob_kp  = pidp->probKaon();
  }
}

int Psip2diPiRhoPi::selectNeutralTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
				   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
					 std::vector<int> & iGam) {
  // loop through neutral tracks
  for(int i=evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) {
   // if (i > m_total_number_of_charged_max) break;

    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i ;
    if(!(*itTrk)->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    
    // TDC window
    if ( !(emcTrk->time() >= m_min_emctime && emcTrk->time() <= m_max_emctime) )
      continue; 

    // Energy threshold
    double abs_costheta(fabs(cos(emcTrk->theta())));
    bool barrel = (abs_costheta < m_costheta_barrel_max); 
    bool endcap = (abs_costheta > m_costheta_endcap_min
		   && abs_costheta < m_costheta_endcap_max);
    double eraw = emcTrk->energy();
    
    if ( !( (barrel && eraw > m_energy_barrel_min)
	    || (endcap && eraw > m_energy_endcap_min)))  continue; 

    // photon isolation: the opening angle between a candidate shower
    // and the closest charged track should be larger than 10 degree 
    CLHEP::Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

    // EMC costheta cut 
    double costhe = cos(emcpos.theta());
    if ( fabs(costhe) >= m_gammaCosCut) continue;
    
    // find the nearest charged track
    double dthe = 200.;
    double dphi = 200.;
    double dang = 200.; 
    for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
      EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
      if(!(*jtTrk)->isExtTrackValid()) continue;
      RecExtTrack *extTrk = (*jtTrk)->extTrack();
      if(extTrk->emcVolumeNumber() == -1) continue;
      CLHEP::Hep3Vector extpos = extTrk->emcPosition();
      double angd = extpos.angle(emcpos);
      double thed = extpos.theta() - emcpos.theta();
      double phid = extpos.deltaPhi(emcpos);
      thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;

      if(fabs(thed) < fabs(dthe)) dthe = thed;
      if(fabs(phid) < fabs(dphi)) dphi = phid;
      if(angd < dang) dang = angd;	    
    }

    if(dang>=200) continue;
    dthe = dthe * 180 / (CLHEP::pi);
    dphi = dphi * 180 / (CLHEP::pi);
    dang = dang * 180 / (CLHEP::pi);
    if (dang < m_photon_iso_angle_min ) continue; 

    iGam.push_back(i);
    //iGam.push_back((*itTrk)->trackId());
  } // end loop neutral tracks     

  m_ngam = iGam.size();

  return iGam.size(); 
}

void Psip2diPiRhoPi::clearNeutralTracks() {
	m_shower_energy->clear();
	m_shower_phi->clear();
	m_shower_theta->clear();
	m_shower_e3->clear();
	m_shower_e5->clear();
	m_shower_status->clear();
	m_shower_numHits->clear();
	m_shower_secondMoment->clear();
	m_shower_latMoment->clear();
	m_shower_a42Moment->clear();
	m_shower_a20Moment->clear();
	m_shower_dE->clear();
}
					

void Psip2diPiRhoPi::storeNeutralTracks(SmartDataPtr<EvtRecEvent> evtRecEvent,
				   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
					 std::vector<int> & iGam) {
	for(int i = 0; i < iGam.size(); i++){
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i] ;
    if(!(*itTrk)->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
		m_shower_energy->push_back(emcTrk->energy());
		m_shower_phi->push_back(emcTrk->phi());
		m_shower_theta->push_back(emcTrk->theta());
		m_shower_e3->push_back(emcTrk->e3x3());
		m_shower_e5->push_back(emcTrk->e5x5());
		m_shower_status->push_back(emcTrk->status());
		m_shower_numHits->push_back(emcTrk->numHits());
		m_shower_secondMoment->push_back(emcTrk->secondMoment());
		m_shower_latMoment->push_back(emcTrk->latMoment());
		m_shower_a42Moment->push_back(emcTrk->a42Moment());
		m_shower_a20Moment->push_back(emcTrk->a20Moment());
		m_shower_dE->push_back(emcTrk->dE());
	}
}

bool Psip2diPiRhoPi::FourC(SmartDataPtr<EvtRecEvent> evtRecEvent,
				   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
					 std::vector<int> & iGam,
					 std::vector<int> & iPGood,
					 std::vector<int> & iMGood,
					 int & m_ig1,
					 int & m_ig2,
					 std::vector<HepLorentzVector> & pip4p,
					 std::vector<HepLorentzVector> & pim4p,
					 std::vector<HepLorentzVector> & gam4p,
					 double & chisq4C) {
	KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
	double chisq = 9999;
	WTrackParameter wvpip1Trk, wvpip2Trk, wvpim1Trk, wvpim2Trk;
  RecMdcKalTrack *Trk_p1 = (*(evtRecTrkCol->begin() + iPGood[0]))->mdcKalTrack();
  RecMdcKalTrack *Trk_p2 = (*(evtRecTrkCol->begin() + iPGood[1]))->mdcKalTrack();
  RecMdcKalTrack *Trk_p3 = (*(evtRecTrkCol->begin() + iMGood[0]))->mdcKalTrack();
  RecMdcKalTrack *Trk_p4 = (*(evtRecTrkCol->begin() + iMGood[1]))->mdcKalTrack();
	wvpip1Trk = WTrackParameter(PION_MASS, Trk_p1->getZHelix(), Trk_p1->getZError());
	wvpip2Trk = WTrackParameter(PION_MASS, Trk_p2->getZHelix(), Trk_p2->getZError());
	wvpim1Trk = WTrackParameter(PION_MASS, Trk_p3->getZHelix(), Trk_p3->getZError());
	wvpim2Trk = WTrackParameter(PION_MASS, Trk_p4->getZHelix(), Trk_p4->getZError());

	for(int i = 0; i < m_ngam-1; i++){
		RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();
		for(int j = i+1; j < m_ngam; j++){
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[j]))->emcShower();

			kmfit->init();
			kmfit->AddTrack(0, wvpip1Trk); 
			kmfit->AddTrack(1, wvpip2Trk); 
			kmfit->AddTrack(2, wvpim1Trk); 
			kmfit->AddTrack(3, wvpim2Trk); 
			kmfit->AddTrack(4, 0.0, g1Trk); 
			kmfit->AddTrack(5, 0.0, g2Trk); 
			kmfit->AddFourMomentum(0, ECMS_4P);
			bool oksq = kmfit->Fit();
			if(oksq) {
			  double chi2 = kmfit->chisq();
			  if(chi2 < chisq) {
			    chisq = chi2;
			    m_ig1 = i;
			    m_ig2 = j;
					pip4p.push_back(kmfit->pfit(0));
					pip4p.push_back(kmfit->pfit(1));
					pim4p.push_back(kmfit->pfit(2));
					pim4p.push_back(kmfit->pfit(3));
					gam4p.push_back(kmfit->pfit(4));
					gam4p.push_back(kmfit->pfit(5));
			  }
			}
    }
  }
	if(chisq > m_chisq4C_cut) return false;
	return true;
}
	

bool Psip2diPiRhoPi::selectPion0(SmartDataPtr<EvtRecEvent> evtRecEvent,
				   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
					 std::vector<HepLorentzVector> & gam4p) {
	double mpi0;
	HepLorentzVector pi0 = gam4p[0] + gam4p[1];
	mpi0 = pi0.m();
	if(fabs(mpi0-PION0_MASS) > m_pi0mass_cut) return false;
	m_pi0_m = mpi0;
	return true;
}

bool Psip2diPiRhoPi::selectJpsi(SmartDataPtr<EvtRecEvent> evtRecEvent,
				   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol,
					 std::vector<HepLorentzVector> & pip4p,
					 std::vector<HepLorentzVector> & pim4p) {
	HepLorentzVector Jpsi4p;
	double mJpsi = 9999;
	double mJpsi0;
	for(int i = 0;i < pip4p.size(); i++){
		for(int j = 0;j < pim4p.size(); j++){
			Jpsi4p = ECMS_4P - pip4p[i] - pim4p[j];
			mJpsi0 = Jpsi4p.m();
			if(fabs(mJpsi0-JPSI_MASS) > m_Jpsimass_cut) continue;
			if(fabs(mJpsi0-JPSI_MASS) > fabs(mJpsi-JPSI_MASS)) continue;
			mJpsi = mJpsi0;
		}
	}
	if(mJpsi == 9999) return false;
	m_Jpsi_m = mJpsi;
	return true;
}

int Psip2diPiRhoPi::psi2s_sample(int runNo){
		//Isample=1,2,3,4,5,6,7 for good data, 
		// normal but with frequent TOF, TRG, DAQ errors
		if(runNo>=25338 && runNo<=25466)
			return 1; 
		//noisy MDC
		else if(runNo>=25467 && runNo<=25688)
			return 0;
		//normal
		else if(runNo>=25689 && runNo<=25750)
			return 2;
			//test run
			else if(runNo>=25751 && runNo<=25783)
				return 0;
				//Malter effect in some MDC cells, HV lowered by 30%, MUC partly disabled
				else if(runNo>=25784 && runNo<=26048)
					return 3;
					//test run
					else if(runNo>=26049 && runNo<=26165)
						return 0; 
						//Malter eff./TOF HV adjust /MUC dead
						else if(runNo>=26166 && runNo<=26274)
							return 0;
							//Malter effect in about 505 MDC inner layers, HV lowered, MUC partly disabled
							else if(runNo>=26275 && runNo<=26523)
								return 4;
								//add CO2, MDC training, machine study
								else if(runNo>=26524 && runNo<=26576)
									return 0;
									//add CO2, MDC training, adjust MDC HV
									else if(runNo>=26577 && runNo<=26632)
										return 0; 
										//CO2 in MDC, normal
										else if(runNo>=26633 && runNo<=26697)
											return 5; 
											//CO2 in MDC, high MDC noise level
											else if(runNo>=26698 && runNo<=26745)
												return 6;
												//test run
												else if(runNo>=26746 && runNo<=26763)
													return 0; 
													//CO2 in MDC, high MDC noise level, full MUC disabled
													else if(runNo>=26764 && runNo<=27090)
														return 7;
														else 
															return 0;
}

int Psip2diPiRhoPi::psi2s_status(int runNo){
		//Istat=0 for bad RUNs, Istat>0 for good RUNs, currently Istat =3 
		int Isample=psi2s_sample(runNo);
		if(Isample==0)
			return 0;
			else if(runNo==25341||runNo==25347||runNo==25348||(runNo>=25350&&runNo<=25352)||runNo==25372||runNo==25372||runNo==25374||runNo==25375||runNo==25382||(runNo>=25388&&runNo<=25391)||(runNo>=25403&&runNo<=25410)||runNo==25412||runNo==25428||(runNo>=25448&&runNo<=25454)||runNo==25456||runNo==25460||runNo==25698||runNo==25699||runNo==25702||runNo==25705||runNo==25711||(runNo>=25715&&runNo<=25719)||runNo==25740||runNo==25789||runNo==25790||runNo==25792||runNo==25793||runNo==25796||runNo==25808||runNo==25810||runNo==25813||runNo==25814||runNo==25823||runNo==25828||runNo==25832||(runNo>=25834&&runNo<=25841)||runNo==25845||runNo==25850||runNo==25852||(runNo>=25864&&runNo<=25872)||(runNo>=25879&&runNo<=25881)||(runNo>=25885&&runNo<=25907)||runNo==25917||runNo==25920||(runNo>=25926&&runNo<= 25932)||runNo==25933||runNo==25941||runNo==25945||(runNo>=25952&&runNo<=25954)||runNo==25973||runNo==25974||(runNo>=26003&&runNo<=26005)||runNo==26012||(runNo>=26014&&runNo<=26019)||runNo==26021||(runNo>=26023&&runNo<=26031)||runNo==26277||(runNo>=26279&&runNo<=26281)||runNo==26283||runNo==26284||runNo==26298||runNo==26301||runNo==26302||runNo==26346||runNo==26366||(runNo>=26379&&runNo<=26388)||runNo==26392||runNo==26397||runNo==26405||(runNo>=26444&&runNo<=26448)||runNo==26465||runNo==26468||runNo==26498||runNo==26507||runNo==26635||runNo==26650||runNo==26651||runNo==26666||runNo==26667||(runNo>=26690&&runNo<=26693)||runNo==26702||runNo==26776||runNo==26777||runNo==26780||(runNo>=26788&&runNo<=26793)||runNo==26797||(runNo>=26803&&runNo<=26811)||runNo==26819||runNo==26826||runNo==26827||(runNo>=26835&&runNo<=26839)||runNo==26842||runNo==26848||runNo==26849||(runNo>=26854&&runNo<=26859)||(runNo>=26876&&runNo<=26974)||runNo==26979||runNo==26984||(runNo>=26986&&runNo<=26988)||runNo==26992||(runNo>=27011&&runNo<=27013)||runNo==27024||(runNo>=27027&&runNo<=27031)||runNo==27070||runNo==27071||runNo==27089)
				return 0;
				else 
					return 3;
}
